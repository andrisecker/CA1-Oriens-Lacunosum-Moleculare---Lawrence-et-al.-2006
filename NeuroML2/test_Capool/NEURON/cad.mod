
NEURON {
    SUFFIX cad
    USEION ca READ cai, cao, ica WRITE cai, ica
    RANGE ica_pump, debugVal0, debugVal1, debugVal2, debugVal3
    GLOBAL vrat0, vrat1, vrat2, vrat3
    :read cai- to init buffer, reads ica- which is the Ca++ influx into outermost shell, reads cao- which is used in pump reaction scheme
    :and writes cai- because it computes Ca++ in the outermost shell, ica- writes pump current here
}

:DEFINE Nannuli 4

UNITS {
    (mol) =    (1)
    (molar) =  (1/liter)
    (mM) =     (millimolar)
    (um) =     (micron)
    (mA) =     (milliamp)
    FARADAY =  (faraday) (10000 coulomb)
    PI = (pi)  (1)   
}

PARAMETER {
    DCa = 0.6   (um2/ms)
    k1buf = 500 (/mM-ms)
    k2buf = 0.5 (/ms)
    k1=1.e10    (um3/s)
    k2=50.e7    (/s)	: k1*50.e-3
    k3=1.e10    (/s)	: k1
    k4=5.e6	    (um3/s)	: k1*5.e-4
	area		(um2)
}

ASSIGNED {
    diam     (um)
    ica      (mA/cm2)
    cai      (mM)
    cao      (mM)
    vrat0
    vrat1
    vrat2
    vrat3
    Kd       (/mM) :dissociation constant for the buffer
    B0       (mM) :initial value of free buffer
    ica_pump (mA/cm2)
    last_ipump (mA/cm2)
    debugVal0 (mM)
    debugVal1 (mM)
    debugVal2 (mM)
    debugVal3 (mM)
}

CONSTANT { volo=1  (liter)} : extracellular volume

STATE {
    ca0       (mM) <1e-10> :ca0 is equivalent to cai
    ca1       (mM) <1e-10>
    ca2       (mM) <1e-10>
    ca3       (mM) <1e-10>
    CaBuffer0 (mM)
    CaBuffer1 (mM)
    CaBuffer2 (mM)
    CaBuffer3 (mM)
    Buffer0   (mM)
    Buffer1   (mM)
    Buffer2   (mM)
    Buffer3   (mM)

    pump      (mol/cm2) <1.e-3>
    pumpca    (mol/cm2) <1.e-15>
}


BREAKPOINT {
    SOLVE state METHOD cnexp

    ica = ica_pump :ensure that the pump current is reckoned in NEURON's calculation of cai
    last_ipump = ica_pump
    
    debugVal0 = ca0
    debugVal1 = ca1 
    debugVal2 = ca2
    debugVal3 = ca3
}


LOCAL factors_done, TotalBuffer, TotalPump

INITIAL {
    : 1st step of initialization -> set up the geometry (by computing vrat[] and frat[])
    if (factors_done == 0) { :this will be executed only ones...
        factors_done = 1
        factors() 
        }
    
    : 2nd step of initialization -> initialize state variables
    Kd = k2buf/k1buf
    TotalBuffer = 1.2
    TotalPump = 0.2

    pump = TotalPump / (1 + 1.e-18*k4*cao/k3)
    pumpca = 2.e-22

    ca0 = cai
    CaBuffer0 = (TotalBuffer*ca0) / (Kd + ca0)
    Buffer0 = TotalBuffer - CaBuffer0
    ca1 = cai
    CaBuffer1 = (TotalBuffer*ca1) / (Kd + ca1)
    Buffer1 = TotalBuffer - CaBuffer1
    ca2 = cai
    CaBuffer2 = (TotalBuffer*ca2) / (Kd + ca2)
    Buffer2 = TotalBuffer - CaBuffer2
    ca3 = cai
    CaBuffer3 = (TotalBuffer*ca3) / (Kd + ca3)
    Buffer3 = TotalBuffer - CaBuffer3

    VERBATIM
    //printf("init cai:%g, 0:%g, 1:%g, 2:%g, 3:%g\n",(cai, ca0, ca1, ca2, ca3));
    ENDVERBATIM

    ica_pump = 0

}

LOCAL frat0, frat1, frat2, frat3: scales the rate constant for model geometry

PROCEDURE factors() {
    :LOCAL r, dr2
    :r = 1/2 :assumes that d = 1 ...
    :dr2 = r/(Nannuli-1)/2
    
    :vrat[0] = 0
    :frat[0] = 2*r
    :FROM i=0 TO Nannuli-2 {
    :    vrat[i] = vrat[i] + PI*(r-dr2/2)*2*dr2 :interior half
    :    r = r - dr2
    :    frat[i+1] = 2*PI*r/(2*dr2) :outer radius of annulus/distance between centers
    :    r = r - dr2
    :    vrat[i+1] = PI*(r+dr2/2)*2*dr2: outer half of annulus
    :}

: ===================================================================================================== 

    :Nannuli = 4, d=1, l=1
    :vrats- scale shell volume (originally: d = 1uM, l = 1mM), shell volume: diam^2 * vrat[i]
    :frats- scale rate constants (frat[i+1] = A[i+1]/dr : surface area between shells i, i+1 / distance between shell centers)
    vrat0 = 0.23998277
    frat0 = 1
    vrat1 = 0.34906585
    frat1 = 15.70796327
    vrat2 = 0.17453293
    frat2 = 9.42477796
    vrat3 = 0.02181662
    frat3 = 3.14159265
   
}

DERIVATIVE state {
    : The STATEs ca[] are intensive variables (concentration, or mass/volume), [...] flux is the time derivative of an extensive variable.
    : This disparity is corrected by specifying STATE volumes with the COMPARTMENT statement...
    : The volume merely multiplies the dSTATE/dt left hand side of the equivalent differential equations,
    : converting it to an extensive quantity and making it consistent with flux terms in units of absolute quantity per time.

    :COMPARTMENT i, diam*diam*vol[i] {ca CaBuffer Buffer}
    :COMPARTMENT (1e10)*area     {pump pumpca}
    :COMPARTMENT (1e15)*volo     {cao}

    : pump
    :~ca[0] + pump <-> pumpca     ((1.e-11)*k1*area, (1.e7)*k2*area)
    :~pumpca       <-> pump + cao ((1.e7)*k3*area, (1.e-11)*k4*area)
    :ica_pump = 2*FARADAY*(f_flux-b_flux) / area

    : all currents except pump
    :~ ca[0] << (-(ica-ica_pump)*PI*diam / (2*FARADAY)): ica is direct Ca++ efflux
    :FROM i=0 TO Nannuli-2 {
    :    ~ ca[i] <-> ca[i+1] (DCa*frat[i+1], DCa*frat[i+1])
    :}
    :dsq = diam*diam
    :FROM i=0 TO Nannuli-1 {
    :    dsqvol = dsq*vrat[i]
    :    ~ ca[i] + Buffer[i] <-> CaBuffer[i] (k1buf*dsqvol, k2buf*dsqvol)
    :}    
    :cai = ca[0]

: =====================================================================================================

    LOCAL dsq, dsqvol0, dsqvol1, dsqvol2, dsqvol3, f_flux, b_flux

    dsq = diam*diam
    dsqvol0 = dsq*vrat0
    dsqvol1 = dsq*vrat1
    dsqvol2 = dsq*vrat2
    dsqvol3 = dsq*vrat3

    ca0' = (- ((ica-last_ipump)*PI*diam) / (2*FARADAY)  - (DCa*frat1)*ca0 + (DCa*frat1)*ca1 - (k1buf*dsqvol0)*ca0*Buffer0 + (k2buf*dsqvol0)*CaBuffer0) - ((1.e-11)*k1*area)*ca0*pump + ((1.e7)*k2*area)*pumpca / dsqvol0
    ca1' = ((DCa*frat1)*ca0 - ((DCa*frat1)+(DCa*frat2))*ca1 + (DCa*frat2)*ca2 - (k1buf*dsqvol1)*ca1*Buffer1 + (k2buf*dsqvol1)*CaBuffer1) / dsqvol1
    ca2' = ((DCa*frat2)*ca1 - ((DCa*frat2)+(DCa*frat3))*ca2 + (DCa*frat3)*ca3 - (k1buf*dsqvol2)*ca2*Buffer2 + (k2buf*dsqvol2)*CaBuffer2) / dsqvol2
    ca3' = ((DCa*frat3)*ca2 - (DCa*frat3)*ca3 - (k1buf*dsqvol3)*ca3*Buffer3 + (k2buf*dsqvol3)*CaBuffer3) / dsqvol3
    
    Buffer0' =   (-(k1buf*dsqvol0)*ca0*Buffer0 + (k2buf*dsqvol0)*CaBuffer0) / dsqvol0
    CaBuffer0' = ((k1buf*dsqvol0)*ca0*Buffer0 - (k2buf*dsqvol0)*CaBuffer0) / dsqvol0
    Buffer1' =   (-(k1buf*dsqvol1)*ca1*Buffer1 + (k2buf*dsqvol1)*CaBuffer1) / dsqvol1
    CaBuffer1' = ((k1buf*dsqvol1)*ca1*Buffer1 - (k2buf*dsqvol1)*CaBuffer1) / dsqvol1
    Buffer2' =   (-(k1buf*dsqvol2)*ca2*Buffer2 + (k2buf*dsqvol2)*CaBuffer2) / dsqvol2
    CaBuffer2' = ((k1buf*dsqvol2)*ca2*Buffer2 - (k2buf*dsqvol2)*CaBuffer2) / dsqvol2
    Buffer3' =   (-(k1buf*dsqvol3)*ca3*Buffer3 + (k2buf*dsqvol3)*CaBuffer3) / dsqvol3
    CaBuffer3' = ((k1buf*dsqvol3)*ca3*Buffer3 - (k2buf*dsqvol3)*CaBuffer3) / dsqvol3
    
    pump' =    (-((1.e-11)*k1*area)*ca0*pump + (((1.e7)*k2*area)+((1.e7)*k3*area))*pumpca - ((1.e-11)*k4*area)*pump*cao/volo) / (1e10)*area
    pumpca' =  (((1.e-11)*k1*area)*ca0*pump - (((1.e7)*k2*area)+((1.e7)*k3*area))*pumpca + ((1.e-11)*k4*area)*pump*cao/volo) / (1e10)*area
    :cao' =     (((1.e7)*k3*area)*pumpca - ((1.e-11)*k4*area)*pump*cao) / (1e15)*volo
    f_flux =   ((1.e7)*k3*area)*pumpca
    b_flux =   ((1.e-11)*k4*area)*pump*cao
    ica_pump = 2*FARADAY*(f_flux-b_flux) / area

    cai = ca0
    
    


}


