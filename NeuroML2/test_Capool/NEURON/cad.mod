
NEURON {
    SUFFIX cad
    USEION ca READ cai, cao, ica WRITE cai, ica
    RANGE ica_pump
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
    FARADAY =  (faraday) (10000 coulomb/mol)
    PI = (pi)  (1)   
}

PARAMETER {
    DCa = 0.6         (um2/ms)
    k1buf = 500       (/mM-ms)
    k2buf = 0.5       (/ms)
    TotalBuffer = 1.2 (mM)
    k1=1e-11    (/mM-ms)
    k2=50e4     (/ms)	
    k3=1e7      (/ms)	
    k4=5e-18	(/mM-ms)
    TotalPump = 0.2   (mol/cm2)
    PCa0 = 2.e-22     (mol/cm2)
	area		      (um2)
}

ASSIGNED {
    diam       (um)
    ica        (mA/cm2)
    cai        (mM)
    cao        (mM)
    vrat0
    vrat1
    vrat2
    vrat3
    Kd         (mM) :dissociation constant for the buffer
    ica_pump   (mA/cm2)
    last_ipump (mA/cm2)
}

:CONSTANT { volo=1  (liter)} : extracellular volume

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

    :ica = ica_pump :ensure that the pump current is reckoned in NEURON's calculation of cai
}


LOCAL factors_done

INITIAL {
    : 1st step of initialization -> set up the geometry (by computing vrat[] and frat[])
    if (factors_done == 0) { :this will be executed only ones...
        factors_done = 1
        factors() 
        }
    
    : 2nd step of initialization -> initialize state variables
    Kd = k2buf/k1buf

    pump = TotalPump / (1 + 1.e3 * k4/k3 * cao)
    pumpca = PCa0

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
    :COMPARTMENT i, diam*diam*vol[i] {ca CaBuffer Buffer}
    :COMPARTMENT (1e10)*area     {pump pumpca}
    :COMPARTMENT (1e15)*volo     {cao}

    : pump
    :~ca[0] + pump <-> pumpca     ((1.e-8)*k1*area, (1.e10)*k2*area)
    :~pumpca       <-> pump + cao ((1.e10)*k3*area, (1.e-8)*k4*area)
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

    LOCAL dsq, dsqvol0, dsqvol1, dsqvol2, dsqvol3, ca0_efl, ca0_dif, ca0_buf, ca0_pump, f_flux, b_flux

    dsq     = diam*diam
    dsqvol0 = dsq*vrat0
    dsqvol1 = dsq*vrat1
    dsqvol2 = dsq*vrat2
    dsqvol3 = dsq*vrat3

    ca0_efl  = - (ica)*PI*diam/dsqvol0 / (2*FARADAY)
    printf("ica: %g, ca_efflux: %g \n",ica, ca0_efl)
    ca0_dif  = - (DCa*frat1/dsqvol0)*ca0 + (DCa*frat1/dsqvol0)*ca1
    ca0_buf  = - k1buf*ca0*Buffer0 + k2buf*CaBuffer0
    ca0_pump = - ((1.e10)*k1*area)*ca0*pump/dsqvol0 + ((1.e10)*k2*area)*pumpca/dsqvol0
    :printf("ica: %g, ica_pump: %g, efflux: %g \n",ica, ica_pump, ca0_efl)

    ca0' = ca0_dif + ca0_buf + ca0_efl: + ca0_pump
    ca1' = (DCa*frat1/dsqvol1)*ca0 - ((DCa*frat1/dsqvol1)+(DCa*frat2/dsqvol1))*ca1 + (DCa*frat2/dsqvol1)*ca2 - k1buf*ca1*Buffer1 + k2buf*CaBuffer1
    ca2' = (DCa*frat2/dsqvol2)*ca1 - ((DCa*frat2/dsqvol2)+(DCa*frat3/dsqvol2))*ca2 + (DCa*frat3/dsqvol2)*ca3 - k1buf*ca2*Buffer2 + k2buf*CaBuffer2
    ca3' = (DCa*frat3/dsqvol3)*ca2 - (DCa*frat3/dsqvol3)*ca3 - k1buf*ca3*Buffer3 + k2buf*CaBuffer3

    Buffer0'   = -k1buf*ca0*Buffer0 + k2buf*CaBuffer0
    CaBuffer0' =  k1buf*ca0*Buffer0 - k2buf*CaBuffer0
    Buffer1'   = -k1buf*ca1*Buffer1 + k2buf*CaBuffer1
    CaBuffer1' =  k1buf*ca1*Buffer1 - k2buf*CaBuffer1
    Buffer2'   = -k1buf*ca2*Buffer2 + k2buf*CaBuffer2
    CaBuffer2' =  k1buf*ca2*Buffer2 - k2buf*CaBuffer2
    Buffer3'   = -k1buf*ca3*Buffer3 + k2buf*CaBuffer3
    CaBuffer3' =  k1buf*ca3*Buffer3 - k2buf*CaBuffer3
    
    :pump'   = (-((1.e10)*k1*area)*ca0*pump + (((1.e10)*k2*area)+((1.e10)*k3*area))*pumpca - ((1.e10)*k4*area)*pump*cao/volo) / (1e10)*area
    :pumpca' =  (((1.e10)*k1*area)*ca0*pump - (((1.e10)*k2*area)+((1.e10)*k3*area))*pumpca + ((1.e10)*k4*area)*pump*cao/volo) / (1e10)*area
    :pump'   = -k1*ca0*pump + (k2+k3)*pumpca - k4*pump*cao
    :pumpca' =  k1*ca0*pump - (k2+k3)*pumpca + k4*pump*cao
    :f_flux =   ((1.e10)*k3*area)*pumpca
    :b_flux =   ((1.e10)*k4*area)*pump*cao
    :ica_pump = 2*FARADAY*(f_flux-b_flux) / area
    :printf("f_flux: %g, b_flux: %g \n", f_flux, b_flux)

    cai = ca0

}


