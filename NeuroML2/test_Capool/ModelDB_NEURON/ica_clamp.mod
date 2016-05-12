
NEURON {
    SUFFIX ica_clamp
    USEION ca WRITE ica
    RANGE ica_clamp
    
}

PARAMETER {

    delay = 50 (ms)
    duration = 100 (ms)
    magnitude = -0.0005 (mA/cm2)
}

ASSIGNED {
    ica_clamp (mA/cm2)
    ica (mA/cm2)
    
}


BREAKPOINT {
    
    if (t <=  delay) {
        ica_clamp = 0
    }
    
    if (t >=  delay  && t <=  duration  +  delay) {
        ica_clamp = magnitude
    }
    
    if (t >=  duration  +  delay) {
        ica_clamp = 0
    }
       
    ica = ica_clamp    
}

