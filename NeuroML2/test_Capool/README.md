### Ca++ dynamics

The model uses a very complex internal Ca++ dynamics, which includes Ca++ diffusion (between concentring rings) buffering and pumping... see more in the original code [ModelDB_NEURON/cad_orig.mod](https://github.com/andrisecker/Lawrence2006-CA1-OLM/blob/master/NeuroML2/test_Capool/ModelDB_NEURON/cad_orig.mod).

    cd ModelDB_NEURON
    python test_Ca_nrn.py

In order to implement it in NeuroML2 we created a simplified(/modified) NEURON version [NEURON/cad.mod](https://github.com/andrisecker/Lawrence2006-CA1-OLM/blob/master/NeuroML2/test_Capool/NEURON/cad.mod) with ODE's (and hard coded variables) instead of the original kinetic schemes (which is producing the same behaviour as the original one).

    cd NEURON
    python test_Ca_nrn.py

Based on the simplified NEURON code we created a NeuroML2 file [../../Capool.channel.nml](https://github.com/andrisecker/Lawrence2006-CA1-OLM/blob/master/NeuroML2/Capool.channel.nml) reprodes the same behaviour without the pump, but has some limitations when used with the Ca++ pump on the membrane. Still in development ...

    cd NeuroML2
    python test_Ca_jnml.py
