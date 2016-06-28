### test Ca++ channels

The driving force of the Ca++ channels ([NEURON/ICaL.mod](https://github.com/andrisecker/Lawrence2006-CA1-OLM/blob/master/NeuroML2/test_Ca_chans/NEURON/ICaL.mod), [NEURON/ICaT.mod](https://github.com/andrisecker/Lawrence2006-CA1-OLM/blob/master/NeuroML2/test_Ca_chans/NEURON/ICaT.mod)) in the model neruon is based on a modified (Jaffe 1994) Goldman–Hodgkin–Katz flux equation [wikipedia](https://en.wikipedia.org/wiki/GHK_flux_equation).

We tested the NeuroML2 implementation of this channels with a simple Ca++ buffer. See more about the Ca++ dynamics of the model [here](https://github.com/andrisecker/Lawrence2006-CA1-OLM/tree/master/NeuroML2/test_Capool).

    cd NEURON
    python test_Ca_nrn.py
    cd ../NeuroML2
    python test_Ca_jnml.py

> Note: One can use the NeuroML2 version by including [NeuroML2/channelDensityGHK2.xml](https://github.com/andrisecker/Lawrence2006-CA1-OLM/blob/master/NeuroML2/test_Ca_chans/NeuroML2/channelDensityGHK2.xml), but in order to convert it to NEURON the development version of jnml is needed [see descp. here](https://github.com/NeuroML/jNeuroML).
