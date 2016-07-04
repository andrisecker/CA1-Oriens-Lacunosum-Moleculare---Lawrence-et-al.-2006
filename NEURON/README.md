### NEURON scripts for model

This folder contains files from the original model scripts [ModelDB:102288](https://senselab.med.yale.edu/modelDB/showModel.cshtml?model=102288) and additional tester scripts (in order to create a .mep file, which against LEMS/NeuroML2 implementation could be tested).

To run the scripts, [install NEURON](https://www.neuron.yale.edu/neuron/download) and run:

    git clone https://github.com/andrisecker/Lawrence2006-CA1-OLM.git  # clone git repository
    cd Lawrence2006-CA1-OLM/NEURON
    nrnivmodl  # compile .mod files
    nrngui test_lawrenceolm.hoc  # runs a simulation (single cell, current clamp) and saves data into lawrenceolm.dat
