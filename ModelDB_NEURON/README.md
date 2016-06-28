### Original NEURON scripts for model

To run the scripts, [install NEURON](https://www.neuron.yale.edu/neuron/download) and run:

    git clone https://github.com/andrisecker/Lawrence2006-CA1-OLM.git # clone git repository
    cd Lawrence2006-CA1-OLM/ModelDB_NEURON
    nrnivmodl  # compile .mod files
    nrngui mosinit.hoc  # runs the simulation

Once the simulation starts you can create subplots of Fig 9 by clicking the appropriately labeled buttons or you can create the whole figure by clicking the create all of fig 9 button (lower right in control box).

> Note: The variable step method was used in this demo, if you would like to run the simulation as in the publication you can run cvode_active(0)!
