from pyneuroml.neuron import export_to_neuroml2
# from pyneuroml.neuron.nrn_export_utils import clear_neuron

export_to_neuroml2("../NEURON/tester.hoc", 
                   "test.cell.nml", 
                   includeBiophysicalProperties=True,
                   known_rev_potentials={'Ih':-32})

