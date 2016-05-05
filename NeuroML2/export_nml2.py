from pyneuroml.neuron import export_to_neuroml2
# from pyneuroml.neuron.nrn_export_utils import clear_neuron

export_to_neuroml2("../NEURON/tester.hoc", 
                   "test.cell.nml", 
                   includeBiophysicalProperties=True,
                   known_rev_potentials={'Ih':-32.9, 'cal':0, 'cat':0,
                                         'kca':-95, 'Ika':-95, 'IM':-95, 'Ikdrf':-95, 'Ikdrs':-95, 'passsd':-70,
                                         'Ikdrfaxon':-95, 'Ikdrsaxon':-95, 'passaxon':-70})

