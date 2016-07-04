# -*- coding: utf-8 -*-
""" 
    Helper method to handle Ca++ dynamics in the pyramidal cell model of Hemond 2008 (implemented by Migliore)
    1) the code deletes numberInternalDivisions from annotation (makes it a valid nml file)
    2) iterates over dendrites and apical_dendrites and calculate the average diam of segments (the depth of Ca++ shell is changed (from hoc level) according to the diameter in the original code)
    Authors: András Ecker, Padraig Gleeson
"""

import neuroml

import neuroml.loaders as loaders
import neuroml.writers as writers

import numpy as np


def helper_morphology(morph_file):   

    doc = loaders.NeuroMLLoader.load(morph_file)
    print("Loaded morphology file from: %s"%morph_file)
    
    dSegs = []
    # iterates over segment groups
    for segGroup in doc.cells[0].morphology.segment_groups:
        # deletes numberInternalDivisions from annotation (will be valid -> it will be possible to visualize it on OSB)
        if segGroup.annotation:
            segGroup.annotation = None
       
        # gets the seg.ids of segmentGroups (dendrite_group)
        if segGroup.neuro_lex_id:
            if segGroup.id[:segGroup.id.rfind('_')] == "dendrite":
                for member in segGroup.members:
                    dSegs.append(member.segments)
    
    
    dSegDiams = []
    # iterates over segments and gets the diameter of the segment (distal.diam)
    for seg in doc.cells[0].morphology.segments:
        d = seg.distal.diameter
        if seg.id in dSegs:
            dSegDiams.append(d)
            

    capool = ''
    capool += '\t<LawrenceCaConcentrationModel id="capool_diam_%s" ion="ca" diam_nml="%gum" DCa="0.6per_ms" TotalBuffer="1.2mM" k1buf="500per_ms" k2buf="0.5per_ms"/>\n\n'%(str(np.mean(dSegDiams)).replace('.', '_'), np.mean(dSegDiams))

    # write out to file... afterwards copy them (manually) to Capool.channel.nml TODO: automate this... 
    capools = open("capools.txt", 'w')
    capools.write(capool)
    capools.close()

    writers.NeuroMLWriter.write(doc, morph_file)
    print("Saved modified file to: " + morph_file)


if __name__ == "__main__":

    helper_morphology("LawrenceOLM.cell.nml")


