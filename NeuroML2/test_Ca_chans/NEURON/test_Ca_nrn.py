#!/usr/bin/ipython -i

from neuron import *
from nrn import *
#from neuron import gui


def create_comp(name = 'soma'):
    
    comp = h.Section('soma')

    comp.insert('cal')
    comp.insert('cat')
    comp.insert('kca')

    comp.insert('cabuff')
    comp.insert('pas')

    # comp.insert('CaClamp')

    comp.nseg = 1
    comp.L = 10
    comp.diam = 1
    
    comp(0.5).g_pas = 3e-4
    comp(0.5).e_pas = -65

    comp(0.5).cal.gcalbar = 2e-5
    comp(0.5).cat.gbar = 2e-5
    comp(0.5).kca.gkbar = 1e-4
    
    #phi will be multiplied by ica _density_
    area = h.area(0.5)
    phi = 3e-3
    print '0.1 * phi * area to be used in lems', 0.1 * phi * area 
    comp(0.5).cabuff.phi = phi

    h.cao0_ca_ion = 2
    h.cai0_ca_ion = 5e-6

    h.celsius = 24

    return comp

    
def plot_timeseries(vdict, varlist):
    from pylab import plot, show, figure, title
    t = vdict['t']
    for n in varlist:
        figure()
        plot(t, vdict[n], label=n)
        title(n)
    
    show()

def create_dumps(section, varlist):
    recordings = {n: h.Vector() for n in varlist}

    for (vn, v) in recordings.iteritems():
        v.record(section(0.5).__getattribute__('_ref_' + vn))
    
    recordings['t'] = h.Vector()
    recordings['t'].record(h._ref_t)
    return recordings 


def dump_to_file(vdict, varlist, fname='nrn_ca.dat'):
    from numpy import savetxt, array

    vnames = ['t'] + varlist
    X = array([vdict[x].to_python() for x in vnames]).T
    savetxt(fname, X)


def run(tstop=10, dt=0.001):
    h.dt = dt
    h.finitialize(-65)
    h.fcurrent()
    h.frecord_init()
    while h.t < tstop:
        h.fadvance()



comp = create_comp('soma')

stim = h.IClamp(0.5, sec=comp)
stim.delay = 4
stim.dur = 6.0
stim.amp = 0.005

varlist = ['v', 'ica_cal', 'ica_cat', 'cai']
ds = create_dumps(comp, varlist)

run(50, 0.001)

plot_timeseries(ds, varlist)
dump_to_file(ds, varlist)
