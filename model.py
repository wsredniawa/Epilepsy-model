# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 12:14:09 2020

@author: Wladek
"""


import os
from neuronpp.cells.cell import Cell
from neuronpp.core.cells.netstim_cell import NetStimCell
from neuronpp.utils.record import Record
from neuronpp.utils.graphs.network_status_graph import NetworkStatusGraph
from neuron import h
from neuronpp.utils.simulation import Simulation
from mayavi import mlab
import pylab as py
import numpy as np 

h.load_file('stdrun.hoc')

def LOT(source, weights):
    cell = Cell(name="lot", compile_paths="mods")
    cell.add_sec("axon", diam=1, l=1, nseg=1)
    cell.insert("pas")
    cell.insert("hh")
    secs = cell.filter_secs("axon")
    cell.insert("kdifl")
    syn1 = cell.add_synapse(source=source, seg=secs(0),
                            mod_name="ExpSyn", netcon_weight=weights[0])
    axon = cell.filter_secs("axon")
    recs = Record(axon(0.5), variables='v')
    return cell, recs, axon

def piriform_cell(source, weights):
    cell = Cell(name="pyramidal")
    cell.load_morpho(filepath="pyr.swc")
    posswc = np.loadtxt('pyr.swc').T
    # cell.add_sec("soma", diam=1, l=1000, nseg=1)
    cell.insert("pas")
    cell.insert("hh")
    cell.insert("kdiff2")
    secs = cell.filter_secs("soma")
    syns = cell.add_synapse(source=source, seg=secs(0), 
                            mod_name="Exp2Syn", netcon_weight=weights)
    soma = cell.filter_secs("soma")
    recs = Record(soma(0.5), variables='v')
    return cell, soma, posswc, recs

def vis_3d_nueron(pir_cell, csd=np.array([0])):
    mlab.figure(bgcolor=(0.2, 0.2, 0.2), size=(1000, 800))
    # mlab.points3d(pos[2], pos[3], pos[4], scale_mode='none', scale_factor=3)
    morph_x = []
    morph_y = []
    morph_z = []
    kos = []
    for i,seg in enumerate(pir_cell.secs): 
        # print(i, seg)
        ko = seg.hoc.psection()['ions']['k']['ko'][0]
        pts3d = seg.hoc.psection()['morphology']['pts3d']
        for point in pts3d:
            morph_x.append(point[0])
            morph_y.append(point[1]) 
            morph_z.append(point[2])
            kos.append(ko)
    fig = mlab.points3d(morph_x,morph_y,morph_z, kos,colormap = 'bwr', 
                        vmax =3.008, vmin=3.005,
                        scale_mode='none', scale_factor=3)
    return fig
    
def virtual_points3d(coords, figure=None, scale_factor=None, color=None, 
    name=None):
    c = np.array(coords)
    source = mlab.pipeline.scalar_scatter(c[:,0], c[:,1], c[:,2],
                                          figure=figure)
    return mlab.pipeline.glyph(source, scale_mode='none', 
                               scale_factor=scale_factor,
                               mode='sphere', figure=figure, color=color, name=name)

if __name__ == '__main__':
    # Create NetStim
    stim = NetStimCell("stim").make_netstim(start=21, number=1000, interval=50)
    # Create population 1
    lot_cell, axon_recs, axon = LOT(stim, [0.1, 0.1])
    # pop1.record()
    pir_cell, soma, swc, recs  = piriform_cell(stim, .1)
    fig = vis_3d_nueron(pir_cell)
    # h.setpointer(axon.hoc._ref_ko, 'kog', soma.hoc.kdiff2)
    h.setpointer(axon(.5).hoc._ref_ko, 'kog', soma(.5).hoc.kdiff2)
    # Create connectivity graph grouped by populations, with weighs and spike rates updated
    # h.tstop = 80
    # h.run()
    # axon_recs.plot()
    # Run
    sim = Simulation(init_v=-70, warmup=2)
    for i in range(1000):
        sim.run(runtime=2)
        kos = []
        for i,seg in enumerate(pir_cell.secs): 
            ko = seg.hoc.psection()['ions']['k']['ko'][0]
            pts3d = seg.hoc.psection()['morphology']['pts3d']
            for point in pts3d: kos.append(ko)
        ms = fig.mlab_source    
        ms.reset(s=np.array(kos))
        print(np.array(kos)[0])
        recs.plot(animate=True, )
        # pir_rec_k.plot(animate=True)
        # pop1.plot(animate=True)
        # pop2.plot(animate=True)
        # pop3.plot(animate=True)