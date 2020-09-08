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
    cell.add_sec("axon", diam=1, l=1000, nseg=10)
    # cell.insert("pas")
    cell.insert("hh")
    cell.insert("nakpumpder")
    cell.insert("kdifl")
    axon = cell.filter_secs("axon")
    for i in axon.hoc: 
        i.nakpumpder.totalpump = 0.04 
        i.kdifl.fhspace = 50
        # i.kdifl.
    syn1 = cell.add_synapse(source=source, seg=axon(0),
                            mod_name="ExpSyn", netcon_weight=weights[0])
    syn1.hoc.tau = 0.5
    return cell, axon, syn1

def piriform_cell(source, weights):
    cell = Cell(name="pyramidal")
    cell.load_morpho(filepath="pyr.swc")
    posswc = np.loadtxt('pyr.swc').T
    # cell.add_sec("soma", diam=1, l=1000, nseg=1)
    cell.insert("pas")
    cell.insert("hh")
    #cell.insert("nakpumpder")
    cell.insert("kdiff2", 'dend[0]')
    secs = cell.filter_secs("soma")
    syns = cell.add_synapse(source=source, seg=secs(0), 
                            mod_name="Exp2Syn", netcon_weight=weights)
    soma = cell.filter_secs("soma")
    return cell, soma, posswc

def vis_3d_nueron(pir_cell, csd=np.array([0])):
    mlab.figure(bgcolor=(0.2, 0.2, 0.2), size=(1000, 800))
    mlab.points3d(swc[2], swc[3], swc[4], scale_mode='none', scale_factor=3)
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
    
def virtual_points3d(coords, figure=None, scale_factor=None, color=None, 
    name=None):
    c = np.array(coords)
    source = mlab.pipeline.scalar_scatter(c[:,0], c[:,1], c[:,2],
                                          figure=figure)
    return mlab.pipeline.glyph(source, scale_mode='none', 
                               scale_factor=scale_factor,
                               mode='sphere', figure=figure, color=color, name=name)

if __name__ == '__main__':
    stim = NetStimCell("stim").make_netstim(start=10, number=1000, interval=200)
    lot_cell, axon, syn1 = LOT(stim, [1])
    pir_cell, soma, swc  = piriform_cell(axon(1), 1)
    dend0 = pir_cell.filter_secs('dend[0]')
    soma = pir_cell.filter_secs('soma')
    h.setpointer(axon(.5).hoc._ref_ko, 'kog', dend0(0.5).hoc.kdiff2)
    # vis_3d_nueron()
    rec_list = [axon(0.5), dend0(0.5), soma(0.5)]
    recs_v, recs_ko = Record(rec_list, variables='v'), Record(rec_list, variables='ko')
    sim = Simulation(init_v=-70, warmup=2, with_neuron_gui=False)
    for i in range(9000):
        sim.run(runtime=25, stepsize=25)
        # for i,seg in enumerate(pir_cell.secs): 
            # ko = seg.hoc.psection()['ions']['k']['ko'][0]
            # pts3d = seg.hoc.psection()['morphology']['pts3d']
            # for point in pts3d: kos.append(ko)
        # ms = fig.mlab_source    
        # ms.reset(s=np.array(kos))
        # print(np.array(kos)[0])
        recs_v.plot(animate=True, y_lim=[-80,60], position='merge')
        # recs_ko.plot(animate=True, y_lim=[0,15], position='merge')
        # pir_rec_k.plot(animate=True)
        # pop1.plot(animate=True)
        # pop2.plot(animate=True)
        # pop3.plot(animate=True)