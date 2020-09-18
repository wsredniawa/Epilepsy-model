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
import pylab as py
import numpy as np 

h.load_file('stdrun.hoc')

def vis_3d_nueron(cell, swc, csd=np.array([0])):
    from mayavi import mlab
    mlab.figure(bgcolor=(0.2, 0.2, 0.2), size=(1000, 800))
    mlab.points3d(swc[2], swc[3], swc[4], scale_mode='none', scale_factor=3)
    mlab.plot3d(swc[2,1:], swc[3,1:], swc[4, 1:], color=(1,1,1),tube_radius=2)
    # morph_x = []
    # morph_y = []
    # morph_z = []
    # kos = []
    # for i,seg in enumerate(cell.secs): 
    #     # print(i, seg)
    #     ko = seg.hoc.psection()['ions']['k']['ko'][0]
    #     pts3d = seg.hoc.psection()['morphology']['pts3d']
    #     for point in pts3d:
    #         morph_x.append(point[0])
    #         morph_y.append(point[1]) 
    #         morph_z.append(point[2])
    #         kos.append(ko)
    # mlab.points3d(morph_x,morph_y,morph_z, kos,colormap = 'bwr', vmax =3.008, vmin=3.005, scale_mode='none', scale_factor=3)

def LOT(source, weights, pos=(0,0,0)):
    cell = Cell(name="lot", compile_paths="mods")
    cell.add_sec("axon", diam=1, l=1000, nseg=10)
    cell.insert("hh")
    cell.insert("nakpumpder")
    cell.insert("kdifl")
    axon = cell.filter_secs("axon")
    for i in axon.hoc: 
        i.nakpumpder.totalpump = 0.05
        i.kdifl.fhspace = 25
    syn1 = cell.add_synapse(source=source, seg=axon(0),
                            mod_name="ExpSyn", netcon_weight=weights[0])
    syn1.hoc.tau = 0.5
    return cell, axon, syn1

def pir_cell(source, weights, pos=(0,0,0)):
    # global src1, sec2
    cell = Cell(name="pyramidal")
    cell.load_morpho(filepath="pir_model.swc")
    posswc = np.loadtxt('pir_model.swc').T
    cell.insert("pas")
    cell.insert("hh")
    cell.insert("kdiff2", 'dend[0]')
    sec1 = cell.filter_secs("dend[0]")
    syn1 = cell.add_synapse(source=source, seg=sec1(0), 
                            mod_name="Exp2Syn", netcon_weight=weights[0])
    src1 = cell.filter_secs('axon')
    sec2 = cell.filter_secs("dend[1]")
    syn2 = cell.add_synapse(source=src1(1), seg=sec2(1), 
                            mod_name="Exp2Syn", netcon_weight=weights[1])
    syns = [syn1, syn2]
    return cell, syns, posswc

def inh_cell(src, seg_names, weights, pos=(0,0,0)):
    cell = Cell(name="Inhibitory")
    cell.load_morpho(filepath="inh_model.swc")
    posswc = np.loadtxt('inh_model.swc').T
    cell.insert("pas")
    cell.insert("hh")
    sec1 = cell.filter_secs(seg_names[0])
    syn1 = cell.add_synapse(source=src, seg=sec1(0), 
                            mod_name="Exp2Syn", netcon_weight=weights[0])
    src1 = cell.filter_secs('axon[0]')
    syn2 = cell.add_synapse(source=src1(1), seg=seg_names[1], 
                            mod_name="Exp2Syn", netcon_weight=weights[1])
    syn2.hoc.e = -90
    return cell, [syn1], posswc

if __name__ == '__main__':
# LOT stimualtion
    stim = NetStimCell("stim").make_netstim(start=10, number=10, interval=200)
    lot_cell, axon, syn1 = LOT(stim, [1])
# Network
    pir1, syns_pir, pos_pir1  = pir_cell(axon(1), [1,1])
    inh1, syns1, pos_inh1  = inh_cell(axon(1), 
                                      ['dend[0]', pir1.filter_secs('dend[1]')(1)], [1,1])
    inh2, syns2, pos_inh2  = inh_cell(pir1.filter_secs('axon')(1), 
                                      ['dend[0]', pir1.filter_secs('soma')(1)], [1,1])
    inh3, syns3, pos_inh3  = inh_cell(pir1.filter_secs('axon')(1), 
                                        ['dend[0]', lot_cell.filter_secs('axon')(0)], [1,1])
# Potassium leak
    dend0 = pir1.filter_secs('dend[0]')
    soma = pir1.filter_secs('soma')
    h.setpointer(axon(.5).hoc._ref_ko, 'kog', dend0(0.5).hoc.kdiff2)
# Visualization
    vis_3d_nueron(pir1, pos_inh1)
    rec_list = [axon(0.5), dend0(0.5)]#, inh3.filter_secs('soma')(0.5)]
    recs_v, recs_ko = Record(rec_list, variables='v'), Record(rec_list, variables='ko')
    # sim = Simulation(init_v=-70, warmup=2, with_neuron_gui=False)
    # for i in range(100):
    #     sim.run(runtime=25, stepsize=25)
    #     recs_v.plot(animate=True, y_lim=[-80,60], position='merge')
    #     recs_ko.plot(animate=True, y_lim=[0,15], position='merge')