# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 12:27:43 2020

@author: Wladek
"""

from neuron import h
from hoc2swc import hoc2swc

hoc2swc("pir.hoc", "pir_model.swc", separate_process=False)