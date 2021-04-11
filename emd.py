#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 10:40:29 2021

@author: carlos
"""


from numpy import genfromtxt
#my_data1 = genfromtxt('eventos/event_1', delimiter=' ')
#my_data2 = genfromtxt('eventos/event_2', delimiter=' ')
import energyflow as ef

import glob
import os
import numpy as np

file_list = sorted(glob.glob(os.path.join(os.getcwd(), "eventos", "event_*")))
hist_list = sorted(glob.glob(os.path.join(os.getcwd(), "histogramas2d", "hist_*")))
event_list = sorted(glob.glob(os.path.join(os.getcwd(), "eventosfin", "ev_*")))


corpus = []
corpus_hist = []
corpus_event= []

for file_path in file_list:
    with open(file_path) as f_input:
        corpus.append(genfromtxt(f_input, delimiter=' '))
        
for file_path in hist_list:
    with open(file_path) as f_input:
        corpus_hist.append(genfromtxt(f_input, delimiter=' '))
        
for file_path in event_list:
    with open(file_path) as f_input:
        corpus_event.append(genfromtxt(f_input, delimiter=' '))

#print(corpus)

res=np.sqrt(ef.emd.emds_wasserstein(events0=corpus, events1=None, R=80, beta=2.0, norm=False, gdim=3, mask=False,
                                                       external_emd_handler=None,
                                                       n_jobs=-1, print_every=0, verbose=0,
                                                       throw_on_error=True, n_iter_max=100000,
                                                       epsilon_large_factor=10000.0,
                                                       epsilon_small_factor=1.0))


np.savetxt("array.txt", res, fmt="%s")