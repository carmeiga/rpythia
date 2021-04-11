#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 10:09:53 2021

@author: carlos
"""

import matplotlib
import matplotlib.pyplot as plt
import ot
import pandas as pd
from sklearn.manifold import MDS
from math import floor
from math import sqrt

# Compute OT distance matrix, with TV regularization
def distance_matrix(A, C, tau=0):
    D = np.zeros((A.shape[1], A.shape[1]))
    for i in range(A.shape[1]):
        D[i, :i] = D[:i, i] = ot.emd2(np.ascontiguousarray(A[:,i]), A[:,:i], C) + tau*C.max()*abs(A[:,[i]] - A[:,:i]).sum(0)
    return D


# Distance of two bins on the 1-D torus
def torus_distance(i, j, n):
    d1 = abs(j - i)
    d2 = n - d1
    return min(d1, d2)



# Generate random distance-like matrix
def random_distance(n, seed=None):
    np.random.seed(seed)
    D = np.random.rand(n, n)
    D += D.T
    np.fill_diagonal(D, 0)
    return D/D.max()



# Hilbert metric between two matrices
def hilbert(D1, D2):
    idx = (np.eye(D1.shape[0]) != 1) # exclude diagonal
    div = np.log(D1[idx]/D2[idx])
    return div.max() - div.min()



# Plot template das miñas
plt.title('Exemplo rapi vs phi')
plt.imshow(corpus_hist[2], cmap='Blues')
plt.colorbar();



# Load dataset
mnist = pd.read_csv('mnist_train_small.csv')


# Subset dataset
n_images = mnist.shape[0]//70
A = np.array(mnist.T, dtype=float)[1:,:n_images]
labels = mnist['6'][:n_images]



# Keep only 0s and 1s
idx = np.array(labels <= 1)
labels = np.array(labels)[idx]
A = A[:,idx]



# Sort by label
idx = np.argsort(labels)
labels = labels[idx]
A = A[:,idx]


# Shape
m, n = A.shape[0], A.shape[1]
print(m, n)





# Display
plt.imshow(A[:,0].reshape((28, 28)), cmap='gray');

mm=len(corpus_hist)
n=15
A_gauss_2D = np.zeros((n*n, mm))
for k in range(mm):
    plane = np.zeros((n, n))
    plane=corpus_hist[k]
    A_gauss_2D[:,k] = plane.reshape(-1)

plt.imshow(A_gauss_2D[:,1].reshape((15, 15)), cmap='gray');

plt.imshow(np.flipud(A_gauss_2D[:,1].reshape((15, 15)).T), cmap='gray');

plt.imshow(corpus_hist[0].reshape(-1).reshape(15,15),cmap='grey')
A=A_gauss_2D #Para probar separar soft de hard
m, n = A.shape[0], A.shape[1]

A += 1e-6
B = A.T
A = A/A.sum(0)
B = B/B.sum(0)





# Generate random ground cost and distance
C = random_distance(m, seed=42)
D = random_distance(n, seed=24)

n_iter =  10
for k in range(n_iter):
  print('Computing D, iteration', k + 1)
  D = distance_matrix(A, C)
  D /= D.max()
  
  print('Computing C, iteration', k + 1)
  C = distance_matrix(B, D)
  C /= C.max()

plt.title('Ground cost')
plt.imshow(C);

np.savetxt("distD.txt", D, fmt="%s")
np.savetxt("distC.txt", C, fmt="%s")



plt.title('Distance matrix')
plt.imshow(D);

embed = MDS(n_components=2, dissimilarity='precomputed').fit_transform(dist)
plt.figure(figsize=(30, 30))
plt.scatter(embed[1:200, 0], embed[1:200, 1], s=60, marker=r'$s$');
plt.scatter(embed[200:400, 0], embed[200:400, 1], s=60, marker=r'$h$');
plt.scatter(embed[400:600, 0], embed[400:600, 1], s=60, marker=r'$w$');
plt.scatter(embed[600:800, 0], embed[600:800, 1], s=60, marker=r'$t$');

C=np.loadtxt('distC.txt'); 

# 0-specific pixels

resfinal = np.zeros((len(event_list), len(event_list)))

for k in range(len(event_list)):
 
    for p in range(len(event_list)):
      
        if k<=p:
            break
        
  
    
        eva=corpus_event[k]
        evb=corpus_event[p]
        
        custos=np.zeros((len(eva),len(evb)))
        
        ja=eva[:,2].astype(int)
        ia=eva[:,1].astype(int)
        
        jb=evb[:,2].astype(int)
        ib=evb[:,1].astype(int)
        
        
        
        for k1 in range(len(eva)):
     
            for p1 in range(len(evb)):
     
                custos[k1,p1]=C[15*ia[k1] + ja[k1]].reshape(15, 15)[ib[p1],jb[p1]]
      
        resfinal[k,p]=sqrt(ef.emd.emd_wasserstein(eva, evb, dists=custos, R=np.amax(custos), beta=2.0, norm=False, gdim=2, mask=False,
                                             return_flow=False, do_timing=False,
                                             n_iter_max=100000000,
                                             epsilon_large_factor=1000000000000000000.0, epsilon_small_factor=1.0))
           
    



dist=(resfinal+resfinal.T)/2

plt.imshow(dist);

np.savetxt("disthaler.txt", dist, fmt="%s")
np.savetxt("distpeyre.txt", D, fmt="%s")



i, j = 12, 1
CC = C[15*i + j].reshape(15, 15)
plt.imshow(CC);

