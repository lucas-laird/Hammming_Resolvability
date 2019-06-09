from sys import *
import os
import time
import helper_funcs as hf
import numpy as np
from itertools import product
import pickle

def generate_hypercube_RSet(k):
    nodes = [''.join(x) for x in product('01',repeat = k)]
    np.random.shuffle(nodes)
    r = nodes.pop()
    R = [r]
    rank = 1
    is_resolving = False
    G,fs,z = hf.hypercube_polys(k)
    Gis = []
    for i,f in enumerate(fs):
        Gi = hf.groebnerbasis(True,G+[f],1,z,order = 'lex')
        Gis.append(Gi)
    A = hf.hypercube_matrix(R)
    lin_fcns = hf.make_linearEqns(A,z)
    for i,Gi in enumerate(Gis):
        Gi = hf.groebnerbasis(True,list(Gi)+lin_fcns,1,z,order = 'lex')
        Gis[i] = Gi
    while is_resolving:
        is_resolving = True
        r = node.pop()
        temp = hf.hypercube_matrix([r])
        A_temp = np.vstack([A,temp])
        new_rank = np.linalg.matrix_rank(A_temp)
        if new_rank > rank:
        	A = A_temp
        	rank = new_rank
	        R.append(node)
	        lin_fcns = hf.make_linearEqns(temp,z)
	        for i,Gi in enumerate(Gis):
	            Gi = hf.groebnerbasis(True,list(Gi)+lin_fcns,1,z,order = 'lex')
	            Gis[i] = Gi
	            if list(Gi) != [1]:
	                is_resolving = False
    return(R)


ks = np.arange(8,17,1)
num_sets = 1
hypercube_sets = {}
hypercube_times = {}
time_taken = 0
for k in ks:
    string = "Starting k = {}, time elapsed = {}\r".format(k,time_taken)
    stdout.write(string)
    stdout.flush()
    res_sets = []
    times = []
    for i in range(num_sets):
        start = time.time()
        R = generate_hypercube_RSet(k)
        end = time.time()
        time_taken += end-start
        times.append(end-start)
        res_sets.append(R)
    key = '{}'.format(k)
    hypercube_times[key] = times
    hypercube_sets[key] = res_sets

with open('Hcube_Rsets.pkl', 'wb') as f:
	pickle.dump(hypercube_sets,f)

with open('Hcube_times.pkl', 'wb') as f:
	pickle.dump(hypercube_times,f)	
print("All Done")

