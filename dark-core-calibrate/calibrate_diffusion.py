# coding: utf-8
import numpy as np
import pandas as pd
import scipy as sp
import subprocess
import os
from scipy import optimize
import pickle
import sys
from functools import lru_cache, wraps
import json

def np_cache(function):
    @lru_cache()
    def cached_wrapper(hashable_array):
        array = np.array(hashable_array)
        return function(array)

    @wraps(function)
    def wrapper(array):
        return cached_wrapper(tuple(array))

    # copy lru_cache attributes over too
    wrapper.cache_info = cached_wrapper.cache_info
    wrapper.cache_clear = cached_wrapper.cache_clear

    return wrapper

np.set_printoptions(precision=10)
np.random.seed(42) # for reproducibility 

Z_X_solar  =  0.02293 #0.02307 

# Directories and save files 
save_dir = "LOGS"
if not os.path.exists(save_dir):
    os.mkdir(save_dir)

age = 4.572
if len(sys.argv) > 2:
    age = sys.argv[2]

X_names = ['Y', 'Z', 'a']
#X = [0.27364444564903373, 0.01871176763464224, 1.7974118486632082]
#X = [0.2736202528288705, 0.01870713324531878, 1.797415951639218]
X = [0.272832246246596, 0.01853349976050075, 1.79923729991001] ## solar 
X = [0.27203765032798743, 0.018539064435792205, 1.8054920767517246] # solar MRG
#X = [0.27114379039643144, 0.018502541737273427, 1.8276307544675756] ## 1e-3
#X = [0.2726686381139602, 0.018531551318191467, 1.8020648672386708] ## 1e-4
#X = [0.2728348527, 0.0185336848, 1.7988029215]
#X = [0.2725708417, 0.018683564,  1.8081185829]
#X = [0.25710512413555603, 0.018345219857917257, 2.1257620067598566] ## 10^-2

if os.path.exists('optimized_params.json'):
    with open('optimized_params.json', 'r') as f:
        X = json.load(f)

X_var   = [0.005,        0.005,           0.01]
bounds  = [(0.25, 0.29), (0.012, 0.022), (1.5, 2.5)]
"""
X_names = [ 'Y',   'a']
X       = [0.26, 1.665]
X_var   = [0.005, 0.01]
bounds  = [(0.25, 0.29), (1.5, 2.5)]
print(X_names)
"""

P = lambda x: x#(x-X)/X_var
R = lambda x: x#x*X_var+X
lower = P(np.array([bound[0] for bound in bounds]))
upper = P(np.array([bound[1] for bound in bounds]))

def get_flags(names, args):
    return ' -'.join([''] + list(map(lambda x,y: x+' '+str(y), names, args)))

@np_cache
def calibrate(theta, fast=True, single=False):
    _theta = R(np.copy(theta))
    print("parameters:", _theta)
    
    Y, Z, alpha = _theta
    #Z = (1-Y) * Z_X_solar / (1+Z_X_solar) # [Fe/H] = 0 = log10(Z/X/0.02293); Z=1-Y-X
    
    if np.any(theta < lower) or np.any(theta > upper):
        print('out of bounds')
        if single: 
            return 1
        return np.ones(len(theta))
    
    bash_cmd = get_flags(X_names, _theta)
    print(bash_cmd)
    
    full = "./mesa.sh"+bash_cmd
    if fast:
        full = full + ' -f'
    print(full)
    with open('temp.txt', 'w') as output:
        process = subprocess.Popen(full.split(), 
            shell=False, stdout=output)
        process.wait(timeout=60000)
    
    # process the results 
    hist_file = os.path.join(save_dir, 'history.data')
    pro_file  = os.path.join(save_dir, 'solar.data')
    if not os.path.isfile(hist_file) or not os.path.isfile(pro_file): 
        print('No output')
        if single: 
            return 1
        return np.ones(len(theta))
    DF = pd.read_table(hist_file, sep='\s+', skiprows=5)
    prof = pd.read_table(pro_file, sep='\s+', skiprows=5)
    
    # correct MESA values 
    R_correct = 6.9566/6.957 #6.9599/6.957
    logR = DF['log_R'].values[-1]
    logR = np.log10((R_correct*10**logR))
    logL = DF['log_L'].values[-1]
    Fe_H = np.log10(prof.z.values[0] / prof.x.values[0] / Z_X_solar)
    print('log R:', logR)
    print('log L:', logL)
    print('[Fe/H]:', Fe_H)
    
    if np.abs(logR) < 1e-7 and np.abs(logL) < 1e-7 and np.abs(Fe_H) < 1e-5:
        print('good enough!')

        # Save optimized parameters
        with open('optimized_params.json', 'w') as f:
            json.dump(_theta.tolist(), f)

        sys.exit()
    
    if single:
        return np.log10(np.abs(logR)) + \
               np.log10(np.abs(logL)) + \
               np.log10(np.abs(Fe_H)) 
    return np.array([logR, logL, Fe_H])

eps = np.finfo(float).eps
result = optimize.root(calibrate, x0=P(np.array(X)))

print(result)

