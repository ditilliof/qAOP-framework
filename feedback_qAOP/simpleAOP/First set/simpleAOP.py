import os
import matplotlib
import numpy as np
import pandas as pd
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.axes as axes
from matplotlib.backends.backend_pdf import PdfPages
#import ipywidgets as widgets
#from IPython.display import display
import seaborn as sns
import random
import timeit
from multiprocessing import Pool
from datetime import datetime


#model name
MODEL = "simpleAOP"
DATETIME = datetime.today().strftime('%Y%m%d_%H%M%S')

#Create directories

pathToOutput = "/home/ditilliof/Phd material/qAOP_framework/feedback qAOP/" + MODEL + "/"

# Create directories if they do not exist
if os.path.isdir(pathToOutput):
    print(pathToOutput + " does already exist. Overwriting previous output...")
else:
    print(pathToOutput + " does not exist. Creating a new directory and writing the output to " + pathToOutput)
    os.mkdir(pathToOutput)



##Define the model 
### NB: in the future you will need to create a separate file for this


def AOPP(ic, t, args):
    
    k, d1, s2, d2, eps= args
    
    KE1, KE2 = ic #initial conditions, ic is the vector [KE1(0), KE2(0)]
    
    dKE1 =  - d1  * KE1 + k*KE2 + eps
    dKE2 = s2*(KE1**n)/(1+(KE1**n)) - d2 * KE2 + eps

    return [dKE1, dKE2]

print("the model is defined")

#we fix values for the parameters
n=3
t = [0,10000]

print("parameters are defined")

#Define argument function of the starmap (multiprocessing)

def solveAOP(k,d1,s2,d2,eps):
    result= odeint(AOPP,[1,1], t, args = ([0, d1, s2, d2, eps],))
    ic = [result[len(t)-1, 0], result[len(t)-1, 1]]
    result = odeint(AOPP,ic, t, args = ([k, d1, s2, d2, eps],))
    return result[1,1]

print("solve function is defined")

ress = 50
t = [0,1000]

d1_set = np.linspace(0,1, ress)
s2_set = np.linspace(3,20, ress)
d2_set = np.linspace(0.7,1, ress)
eps_set = np.linspace(0, 1, ress)
k = 0.4
print("sets are defined")

start = timeit.default_timer()


#I need 87 seconds to define 25^6 items, 98 hours to define 100^6
items = []
for d1 in d1_set:
    for s2 in s2_set:
        for d2 in d2_set:
            for eps in eps_set:
                 items.append((k,d1,s2,d2,eps))
                        
print("items are defined")
#print("items: {}".format(items))
stop = timeit.default_timer()
  
NCPU = 56
pool = Pool(NCPU)
start = timeit.default_timer()

KE2_vec = pool.starmap(solveAOP, items)
KE2_max = max(KE2_vec)
AO_sc = KE2_vec/KE2_max
stop = timeit.default_timer()

#print(KE2_vec)
print("The time of simulation for the parallel computation is {} seconds".format(stop - start))

plt.hist(AO_sc, bins = 50, range = (0,1), density = True)
plt.gca().set(title='AO score distribution, k = 0.4', xlabel = "AO score")
plt.savefig(pathToOutput + DATETIME + "_FDT_" + "_AOdistribution.pdf")



