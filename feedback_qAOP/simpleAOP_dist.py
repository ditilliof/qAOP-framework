import os
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import random
import timeit
from multiprocessing import Pool
from datetime import datetime

# Model name
MODEL = "simpleAOP"
DATETIME = datetime.today().strftime('%Y%m%d_%H%M%S')

# Create directories
pathToOutput = f"./{MODEL}/"

# Create directories if they do not exist
if os.path.isdir(pathToOutput):
    print(f"{pathToOutput} does already exist. Overwriting previous output...")
else:
    print(f"{pathToOutput} does not exist. Creating a new directory and writing the output to {pathToOutput}")
    os.mkdir(pathToOutput)

# Define the model
def AOPP(ic, t, args):
    k, d1, s2, d2, eps = args
    KE1, KE2 = ic  # initial conditions, ic is the vector [KE1(0), KE2(0)]
    
    dKE1 = -d1 * KE1 + k * KE2 + eps
    dKE2 = s2 * (KE1**n) / (1 + (KE1**n)) - d2 * KE2 + eps

    return [dKE1, dKE2]

print("The model is defined")

# Fix values for the parameters
n = 3
t = np.linspace(0, 1000, 100)  # Properly define the time array

print("Parameters are defined")

# Define argument function of the starmap (multiprocessing)
def solveAOP(k, d1, s2, d2, eps):
    result = odeint(AOPP, [1, 1], t, args=([0, d1, s2, d2, eps],))
    ic = [result[-1, 0], result[-1, 1]]
    result = odeint(AOPP, ic, t, args=([k, d1, s2, d2, eps],))
    return result[1, 1]

print("Solve function is defined")



# Means of the distributions
mean_d1 = 1
mean_s2 = 4
mean_d2 = 1
mean_eps = 0.2

# Standard deviations of the distributions (example values, you can adjust)
sd_d1 = 0.1
sd_s2 = 0.3
sd_d2 = 0.1
sd_eps = 0.01  # Assuming a small noise level

k = 10
print("Means and standard deviations for distributions are defined")

start = timeit.default_timer()

# Number of iterations
num_iterations = 1000000

# Define items for multiprocessing by sampling from distributions
items = [(k,
          np.random.normal(mean_d1, sd_d1),
          np.random.normal(mean_s2, sd_s2),
          np.random.normal(mean_d2, sd_d2),
          np.random.normal(mean_eps, sd_eps)) for _ in range(num_iterations)]

print("Items are defined")
stop = timeit.default_timer()

NCPU = 56
pool = Pool(NCPU)
start = timeit.default_timer()

KE2_vec = pool.starmap(solveAOP, items)
KE2_max = max(KE2_vec)
AO_sc = KE2_vec / KE2_max
stop = timeit.default_timer()

print(f"The time of simulation for the parallel computation is {stop - start} seconds")

# Plotting the histogram
plt.hist(AO_sc, bins=500, range=(0, 1), density=True)
plt.gca().set(title='AO score distribution, k = 10', xlabel='AO score')
plt.savefig(f"{pathToOutput}{DATETIME}_FDT_AOdistribution.pdf")
plt.show()
