import os
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import timeit
from multiprocessing import Pool
from datetime import datetime

# Model name
MODEL = "simpleAOP"
DATETIME = datetime.today().strftime('%Y%m%d_%H%M%S')

# Create directories
base_output_path = f"./{MODEL}/"
timestamped_output_path = os.path.join(base_output_path, DATETIME)
if not os.path.exists(timestamped_output_path):
    os.makedirs(timestamped_output_path)

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

# Means and standard deviations for the distributions
mean_d1 = 1
mean_s2 = 4
mean_d2 = 1
mean_eps = 0.2

# Standard deviations of the distributions (example values, you can adjust)
sd_d1 = 0.1
sd_s2 = 0.3
sd_d2 = 0.1
sd_eps = 0.01  # Assuming a small noise level

# k_values for the distributions and additional k_values for the last panel
distribution_k_values = [0.5, 0.625, 0.75, 0.875, 1]
additional_k_values = np.linspace(0.1, 1, 20).tolist()  # More k values for smoother line
all_k_values = sorted(distribution_k_values + additional_k_values)
num_iterations = 1000000

start = timeit.default_timer()

# Define items for multiprocessing by sampling from distributions
def generate_items(k):
    return [(k,
             np.random.normal(mean_d1, sd_d1),
             np.random.normal(mean_s2, sd_s2),
             np.random.normal(mean_d2, sd_d2),
             np.random.normal(mean_eps, sd_eps)) for _ in range(num_iterations)]

print("Items are defined")
stop = timeit.default_timer()

NCPU = 56
pool = Pool(NCPU)

# Function to run the simulation and normalize results
def run_simulation(k_values):
    all_results = {}
    max_ke2 = 0

    # First, run the simulation for k=1 to determine max KE2
    print("Running simulation for k = 1 to determine max KE2")
    items = generate_items(1)
    KE2_vec = pool.starmap(solveAOP, items)
    max_ke2 = max(KE2_vec)
    all_results[1] = [ke2 / max_ke2 for ke2 in KE2_vec]

    # Now run for all other k values and normalize with max_ke2 from k=1
    for k in k_values:
        if k == 1:
            continue  # Skip k=1 as it is already processed
        print(f"Running simulation for k = {k}")
        items = generate_items(k)
        KE2_vec = pool.starmap(solveAOP, items)
        AO_sc = [ke2 / max_ke2 for ke2 in KE2_vec]
        all_results[k] = AO_sc

    return all_results

start = timeit.default_timer()
all_results = run_simulation(all_k_values)
stop = timeit.default_timer()

print(f"The time of simulation for the parallel computation is {stop - start} seconds")

# Plotting the histogram for each k in distribution_k_values
for k in distribution_k_values:
    AO_sc = all_results[k]
    plt.figure(figsize=(6, 4))
    plt.hist(AO_sc, bins=500, range=(0, 1), density=True)
    plt.title(f'AO score distribution, k = {k}', fontsize=18)
    plt.xlabel('AO score', fontsize=20)
    plt.ylabel('Density', fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=18)
    plt.grid(False)
    plt.savefig(f"{timestamped_output_path}/FDT_AOdistribution_k{k}.pdf")
    plt.close()

print("Histograms saved for all k values")

# Calculate percentage of AO_sc values above the threshold for each k
threshold = 0.3
percentages = {}

for k, AO_sc in all_results.items():
    count_above_threshold = sum(1 for x in AO_sc if x > threshold)
    percentage_above_threshold = (count_above_threshold / len(AO_sc)) * 100
    percentages[k] = percentage_above_threshold

# Sort the percentages by k value
sorted_k_values = sorted(percentages.keys())
sorted_percentages = [percentages[k] for k in sorted_k_values]

# Create composite plot
fig, axs = plt.subplots(2, 3, figsize=(18, 10))

# Plot AO score distributions
for i, k in enumerate(distribution_k_values):
    AO_sc = all_results[k]
    ax = axs[i // 3, i % 3]
    ax.hist(AO_sc, bins=500, range=(0, 1), density=True)
    ax.set_title(f'AO score distribution, k = {k}', fontsize=18)
    ax.set_xlabel('AO score', fontsize=20)
    ax.set_ylabel('Density', fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=18)
    ax.grid(False)

# Plot the percentage of AO_sc values above the threshold vs k
ax = axs[1, 2]
ax.plot(sorted_k_values, sorted_percentages, linewidth=2)
ax.set_xlabel('k (a.u)', fontsize=20)
ax.set_ylabel(f'Percentage of AO_sc > {threshold}', fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=18)
ax.grid(False)

# Adjust layout
plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1, wspace=0.4, hspace=0.6)
plt.savefig(f"{timestamped_output_path}/composite_plot_adjusted_titles.pdf")
plt.show()

print(f"Composite plot saved as {timestamped_output_path}/composite_plot_adjusted_titles.pdf")
