import os
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for headless mode
import matplotlib.pyplot as plt

def extract_time_from_file(filename):
    """Extract times from a file where lines contain elapsed time values ending in 'seconds'."""
    print(f"Extracting times from {filename}...")
    times = []
    with open(filename, 'r') as file:
        for line in file:
            print(f"Checking line: {line.strip()}")  # Debug output
            match = re.search(r'([-+]?[0-9]*\.?[0-9]+)\s*seconds$', line)
            if match:
                times.append(float(match.group(1)))
    print(f"Extracted {len(times)} time values from {filename}.")
    return np.array(times)

def process_files_in_directory(directory):
    """Process all .out files in the directory, extracting times and computing averages and sums."""
    print(f"Processing .out files in {directory}...")
    time_data_mean = {}
    time_data_sum = {}
    processor_counts = []
    
    for filename in os.listdir(directory):
        if filename.endswith(".out"):
            match = re.search(r'(\d+)\.out$', filename)
            if match:
                num_processors = int(match.group(1))
                print(f"Processing {filename} for {num_processors} processors...")
                processor_counts.append(num_processors)
                
                filepath = os.path.join(directory, filename)
                times = extract_time_from_file(filepath)
                
                if times.size > 0:
                    time_data_mean[num_processors] = np.mean(times)
                    time_data_sum[num_processors] = np.sum(times)
                    print(f"Computed mean: {time_data_mean[num_processors]:.5f}, sum: {time_data_sum[num_processors]:.5f} for {num_processors} processors.")
    
    return time_data_mean, time_data_sum

def plot_performance(time_data_mean, time_data_sum):
    """Plot the averaged time, cumulative time, speedup, and efficiency vs. number of processors."""
    print("Generating performance plots...")
    sorted_processors = sorted(time_data_mean.keys())
    avg_times = [time_data_mean[p] for p in sorted_processors]
    cum_times = [time_data_sum[p] for p in sorted_processors]
    
    # Identify T(1) using the smallest processor count available
    min_processors = min(sorted_processors)
    T1 = time_data_mean[min_processors] * min_processors
    
    speedup = [T1 / time_data_mean[p] for p in sorted_processors]
    efficiency = [speedup[i] / sorted_processors[i] for i in range(len(sorted_processors))]
    
    # Plot Averaged Time
    plt.figure(figsize=(8,6))
    plt.plot(sorted_processors, avg_times, marker='o', linestyle='-', color='b')
    plt.xscale('linear')
    plt.yscale('log')
    plt.xlabel("Processors")
    plt.ylabel("Time (s)")
    plt.title("Averaged Walking Search Performance (for 1e4 particles)")
    plt.grid(True, which="both", linestyle='--', linewidth=0.5)
    plt.savefig("performance_plot_avg.png", dpi=300)
    print("Saved averaged performance plot as performance_plot_avg.png")
    
    # Plot Cumulative Time
    plt.figure(figsize=(8,6))
    plt.plot(sorted_processors, cum_times, marker='o', linestyle='-', color='purple')
    plt.xscale('linear')
    plt.yscale('log')
    plt.xlabel("Processors")
    plt.ylabel("Time (s)")
    plt.title("Cumulative Walking Search Performance (for 1e4 particles)")
    plt.grid(True, which="both", linestyle='--', linewidth=0.5)
    plt.savefig("performance_plot_cum.png", dpi=300)
    print("Saved cumulative performance plot as performance_plot_cum.png")
    
    # Plot Speedup
    plt.figure(figsize=(8,6))
    plt.plot(sorted_processors, speedup, marker='o', linestyle='-', color='g')
    plt.xscale('linear')
    plt.xlabel("Processors")
    plt.ylabel("Speedup")
    plt.title("Speedup of Walking Search Performance")
    plt.grid(True, linestyle='--', linewidth=0.5)
    plt.savefig("performance_plot_speedup.png", dpi=300)
    print("Saved speedup plot as performance_plot_speedup.png")
    
    # Plot Efficiency
    plt.figure(figsize=(8,6))
    plt.plot(sorted_processors, efficiency, marker='o', linestyle='-', color='r')
    plt.xscale('linear')
    plt.xlabel("Processors")
    plt.ylabel("Efficiency")
    plt.title("Efficiency of Walking Search Performance")
    plt.grid(True, linestyle='--', linewidth=0.5)
    plt.savefig("performance_plot_efficiency.png", dpi=300)
    print("Saved efficiency plot as performance_plot_efficiency.png")

def main():
    directory = os.getcwd()  # Use the current directory
    print("Starting main process...")
    time_data_mean, time_data_sum = process_files_in_directory(directory)
    if time_data_mean:
        plot_performance(time_data_mean, time_data_sum)
    else:
        print("No valid data found.")
    print("Process completed.")

if __name__ == "__main__":
    main()
