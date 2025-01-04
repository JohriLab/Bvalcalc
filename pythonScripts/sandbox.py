from concurrent.futures import ProcessPoolExecutor
import os

# Function to be run in parallel
def compute(n):
    return n * n

if __name__ == '__main__':
    # Define your input data
    data = range(1,11)

    # Use a number of workers suitable for the system
    num_workers = os.cpu_count()  # Dynamically determine available cores

    # Use ProcessPoolExecutor for cross-platform parallelism
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        results = list(executor.map(compute, data))

    print("Results:", results)
