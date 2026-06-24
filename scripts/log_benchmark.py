import csv
import os
from datetime import datetime

OUTPUT_DIR = "data/benchmark_results"
OUTPUT_FILE = os.path.join(OUTPUT_DIR, "benchmark_results.csv")
HEADERS = ["timestamp", "grid_size", "kernel_name", "time_ms", "occupancy", "memory_bw"]

def log(grid_size, kernel_name, time_ms, occupancy, memory_bw):
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    file_exists = os.path.isfile(OUTPUT_FILE)
    with open(OUTPUT_FILE, "a", newline="") as f:
        writer = csv.writer(f)
        if not file_exists:
            writer.writerow(HEADERS)
        writer.writerow([datetime.now().isoformat(), grid_size, kernel_name, time_ms, occupancy, memory_bw])

#test 
if __name__ == "__main__":
    log(1000, "x_sweep_baseline", 3.2, 0.75, 210.5)