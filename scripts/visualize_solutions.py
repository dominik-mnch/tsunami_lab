#!/usr/bin/env python3
"""
Script to visualize tsunami lab solutions as an animated GIF.
Reads all solution_*.csv files and creates a plot of height vs x coordinate
for each time step, then combines them into an animated GIF.
"""

import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import imageio
from pathlib import Path
import numpy as np

# Get the project root directory (parent of the scripts directory)
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def get_solution_files(solutions_dir=None):
    """
    Get all solution files sorted by number.
    """
    if solutions_dir is None:
        solutions_dir = os.path.join(PROJECT_ROOT, "solutions")
    
    print(f"Looking for solution files in: {solutions_dir}")
    print(f"Directory exists: {os.path.exists(solutions_dir)}")
    
    if os.path.exists(solutions_dir):
        print(f"Directory contents: {os.listdir(solutions_dir)}")
    
    pattern = os.path.join(solutions_dir, "solution_*.csv")
    print(f"Using glob pattern: {pattern}")
    
    files = glob.glob(pattern)
    print(f"Found {len(files)} solution files")
    
    # Sort by the numeric part in the filename
    files.sort(key=lambda x: int(Path(x).stem.split('_')[1]))
    
    return files

def read_solution_file(filepath):
    """
    Read a solution CSV file and return x, height, bathymetry data, and simulation time.
    """
    # Read CSV data normally
    df = pd.read_csv(filepath)

    # Sort by x coordinate to ensure proper plotting
    df = df.sort_values('x')

    # Extract simulation time (same for all rows)
    sim_time = None
    if 'time' in df.columns:
        sim_time = df['time'].iloc[0]

    if 'bathymetry' in df.columns:
        bathymetry = df['bathymetry'].values
    else:
        # Backward compatibility for older solution files.
        bathymetry = np.zeros_like(df['height'].values)

    return df['x'].values, df['height'].values, bathymetry, sim_time

def create_plot_image(x_data, height_data, bathymetry_data, time_step, sim_time, output_path):
    """
    Create a plot of water column, water surface, and seafloor vs x and save as a PNG image.
    """
    fig, ax = plt.subplots(figsize=(12, 6))
    surface_data = height_data + bathymetry_data
    
    ax.plot(x_data, surface_data, color='tab:blue', linewidth=2, label='Water Surface')
    ax.plot(x_data, bathymetry_data, color='saddlebrown', linewidth=2, label='Seafloor')
    ax.fill_between(x_data, bathymetry_data, surface_data, color='tab:blue', alpha=0.3)
    ax.fill_between(x_data, bathymetry_data, np.min(bathymetry_data) - 0.05 * max(1.0, np.max(height_data)), color='peru', alpha=0.2)
    
    ax.set_xlabel('X Coordinate (m)', fontsize=12)
    ax.set_ylabel('Height (m)', fontsize=12)
    if sim_time is not None:
        ax.set_title(f'Water Surface and Seafloor at t = {sim_time:.2f} s', fontsize=14)
    else:
        ax.set_title(f'Water Surface and Seafloor at Time Step {time_step}', fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right')
    
    # Set y-axis to show both the seafloor and the water surface.
    y_min = np.min(bathymetry_data)
    y_max = np.max(surface_data)
    y_pad = max(1e-6, 0.1 * (y_max - y_min if y_max > y_min else 1.0))
    ax.set_ylim(bottom=y_min - y_pad, top=y_max + y_pad)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=100, bbox_inches='tight')
    plt.close(fig)
    
    return output_path

def create_gif(solutions_dir=None, output_gif="tsunami_animation.gif", temp_dir=".temp_frames"):
    """
    Main function to create an animated GIF from solution files.
    """
    if solutions_dir is None:
      solutions_dir = os.path.join(PROJECT_ROOT, "solutions")
    # Create temporary directory for frames
    os.makedirs(temp_dir, exist_ok=True)
    
    # Get all solution files
    solution_files = get_solution_files(solutions_dir)
    
    if not solution_files:
        print(f"No solution files found in {solutions_dir}")
        return
    
    print(f"Found {len(solution_files)} solution files")
    
    # Generate plots for each solution
    image_files = []
    for idx, filepath in enumerate(solution_files):
        print(f"Processing {os.path.basename(filepath)} ({idx + 1}/{len(solution_files)})")
        
        x_data, height_data, bathymetry_data, sim_time = read_solution_file(filepath)
        temp_image_path = os.path.join(temp_dir, f"frame_{idx:03d}.png")
        create_plot_image(x_data, height_data, bathymetry_data, idx, sim_time, temp_image_path)
        image_files.append(temp_image_path)
    
    # Create animated GIF
    print(f"Creating animated GIF: {output_gif}")
    images = [imageio.imread(img) for img in image_files]
    imageio.mimsave(output_gif, images, duration=0.1)
    
    # Clean up temporary frames
    import shutil
    shutil.rmtree(temp_dir)
    
    print(f"Animation saved to {output_gif}")

if __name__ == "__main__":
    create_gif()
