#!/usr/bin/env python3
"""
Script to read station CSV files and create plots for water height, momentum_x, and momentum_y over time.
"""

import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path


def find_first_arrival(times, heights, margin=0.05):
    """
    Return the time at which water height first exceeds initial value + margin.

    Parameters:
    -----------
    times : array-like
        Time values
    heights : array-like
        Water height values
    margin : float
        Threshold above initial water height to consider as arrival

    Returns the arrival time, or None if the threshold is never exceeded.
    """
    if len(heights) == 0:
        return None
    initial = heights.iloc[0]
    threshold = initial + margin
    mask = heights > threshold
    if not mask.any():
        return None
    return times[mask.idxmax()]


def plot_station_data(output_dir="./stations/output", save_dir="./stations/plots", arrival_margin=0.05):
    """
    Read all CSV files in the output directory and create plots for each station.
    
    Parameters:
    -----------
    output_dir : str
        Directory containing the station CSV files
    save_dir : str
        Directory where plots will be saved
    arrival_margin : float
        Minimum increase above initial water height (in metres) to detect tsunami arrival
    """
    
    # Create save directory if it doesn't exist
    Path(save_dir).mkdir(parents=True, exist_ok=True)
    
    # Find all CSV files in the output directory
    csv_files = glob.glob(os.path.join(output_dir, "*.csv"))
    
    if not csv_files:
        print(f"No CSV files found in {output_dir}")
        return
    
    print(f"Found {len(csv_files)} station files. Creating plots...")
    
    for csv_file in csv_files:
        station_name = Path(csv_file).stem
        
        try:
            # Read the CSV file
            df = pd.read_csv(csv_file)
            
            # Check if required columns exist
            required_cols = ['time', 'water_height', 'momentum_x', 'momentum_y']
            if not all(col in df.columns for col in required_cols):
                print(f"Warning: {station_name} missing required columns. Skipping...")
                continue
            
            # Create figure with single plot
            fig, ax = plt.subplots(figsize=(12, 6))
            fig.suptitle(f'Station: {station_name}', fontsize=16, fontweight='bold')
            
            # Plot all three variables on the same axes
            ax.plot(df['time'], df['water_height'], 'b-', linewidth=2, label='Water Height (m)')
            ax.plot(df['time'], df['momentum_x'], 'r-', linewidth=2, label='Momentum X (m²/s)')
            ax.plot(df['time'], df['momentum_y'], 'g-', linewidth=2, label='Momentum Y (m²/s)')

            arrival_time = find_first_arrival(df['time'], df['water_height'], margin=arrival_margin)
            if arrival_time is not None:
                ax.axvline(x=arrival_time, color='purple', linestyle='--', linewidth=1.5,
                           label=f'First arrival (t={arrival_time:.1f} s)')

            ax.set_xlabel('Time (s)', fontsize=12)
            ax.set_ylabel('Value', fontsize=12)
            ax.grid(True, alpha=0.3)
            ax.legend(loc='best', fontsize=11)
            ax.set_title('Water Height and Momenta over Time', fontsize=12)
            
            plt.tight_layout()
            
            # Save the plot
            output_path = os.path.join(save_dir, f"{station_name}_plot.png")
            plt.savefig(output_path, dpi=150, bbox_inches='tight')
            print(f"✓ Created plot for {station_name}: {output_path}")
            
            plt.close()
            
        except Exception as e:
            print(f"Error processing {station_name}: {e}")
            continue
    
    print(f"\nAll plots saved to {save_dir}")


if __name__ == "__main__":
    plot_station_data()
