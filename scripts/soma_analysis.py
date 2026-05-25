import numpy as np
import csv

data = []
with open('/mnt/g/VSCodeProjekte/tsunami_lab/data/soma_epicenter_tohoku_2011_bathymetry.csv') as f:
    reader = csv.DictReader(f)
    for row in reader:
        data.append((float(row['z']), float(row['Points:0']), float(row['Points:1'])))

z = np.array([d[0] for d in data])
x = np.array([d[1] for d in data])
y = np.array([d[2] for d in data])

print('=== Cut overview ===')
print(f'Points: {len(z)}')
print(f'First: z={z[0]:.1f} m  x={x[0]:.0f} m  y={y[0]:.0f} m  (inland)')
print(f'Last:  z={z[-1]:.1f} m  x={x[-1]:.0f} m  y={y[-1]:.0f} m  (deep ocean)')

# --- Coast transition (z changes sign land -> ocean) ---
trans = np.where(np.diff(np.sign(z)))[0]
assert len(trans) >= 1, 'no coast crossing found'
c = trans[0]  # first land-to-ocean transition = Soma coast
print(f'\n=== Coast (Soma shoreline) ===')
print(f'Index {c}: z goes {z[c]:.1f} -> {z[c+1]:.1f} m')
print(f'Coordinates: x={x[c]:.0f} m, y={y[c]:.0f} m')

# --- Epicenter: origin (0,0) of the simulation coordinate system ---
# The simulation domain is centred so that the epicenter is at (0,0).
# Soma is ~128 km west and ~55 km south => epicenter at (+128000, +55000)
# relative to Soma, i.e. near (0,0) in the global system.
dist_to_origin = np.sqrt(x**2 + y**2)
epic_idx = int(np.argmin(dist_to_origin))
print(f'\n=== Epicenter (closest point to origin) ===')
print(f'Index {epic_idx}: z={z[epic_idx]:.1f} m  x={x[epic_idx]:.0f} m  y={y[epic_idx]:.0f} m')
print(f'Distance from origin: {dist_to_origin[epic_idx]/1000:.1f} km')

# --- Segment lengths ---
ds = np.sqrt(np.diff(x)**2 + np.diff(y)**2)

# --- Path length: epicenter -> Soma coast ---
i_start = min(c, epic_idx)
i_end   = max(c, epic_idx)
path_length = ds[i_start:i_end].sum()
print(f'\n=== Path epicenter -> Soma coast ===')
print(f'Path length: {path_length/1000:.1f} km')

# --- Travel time using water depth at the epicenter ---
g = 9.81
h_epic = -z[epic_idx]          # water depth at the epicenter (positive)
lam    = (g * h_epic) ** 0.5   # wave speed at the epicenter
travel_time = path_length / lam

print(f'\n=== Result ===')
print(f'Epicenter depth: {h_epic:.1f} m')
print(f'Wave speed:      {lam:.1f} m/s')
print(f'Travel time:     {travel_time:.0f} s  =  {travel_time/60:.1f} min')
