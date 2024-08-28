import numpy as np
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt
from tqdm import tqdm

def calculate_free_energy_3d(data, bins=100, temperature=300):
    kB = 1.380649e-23  # Boltzmann constant in J/K
    kB_kcal_per_mol_per_K = kB * 2.39006e-4 * 6.022e23  # Convert to kcal/mol/K
    
    # Calculate the histogram and free energy
    histogram, edges = np.histogramdd(data, bins=bins, density=True)
    print(histogram)
    free_energy = -kB_kcal_per_mol_per_K * temperature * np.log(histogram + 1e-10)
    free_energy -= np.min(free_energy)  # Normalize to set min to 0
    
    return free_energy, edges

def interpolate_free_energy(data, edges, free_energy):
    #print(np.shape(edges[::-1, ::-1]))
    x_edges, y_edges, z_edges = edges#[::-1, ::-1]
    
    # Clip the data to avoid out-of-bounds error
    data[:, 0] = np.clip(data[:, 0], x_edges[0], x_edges[-2])
    data[:, 1] = np.clip(data[:, 1], y_edges[0], y_edges[-2])
    data[:, 2] = np.clip(data[:, 2], z_edges[0], z_edges[-2])
    
    interpolator = RegularGridInterpolator((x_edges[:-1], y_edges[:-1], z_edges[:-1]), free_energy)
    free_energy_values = interpolator(data)
    return free_energy_values

def show_scatter(data: np.ndarray, color: np.ndarray, output: str):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.set_facecolor("white")
    ff = ax.scatter(data[:, 0], data[:, 1], data[:, 2], c=color, vmin=0, vmax=2.5, cmap='inferno')
    #ax.view_init(90, 90, 90)
    plt.colorbar(ff)
    plt.savefig(output, dpi=300)
    plt.close()

def main(dihedrals, e_rmsd, ZPrj4, z, out_image):
    
    # Example usage

    for i in range(dihedrals.shape[1]):
        dihedral_test = dihedrals[:,i,0]
    
    plt.plot(dihedral_test[::100])
    plt.savefig(f'{out_image}/cos_dihedral_test.png', dpi=300)
    plt.close()
    
    print(np.shape(ZPrj4[:, :3]))
    free_energy, edges = calculate_free_energy_3d(ZPrj4[:, 0:3], bins=3)
    print(free_energy)
    print(edges)
    free_energy_values = interpolate_free_energy(ZPrj4[:, 0:3], edges, free_energy)
    # Mask out values above the threshold by setting them to NaN
    threshold = 2.5
    masked_free = np.where(free_energy_values > threshold, np.nan, free_energy_values)
    show_scatter(ZPrj4[:, 0:3], masked_free, f'{out_image}/ZPrj4.png')
    
    #z = np.load(zauto)#'run-4/selection_runs/selection-0/z.npy')
    free_energy, edges = calculate_free_energy_3d(z[:, 0:3], bins=3)
    print(free_energy)
    print(edges)
    free_energy_values = interpolate_free_energy(z[:, 0:3], edges, free_energy)

    show_scatter(z[:, 0:3], free_energy_values, f'{out_image}/z_free.png')

dih_data = np.load('output_total/selection_runs/selection-0/dihedrals.0.npy')
ermsd = np.load('output_total/e_rmsd.0.npy')
print(np.shape(ermsd))
for rep in range(1, 5):
    dih_data = np.concatenate((dih_data, np.load(f'output_total/selection_runs/selection-0/dihedrals.{rep}.npy')))
    #ermsd = np.concatenate((ermsd, np.load(f'output_total/e_rmsd.{rep}.npy')), axis=1)
    #print(np.shape(ermsd))

z = np.load('output_total/selection_runs/selection-0/z.0.npy')
sd4 = np.load('output_total/selection_runs/selection-0/sd4_projection.0.npy')
print(np.shape(z))
print(np.shape(sd4))

main(dih_data, ermsd, sd4, z, 'images')
