import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.interpolate import RegularGridInterpolator

def calculate_free_energy_nd(data, bins=100, temperature=300):
    kB = 1.380649e-23  # Boltzmann constant in J/K
    kB_kcal_per_mol_per_K = kB * 2.39006e-4 * 6.022e23  # Convert to kcal/mol/K

    # Estimate the probability density in N-dimensional space
    histogram, edges = np.histogramdd(data, bins=bins, density=True)
    
    # Calculate free energy in kcal/mol
    free_energy = -kB_kcal_per_mol_per_K * temperature * np.log(histogram + 1e-10)
    free_energy -= np.min(free_energy)  # Normalize to set min to 0
    
    return free_energy, edges

def project_free_energy_to_3d(free_energy, edges, dims=(0, 1, 2)):
    # Integrate/sum over the dimensions not selected for the projection
    for dim in tqdm(range(free_energy.ndim)):
        if dim not in dims:
            try:
                free_energy = np.sum(free_energy, axis=dim)
            except:
                continue
    # Get the 3D edges for the selected dimensions
    edges_3d = [edges[dim] for dim in dims]
    
    return free_energy, edges_3d

def interpolate_free_energy(data, edges, free_energy, dims=(0, 1, 2)):
    # Interpolating the 3D projected free energy
    x_edges, y_edges, z_edges = edges
    
    data_3d = data[:, dims]
    
    interpolator = RegularGridInterpolator((x_edges[:-1], y_edges[:-1], z_edges[:-1]), free_energy)
    
    free_energy_values = interpolator(data_3d)
    
    return free_energy_values

def show_scatter(data: np.ndarray, color: np.ndarray, output: str):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.set_facecolor("white")
    ff = ax.scatter(data[:, 0], data[:, 1], data[:, 2], c=color)
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
    free_energy, edges = calculate_free_energy_nd(ZPrj4, bins=10)
    free_energy, edges = project_free_energy_to_3d(free_energy, edges, dims=(0, 1, 2))
    print(free_energy)
    print(edges)
    free_energy_values = interpolate_free_energy(ZPrj4, edges, free_energy)
    show_scatter(ZPrj4[:, :3], free_energy_values, f'{out_image}/ZPrj4.png')
    
    #z = np.load(zauto)#'run-4/selection_runs/selection-0/z.npy')
    free_energy, edges = calculate_free_energy_3d(z[:, :3], bins=3)
    print(free_energy)
    print(edges)
    free_energy_values = interpolate_free_energy(z[:, :3], edges, free_energy)

    show_scatter(z[:, :3], free_energy_values, f'{out_image}/z_free.png')

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
