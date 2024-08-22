import itertools
from pathlib import Path
from typing import TYPE_CHECKING, Dict, List

import MDAnalysis as mda
import numpy as np
from tqdm import tqdm
from mdlearn.data.preprocess.align import iterative_means_align
from mdlearn.data.preprocess.decorrelation.spatial import SD2, SD4
from mdlearn.nn.models.ae.linear import LinearAETrainer
from mdlearn.utils import PathLike
import argparse

if TYPE_CHECKING:
    import numpy.typing as npt

def calc_dihedral(p0, p1, p2, p3):
    """Praxeolitic formula
    1 sqrt, 1 cross product"""

    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.arctan2(y, x)#np.degrees(np.arctan2(y, x))

def _calc_dihedral(u1, u2, u3, u4):
    """ Calculate dihedral angle method. From bioPython.PDB
    (adapted to np.array)
    Calculate the dihedral angle between 4 vectors
    representing 4 connected points. The angle is in
    [-pi, pi].
    """

    a1 = u2 - u1
    a2 = u3 - u2
    a3 = u4 - u3

    v1 = np.cross(a1, a2)
    v1 = v1 / (v1 * v1).sum(-1)**0.5
    v2 = np.cross(a2, a3)
    v2 = v2 / (v2 * v2).sum(-1)**0.5
    porm = np.sign((v1 * a3).sum(-1))
    rad = np.arccos((v1*v2).sum(-1) / ((v1**2).sum(-1) * (v2**2).sum(-1))**0.5)
    if not porm == 0:
        rad = rad * porm

    return rad

def calc_phi_psis (coords, names, resids):
    import math
    #print(resids)
    uniq_resids = set(np.array(resids)[0])
    #print(len(uniq_resids))
    #print(coords.shape[0])
    dihedral_data = np.zeros((coords.shape[0], 4, len(uniq_resids)-2))
    #psi_arr = np.zeros((coords.shape[0], 4, len(uniq_resids)-2))
    for t, coord_t in tqdm(enumerate(coords)):
        coord_t = coord_t.T
        names_t = names[t]
        for cit in range (4, len(coord_t)-4):
           if cit%4==0:
                '''
                Determining the atoms N, CA, C, Npositive
                '''
                c_neg_atom = coord_t[cit-2]
                n_atom = coord_t[cit]
                calph_atom = coord_t[cit+1]
                c_atom = coord_t[cit+2]
                n_plus_atom = coord_t[cit+4]
                #print(f"Getting phi, psi with: \n \
                #        {names_t[cit-2]}\n \
                #        {names_t[cit]}\n \
                #        {names_t[cit+1]}\n \
                #        {names_t[cit+2]}\n \
                #        {names_t[cit+4]}")

                '''
                Plug into calc_dihedral to get dihedral at int(cit/4) position
                '''
                phi_ang = calc_dihedral(c_neg_atom, n_atom, calph_atom, c_atom)
                psi_ang = calc_dihedral(n_atom, calph_atom, c_atom, n_plus_atom)
                
                dihedral_data[t][0][int(cit/4)-1] = math.cos(phi_ang)
                dihedral_data[t][1][int(cit/4)-1] = math.sin(phi_ang)
                dihedral_data[t][2][int(cit/4)-1] = math.cos(psi_ang)
                dihedral_data[t][3][int(cit/4)-1] = math.sin(psi_ang)

                #phi_arr[t][int(cit/4)-1] = calc_dihedral(c_neg_atom, n_atom, calph_atom, c_atom)
                #psi_arr[t][int(cit/4)-1] = calc_dihedral(n_atom, calph_atom, c_atom, n_plus_atom)
                #print(phi_arr)
                #print(psi_arr)
    print(dihedral_data.shape)
    return dihedral_data#phi_arr, psi_arr


def parse_position_groups(
    pdb_file: PathLike, traj_file: PathLike, selections: List[str]
) -> Dict[str, "npt.ArrayLike"]:
    u = mda.Universe(str(pdb_file), str(traj_file))
    atom_groups = {sel: u.select_atoms(sel) for sel in selections}
    #print(atom_groups)
    position_groups = {sel: [] for sel in selections}
    name_groups = {sel: [] for sel in selections}
    resid_groups = {sel: [] for sel in selections}
    for _ in u.trajectory:#[0:100]:#40000]:
        for sel in selections:
            names = atom_groups[sel].names
            resids = atom_groups[sel].resids
            positions = atom_groups[sel].positions.transpose()
            position_groups[sel].append(positions)
            name_groups[sel].append(names)
            resid_groups[sel].append(resids)
    # Convert from list to np.array
    position_groups = {k: np.array(v) for k, v in tqdm(position_groups.items())}
    name_groups = {k: np.array(v) for k, v in tqdm(name_groups.items())}
    resid_groups = {k: np.array(v) for k, v in tqdm(resid_groups.items())}
    return position_groups, name_groups, resid_groups


def main(pdb_file,
         traj_file,
         selection_phrase,
         outpath,
        ):

    selections = [selection_phrase] ## one off on both ends for dihedral calcs
    alignment_workers = 10
    num_sd4_subspaces = 10
    latent_dim = 8
    autoencoder_params = dict(
        input_dim=num_sd4_subspaces,
        latent_dim=latent_dim,
        hidden_neurons=[32, 16],
        epochs=500,
        checkpoint_log_every=500,
        plot_method=None,
        device="gpu",
    )
    output_path = Path(outpath)
    # Autoencoder outputs for each selection
    selection_run_path = output_path / "selection_runs"
    # Autoencoder output for the aggregate selection
    aggregate_run_path = output_path / "aggregate_run"

    try:
        output_path.mkdir()
    except:
        pass

    try:
        selection_run_path.mkdir()
    except:
        pass
    
    try:
        aggregate_run_path.mkdir()
    except:
        pass

    # Parse atom group positions from trajectory file
    position_groups, name_groups, resid_groups = parse_position_groups(pdb_file, traj_file, selections)
    all_positions = np.concatenate(list(position_groups.values()), axis=2)
    name_vals = name_groups[selections[0]]
    resid_vals = resid_groups[selections[0]]

    # Run alignment of all positions and save results
    _, _, e_rmsd, aligned_coords = iterative_means_align(
        all_positions, num_workers=alignment_workers
    )
    np.save(output_path / "e_rmsd.npy", e_rmsd)
    np.save(output_path / "all_positions.npy", all_positions)
    np.save(output_path / "aligned_coords.npy", aligned_coords)
                                                                                                           
    group_lengths = [positions.shape[2] for positions in position_groups.values()]
    start_ind = 0
    aligned_position_groups = {}
    for sel, length in zip(selections, group_lengths):
        aligned_position_groups[sel] = position_groups[sel][
            :, :, start_ind : start_ind + length
        ]
        start_ind += length
    # Compute and store the latent space (num_frames, latent_dim) for each selection

    latent_spaces = {}
    for i, (selection, positions) in enumerate(aligned_position_groups.items()):
        dihedral_data = calc_phi_psis (positions, name_vals, resid_vals)
        selection_path = selection_run_path / f"selection-{i}"
        selection_path.mkdir()
        # Log the selection string to document the model directories
        with open(selection_path / "selection.txt", "w") as f:
            f.write(selection)
                                                                                                           
        np.save(selection_path / "dihedrals.npy", dihedral_data)
        num_frames, _, num_atoms = dihedral_data.shape
        # Run SD2 and save results
        Y, S, B, U = SD2(dihedral_data.reshape(num_frames, num_atoms * 4), m=num_atoms * 4)
        np.save(selection_path / "Y.npy", Y)
        np.save(selection_path / "S.npy", S)
        np.save(selection_path / "B.npy", B)
        np.save(selection_path / "U.npy", U)
        # Run SD4 and save results
        W = SD4(
            Y[0:num_sd4_subspaces, :],
            m=num_sd4_subspaces,
            U=U[0:num_sd4_subspaces, :],
        )
        coords = np.reshape(dihedral_data, (num_frames, 4 * num_atoms)).T
        avg_coords = np.mean(coords, 1)
        tmp = np.reshape(np.tile(avg_coords, num_frames), (num_frames, 4 * num_atoms)).T
        devs = coords#- tmp
        sd4_projection: "npt.ArrayLike" = W.dot(devs).T
        np.save(selection_path / "W.npy", W)
        np.save(selection_path / "sd4_projection.npy", sd4_projection)
        # At this point, the SD4 projection results in an array with
        # shape (num_frames, num_sd4_subspaces)
        # Run autoencoder and save results
        autoencoder = LinearAETrainer(**autoencoder_params)
        autoencoder.fit(sd4_projection, output_path=selection_path / "autoencoder")
        z, _ = autoencoder.predict(sd4_projection)
        latent_spaces[selection] = z
        np.save(selection_path / "z.npy", z)
    return
   


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('-p',
                        '--pdbinp',
                        type=str,
                        help='input pdb with whole system')

    parser.add_argument('-t',
                        '--trajinp',
                        type=str,
                        help='trajectory file (dcd format)')

    parser.add_argument('-s',
                        '--selphrase',
                        type=str,
                        help='output pdb')

    parser.add_argument('-o',
                        '--out',
                        type=str,
                        help='output traj (dcd)')

    args = parser.parse_args()


    main(args.pdbinp,
         args.trajinp,
         args.selphrase,
         args.out,
        )


