import MDAnalysis as mda
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-p',
                    '--pdb_inp',
                    type=str,
                    help='input pdb with whole system')

parser.add_argument('-t',
                    '--traj_inp',
                    type=str,
                    help='trajectory file (dcd format)')

parser.add_argument('-op',
                    '--outpdb',
                    type=str,
                    help='output pdb')
                                             
parser.add_argument('-ot',
                    '--outtraj',
                    type=str,
                    help='output traj (dcd)')

args = parser.parse_args()


# Load the universe with your topology and trajectory files
u = mda.Universe(args.pdb_inp, args.traj_inp)

# Select only the protein
protein = u.select_atoms('protein')

# Write the protein-only PDB file
with mda.Writer(args.outpdb, protein.n_atoms) as pdb_out:
    pdb_out.write(protein)

# Write the protein-only trajectory as a DCD file
with mda.Writer(args.outtraj, protein.n_atoms) as dcd_out:
    for ts in u.trajectory:
        dcd_out.write(protein)

