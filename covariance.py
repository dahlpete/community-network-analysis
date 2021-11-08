import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align

def covariance(u,structure_file,selection='all'):
	u.trajectory[0]
	reference = u.select_atoms(selection)
	reference.write('ref_frame.pdb')
	ref_u = mda.Universe(structure_file,'ref_frame.pdb')
	reference_pos = reference.positions
	atoms = u.select_atoms(selection)
	align_sel = 'backbone or name FE CHA CHB CHC CHD C1A C2A C3A C4A C1B C2B C3B C4B C1C C2C C3C C4C C1D C2D C3D C4D NA NB NC ND'
	alignment = align.AlignTraj(u,ref_u,select=align_sel,filename='rmsfit.dcd')
	alignment.run()

	u = mda.Universe(structure_file,'rmsfit.dcd')

	natoms = len(atoms)
	nframes = len(u.trajectory)

	
	# initialize variables
	mean_struc = np.zeros(natoms*3)
	covariance  = np.zeros([natoms*3,natoms*3])

	# Compute the mean atomic positions
	for i, ts in enumerate(u.trajectory[:]):
		
		atoms = u.select_atoms(selection)
		mean_struc += atoms.positions.ravel()
	mean_struc /= nframes


	# Compute the covariance matrix
	for i, ts in enumerate(u.trajectory[:]):
		x = atoms.positions.ravel()
		x -= mean_struc
		covariance += np.dot(x[:,np.newaxis], x[:, np.newaxis].T)

	covariance /= nframes-1

	return covariance
	

	
