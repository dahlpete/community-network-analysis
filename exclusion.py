# This function identifies the set of atoms to exclude from the list of interacting atoms
# Default exclusion condition: > 5 Angstroms away for >= 70 % of simultion
import numpy as np
import csv
import time
import MDAnalysis as mda
from MDAnalysis.analysis import distances as dst

def exclusion(u,numframes,dist_thresh=5.0,selection='all'):
	start = time.time()
	u.trajectory[0]
	atom_selection = u.select_atoms(selection)
	natoms = len(atom_selection)
	#indices = atom_selection._ix

	fileOUT = open('notfive.txt','w')
	#dist_array = np.zeros([natoms,natoms,numframes])
	dist_array = np.zeros([natoms,natoms])
	count_array = np.zeros([natoms,natoms])
	temp = np.zeros([natoms,natoms])
	#mem_usg = dist_array.nbytes * 1e-9
	#print('memory usage: %s GB' % mem_usg)
	print('starting traj loop')
	for i,ts in enumerate(u.trajectory):
		atom_coords = atom_selection.positions
		dist_array = dst.distance_array(atom_coords,atom_coords)
		gt5 = dist_array > dist_thresh
		temp[gt5] = 1.0
		count_array += temp
	count_array /= numframes	
	
	#mem_usg = dist_array.nbytes * 1e-9
	#print('memory usage: %s GB' % mem_usg)
	print('traj loop done')
	#gt5 = dist_array > dist_thresh
	#count_array = np.zeros([natoms,natoms,numframes])
	#count_array[gt5] = 1.0
	#count_sum = np.sum(count_array,axis=2)
	#count_sum = count_sum / numframes
	row_idx,col_idx = np.where(count_array >= 0.7)
	exclude_coordinates = list(zip(row_idx,col_idx))	
	exclude_coordinates = np.array(exclude_coordinates)
	for i in range(natoms):
		coi = np.where(exclude_coordinates[:,0] == i)
		values = exclude_coordinates[coi,1]
		fileOUT.writelines('%s ' % item for item in values[0])
		fileOUT.write('\n')
	
	end = time.time()
	total_time = end - start
	print('Time to write exclusion list: %.2f seconds' % total_time)
