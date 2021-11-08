# This script computes the mutual information (MI) from the covariance between atomic displacements

import numpy as np
import time
from numpy import linalg as LA
import csv
import exclusion as ex
import MDAnalysis as mda
import MDAnalysis.analysis.encore.covariance as covar
import MDAnalysis.analysis.pca as pca
import covariance as covar

#_______________________________________
# INPUTS
structure_file = '/gpfs/scratch60/fas/batista/pd455/omcs_structure/2chains/fully_ox/CNA/matrix_calculation/extract_protein/myprot_reference_ox.psf'
dcd_trajectory = '/gpfs/scratch60/fas/batista/pd455/omcs_structure/2chains/fully_ox/CNA/matrix_calculation/extract_protein/trunc_prot_only_traj_ox.dcd'
#______________________________________


def MI(covar_matrix,natoms):
	# This function computes the linear mutual information matrix
	print('computing linear mutual information matrix\n')
	twobodycorr = np.zeros([natoms,natoms])
	for i in range(0,natoms):
		temp = mat_det(covar_matrix,i,i,3)
		for j in range(0,i):
			temp1 = mat_det(covar_matrix,j,j,3)
			temp2 = mat_det(covar_matrix,i,j,6)
			temp3 = temp*temp1/temp2
			temp3 = np.log(temp3)
			twobodycorr[i,j]=temp3
			twobodycorr[j,i] = twobodycorr[i,j]
	twobodycorr = 0.5 * twobodycorr

	return twobodycorr


def gen_corr(LMI_matrix,natoms,cen=None):
	# This function converts the linear mutual information into the generalized
	# correlation coefficient
	print('computing general correlation matrix\n')
	fileOUT = open('g_corr_list.txt','w')
	matOUT = open('g_corr_mat.txt','w')

	for j in range(0,natoms):
		for i in range(0,j):
			LMI_matrix[i,j] = 1.0 - np.exp((-2.0/3.0) * LMI_matrix[i,j])
			LMI_matrix[i,j] = np.sqrt(LMI_matrix[i,j])
			LMI_matrix[j,i] = LMI_matrix[i,j]
			print('%s	%s	%.3f' % (i,j,LMI_matrix[i,j]),file=fileOUT)
		LMI_matrix[j,j] = 1.0

	# Write general correlation matrix
	print("%s\n" % natoms,file=matOUT)
	with matOUT as f:
		csv_writer = csv.writer(f,delimiter = ' ')
		csv_writer.writerows(LMI_matrix)

	matOUT.close()

	if cen:
		centrality(LMI_matrix,natoms)
		
	
			
def mat_det(covar_matrix,ii,jj,n):	
	matrix = np.zeros([n,n])
	l = 1
	if n == 6:
		for i in range(0,3):
			for j in range(0,3):
				matrix[i,j]=covar_matrix[ii*3+i,ii*3+j]
				matrix[i,j+3]=covar_matrix[ii*3+i,jj*3+j]
				matrix[i+3,j]=covar_matrix[jj*3+i,ii*3+j]
				matrix[i+3,j+3]=covar_matrix[jj*3+i,jj*3+j]
	elif n == 3:
		for i in range(0,3):
			for j in range(0,3):
				matrix[i,j]=covar_matrix[ii*3+i,jj*3+j]
	
	determinant = LA.det(matrix)
	return determinant


def centrality(g_corr_mat,natoms):
	cenOUT = open('eigvec_centrality.txt','w')
	
	[evalues,evectors] = LA.eig(g_corr_mat)
	idx = evalues.argsort()[::-1]
	evalues = evalues[idx]
	print(evalues[0])
	evectors = abs(evectors[:,idx])
	print(evectors[0:10,0])
	
	for i in range(natoms):
		sum_array = 0
		for j in range(natoms):
			sum_array += g_corr_mat[i,j]*evectors[j,0]
		cen_value = (1.0/evalues[0]) * sum_array
		print("%s       %s" % (i,cen_value), file=cenOUT)	
	

start = time.time()

## Use MDAnalysis to compute the covariance matrix from your trajectory
u = mda.Universe(structure_file,dcd_trajectory)
nframes = len(u.trajectory)

# compute the covariance matrix
print('computing covariance matrix\n')
selection = 'all'
natoms = len(u.select_atoms(selection))
#print(natoms)
covar_matrix = covar.covariance(u,structure_file,selection)

# compute the mutual information matrix from the covariance matrix
LMI = MI(covar_matrix,natoms)

# compute and output the general correlation matrix
gen_corr(LMI,natoms,cen=True)

# compute the exclusion list
ex.exclusion(u,nframes,selection=selection)

end = time.time()
print('Done')
print('WallClock: %.6f' % (end-start))
