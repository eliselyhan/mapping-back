#used tabs instead of spaces

import pandas as pd
import numpy as np
import starfile
from ipdm import *
from filter_distances import *
from eulerangles import euler2matrix
from polysomes_from_distance_matrix import *
from distance_matrix import *
from disjoint_set import DisjointSet
import argparse
import os

parser = argparse.ArgumentParser()

parser.add_argument("-f", "--filename", type=str, required=True, help="star file with the ribosome data")
parser.add_argument("-c", "--center_mass", nargs=3, type=float, required=True, help="x,y,z coordinates of the subtomogram-averaged ribosome's center of mass")
parser.add_argument("-e", "--entry", nargs=3, type=float, required=True, help="x,y,z coordinates of the ribosome's entry point")
parser.add_argument("-x", "--exit", nargs=3, type=float, required=True, help="x,y,z coordinates of the ribosome's exit point")

args = parser.parse_args()

filename = args.filename
center_mass = args.center_mass
entry = args.entry
exit = args.exit

# load data with duplicate ribosomes removed
df_all = starfile.read(filename)
df_particles = df_all['particles']

# names of tomograms
tomoNames = pd.unique(df_particles['rlnMicrographName'])
#print("tomoNames: {}".format(tomoNames))

start_index = 0;

# lines in dataframe that lie in polysomes
df_all_poly = pd.DataFrame();
#df_all_distinct = pd.DataFrame();

# angstroms per pixel
# we are doing everything in angstroms now
angPix = 4;

# center of box
box_shift = np.array([35,35,35])*angPix;

# center of mass of ribosome
cm_avg = np.array(center_mass);

# iterate through each tomogram
for tomogram in tomoNames:
	df = df_particles[df_particles['rlnMicrographName'] == tomogram ]
	x_cm = df['rlnCoordinateX'].values
	y_cm = df['rlnCoordinateY'].values
	z_cm = df['rlnCoordinateZ'].values

	# n x 3 matrix of center of mass coordinates of ribosomes in the overall tomogram
	cm = np.column_stack((x_cm, y_cm, z_cm))*angPix
	n = len(x_cm)
	#removed = []
	'''distance_matrix = get_distances(cm, cm);
	duplicates_indices = filter_duplicates(distance_matrix, 2.0);

	kept =[]
	if (len(duplicates_indices[0]) > 0):
		ds = DisjointSet()
		rows = set(duplicates_indices[0])
		columns = set(duplicates_indices[1])
		for row in rows: 
			ds.find(row);
		for col in columns: 
			ds.find(col);
		for i in range(len(duplicates_indices[0])):
			ds.union(duplicates_indices[0][i], duplicates_indices[1][i]);

		dup_sets = list(ds.itersets())
		total = 0
		for s in dup_sets:
			total += len(s)
			counter = 0;
			for val in s:
				if counter < len(s)-1:
					removed.append(val) 
					counter +=1
				else:
					kept.append(val)

		df = df.iloc[kept,:]
		df_all_distinct = df_all_distinct.append(df)
	else:
		df_all_distinct = df_all_distinct.append(df)
	starfile.write(df_all_distinct, 'run_data_distinct.star', overwrite =True)'''
	end_index = start_index + n;
	#exit_shift = np.array([327.52, 260.55, 304.01]) - box_shift
	#entry_shift = np.array([228.33, 317.68, 269.50]) - box_shift

	# make center of the box the new origin
	exit_angst = np.array(exit) - box_shift
	entry_angst = np.array(entry) - box_shift
	
	exit_adjusted = np.zeros((n,3))
	entry_adjusted = np.zeros((n,3))

	count = 0;

	for i in range(start_index, end_index):
		
		#test: adding shiftAngst
		shiftAngst = df.loc[i, ['rlnOriginXAngst', 'rlnOriginYAngst', 'rlnOriginZAngst']].to_numpy()
		cm[count, :] = cm[count, :] - shiftAngst;

		eulers_relion = df.loc[i, ['rlnAngleRot', 'rlnAngleTilt', 'rlnAnglePsi']].tolist()
		#print(eulers_relion)
		eulers_radians = [np.radians(x) for x in eulers_relion]
		phi, theta, psi = eulers_radians
		rot_phi = np.array([[np.cos(phi), -np.sin(phi), 0],
							[np.sin(phi), np.cos(phi), 0],
							[0, 0, 1]])
		rot_theta = np.array([[np.cos(theta), 0, np.sin(theta)],
								[0, 1, 0],
								[-np.sin(theta),0,np.cos(theta)]])
		rot_psi = np.array([[np.cos(psi), -np.sin(psi), 0],
							[np.sin(psi), np.cos(psi), 0],
							[0, 0, 1]])
		#rotm = euler2matrix(eulers_relion, axes='zyz', intrinsic=True, right_handed_rotation=True)

		rotm = np.dot(rot_phi, np.dot(rot_theta, rot_psi))
		rotm = rotm.T;
		exit_rot = np.dot(rotm, exit_angst) + box_shift
		entry_rot = np.dot(rotm, entry_angst) + box_shift

		cm_rot = cm_avg - box_shift;
		cm_rot = np.dot(rotm, cm_rot)
		cm_rot = cm_rot + box_shift

		# vector from rotated center of mass to the entry/exit
		exit_shift = exit_rot - cm_rot
		entry_shift = entry_rot - cm_rot
		exit_adjusted[count,:] = cm[count,:] + exit_shift
		entry_adjusted[count,:] = cm[count,:] + entry_shift
		#print('cm{}, exit{}, entry{}'.format(cm[count,:], exit_adjusted[count,:], entry_adjusted[count,:]));

		count += 1; 
	

	adjusted_coordinates = np.column_stack((entry_adjusted, exit_adjusted))

	# save entry and exit coordinates to csv
	dir_name = filename.rpartition(".")[0].split("/")[-1] + "_tomogram_entry_exit/"
	#print(dir_name)
	os.makedirs(dir_name, exist_ok=True)
	np.savetxt(dir_name + tomogram + ".csv", adjusted_coordinates, delimiter=",")

	start_index = end_index

	entry_x = entry_adjusted[:,0]
	entry_y = entry_adjusted[:,1]
	entry_z = entry_adjusted[:,2]
	exit_x = exit_adjusted[:,0]
	exit_y = exit_adjusted[:,1]
	exit_z = exit_adjusted[:,2]
	

	x_distances = ipdm_1d(entry_x, exit_x)
	y_distances = ipdm_1d(entry_y, exit_y)
	z_distances = ipdm_1d(entry_z, exit_z)

	ipdm_matrix = ipdm_3d(x_distances, y_distances, z_distances)
	filtered_matrix = find_smaller_distances(ipdm_matrix)
	polysome_matrix = apply_threshold(filtered_matrix)
	polysomes = find_polysomes(polysome_matrix)
	indices = get_ribosomes_in_polysomes(polysomes)

	if len(indices) > 0:
		#print(tomogram, indices)
		df_polysomes = df.iloc[indices, :]
		df_all_poly = df_all_poly.append(df_polysomes)
	#print(df_all_poly)

#starfile.write(df_all_poly, 'run_data_polysomes.star',overwrite=True)
