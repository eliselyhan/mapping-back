import pandas as pd
import numpy as np
import starfile
from distance_matrix import *
from disjoint_set import DisjointSet
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-f", "--filename", type=str, required=True, help="star filename")
parser.add_argument("-t", "--threshold", type=float, required=True, help="number of pixels threshold for removing duplicates (e.g., 30)")

args = parser.parse_args()

filename = args.filename
threshold = args.threshold

df_all = starfile.read(filename)
df_particles = df_all['particles']
df_all_distinct = pd.DataFrame()

tomoNames = pd.unique(df_particles['rlnMicrographName'])

for tomogram in tomoNames:
	df = df_particles[df_particles['rlnMicrographName'] == tomogram ]
	x_cm = df['rlnCoordinateX'].values
	y_cm = df['rlnCoordinateY'].values
	z_cm = df['rlnCoordinateZ'].values

	# n x 3 matrix of center of mass coordinates
	cm = np.column_stack((x_cm, y_cm, z_cm))
	n = len(x_cm)
	removed = []

	distance_matrix = get_distances(cm, cm);
	duplicates_indices = filter_duplicates(distance_matrix, threshold);

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

starfile.write(df_all_distinct, filename.split(".")[0] +'_distinct.star', overwrite =True)
