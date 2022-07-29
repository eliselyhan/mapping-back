import argparse
import pandas as pd
import numpy as np
import os
from ipdm import *
from filter_distances import *
from polysomes_from_distance_matrix import *
from distance_matrix import *

parser = argparse.ArgumentParser()

parser.add_argument("-c", "--coordinate_dir", required=True, type=str, help="relative path containing csv files with entry and exit coordinates (MUST INCLUDE / AT THE END)")
parser.add_argument("-t", "--threshold", required=True, type=float, help="threshold for ribosomes considered to be in polysomes, in angstroms")
parser.add_argument("-o", "--output", required=True, type=str, help="output path of the csv file containing the polysome indices for each ribosome")

args = parser.parse_args()

coordinate_dir = args.coordinate_dir
threshold = args.threshold
output = args.output

csv_files = os.listdir(coordinate_dir)

monosome_indices = {}
dimer_indices = {}
polychain_indices = {}

#print("csv_files: {}".format(csv_files))

for csv in csv_files:
	df = pd.read_csv(coordinate_dir + csv, header=None)
	print(df)
	coordinates = df.to_numpy()
	tomogram = csv.split(".")[0]
	

	entry_x = coordinates[:,0]
	entry_y = coordinates[:,1]
	entry_z = coordinates[:,2]
	exit_x = coordinates[:,3]
	exit_y = coordinates[:,4]
	exit_z = coordinates[:,5]


	x_distances = ipdm_1d(entry_x, exit_x)
	y_distances = ipdm_1d(entry_y, exit_y)
	z_distances = ipdm_1d(entry_z, exit_z)

	ipdm_matrix = ipdm_3d(x_distances, y_distances, z_distances)

	filtered_matrix = find_smaller_distances(ipdm_matrix)
	polysome_matrix = apply_threshold(filtered_matrix, threshold)
	polysomes = find_polysomes(polysome_matrix)

	monosome_indices[tomogram] = get_ribosomes_in_polysomes(polysomes)[0]
	dimer_indices[tomogram] = get_ribosomes_in_polysomes(polysomes)[1]
	polychain_indices[tomogram] = get_ribosomes_in_polysomes(polysomes)[2]
	
	print(tomogram)
	print("monosomes: {}; dimers: {}; polychains: {}".format(monosome_indices, dimer_indices, polychain_indices))
	print("===========================================================")

monosomes_df = pd.DataFrame.from_dict(monosome_indices, orient='index')
dimers_df = pd.DataFrame.from_dict(dimer_indices, orient='index')
polychains_df = pd.DataFrame.from_dict(polychain_indices, orient='index')

monosomes_df.to_csv(output+'_monosomes.csv')
dimers_df.to_csv(output+'_dimers.csv')
polychains_df.to_csv(output+'_polychains.csv')
