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
parser.add_argument("-t", "--threshold", required=True, type=str, help="threshold for ribosomes considered to be in polysomes, in angstroms")

args = parser.parse_args()

coordinate_dir = args.coordinate_dir
threshold = args.threshold

csv_files = os.listdir(coordinate_dir)

for csv in csv_files:
	df = pd.read_csv(coordinate_dir + csv)
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
	indices = get_ribosomes_in_polysomes(polysomes)

	print(tomogram)
	print("Ribosomes in polysomes: {}".format(indices))
	print("===========================================================")
