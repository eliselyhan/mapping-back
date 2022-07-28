import pandas as pd
import starfile
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-r", "--ribosome_filename", type=str, required=True, help="star file with the ribosome data")
parser.add_argument("-p", "--polysome_filename", type=str, required=True, help="csv file with the ribosome in polysome indices for each tomogram")
parser.add_argument("-o", "--output", type=str, required=True, help="the filename that you write the your output star file to (relative path)")

args = parser.parse_args()

ribosome_filename = args.ribosome_filename
polysome_filename = args.polysome_filename
output = args.output

df_all = starfile.read(ribosome_filename)
df_particles = df_all['particles']
df_all_poly = pd.DataFrame()

# TODO: get tomogram name from polysome csv file

polysomes = pd.read_csv(polysome_filename)
polysomes.set_index("Unnamed: 0", inplace=True)

#print(polysomes.loc["20220106_5991-2_L3_ts002"])
#print(polysomes.head())

tomograms = list(polysomes.index)
#print(tomograms)


for tomogram in tomograms:
	indices = list(polysomes.loc[tomogram])
	indices = [x for x in indices if pd.isnull(x) == False]


	tomogram = tomogram + ".mrc.tomostar"
	df = df_particles[df_particles['rlnMicrographName'] == tomogram]
	df_polysomes = df.iloc[indices, :]
	df_all_poly = df_all_poly.append(df_polysomes)

print(df_all_poly)
starfile.write(df_all_poly, output, overwrite=True)
