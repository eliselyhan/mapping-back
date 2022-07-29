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

df_all_mono = pd.DataFrame()
df_all_dimer = pd.DataFrame()
df_all_poly = pd.DataFrame()


# TODO: get tomogram name from polysome csv file

monosomes = pd.read_csv(polysome_filename+'_monosomes.csv')
dimers = pd.read_csv(polysome_filename+'_dimers.csv')
polychains = pd.read_csv(polysome_filename+'_polychains.csv')

monosomes.set_index("Unnamed: 0", inplace=True)
dimers.set_index("Unnamed: 0", inplace=True)
polychains.set_index("Unnamed: 0", inplace=True)

#print(polysomes.loc["20220106_5991-2_L3_ts002"])
#print(polysomes.head())

df_list = [monosomes, dimers, polychains];
df_names = ['monosomes', 'dimers', 'polychains']
df_all_list = [df_all_mono, df_all_dimer, df_all_poly]

for i in range(3):
	df_class = df_list[i]

	tomograms = list(df_class.index)
	for tomogram in tomograms:
		indices = list(df_class.loc[tomogram])

		#mono_indices = list(monosomes.loc[tomogram])
		#dimer_indices = list(dimers.loc[tomogram])
		#poly_indices = list(polychains.loc[tomogram])

		indices = [x for x in indices if pd.isnull(x) == False]
		#mono_indices = [x for x in mono_indices if pd.isnull(x) == False]
		#dimer_indices = [x for x in dimer_indices if pd.isnull(x) == False]
		#poly_indices = [x for x in poly_indices if pd.isnull(x) == False]

		tomogram = tomogram + ".mrc.tomostar"

		df = df_particles[df_particles['rlnMicrographName'] == tomogram]
		
		df_filtered = df.iloc[indices, :]

		#df_monosomes = df.iloc[mono_indices, :]
		#df_dimers = df.iloc[dimer_indices, :]
		#df_polychains = df.iloc[poly_indices, :]

		df_all_list[i]  = df_all_list[i].append(df_filtered)
		print(df_all_list[i])
		#df_all_mono = df_all_mono.append(df_monosomes)
		#df_all_dimer= df_all_mono.append(df_dimers)
		#df_all_poly = df_all_mono.append(df_polychains)

df_all_mono, df_all_dimer, df_all_poly = df_all_list

print('mono:', df_all_mono)
print('dimer:',df_all_dimer)
print('poly:', df_all_poly)

starfile.write(df_all_mono, output+'_monosomes.star', overwrite=True)
starfile.write(df_all_dimer, output+'_dimers.star', overwrite=True)
starfile.write(df_all_poly, output+'_polychains.star', overwrite=True)
