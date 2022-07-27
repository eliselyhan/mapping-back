# Polysome Analysis
## ipdm.py
* Generates interpoint distance matrix betweeen the entries and exits for all ribosomes
## filter_distances.py
* Uses the interpoint distance matrix from ipdm.py
* Keeps the smaller value between the index pairs (i,j) and (j,i)
* Remove all entries above a certain distance threshold
## polysomes_from_distance_matrix.py
* Uses the filtered distance matrix from filter_distances.py
* Produces a list of sets; each set corresponds to a polysome
