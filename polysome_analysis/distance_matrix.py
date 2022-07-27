from scipy.spatial import distance_matrix
import numpy as np

def get_distances(x: np.ndarray, y: np.ndarray):
	mat = distance_matrix(x,y);
	np.fill_diagonal(mat, 10**6);
	return mat

def filter_duplicates(d: np.ndarray, threshold: float):
    n = d.shape[0]
    filtered = d < threshold
    indicies = np.where(filtered)

    return indicies
