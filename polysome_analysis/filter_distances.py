import numpy as np
import numpy.ma as ma

def find_smaller_distances(ipdm: np.ndarray) -> np.ndarray:
    ipdm_T = ipdm.T
    upper_triangle_mask = (ipdm > ipdm_T)
    lower_triangle_mask = np.logical_not(upper_triangle_mask.T)
    
    upper_triangle_mask = upper_triangle_mask.astype(int)
    lower_triangle_mask = lower_triangle_mask.astype(int)

    greater_distances_mask = upper_triangle_mask + lower_triangle_mask
    filtered = ma.array(ipdm, mask=greater_distances_mask).filled(10**6)
    np.fill_diagonal(filtered, 10**6)
    return filtered

def apply_threshold(greater_ipdm: np.ndarray, threshold: float=70) -> np.ndarray:
    poly_mask = greater_ipdm > threshold
    filtered_ipdm = ma.array(greater_ipdm, mask=poly_mask).filled(0)

    return filtered_ipdm 

