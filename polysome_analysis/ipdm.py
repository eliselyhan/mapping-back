import numpy as np

def ipdm_1d(entry: np.ndarray, exit: np.ndarray) -> np.ndarray:
    a = entry
    b = exit

    upper_triangle = np.subtract(b[:,None], a)

    ipdm = np.abs(upper_triangle)
    
    return ipdm

def ipdm_3d(x_distances: np.ndarray, y_distances: np.ndarray, z_distances: np.ndarray) -> np.ndarray:
    x_square = np.square(x_distances)
    y_square = np.square(y_distances)
    z_square = np.square(z_distances)

    ipdm_3d = np.sqrt(x_square + y_square + z_square)

    return ipdm_3d

if __name__ == "__main__":
    a = np.random.randint(0,10,5)
    b = np.random.randint(0,10,5)
    print(a)
    print(b)
    print(ipdm_1d(a,b))
