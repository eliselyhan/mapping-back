import numpy as np
from disjoint_set import DisjointSet

def find_polysomes(filterd_ipdm: np.ndarray) -> list:
    n = filterd_ipdm.shape[0]
    
    row_indices = np.nonzero(filterd_ipdm)[0]
    column_indices = np.nonzero(filterd_ipdm)[1]

    nonzero_indices = np.hstack((row_indices[:,None], column_indices[:,None]))
    m = nonzero_indices.shape[0]

    # use disjoint sets to track polysomes
    ds = DisjointSet()
    for i in range(n):
        ds.find(i)

    for i in range(m):
        x = nonzero_indices[i,0]
        y = nonzero_indices[i,1]
        ds.union(x,y)
            
    return list(ds.itersets())

def get_ribosomes_in_polysomes(polysomes):
    ribosomes_in_polysomes = []
    for set in polysomes:
        if len(set) > 1:
            for index in set:
                ribosomes_in_polysomes.append(index)
    
    return ribosomes_in_polysomes

if __name__ == "__main__":
    
    x = np.array([[0,1,0,0,0],
                [0,0,0,0,0],
                [0,0,0,0,0],
                [0,0,0,0,1],
                [0,0,0,0,0]])

    polysomes = find_polysomes(x)
    print(polysomes)
    print(get_ribosomes_in_polysomes(polysomes))
