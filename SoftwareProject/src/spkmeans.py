import pandas as pd 
import numpy as np
import sys
from cmath import inf
import mykmeanssp as mk


def main():
    try:

        data = processData(sys.argv)
        if (data == None):
            print("Invalid Input")
            return 0
        k = data[0]
        goal = data[1]
        N = data[2]
        dim = data[3]
        mat = data[4]
        mat_copy = mat.tolist()
        MAX_ITER = 300

        if goal == "spk":
            T = mk.spkFit(mat_copy, N, dim, k)
            K = len(T[0])
            
            if T == None:
                raise Exception

            list_of_centroids, list_of_indices = kmeansPP(T, N, K, K)
            if list_of_centroids == None or list_of_indices == None:
                raise Exception

            ord_mat = fixMat(T, list_of_indices)
            final_centroids = mk.kmeansFit(ord_mat, N, K, K, MAX_ITER)
            if final_centroids == None:
                raise Exception

            # Printing the indices chosen by the Kmeans++ algorithm 
            print(','.join(str(item) for item in list_of_indices))

            # Printing the calculated final centroids from the Kmeans algorithm
            for i in range(K):
                print(','.join(str("%.4f"%item) for item in final_centroids[i]))


        elif goal == "wam":
            weightedMatrix = mk.weightedFit(mat_copy, N, dim)
            if weightedMatrix == None:
                raise Exception
            #Printing the weighted matrix
            for i in range(N):
                print(','.join(str("%.4f"%item) for item in weightedMatrix[i]))

        elif goal == "ddg":
            diagonalMatrix = mk.diagonalFit(mat_copy, N, dim)

            if diagonalMatrix == None:
                raise Exception
            
            #Printing the diagonal matrix
            for i in range(N):
                print(','.join(str("%.4f"%item) for item in diagonalMatrix[i]))


        elif goal == "lnorm":
            lnormMatrix = mk.lnormFit(mat_copy, N, dim)

            if lnormMatrix == None:
               raise Exception
            
            #Printing the lnorm matrix
            for i in range(N):
                print(','.join(str("%.4f"%item) for item in lnormMatrix[i]))

        elif goal == "jacobi":
            eiganVectorsMatrix, eiganValues = mk.jacobiFit(mat_copy, N)   
            if eiganVectorsMatrix == None or eiganValues == None:
                raise Exception
            #Printing the eiganvalues
            print(','.join(str("%.4f"%item) for item in eiganValues))

            #Printing the corresponding eigenvectors (as columns)
            for i in range(N):
                print(','.join(str("%.4f"%item) for item in eiganVectorsMatrix[i]))

    except:                
        print("An Error Has Occurred")
        return

"""
    processData validate and orginaize user input
    Parameters
    ----------
    arr: list
        User arguments from input
    Returns
    -------
    "": list
        list that contains:
        k: int
            Number of required clusters
        goal: enum
            The excecution goal from INPUT_GOAL
        N: int
            Number of vectors
        dim: int
            Dimention of vectors
        vectors: NumPy array
            NumPy array of all vectors from user file input path
            *Except Jacobi, in this case the array represent symmetric matrix
"""
def processData(arr):

    INPUT_GOAL = ["spk", "wam", "ddg", "lnorm", "jacobi"]
    if (len(arr) != 4):
        return None
    
    try:
        k = int(arr[1])
        goal = arr[2]
        file_name = arr[3]
        if goal not in INPUT_GOAL:
            raise Exception

        file = open(file_name, 'rb')
        vectors = pd.read_csv(file, header=None)
        vectors = vectors.to_numpy()

        N = vectors.shape[0]
        dim = vectors.shape[1]

        if (k < 0 or k > N):
            raise Exception
        
    except:
        return None

    return [k, goal , N, dim, vectors] 

"""
    kmeanPP excecute Kmeans++ algorithm
    Parameters
    ----------
    T: list 
        List of vectors to excecuted kmeans algorithm
    N: int
        Number of vectors
    k: int
        Number of clusters
    dim: int
        Dimention of vectors in T
    Returns
    -------
    "": tuple
         list_centroids: list
            List of vectros to be selected as intial centroids 
        list_of_indices: list
            List of the indices of the selected intial centroids 
"""
def kmeansPP(T, N, K, dim):
    np.random.seed(0)
    i = 0
    random_ind = np.random.choice(N)
    list_of_indices = [random_ind]
    list_of_centroids = [T[random_ind]]
    dist_array = np.zeros(N) 
    prob_Array = np.zeros(N)

    while (len(list_of_indices) < K):
        for r in range(N):
            low = float('inf')
            for j in range(i+1):
                dist = np.square(np.linalg.norm(np.asarray(T[r]) 
                - np.asarray(list_of_centroids[j])))

                if dist < low:
                    low = dist
            dist_array[r] = low

        prob_array = dist_array / np.sum(dist_array) 
        i += 1
        random_ind = np.random.choice(N, p = prob_array)

        list_of_indices.append(random_ind)
        list_of_centroids.append(T[random_ind])
    
    return (list_of_centroids, list_of_indices)

"""
    fixmat reorganize matrix such that the centroids will be first
    Parameters
    ----------
    mat: list
        List of vectors
    list_of_indices: list
        List of the indices of the selected intial centroids 
    Returns
    -------
    result: list
        New list that contains organized vectors
"""
def fixMat(mat, list_of_indices):
    result = []

    for ind in list_of_indices:
        result.append(mat[ind])

    for i in range(len(mat)):
        if i not in list_of_indices:
            result.append(mat[i])

    return result

main()