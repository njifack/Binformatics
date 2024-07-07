import numpy as np
import random
from Bio import Phylo
import io
distanceMatrix = np.array([[0, 12,12,13,15,15],
                           [12, 0, 2, 6, 8, 8],
                           [12, 2, 0, 6, 9, 9],
                           [13, 6, 6, 0, 8, 8],
                           [15, 8, 9, 8, 0, 4],
                           [15, 8, 9, 8, 4, 0]])
speciesList = ["M_Spacii", "T_Pain", "G_Unit", "Q_Doba", "R_Mani", "A_Finch"]


def UPGMA(dM, sp):
    while len(dM) > 1:
        leastRow, leastCol = findSmallest(dM)
        dM = updateMatrix(dM, leastRow, leastCol)
        sp = updateSpecies(sp, leastRow, leastCol)
        print(sp[-1])
    # After clustering is completed, generate and display the tree
    tree = generateTree(sp[-1])

    Phylo.draw(tree)


def findSmallest(dM):
    smallest_val = np.min(dM[np.nonzero(dM)])
    if smallest_val == 0:
        return None, None
    indices = np.where(dM == smallest_val)
    rand_index = random.choice(range(len(indices[0])))
    row = indices[0][rand_index]
    col = indices[1][rand_index]
    return row, col


def updateSpecies(sp, r, c):
    distance = distanceMatrix[r][c] / 2
    new_species = f"({sp[r]}, {sp[c]}:{distance})"
    sp[r] = new_species
    del sp[c]
    return sp


def updateMatrix(dM, row, col):
    n = len(dM)
    new_row = np.zeros(n - 1)
    new_matrix = np.zeros((n - 1, n - 1))

    new_row_idx = 0
    for i in range(n):
        if i != row and i != col:
            new_col_idx = 0
            for j in range(n):
                if j != row and j != col:
                    new_matrix[new_row_idx, new_col_idx] = dM[i, j]
                    new_col_idx += 1
            new_row[new_row_idx] = (dM[row, i] + dM[col, i]) / 2
            new_row_idx += 1

    new_matrix[-1, :] = new_row
    new_matrix[:, -1] = new_row

    return new_matrix


def generateTree(tree_str):
    return Phylo.read(io.StringIO(tree_str), "newick")


UPGMA(np.array(distanceMatrix), speciesList)
