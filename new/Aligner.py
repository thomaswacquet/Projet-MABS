import numpy as np 

class Alignment:
    def __init__(self, alignedU, alignedV, startCoords, endCoords, mat):
        self.startCoords = startCoords
        self.endCoords = endCoords
        self.mat = mat
        self.alignedU = alignedU
        self.alignedV = alignedV

    def Score(self):
        return self.mat[self.endCoords[0]][self.endCoords[1]]


class Aligner:
    def __init__(self, U, V, cost):
        self.__U = U
        self.__V = V
        self.__cost = cost

    def DPmatrix(self):
        res = np.zeros((len(self.__V)+1, len(self.__U)+1), dtype=np.int8)

        for idxV in range(len(self.__V)):
            for idxU in range(len(self.__U)):
                res[idxV+1, idxU+1] = max(
                    res[idxV, idxU] + self.__cost(self.__U[idxU], self.__V[idxV]),
                    res[idxV, idxU+1] + self.__cost(self.__U[idxU], ""),
                    res[idxV+1, idxU] + self.__cost("", self.__V[idxV]),
                    0
                )
        return res

    def SimilarityScore(self):
        mat = self.DPmatrix()
        row, col = self.SimilarityScoreCoords(mat)
        return mat[row][col]

    def __SimilarityScoreCoords(self, mat):
        max = 0
        scoreCoords = (0, 0)
        for row in range(len(mat)):
            for col in range(len(mat[0])):
                if mat[row][col] > max:
                    max = mat[row][col]
                    scoreCoords = (row, col)
        return scoreCoords
    
    def FindAlignment(self):
        mat = self.DPmatrix()
        endCoords = self.__SimilarityScoreCoords(mat)

        idxRow = endCoords[0]
        idxCol = endCoords[1]
        alignedU = ""
        alignedV = ""

        # stocke les coordonnées de départ (correspondant à la valeur 0 dans la matrice)
        startCoords = endCoords

        while mat[idxRow][idxCol] != 0:
            startCoords = (idxRow, idxCol)

            # trouver quelle case est utilisée pour calculer la case courante
            if mat[idxRow-1][idxCol-1] + self.__cost(self.__V[idxRow-1], self.__U[idxCol-1]) == mat[idxRow][idxCol]:
                # case gauche-haut
                alignedU += self.__U[idxCol-1]
                alignedV += self.__V[idxRow-1]
                idxRow -= 1
                idxCol -= 1
            elif mat[idxRow-1][idxCol] + self.__cost(self.__V[idxRow-1], "") == mat[idxRow][idxCol]:
                # case haut
                alignedU += "-"
                alignedV += self.__V[idxRow-1]
                idxRow -= 1
            else:
                # case gauche
                alignedV += "-"
                alignedU += self.__U[idxCol-1]
                idxCol -= 1

        # print(alignedU[::-1])
        # print(alignedV[::-1])
        alignment = Alignment(alignedU[::-1], alignedV[::-1], startCoords, endCoords, mat)
        return alignment
#        return ((alignedU[::-1], alignedV[::-1]), startCoords)
    
    def DPmatrix_v2(self):
        score = 0
        ligne1 = np.zeros(len(self.__U)+1, dtype=np.int8)
        ligne2 = np.zeros(len(self.__U)+1, dtype=np.int8)

        for idxV in range(len(self.__V)):
            for idxU in range(len(self.__U)):
                # calcule le max entre la case de gauche, du haut, et de gauche-haut sur la diagonale et 0
                ligne2[idxU+1] = max(
                    ligne1[idxU] + self.__cost(self.__U[idxU], self.__V[idxV]),
                    ligne1[idxU+1] + self.__cost(self.__U[idxU], ""),
                    ligne2[idxU] + self.__cost("", self.__V[idxV]),
                    0
                )
                if ligne2[idxU+1] > score:
                    score = ligne2[idxU+1]

            ligne1 = ligne2
            ligne2 = np.zeros(len(self.__U)+1, dtype=np.int8)
        return score
    
    def Score(self):
        return self.FindAlignment().Score()
    
