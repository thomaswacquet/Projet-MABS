import numpy as np 
from Alignment import Alignment

class Aligner:
    """
    Classe contenant les méthodes de calcul de l'alignment
    """
    def __init__(self, U, V, cost):
        self.__U = U
        self.__V = V
        self.__cost = cost

    def DPmatrix(self):
        """
        Calcule la matrice d'alignement pour U et V

        Returns
        -------
        res : np.array
            matrice d'alignement
        """
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
        """
        Retourne le score à partir de la matrice
        """
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
    
    def __Traceback(self, mat):
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
                # en cas de match/mismatch
                alignedU += self.__U[idxCol-1]
                alignedV += self.__V[idxRow-1]
                idxRow -= 1
                idxCol -= 1
            elif mat[idxRow-1][idxCol] + self.__cost(self.__V[idxRow-1], "") == mat[idxRow][idxCol]:
                # en cas d'insertion
                alignedU += "-"
                alignedV += self.__V[idxRow-1]
                idxRow -= 1
            else:
                # en cas de délétion
                alignedV += "-"
                alignedU += self.__U[idxCol-1]
                idxCol -= 1

        return Alignment(alignedU[::-1], alignedV[::-1], startCoords, endCoords, mat)

    def FindAlignment(self):
        """
        Retourne l'alignement entre U et V en calculant la matrice puis en calculant le traceback

        Returns
        -------
        alignment : Alignment
        """
        mat = self.DPmatrix()
        return self.__Traceback(mat)
    
    def DPmatrix_v2(self):
        """
        Calcule le score pour U et V en calculant les valeurs de la matrice en stockant 2 lignes successives de la matrice

        Returns
        -------
        score : int
            score d'alignement
        """
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
        """
        Calcule le score

        Returns
        -------
        score : int
        """
        return self.DPmatrix_v2()