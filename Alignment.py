class Alignment:
    """Classe qui contient un alignement (valeurs de U et de V alignés avec leurs coordonnés ainsi que la matrice d'alignement)"""
    def __init__(self, alignedU, alignedV, startCoords, endCoords, mat):
        self.startCoords = startCoords
        self.endCoords = endCoords
        self.mat = mat
        self.alignedU = alignedU
        self.alignedV = alignedV

    def Score(self):
        return self.mat[self.endCoords[0]][self.endCoords[1]]