class Alignment:
    """Classe contenant un alignement (valeurs de U et de V alignées avec leurs coordonnées 
        ainsi que la matrice d'alignement)"""
    def __init__(self, alignedU, alignedV, startCoords, endCoords, mat):
        self.startCoords = startCoords
        self.endCoords = endCoords
        self.mat = mat
        self.alignedU = alignedU
        self.alignedV = alignedV

    def Score(self):
        """Retourne le score de l'alignement entre U et V"""
        return self.mat[self.endCoords[0]][self.endCoords[1]]