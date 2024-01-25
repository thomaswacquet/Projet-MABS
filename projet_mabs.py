# oral semaine du 11 mars ?

# costs values :
    # identity : +2
    # substitution : -1
    # indel : -2

import sys
import re
if 'Cost' not in sys.modules:
    from Cost import Cost

def ParseFasta(Fi):
    """Lecture d'un fichier au format fasta

        Returns
        -------
        read_file : List
            Liste avec les séquences et leur identifiant
    """
    read_file = ''
    with open(Fi) as fileIn: 
        for ligneLue in fileIn: 
            if ligneLue.find('>') != -1 and ligneLue[0] in ['A', 'C', 'G', 'T', '>'] :
                if read_file != '' :
                    read_file += '\n'
                read_file += ligneLue
            elif not ligneLue[0] in ['A', 'C', 'G', 'T', '>'] :
                raise TypeError
            else :
                ligneLue = ligneLue.rstrip()
                read_file += ligneLue 
    read_file = read_file.split('\n')
    return read_file

def ParseConfig(Fi):
    """Lecture d'un fichier au format fasta

        Returns
        -------
        read_file : List
            Liste avec les séquences et leur identifiant
    """
    liste_keys = []
    liste_values = []
    i = 0
    with open(Fi) as fileIn: 
        for ligneLue in fileIn: 
            score = re.findall(r'[-+]\d+', ligneLue)
            keys = re.findall(r'\b[a-zA-Z]+\b', ligneLue)
            if len(score) >= 1 :
                liste_values.append(score[0])
            if len(keys) >= 1 :
                liste_keys.append(keys[0])
    liste_keys.pop(0)
    dictionary = dict(zip(liste_keys, liste_values))
    return dictionary

def DPmatrix(U, V, cost):
    """calcule la matrice de programmation dynmaique de
    l'algorithme de Smith et Waterman pour les séquences U et V"""
    i = 0
    j = 0
    liste = []
    matrice = []
    while i < len(V):
        while j < len(U) :
            cout = Cost(cost, ['A', 'T', 'C', 'G'])
            a=
            c = cout(U[j], V[i])
            liste.append(c)
            j += 1
        matrice.append(liste)
        liste = []
        j = 0
        i += 1
    return matrice


if __name__ == '__main__' :
    seqs = ParseFasta('sequences_2.fa')
    score = ParseConfig('score.txt')
    print(DPmatrix(seqs[1], seqs[3], score))




