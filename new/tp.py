from calendar import c
import sys #à retirer ?
from IOAlignment import IOAlignment # modifier le parser de score pour qu'il prenne l'alphabet

from Aligner import Alignment
from typing import Type
from Cost import Cost
import numpy as np
import random
import time #à retirer ?
from Aligner import Aligner
import argparse

# 20 min d'oral => 10 min de présentations + 10 min de questions


# afficher matrice
# a = "    " + " ".join(list(U))
# print(a)
# for row in range(len(mat)):
# 	if row != 0:
# 		print(V[row-1] + " ", end="")
# 	else:
# 		print("  ", end="")
# 	for col in range(len(mat[0])):
# 		print(str(mat[row][col]) + " ", end="")
# 	print()
# print("  ")

def parse():
    """
    Fonction pour créer un 'parser' pour les options sur le terminal
    """
    parser = argparse.ArgumentParser(description="Parse arguments of #nom du fichier#")
    parser.add_argument('-f', dest='file', type = str,  required=True,
                        help="fichier à parser contenant les séquences")
    parser.add_argument('-s', dest='seuil', required=False, type = int,
                        help="seuil minimal pour le score d'alignement | par défaut s = 5", default=5) 
    parser.add_argument('-n', dest='nb_seq', required=False, type = int,
                        help="nombre de séquences aléatoires à construire | par défaut n = 15",  default=15) 
    parser.add_argument('-c_abc', dest='cost_alphabet', required=True, type =str,
                        help="alphabet des séquences") # à mettre dans le fichier des scores 
    parser.add_argument('-c_scr', dest='cost_scores', required=True, type = str,
                        help="fichier contenant les scores")
    parser.add_argument('-l', dest='length', required=False,type=int, 
                        help='nombre de bases à afficher par ligne | par défaut l = 10',  default=10) 
    return parser.parse_args()


def DisplayAlign(fi, s, n, cost, l):
	"""Affichage des séquences, score de l'alignement et alignement des 
	séquences données dans le fichier fasta en entrée
	
    Parameters:
    --------
    fi : File                         
        fichier fasta à parser      
    s : Int
        seuil minimal pour le score d'alignement
    n : Int
        nombre de séquences aléatoires à construire
	cost : Cost
        objet de la classe Cost contenant les scores et l'alphabet
    l : Int
        nombre de bases à afficher par ligne

		
	Returns
    -------
    score : Int
        score de l'alignement
	"""
	# ParseFasta parse les séquences en dictionnaire 
	# sous la forme {label : séquence}
	sequences = fi
	seqVals = list(sequences.values())
	U = seqVals[0]
	V = seqVals[1]
	seqLabels = list(sequences.keys())
	labelU = seqLabels[0]
	labelV = seqLabels[1]

	# Calcule de la matrice et du score
	aligner = Aligner(U, V, cost)
	alignment = aligner.FindAlignment()

	score  = alignment.Score()

	# Affichage du résultat
	print(">" + labelU)
	print(U)
	print(">" + labelV)
	print(V + "\n")
	print("Alignment score: " + str(score) + "\n")
	print("Alignment:\n")

	nbChars = max(len(str(alignment.startCoords[0])), len(str(alignment.startCoords[1])))

	for i in range(0, len(alignment.alignedU), l):
		# Détermine les positions de départ et de fin
		Ustart = i + alignment.startCoords[1]
		Vstart = i + alignment.startCoords[0]
		Uend = Ustart + min(i + l, alignment.endCoords[1])
		Vend = Vstart + min(i + l, alignment.endCoords[0])
		Upart = alignment.alignedU[i:i+l]
		Vpart = alignment.alignedV[i:i+l]

		print(f"{labelU}   {str(Ustart).rjust(nbChars)}  {Upart}  {Uend}")
		display = " " * (len(labelU) + 3 + nbChars + 2)
		for j in range(len(Upart)):
			if Upart[j] == Vpart[j]:
				display += "|"
			else:
				display += " "
		print(display)
		print(f"{labelV}   {str(Vstart).rjust(nbChars)}  {Vpart}  {Vend}")
		print()


	DisplayDist(U, V, s, n, cost, score)
	return score

def Significance(U, V, n, cost):
	"""Création de n séquences aléatoires et retourne les scores d'alignement
	
	Parameters:
    --------
    U : File                         
        séquence U
    V : Int
        séquence V
    n : Int
        nombre de séquences aléatoires à construire
	cost : Cost
        objet de la classe Cost contenant les scores et l'alphabet
		
	Returns
    -------
    scores : List
        scores des alignements
	"""
	a = U.count("A")
	g = U.count("G")
	t = U.count("T")
	c = U.count("C")
	scores = []

	for i in range(n):
		# crée la chaîne aléatoire
		rnd = list(a * "A" + g * "G" + c * "C" + t * "T")
		random.shuffle(rnd)
		seq = ''.join(rnd)
		aligner = Aligner(seq, V, cost)

		score = aligner.Score()
		scores.append(score)
	return scores

def DisplayDist(U, V, s, n, cost, alignmentScore):
	"""Affichage de la distribution des scores calculée par Significance
	
	Parameters:
    --------
    U : File                         
        séquence U
    V : Int
        séquence V
    n : Int
        nombre de séquences aléatoires à construire
	cost : Cost
        objet de la classe Cost contenant les scores et l'alphabet
	alignmentScore : Int
        score de l'alignement
	"""
	scores = Significance(U, V, n, cost)

	count = 0
	for score in scores:
		if score > s:
			count += 1
	print(f"Significance of the score: {count / n * 100}% ({n} sequences)")

	print(" score       #")
	lenScore = len("score")+1
	for i in range(max(scores)+1):
		# Espacer les scores afin qu'ils soient alignés sur "score"
		count = scores.count(i)
		equals = int(count/n*100) * "="

		offset = 0

		# Si le score courant est égal au score d'alignement, ajouter une *
		if i == alignmentScore:
			print("*", end="")
			offset = -1
		print(f"{str(i).rjust(lenScore + offset)}     {str(count).rjust(3)}: {equals}")


# ------------------------------ méthode de calcul du score sans matrice complète (Q9) --------------------


# def Significance_v2(U, V, s, n, cost):
# 	a = U.count("A")
# 	g = U.count("G")
# 	t = U.count("T")
# 	c = U.count("C")
# 	scores = []
# 	#count = 0
# 	for i in range(n):
# 		rnd = list(a * "A" + g * "G" + c * "C" + t * "T")
# 		random.shuffle(rnd)
# 		seq = ''.join(rnd)

# 		aligner = Aligner()
# 		score = DPmatrix_v2(seq, V, cost)
# #		score = SimilarityScore(mat)
# 		#if score >= s:
# 		scores.append(score)
# #			count += 1
# #			count += 1
# 	return scores

# def DisplayDist_v2(U, V, s, n, cost, scoreAlign):
# 	scores = Significance_v2(U, V, s, n, cost)

# 	count = 0
# 	for score in scores:
# 		if score > s:
# 			count += 1
# 	print(f"Significance of the score: {count / n * 100}% ({n} sequences)")

# 	print(" score       #")
# 	lenScore = len("score")
# 	for i in range(max(scores)+1):
# 		# ajoute des espaces à gauche de tel sorte que les scores soient alignés sur "score"
# 		count = scores.count(i)
# 		equals = int(count/n*100) * "="

# 		offset = 0

# 		# si le score courant est égal au score d'alignement, mettre un *
# 		if i == scoreAlign:
# 			print("*", end="")
# 			offset = -1

# 		print(f"{str(i).rjust(lenScore + offset)}     {str(count).rjust(3)}: {equals}")


if __name__ == "__main__": 
    OPTIONS = parse()
    file = IOAlignment(OPTIONS.file)
    file = file.ParseFasta()
    alphabet = OPTIONS.cost_alphabet
    scores = IOAlignment(OPTIONS.cost_scores)
    scores = scores.ParseConfig()
    #string_to_break = OPTIONS.cost_scores
    #scores = string_to_break.split()
    # ide = int(scores[0])
    # sbt = int(scores[1])
    # ind = int(scores[2])
    # scores =  {"identity": ide, "substitution": sbt, "indel": ind}
    #scores = {"identity": 2, "substitution": -1, "indel": -2}
    # cost = Cost(scores, "ATCG")
    #U = "TGTTACGG"
    #V = "GGTTGACTA"
    # l = 10

    score = DisplayAlign(file, OPTIONS.seuil,  OPTIONS.nb_seq, Cost(scores, alphabet), OPTIONS.length)
    #score = DisplayAlign("sequences_2.fa", 5, 30, cost, l)


# for n in [5000, 10000, 100000]:
# 	time1_start = time.time()
# 	DisplayDist(U, V, 7, n, cost, score)
# 	time1_end = time.time()

# 	duration1 = time1_end - time1_start

# 	time2_start = time.time()
# 	DisplayDist_v2(U, V, 7, n, cost, score)
# 	time2_end = time.time()

# 	duration2 = time2_end - time2_start
# 	string += f"{n} {duration1} {duration2}\n"
# print(string)

#print(DPmatrix_v2(U, V, cost))
#print("---------------")
#print(DPmatrix(U, V, cost))
#print(max)
#print(r)