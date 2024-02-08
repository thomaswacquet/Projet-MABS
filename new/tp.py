from calendar import c
import sys
from IOAlignment import IOAlignment

from Aligner import Alignment
from typing import Type
from Cost import Cost
import numpy as np
import random
import time
from Aligner import Aligner

# 15 min de présentation


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



def DisplayAlign(fi, s, n, cost, l):
	# parse le fichier fasta et trouve les séquences U et V et leurs label
	sequences = IOAlignment.ParseFasta(fi)
	seqVals = list(sequences.values())
	U = seqVals[0]
	V = seqVals[1]
	seqLabels = list(sequences.keys())
	labelU = seqLabels[0]
	labelV = seqLabels[1]

	# calcule la matrice et le score
	aligner = Aligner(U, V, cost)
	alignment = aligner.FindAlignment()

	score  = alignment.Score()

	# affiche le résultat
	print(">" + labelU)
	print(U)
	print(">" + labelV)
	print(V + "\n")
	print("Alignment score: " + str(score) + "\n")
	print("Alignment:\n")

	nbChars = max(len(str(alignment.startCoords[0])), len(str(alignment.startCoords[1])))

	for i in range(0, len(alignment.alignedU), l):
		# détermine les positions de départ et de fin
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
	scores = Significance(U, V, n, cost)

	count = 0
	for score in scores:
		if score > s:
			count += 1
	print(f"Significance of the score: {count / n * 100}% ({n} sequences)")

	print(" score       #")
	lenScore = len("score")+1
	for i in range(max(scores)+1):
		# ajoute des espaces à gauche de tel sorte que les scores soient alignés sur "score"
		count = scores.count(i)
		equals = int(count/n*100) * "="

		offset = 0

		# si le score courant est égal au score d'alignement, mettre un *
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



scores = {"identity": 2, "substitution": -1, "indel": -2}
cost = Cost(scores, "ATCG")
U = "TGTTACGG"
V = "GGTTGACTA"
l = 10

score = DisplayAlign("sequences_2.fa", 5, 30, cost, l)

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