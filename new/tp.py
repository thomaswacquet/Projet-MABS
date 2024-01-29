from calendar import c
import re
import sys
from typing import Type
from Cost import Cost
import numpy as np
import random

def ParseFasta(Fi):
	f = open(Fi, "r")
	currentSeq = ""
	res = dict()
	for l in f:
		l = l.rstrip()
		if l[0] == ">":
			currentSeq = l[1:]
			res[currentSeq] = ""
		else:
			if currentSeq == "":
				raise TypeError
			if not re.match("[ACTG]*", "ATGC"):
				raise TypeError
			res[currentSeq] += l

	return res

def ParseConfig(Fi):
	f = open(Fi, "r")
	res = dict()

	for l in f:
		l = l.rstrip()
		if l != "" and not re.match("^ *#.*", l):
			m = re.match("([a-z]+) *: *([+-]\\d+)", l)
			if not m:
				raise TypeError
			if m.group(1) == "identity" or m.group(1) == "substitution" or m.group(1) == "indel":
				res[m.group(1)] = int(m.group(2))
			else:
				raise TypeError
	return res

def DPmatrix(U, V, cost):
	res = np.zeros((len(V)+1, len(U)+1), dtype=np.int8)

	for idxV in range(len(V)):
		for idxU in range(len(U)):
			res[idxV+1, idxU+1] = max(
				res[idxV, idxU] + cost(U[idxU], V[idxV]),
				res[idxV, idxU+1] + cost(U[idxU], ""),
				res[idxV+1, idxU] + cost("", V[idxV]),
				0
			)
	return res



def SimilarityScore(mat):
	row, col = SimilarityScoreCoords(mat)
	return mat[row][col]

def SimilarityScoreCoords(mat):
	max = 0
	scoreCoords = (0, 0)
	for row in range(len(mat)):
		for col in range(len(mat[0])):
			if mat[row][col] > max:
				max = mat[row][col]
				scoreCoords = (row, col)
	return scoreCoords


def FindAlignment(mat, scoreCoords, U, V, cost):
	idxRow = scoreCoords[0]
	idxCol = scoreCoords[1]
	alignedU = ""
	alignedV = ""

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

	last = scoreCoords

	while mat[idxRow][idxCol] != 0:
		last = (idxRow, idxCol)

		# trouver quelle case est utilisée pour calculer la case courante
		if mat[idxRow-1][idxCol-1] + cost(V[idxRow-1], U[idxCol-1]) == mat[idxRow][idxCol]:
			# case gauche-haut
			alignedU += U[idxCol-1]
			alignedV += V[idxRow-1]
			idxRow -= 1
			idxCol -= 1
		elif mat[idxRow-1][idxCol] + cost(V[idxRow-1], "") == mat[idxRow][idxCol]:
			# case haut
			alignedU += "-"
			alignedV += V[idxRow-1]
			idxRow -= 1
		else:
			# case gauche
			alignedV += "-"
			alignedU += U[idxCol-1]
			idxCol -= 1

	# retourne tuple car plus optimisé que liste
	return ((alignedU[::-1], alignedV[::-1]), last)

def Align(fi, cost):
	sequences = ParseFasta(fi)
	seqVals = list(sequences.values())
	U = seqVals[0]
	V = seqVals[1]
	seqLabels = list(sequences.keys())
	labelU = seqLabels[0]
	labelV = seqLabels[1]


	mat = DPmatrix(U, V, cost)
	endCoords = SimilarityScoreCoords(mat)
	score  = mat[endCoords[0]][endCoords[1]]


	print(">" + labelU)
	print(U)
	print(">" + labelV)
	print(V + "\n")
	print("Alignment score: " + str(score) + "\n")
	print("Alignment:\n")

	(alignedU, alignedV), startCoords = FindAlignment(mat, endCoords, U, V, cost)

	print(f"{labelU}\t{startCoords[1]}  {alignedU}  {endCoords[1]}")
	display = " " * len(labelU) + "\t   " 
	for i in range(len(alignedU)):
		if alignedU[i] == alignedV[i]:
			display += "|"
		else:
			display += " "
	print(display)
	print(f"{labelV}\t{startCoords[0]}  {alignedV}  {endCoords[0]}")



def Significance(V, U, s, n, cost):
	a = U.count("A")
	g = U.count("G")
	t = U.count("T")
	c = U.count("C")
	scores = []
	count = 0
	for i in range(n):
		rnd = list(a * "A" + g * "G" + c * "C" + t * "T")
		random.shuffle(rnd)
		seq = ''.join(rnd)

		mat = DPmatrix(seq, V, cost)
		score = SimilarityScore(mat)
		scores.append(score)
#		if score >= s:
#			count += 1
	return scores

def DisplayDist(V, U, s, n, cost):
	scores = Significance(V, U, s, n, cost)
	print(scores)
	for i in range(len(scores)):
		print(str(i) + ": " + "=" * scores.count(i))



scores = {"identity": 2, "substitution": -1, "indel": -2}
cost = Cost(scores, "ATCG")
#U = "AGTTTTCAG"
#V = "CTCATT"
U = "TGTTACGG"
V = "GGTTGACTA"
#mat = DPmatrix(U, V, cost)
#scoreCoords = SimilarityScore(mat)
#FindAlignment(mat, scoreCoords, U, V, cost)
#Align("sequences.fa", cost)
#print(max)

DisplayDist("A", "AGC", 2, 5, cost)
#print(r)