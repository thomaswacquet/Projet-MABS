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
			currentSeq = l[1]
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
	max = 0
	scoreCoords = (0, 0)
	for row in range(len(mat)):
		for col in range(len(mat[0])):
			if mat[row][col] > max:
				max = mat[row][col]
				scoreCoords = (row, col)
	return scoreCoords

def Significance(V,U,s,n,cost):
	a = U.count("A")
	g = U.count("G")
	t = U.count("T")
	c = U.count("C")
	for i in range(n):
		s = a * "A" + g * "G" + c * "C" + t * "T"
		rnd = list(s)
		random.shuffle(rnd)
		seq = ''.join(rnd)

def FindAlignment(mat, scoreCoords, U, V, cost):
	idxRow = scoreCoords[0]
	idxCol = scoreCoords[1]
	alignedU = ""
	alignedV = ""

	
	alignedU = U[idxCol-1]
	alignedV = V[idxRow-1]

	idxCol -= 1
	idxRow -= 1

	while idxRow > 0 and idxCol > 0:
		idxU = idxCol + 1
		idxV = idxRow + 1

		#if (mat[idxRow-1][idxCol]):


		return
		maxVal = max(
			mat[idxRow-1][idxCol],
			mat[idxRow][idxCol-1],
			mat[idxRow-1][idxCol-1]
		)
		#print(mat[idxRow][idxCol])
		print(maxVal)

		if (maxVal == mat[idxRow-1][idxCol]):
			alignedU += "-"
			alignedV += V[idxRow-1]
		elif (maxVal == mat[idxRow][idxCol-1]):
			alignedU += U[idxCol-1]
			alignedV += "-"
		else:
			alignedU += U[idxCol-1]
			alignedV += V[idxRow-1]

		idxCol -= 1
		idxRow -= 1

	print(alignedU)
	print(alignedV)

	print("    " + " ".join(list(U)))
	print(mat)


scores = {"identity": 2, "substitution": -1, "indel": -2}
cost = Cost(scores, "ATCG")
U = "TGTTACGG"
V = "GGTTGACTA"
mat = DPmatrix(U, V, cost)
scoreCoords = SimilarityScore(mat)
FindAlignment(mat, scoreCoords, U, V, cost)
#print(max)

#Significance("A", "AGC", "A", 5, "A")