from IOAlignment import IOAlignment
from CostIUPAC import CostIUPAC
import random
from Aligner import Aligner
import argparse
import time

def DisplayMatrice(mat, U, V):
	"""
	Affiche la matrice mat par rapport aux séquences U et V
	"""
	# afficher matrice
	a = "    " + " ".join(list(U))
	print(a)
	for row in range(len(mat)):
		if row != 0:
			print(V[row-1] + " ", end="")
		else:
			print("  ", end="")
		for col in range(len(mat[0])):
			print(str(mat[row][col]) + " ", end="")
		print()
	print("  ")

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
	parser.add_argument('-c', dest='config_file', required=False, type = str,
						help="fichier contenant la configuration")
	parser.add_argument('-l', dest='length', required=False,type=int, 
						help='nombre de bases à afficher par ligne | par défaut l = 10',  default=10) 
	# display_matrix contient True si -M est spécifié
	parser.add_argument('-M', dest='display_matrix', required=False, 
						help='afficher la matrice', action='store_true')

	return parser.parse_args()


def DisplayAlign(fi, s, n, cost, l, displayMatrix):
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
	sequences = IOAlignment(fi).ParseFasta()
	seqVals = list(sequences.values())
	U = seqVals[0]
	V = seqVals[1]
	seqLabels = list(sequences.keys())
	labelU = seqLabels[0]
	labelV = seqLabels[1]

	# Calcule de la matrice et du score
	aligner = Aligner(U, V, cost)
	alignment = aligner.FindAlignment()
	if displayMatrix:
		DisplayMatrice(alignment.mat, U, V)

	score  = alignment.Score()

	# Affichage du résultat
	print(">" + labelU)
	print(U)
	print(">" + labelV)
	print(V + "\n")
	print("Alignment score: " + str(score) + "\n")
	print("Alignment:\n")


	for i in range(0, len(alignment.alignedU), l):
		# Détermine les positions de départ et de fin
		Ustart = i + alignment.startCoords[1]
		Vstart = i + alignment.startCoords[0]
		Uend = Ustart + min(i + l, alignment.endCoords[1])
		Vend = Vstart + min(i + l, alignment.endCoords[0])
		Upart = alignment.alignedU[i:i+l]
		Vpart = alignment.alignedV[i:i+l]
		maxLabelLength = max(len(labelU), len(labelV))

		nbChars = max(len(str(Ustart)), len(str(Vstart)))

		strU = f"{labelU.ljust(maxLabelLength)}   {str(Ustart).rjust(nbChars)}  "
		strV = f"{labelV.ljust(maxLabelLength)}   {str(Vstart).rjust(nbChars)}  "

		print(f"{strU}{Upart}  {Uend}")

		display = " " * (len(strU))
		for j in range(len(Upart)):
			if Upart[j] == Vpart[j]:
				display += "|"
			else:
				display += " "
		print(display)

		print(f"{strV}{Vpart}  {Vend}")
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

if __name__ == "__main__": 
	OPTIONS = parse()

	if OPTIONS.config_file != None:
		config = IOAlignment(OPTIONS.config_file).ParseConfig()
	else:
		# Création de valeur par défaut pour la configuration
		# si l'utilisateur n'en choisit pas
		config = dict()
		config["combs"] = {}
		config["substitution"] = -1
		config["indel"] = -2
		config["identity"] = 2
		config["alphabet"] = "ACGT"

	# time1_start = time.time()
	score = DisplayAlign(OPTIONS.file, OPTIONS.seuil,  OPTIONS.nb_seq, CostIUPAC(config, config["alphabet"], config["combs"]), OPTIONS.length, OPTIONS.display_matrix)
	# time1_end = time.time()

	# print(f"Temps d'exécution : {time1_end - time1_start}")