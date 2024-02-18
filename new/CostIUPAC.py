import re

class CostIUPAC:
	def __init__(self, costs, alphabetFile):
		self.__indel = costs["indel"]
		self.__subst = costs["substitution"]
		self.__identity = costs["identity"]
		self.__alphabet, self.__combs = self.__readCombinations(alphabetFile)

	def __call__(self, c1, c2):
		return self.__compare(c1, c2)

	def __compare(self, c1, c2):
		# indels
		if c1 == "" or c2 == "":
			return self.__indel

		# si une des lettres n'est pas dans l'alphabet renvoyer une erreur
		if c1 in self.__alphabet:
			if c2 in self.__alphabet:
				# si les 2 lettres de l'alphabet sont égaux
				if c1 == c2:
					return self.__identity
				else: # 2 lettres de l'alphabets non égales
					return self.__subst

			# si c2 n'est ni dans l'alphabet ni une lettre IUPAC
			if not c2 in self.__combs:
				raise ValueError

			# c2 est une lettre IUPAC
			if c1 in self.__combs[c2]:
				return 0
			# intersection vide
			return self.__subst

		# c1 est une lettre IUPAC mais c2 est dans l'alphabet
		if c2 in self.__alphabet:
			# c2 est une lettre IUPAC
			if c2 in self.__combs[c1]:
				return 0
			# intersection vide
			return self.__subst

		# c1 comme c2 sont des lettres IUPAC

		for combC1 in self.__combs[c1]:
			if combC1 in self.__combs[c2]:
				return 0
		return self.__subst

	def __readCombinations(self, file):
		f = open(file, "r")
		combs = {}

		# la première ligne correspond à l'alphabet
		alphabet = f.readline().rstrip()

		for l in f:
			l = l.rstrip()

			# essaie de matcher <lettre>\t<lettres séparés par espace>
			# et stocke dans un dictionnaire les combinaisons pour chaque lettre 
			m = re.match("^([A-Z])\t([A-Z](?: [A-Z])*)+", l)
			if m:
				letters = m.group(2).split(" ")
				for letter in letters:
					if not letter in alphabet:
						raise TypeError
				combs[m.group(1)] = letters
			else:
				raise TypeError

		return (alphabet, combs)