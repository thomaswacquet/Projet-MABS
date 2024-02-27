import re

class CostIUPAC:
	"""
	Classe qui permet de comparer 2 caractères à partir des scores et alphabets passés dans le constructeur
	"""

	def __init__(self, costs, alphabet, combs):
		self.__indel = costs["indel"]
		self.__subst = costs["substitution"]
		self.__identity = costs["identity"]
		self.__alphabet = alphabet
		self.__combs = combs

	def __call__(self, c1, c2):
		return self.__compare(c1, c2)

	def __compare(self, c1, c2):
		# en cas d'indels
		if c1 == "" or c2 == "":
			return self.__indel

		# Si c1 et c2 sont tous deux dans l'alphabet
		if c1 in self.__alphabet:
			if c2 in self.__alphabet:
				# en cas de match
				if c1 == c2:
					return self.__identity
				else: # en cas de mismatch
					return self.__subst

			# Si c2 n'est ni dans l'alphabet ni une lettre IUPAC
			if not c2 in self.__combs:
				raise ValueError(f"{c2} n'est pas dans l'alphabet")

			# c2 est une lettre IUPAC mais pas c1
			if c1 in self.__combs[c2]:
				return 0
			# en cas d'intersection vide
			return self.__subst

		# c1 est une lettre IUPAC mais c2 n'est pas dans l'alphabet
		if c2 in self.__alphabet:
			# c2 est une lettre IUPAC

			if not c1 in self.__combs:
				raise ValueError(f"{c1} n'est pas dans l'alphabet")

			if c2 in self.__combs[c1]:
				return 0
			# en cas d'intersection vide
			return self.__subst

		# Si c1 comme c2 sont des lettres IUPAC
		for combC1 in self.__combs[c1]:
			if combC1 in self.__combs[c2]:
				return 0
		return self.__subst

	