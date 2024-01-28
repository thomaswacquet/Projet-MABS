class Cost:
	def __init__(self, costs, alphabet):
		self.__indel = costs["indel"]
		self.__subst = costs["substitution"]
		self.__identity = costs["identity"]
		self.__alphabet = alphabet

	def __call__(self, c1, c2):
		if not (c1 == "" or c1 in self.__alphabet):
			raise ValueError
		if not (c2 == "" or c2 in self.__alphabet):
			raise ValueError
		return self.__compare(c1, c2)

	def __compare(self, c1, c2):
		if c1 == c2:
			return self.__identity
		elif c1 == "" or c2 == "":
			return self.__indel
		return self.__subst
	
