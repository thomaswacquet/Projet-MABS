import re

class IOAlignment:
    """Classe pour le parsing de fichier fasta et de configuration"""
    def __init__(self, Fi):
        self.__fi = Fi

    def ParseFasta(self):
        """Parse un fichier fasta et retourne un dictionnaire contenant les fasta par labels

        

        Returns:
	    -------
	    res : dict
            dictionnaire des séquences issus du fichier fasta avec comme clés les labels
        """
        f = open(self.__fi, "r")
        currentSeq = ""
        res = dict()

        for l in f:
            l = l.rstrip()

            if l[0] == ">":
                # si la ligne commence par '>', stocker le label
                currentSeq = l[1:]
                res[currentSeq] = ""
            else:
                # sinon ajouter la séquence ADN au label correspondant
                if currentSeq == "":
                    raise TypeError
                if not re.match("[ACTG]*", "ATGC"):
                    raise TypeError
                res[currentSeq] += l

        return res

    def ParseConfig(self):
        """Parse un fichier config et retourne la config (valeurs de couts, alphabet)

        Returns
	    -------
	    res : dict
            dictionnaire des couts et de l'alphabet
        """
        f = open(self.__fi, "r")
        res = dict()
        alphabet = "ATCG"
        combs = dict()

        for l in f:
            l = l.rstrip()

            # si la ligne n'est ni vide, ni un commentaire (commençant par #)
            if l != "" and not re.match("^ *#.*", l):
                m = re.match("([a-z]+) *: *([+-]\\d+)", l)
                if not m:
                    m = re.match("^([A-Z])\t([A-Z](?: [A-Z])*)+", l)
                    if m:
                        letters = m.group(2).split(" ")
                        for letter in letters:
                            if not letter in alphabet:
                                raise TypeError("La lettre est non présente dans l'alphabet")
                        combs[m.group(1)] = letters
                    else:
                        raise TypeError("Erreur de syntaxe du fichier")
                elif m.group(1) == "identity" or m.group(1) == "substitution" or m.group(1) == "indel":
                    res[m.group(1)] = int(m.group(2))
                else:
                    raise TypeError("Erreur de syntaxe du fichier")
        res["alphabet"] = alphabet
        res["combs"] = combs
        return res