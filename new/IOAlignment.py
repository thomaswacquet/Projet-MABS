import re

class IOAlignment:
    def ParseFasta(Fi):
        f = open(Fi, "r")
        currentSeq = ""
        res = dict()

        for l in f:
            l = l.rstrip()

            if l[0] == ">":
                # si la ligne commence par >, stocke le label
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

    def ParseConfig(Fi):
        f = open(Fi, "r")
        res = dict()

        for l in f:
            l = l.rstrip()

            # si la ligne n'est ni vide, ni un commentaire (commençant par #)
            if l != "" and not re.match("^ *#.*", l):

                m = re.match("([a-z]+) *: *([+-]\\d+)", l)
                if not m:
                    raise TypeError
                if m.group(1) == "identity" or m.group(1) == "substitution" or m.group(1) == "indel":
                    res[m.group(1)] = int(m.group(2))
                else:
                    raise TypeError
        return res