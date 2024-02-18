from CostIUPAC import CostIUPAC
from tp import DisplayAlign

c = CostIUPAC({"indel": -2, "identity": 1, "substitution": -1}, "iupac.txt")
score = DisplayAlign("sequences.fa", 4,  20,  c, 20)

