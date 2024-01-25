import sys
if 'ParseConfig' not in sys.modules:
    from projet_mabs import ParseFasta, ParseConfig

class Cost:
    """Calcule le coût de l'opération dans la matrice de programmation"""
    def __init__(self,  dict_opération, alphabet):
        """
        Parameters:
        --------
        dict_exons : Dict
            les éléments de l'exon
        """
        self.__identité = dict_opération['identity']

        self.__subs = dict_opération['substitution']

        self.__indel = dict_opération['indel']

        self.__alphabet = alphabet
    
    def __compare(self, carac_1, carac_2):
        """compare deux caractères et retourne le coût de l'opération"""
        if carac_1 == carac_2 :
            cost = self.__identité
        elif carac_1 == '' or carac_2 == '' :
            cost = self.__indel
        else :
            cost = self.__subs
        return cost

    def __call__(self, carac_1, carac_2):
        if carac_1 not in self.__alphabet or carac_2 not in self.__alphabet :
            raise ValueError
        return self.__compare(carac_1, carac_2)

if __name__ == '__main__' :
    test = Cost(ParseConfig('score.txt'), ['A', 'C', 'G', 'T'])
    t = test('A','A')
    print(t)
