import MATHIEU_Theo_module.MATHIEU_Theo_NW_fct as fct
import MATHIEU_Theo_module.MATHIEU_Theo_custom as cst
import MATHIEU_Theo_module.MATHIEU_Theo_gestion_print as gpt
import sys
# Paramètres par défaut
seqA = "atggcgt"
seqB = "atgagt"
liste_symbole = ['|', ':', ':', '!', ' ']
liste_score = [2, 1, 1, -1, -10, -1]
save = False
file_name = 'save.txt'
# Customisation type alignement True = génomique , type algorithme True = Needleman
seqA, seqB, liste_score, liste_symbole, type_alignement, type_algorithme, save = \
    cst.custom_input(seqA, seqB, liste_score, liste_symbole, save)
# Alignement des sequences
liste_dico, mat_max_traceback = fct.matrix(seqA, seqB, liste_score, liste_symbole, type_alignement, type_algorithme)
# Print des résultats
gpt.print_final(liste_dico, seqA, seqB, mat_max_traceback, type_alignement, type_algorithme, liste_score, False)
if save is True:
    sys.stdout = open(file_name, 'w')
    gpt.print_final(liste_dico, seqA, seqB, mat_max_traceback, type_alignement, type_algorithme, liste_score, save)
    sys.stdout.close()
