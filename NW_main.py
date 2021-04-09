import module.NW_fct as fct
import module.custom as cst
import module.gestion_print as gpt
# Paramètres par défaut
seqA = "tgttacgg"
seqB = "ggttgacta"
liste_symbole = ['|', ':', ':', '!', ' ']
liste_score = [2, 1, 1, -1, -10, -1]
# Customisation type alignement True = génomique , type algorithme True = Needleman
seqA, seqB, liste_score, liste_symbole, type_alignement, type_algorithme = \
    cst.custom_input(seqA, seqB, liste_score, liste_symbole)
# Alignement des sequences
liste_dico, mat_max_traceback = fct.matrix(seqA, seqB, liste_score, liste_symbole, type_alignement, type_algorithme)
# Print des résultats
gpt.print_final(liste_dico, seqA, seqB, mat_max_traceback, type_alignement, type_algorithme)
