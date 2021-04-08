import module.NW_fct as fct
import module.custom as cst
import module.gestion_print as gpt
# Paramètres par défaut
seqA = "tgttacgg"
seqB = "ggttgacta"
liste_symbole = ['|', ':', ':', ':', ' ']
liste_score = [3, -3, -3, -3, -2, -2]
# Customisation type alignement True = génomique , type algorithme True = Needleman
seqA, seqB, liste_score, liste_symbole, type_alignement, type_algorithme = \
    cst.custom_input(seqA, seqB, liste_score, liste_symbole)
# Alignement des sequences
dico_x_aligne, mat_max_traceback = fct.matrix(seqA, seqB, liste_score, liste_symbole, type_alignement, type_algorithme)
# Print des résultats
gpt.print_final(dico_x_aligne, seqA, seqB, mat_max_traceback, type_alignement)

