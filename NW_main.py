import NW_fct as fct
import custom as cst
import gestion_print as gpt
# Parametre par defaut
seqA = "gcatgcu"
seqB = "gattaca"
liste_symbole = ['|', ':', ':', '!', ' ']
liste_score = [1, -1, -1, -1, -1, -1]
# Customisation type alignement True = genomique , type algotithme True = Needleman
seqA, seqB, liste_score, liste_symbole, type_alignement, type_algorithme = \
    cst.custom_input(seqA, seqB, liste_score, liste_symbole)
# alignement des sequences
dico_x_aligne, mat_max_traceback = fct.matrix(seqA, seqB, liste_score, liste_symbole, type_alignement, type_algorithme)
# Print des resultat
gpt.print_final(dico_x_aligne, seqA, seqB, mat_max_traceback)
