import module.NW_fct as fct
import module.custom as cst
import module.gestion_print as gpt
# Parametre par defaut
seqA = "tgttacgg"
seqB = "ggttgacta"
liste_symbole = ['|', ':', ':', '!', ' ']
liste_score = [3, -3, -3, -3, -2, -2]
# Customisation type alignement True = genomique , type algotithme True = Needleman
seqA, seqB, liste_score, liste_symbole, type_alignement, type_algorithme = \
    cst.custom_input(seqA, seqB, liste_score, liste_symbole)
# alignement des sequences
dico_x_aligne, mat_max_traceback = fct.matrix(seqA, seqB, liste_score, liste_symbole, type_alignement, type_algorithme)
# Print des resultat
gpt.print_final(dico_x_aligne, seqA, seqB, mat_max_traceback)

