import NW_fct as fct
import custom as cst
import gestion_print as gpt
seqA = "atggcgt"
seqB = "atgagt"
liste_symbole = ['|', ':', ':', '!', ' ']
liste_score = [2, 1, 1, -1, -10, -1]
seqA, seqB, liste_score, liste_symbole, type_alignement, algo_type = \
    cst.custom_input(seqA, seqB, liste_score, liste_symbole)
dico_x_aligne, mat_max_traceback = fct.matrix(seqA, seqB, liste_score, liste_symbole, type_alignement, algo_type)
# gpt.print_propre_list_test(mat_max_traceback, seqA, seqB)
gpt.print_final(dico_x_aligne, seqA, seqB)