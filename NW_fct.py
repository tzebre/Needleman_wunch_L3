# import des different fichier de fonction
import matrix_prot_dic as mpd
import gestion_matrice as gm


# initialisation des premieres ligne et colone de la matrice
def matrice_initialise(lenA, lenB, liste_score, type_algorithme):
    score_mat = gm.matrix_zero_list(lenA + 1, lenB + 1)
    traceback_mat = gm.matrix_str(lenA + 1, lenB + 1)
    traceback_mat[0][0] = ' '
    if type_algorithme is True:
        score_mat[0][1][0] = liste_score[4]
        score_mat[0][1][1] = 1
        traceback_mat[0][1] = '→'
        score_mat[1][0][0] = liste_score[4]
        score_mat[1][0][1] = 1
        traceback_mat[1][0] = '↓'
    else:
        max_val = max(liste_score[4],0)
        score_mat[0][1][0] = max_val
        score_mat[0][1][1] = 1
        score_mat[1][0][0] = max_val
        score_mat[1][0][1] = 1
        if max_val == liste_score[4]:
            traceback_mat[0][1] = '→'
            traceback_mat[1][0] = '↓'
        else:
            traceback_mat[0][1] = ' '
            traceback_mat[1][0] = ' '
    for i in range(2, lenA + 1):
        if type_algorithme is True:
            score_mat[0][i][0] = score_mat[0][i - 1][0] + liste_score[5]
            traceback_mat[0][i] = '→'
        else:
            max_val = max((score_mat[0][i - 1][0] + liste_score[5]), 0)
            if max_val == (score_mat[0][i - 1][0] + liste_score[5]):
                traceback_mat[0][i] = '→'
            else:
                traceback_mat[0][i] = ' '
        score_mat[0][i][1] = 1
    for j in range(2, lenB + 1):
        if type_algorithme is True:
            score_mat[j][0][0] = score_mat[j - 1][0][0] + liste_score[5]
            traceback_mat[j][0] = '↓'
        else:
            max_val = max((score_mat[j - 1][0][0] + liste_score[5]), 0)
            if max_val == (score_mat[j - 1][0][0] + liste_score[5]):
                traceback_mat[0][j] = '↓'
            else:
                traceback_mat[j][0] = ' '
        score_mat[j][0][1] = 1
    return score_mat, traceback_mat


# creation d'un indice dans un dictionaire qui regroupe les alignement
def creation_alignement(lenA, lenB, nb_alignement, liste_score, type_algorithme):
    dico_x_aligne = {}
    dico_aligne = {"seqA aligne": "", "seqB aligne": "", "score final": int, "score": list, "seq symbole": "",
                   "matrice score": [], "liste traceback": []}
    dico_x_aligne[str(nb_alignement)] = dico_aligne
    dico_x_aligne[str(nb_alignement)]["matrice score"], traceback_mat = matrice_initialise(lenA, lenB, liste_score, type_algorithme)
    return dico_x_aligne, traceback_mat


def traverse_recursive(matrice_traces, col, row, liste_des_traces, trace):
    if matrice_traces[row][col][0] not in "↘↓→":
        if len(liste_des_traces) >= 1:
            liste_des_traces.append(trace[1:])
        else:
            liste_des_traces.append(trace)
        trace = ''
    for symbole in matrice_traces[row][col]:
        if symbole == '↘':
            trace += '↘'
            liste_des_traces = traverse_recursive(matrice_traces, col - 1, row - 1, liste_des_traces, trace)
        if symbole == '↓':
            trace += '↓'
            liste_des_traces = traverse_recursive(matrice_traces, col, row - 1, liste_des_traces, trace)
        if symbole == '→':
            trace += '→'
            liste_des_traces = traverse_recursive(matrice_traces, col - 1, row, liste_des_traces, trace)
    return liste_des_traces


# compare deux nucleotique a aligner et retourne le score associé
def symbole_compare(a, b, score_list):
    pu = ['A', 'G']
    py = ['C', 'T', 'U']
    if a == '-' or b == '-':
        return score_list[4]
    if a == b:
        return score_list[0]
    else:
        if a in pu and b in pu:
            return score_list[1]
        elif a in py and b in py:
            return score_list[2]
        else:
            return score_list[3]


# rempli la matrice de score et la matrice avec toute les direction posssible ordre d'ajout si 3 possibilité → | ↘ | ↓
def rempli_symbole(row, col, diag, down, right, matrice_score, mat_max, nw):
    if matrice_score[row][col][0] == right:
        mat_max[row][col] += '→'
        if right == diag:
            mat_max[row][col] += '↘'
        if right == down:
            mat_max[row][col] += '↓'
    elif matrice_score[row][col][0] == diag:
        mat_max[row][col] += '↘'
        if diag == down:
            mat_max[row][col] += '↓'
    elif matrice_score[row][col][0] == down:
        mat_max[row][col] += '↓'
    elif matrice_score[row][col][0]:
        if nw is True:
            mat_max[row][col] += '↓'
        else:
            mat_max[row][col] += '↓'
    # si on a creer un gap ou mets 1 comme deuxiemme indice de liste sinon 0
    if max(diag, down, right) == down or max(diag, down, right) == right:
        matrice_score[row][col][1] = 1
    else:
        matrice_score[row][col][1] = 0
    return matrice_score, mat_max


# modif ya pas longtemps a verifier avec d'autre test
def rempli_score(lenA, lenB, seqA, seqB, matrice_score, traceback_mat, liste_score, alignement, nw):
    if alignement is True:
        dico_score, liste_char = mpd.custom_dic_genomique(liste_score)
    else:
        dico_score, liste_char = mpd.blosum62()
    for row in range(1, lenB + 1):
        for col in range(1, lenA + 1):
            AA = mpd.sort_char(seqA[col - 1], seqB[row - 1], liste_char)
            down = int
            right = int
            diag = matrice_score[row - 1][col - 1][0] + int(dico_score[AA[0]][AA[1]])
            if matrice_score[row - 1][col][1] == 1:  # precedent superieur est un gap
                down = matrice_score[row - 1][col][0] + liste_score[5]  # decalage en dessous extention gap
            elif matrice_score[row - 1][col][1] == 0:  # precedement non gap
                down = matrice_score[row - 1][col][0] + liste_score[4]  # decalage en dessous ouverture gap
            else:
                print("erreur down")
            if matrice_score[row][col - 1][1] == 1:  # precedent gauche est un gap
                right = matrice_score[row][col - 1][0] + liste_score[5]  # decalage a droite extention gap
            elif matrice_score[row][col - 1][1] == 0:  # precedent non gap
                right = matrice_score[row][col - 1][0] + liste_score[4]  # decalage a droite ouverture gap
            else:
                print("erreur right")
            # bug sur max ca sort pas 1 des 3 seul solution trouvé faire symbole down si ruien sinon sa fait vide
            if nw is True:
                max_val = max([diag, down, right])
            else:
                max_val = max([diag, down, right, 0])
            matrice_score[row][col][0] = max_val
            matrice_score, traceback_mat = rempli_symbole(row, col, diag, down, right,
                                                          matrice_score, traceback_mat, nw)
    for i in traceback_mat:
        print(i)
    return matrice_score, traceback_mat


def nw_aligne(seqA, seqB, trace):
    lenA = len(seqA)
    lenB = len(seqB)
    i = lenB
    j = lenA
    t = 0
    align1 = ""
    align2 = ""
    while i > 0 and j > 0:
        s_current = trace[t]
        if s_current == '↘':
            align1 += seqA[j - 1]
            align2 += seqB[i - 1]
            i -= 1
            j -= 1
        elif s_current == '↓':
            align1 += '-'
            align2 += seqB[i - 1]
            i -= 1
        elif s_current == '→':
            align1 += seqA[j - 1]
            align2 += '-'
            j -= 1
        else:
            print("erreur alignement")
            break
        t += 1
    while j > 0:
        align1 += seqA[j - 1]
        align2 += '-'
        j -= 1
        t += 1
    while i > 0:
        align1 += '-'
        align2 += seqB[i - 1]
        i -= 1
        t += 1
    align1 = align1[::-1]
    align2 = align2[::-1]
    return align1, align2


def sw_aligne(seqA, seqB, i, j, trace):
    t = 0
    align1 = ""
    align2 = ""
    while i > 0 and j > 0:
        if t == len(trace):
            break
        s_current = trace[t]
        if s_current == '↘':
            align1 += seqA[j - 1]
            align2 += seqB[i - 1]
            i -= 1
            j -= 1
        elif s_current == '↓':
            align1 += '-'
            align2 += seqB[i - 1]
            i -= 1
        elif s_current == '→':
            align1 += seqA[j - 1]
            align2 += '-'
            j -= 1
        else:
            break
        t += 1
    while j > 0:
        if t == len(trace):
            break
        align1 += seqA[j - 1]
        align2 += '-'
        j -= 1
        t += 1
    while i > 0:
        if t == len(trace):
            break
        align1 += '-'
        align2 += seqB[i - 1]
        i -= 1
        t += 1
    align1 = align1[::-1]
    align2 = align2[::-1]
    return align1, align2


# retourne un string composé des symbole pour l'alignement
def symbole_alignement_fct(align1, align2, liste_symbole):
    symbole = ""
    for i in range(0, len(align1)):
        comparaison = symbole_compare(align1[i], align2[i], [0, 1, 2, 3, 4])
        symbole += liste_symbole[comparaison]
    return symbole


# calcul des core de l'aligenment
def calcul_score_match_gap_mismatch(liste_symbole, symbole_alignement):
    nb_match = 0
    nb_mismatch_pu = 0
    nb_mismatch_py = 0
    nb_mismatch = 0
    nb_gap_ext = 0
    nb_ouverture_gap = 0
    for i, symbole in enumerate(symbole_alignement):
        if symbole == liste_symbole[0]:
            nb_match += 1
        elif symbole == liste_symbole[1]:
            nb_mismatch_pu += 1
        elif symbole == liste_symbole[2]:
            nb_mismatch_py += 1
        elif symbole == liste_symbole[3]:
            nb_mismatch += 1
        elif symbole == liste_symbole[4] and symbole_alignement[i - 1] != liste_symbole[4]:
            nb_ouverture_gap += 1
        else:
            nb_gap_ext += 1
    score = [nb_match, nb_mismatch_pu, nb_mismatch_py, nb_mismatch, nb_ouverture_gap, nb_gap_ext]
    return score


def liste_depart_max(mat_score, seqA, seqB):
    co_score_max = []
    score_maximal = 0
    row = 0
    col = 0
    while row < len(seqB) + 1:
        while col < len(seqA) + 1:
            if mat_score[row][col][0] == score_maximal:
                co_score_max.append([])
                x = len(co_score_max)
                co_score_max[x - 1].append(row)
                co_score_max[x - 1].append(col)
            if int(mat_score[row][col][0]) > score_maximal:
                co_score_max = []
                score_maximal = int(mat_score[row][col][0])
                co_score_max.append([row, col])
            col += 1
        col = 0
        row += 1
    return co_score_max, score_maximal


# creation d'alignement selon l'algorithme de needleman et wunsh
def matrix(seqA, seqB, liste_score, liste_symbole, type_alignement, type_algorithme):
    lenA = len(seqA)
    lenB = len(seqB)
    nb_alignement = 0
    liste_des_traces = []
    liste_traceback = []
    depart_max = []
    trace = ""
    seqA_align = ""
    seqB_align = ""
    dico_x_aligne, traceback_mat = creation_alignement(lenA, lenB, nb_alignement, liste_score, type_algorithme)
    score_mat = dico_x_aligne[str(nb_alignement)]["matrice score"]
    print(type_alignement)
    score_mat, traceback_mat = rempli_score(lenA, lenB, seqA, seqB, score_mat,
                                            traceback_mat, liste_score, type_alignement, type_algorithme)
    if type_algorithme is True:
        liste_traceback = traverse_recursive(traceback_mat, lenA, lenB, liste_des_traces, trace)
    else:
        depart_max, score_max = liste_depart_max(score_mat, seqA, seqB)
        for i, dep in enumerate(depart_max):
            if i == 0:
                liste_traceback = [traverse_recursive(traceback_mat, dep[1], dep[0], liste_des_traces, trace)]
            else:
                liste_traceback.append(traverse_recursive(traceback_mat, dep[1], dep[0], liste_des_traces, trace))
            liste_des_traces = []
    for i, traceback in enumerate(liste_traceback):
        dico_x_aligne[str(nb_alignement)] = dico_x_aligne['0'].copy()
        dico_x_aligne[str(nb_alignement)]["liste traceback"] = traceback
        if type_algorithme is True:
            seqA_align, seqB_align = nw_aligne(seqA, seqB, traceback)
            dico_x_aligne[str(nb_alignement)]["score final"] = \
                dico_x_aligne[str(nb_alignement)]["matrice score"][lenB][lenA][0]
        else:
            dico_x_aligne[str(nb_alignement)]["score final"] = score_max
            for trace in traceback:
                seqA_align, seqB_align = sw_aligne(seqA, seqB, depart_max[i][0], depart_max[i][1], trace)
        dico_x_aligne[str(nb_alignement)]["seqA aligne"] = seqA_align
        dico_x_aligne[str(nb_alignement)]["seqB aligne"] = seqB_align
        seq_symbole = symbole_alignement_fct(seqA_align, seqB_align, liste_symbole)
        dico_x_aligne[str(nb_alignement)]["seq symbole"] = seq_symbole
        dico_x_aligne[str(nb_alignement)]["score"] = calcul_score_match_gap_mismatch(liste_symbole, seq_symbole)
        nb_alignement += 1
    return dico_x_aligne, traceback_mat
