# import des different fichier de fonction
import matrix_prot_dic as mpd
import gestion_matrice as gm
blo62, liste_char_blo62 = mpd.BLOSUM62()


# initialisation des premieres ligne et colone de la matrice
def matrice_initialise(lenA, lenB, liste_score):
    score_mat = gm.matrix_zero_list(lenA + 1, lenB + 1)
    traceback_liste_mat = gm.matrix_str(lenA + 1, lenB + 1)
    traceback_liste_mat[0][0] = ' '
    score_mat[0][1][0] = liste_score[4]
    score_mat[0][1][1] = 1
    traceback_liste_mat[0][1] ='→'
    score_mat[1][0][0] = liste_score[4]
    score_mat[1][0][1] = 1
    traceback_liste_mat[1][0] = '↓'
    for i in range(2, lenA+1):
        score_mat[0][i][0] = score_mat[0][i - 1][0] + liste_score[5]
        score_mat[0][i][1] = 1
        traceback_liste_mat[0][i] = '→'
    for j in range(2, lenB+1):
        score_mat[j][0][0] = score_mat[j - 1][0][0] + liste_score[5]
        score_mat[j][0][1] = 1
        traceback_liste_mat[j][0] = '↓'
    return score_mat, traceback_liste_mat


# creation d'un indice dans un dictionaire qui regroupe les alignement
def creation_alignement(lenA, lenB, nb_alignement, liste_score):
    dico_x_aligne = {}
    dico_aligne = {"seqA aligne": "", "seqB aligne": "", "score final": int, "score": list, "seq symbole": "",
                   "matrice score": [], "liste traceback": []}
    dico_x_aligne[str(nb_alignement)] = dico_aligne
    dico_x_aligne[str(nb_alignement)]["matrice score"], traceback_liste_mat = matrice_initialise(lenA,
                                                                                                 lenB, liste_score)
    return dico_x_aligne, traceback_liste_mat


def traverse_recursive(matrice_traces, col, row, liste_des_traces, trace):
    if matrice_traces[row][col][0] not in "↘↓→":
        if len(liste_des_traces) >= 1:
            liste_des_traces.append(trace[1:])
        else:
            liste_des_traces.append(trace)
        trace =''
    for symbole in matrice_traces[row][col]:
        if symbole == '↘':
            trace += '↘'
            liste_des_traces = traverse_recursive(matrice_traces,col-1,row-1,liste_des_traces,trace)
        if symbole == '↓':
            trace += '↓'
            liste_des_traces = traverse_recursive(matrice_traces,col,row-1,liste_des_traces,trace)
        if symbole == '→' :
            trace+='→'
            liste_des_traces = traverse_recursive(matrice_traces,col-1,row,liste_des_traces,trace)
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
def rempli_symbole(i, j, diag, down, right, matrice_score, mat_max):
    if matrice_score[i][j][0] == right:
        mat_max[i][j] += '→'
        if right == diag:
            mat_max[i][j] += '↘'
        if right == down:
            mat_max[i][j] += '↓'
    elif matrice_score[i][j][0] == diag:
        mat_max[i][j] += '↘'
        if diag == down:
            mat_max[i][j] += '↓'
    else:
        mat_max[i][j] += '↓'
    # si on a creer un gap ou mets 1 comme deuxiemme indice de liste sinon 0
    if max(diag, down, right) == down or max(diag, down, right) == right:
        matrice_score[i][j][1] = 1
    else:
        matrice_score[i][j][1] = 0
    return matrice_score, mat_max

# modif ya pas longtemps a verifier avec d'autre test
def rempli_score(lenA, lenB, seqA, seqB, matrice_score, mat_max, liste_score, alignement, nw):
    if alignement is True:
        dico_score, liste_char = mpd.custom_dic_genomique(liste_score)
    else:
        dico_score, liste_char = mpd.BLOSUM62()
    for i in range(1, lenB + 1):
        for j in range(1, lenA+1):
            AA = mpd.sort_char(seqA[j-1], seqB[i-1], liste_char)
            down = int
            right = int
            diag = matrice_score[i - 1][j - 1][0] + int(dico_score[AA[0]][AA[1]])
            if matrice_score[i - 1][j][1] == 1:  # precedent superieur est un gap
                down = matrice_score[i - 1][j][0] + liste_score[5]  # decalage en dessous extention gap
            elif matrice_score[i - 1][j][1] == 0:  # precedement non gap
                down = matrice_score[i - 1][j][0] + liste_score[4]  # decalage en dessous ouverture gap
            else:
                print("erreur down")
            if matrice_score[i][j - 1][1] == 1:  # precedent gauche est un gap
                right = matrice_score[i][j - 1][0] + liste_score[5]  # decalage a droite extention gap
            elif matrice_score[i][j - 1][1] == 0:  # precedent non gap
                right = matrice_score[i][j - 1][0] + liste_score[4]  # decalage a droite ouverture gap
            else:
                print("erreur right")
            if nw is True:
                max_val = max(diag, down, right)
            else:
                max_val = max(diag, down, right, 0)
            matrice_score[i][j][0] = max_val
            matrice_score, mat_max = rempli_symbole(i, j, diag, down, right, matrice_score, mat_max)
    return matrice_score, mat_max


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
        t+=1
    while i > 0:
        align1 += '-'
        align2 += seqB[i - 1]
        i -= 1
        t+=1
    align1 = align1[::-1]
    align2 = align2[::-1]
    return align1, align2

# marche pas 
def sw_aligne(seqA, seqB, i, j, trace):
    lenA = len(seqA)
    lenB = len(seqB)
    score_final = 0
    while i >= 0 and j >= 0:
        s_current = matrice_traceback[i][j]
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
    while j > 0:
        align1 += seqA[j - 1]
        align2 += '-'
        j -= 1
    while i > 0:
        align1 += '-'
        align2 += seqB[i - 1]
        i -= 1
    align1 = align1[::-1]
    align2 = align2[::-1]
    return align1, align2, score_final


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


def liste_depart_max(mat_score,seqA,seqB):
    depart_max = []
    max_score = 0
    i = 0
    j = 0
    while i < len(seqB)+1:
        while j < len(seqA)+1:
            if mat_score[i][j][0] == max_score:
                depart_max.append([])
                x = len(depart_max)
                depart_max[x-1].append(i)
                depart_max[x-1].append(j)
            if int(mat_score[i][j][0]) > max_score:
                depart_max = []
                max_score = int(mat_score[i][j][0])
                depart_max.append([i, j])
            j += 1
        j = 0
        i += 1
    return depart_max


# creation d'alignement selon l'algorithme de needleman et wunsh
def matrix(seqA, seqB, liste_score, liste_symbole, type_alignement, algo_type):
    lenA = len(seqA)
    lenB = len(seqB)
    nb_alignement = 0
    liste_des_traces = []
    trace = ""
    dico_x_aligne, traceback_liste_mat = creation_alignement(lenA, lenB, nb_alignement, liste_score)
    score = dico_x_aligne[str(nb_alignement)]["matrice score"]
    score, mat_max_traceback = rempli_score(lenA, lenB, seqA, seqB,score,
                                                traceback_liste_mat, liste_score, type_alignement, algo_type)
    if algo_type is True:
        list_traceback = traverse_recursive(mat_max_traceback, lenA, lenB, liste_des_traces, trace)
    else:
        depart_max = liste_depart_max(score,seqA,seqB)
        for i,dep in enumerate(depart_max):
            list_traceback = []
            list_traceback.append([])
            liste_des_traces= []
            print(traverse_recursive(mat_max_traceback, dep[1], dep[0], liste_des_traces, trace))
            list_traceback[i-1] = (traverse_recursive(mat_max_traceback, dep[1], dep[0], liste_des_traces, trace)[0])
    dico_x_aligne[str(nb_alignement)]["matrice score"] = score
    for traceback2 in list_traceback:
        dico_x_aligne[str(nb_alignement)] = dico_x_aligne['0'].copy()
        dico_x_aligne[str(nb_alignement)]["liste traceback"] = traceback2
        if algo_type is True:
            seqA_aligne, seqB_align = nw_aligne(seqA, seqB, traceback2)
            dico_x_aligne[str(nb_alignement)]["score final"] = \
                dico_x_aligne[str(nb_alignement)]["matrice score"][lenB][lenA][0]
        else:
            seqA_aligne, seqB_align, dico_x_aligne[str(nb_alignement)]["score final"] = \
                sw_aligne(seqA, seqB,traceback2)
        dico_x_aligne[str(nb_alignement)]["seqA aligne"] = seqA_aligne
        dico_x_aligne[str(nb_alignement)]["seqB aligne"] = seqB_align
        seq_symbole = symbole_alignement_fct(seqA_aligne, seqB_align, liste_symbole)
        dico_x_aligne[str(nb_alignement)]["seq symbole"] = seq_symbole
        dico_x_aligne[str(nb_alignement)]["score"] = calcul_score_match_gap_mismatch(liste_symbole, seq_symbole)
        nb_alignement += 1
    return dico_x_aligne, mat_max_traceback
