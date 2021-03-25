# import des different fichier de fonction
import matrix_prot_dic as mpd
import gestion_matrice as gm
blo62, liste_char_blo62 = mpd.BLOSUM62()


# initialisation des premieres ligne et colone de la matrice
def matrice_initialise(lenA, lenB, liste_score):
    score_mat = gm.matrix_zero_list(lenA + 1, lenB + 1)
    traceback_liste_mat = gm.matrix_list_list(lenA + 1, lenB + 1)
    traceback_liste_mat[0][0].append(' ')
    score_mat[0][1][0] = liste_score[4]
    score_mat[0][1][1] = 1
    traceback_liste_mat[0][1].append('→')
    score_mat[1][0][0] = liste_score[4]
    score_mat[1][0][1] = 1
    traceback_liste_mat[1][0].append('↓')
    for i in range(2, lenA+1):
        score_mat[0][i][0] = score_mat[0][i - 1][0] + liste_score[5]
        score_mat[0][i][1] = 1
        traceback_liste_mat[0][i].append('→')
    for j in range(2, lenB+1):
        score_mat[j][0][0] = score_mat[j - 1][0][0] + liste_score[5]
        score_mat[j][0][1] = 1
        traceback_liste_mat[j][0].append('↓')
    return score_mat, traceback_liste_mat


# creation d'un indice dans un dictionaire qui regroupe les alignement
def creation_alignement(lenA, lenB, nb_alignement, liste_score):
    dico_x_aligne = {}
    dico_aligne = {"seqA aligne": "", "seqB aligne": "", "score final": int, "score": list, "seq symbole": "",
                   "matrice score": [], "matrice traceback": []}
    dico_x_aligne[str(nb_alignement)] = dico_aligne
    dico_x_aligne[str(nb_alignement)]["matrice score"], traceback_liste_mat = matrice_initialise(lenA,
                                                                                                 lenB, liste_score)
    return dico_x_aligne, traceback_liste_mat


# separe en plsuieur matrice la matrice de traceback avec toute les directions possible
def sep_mat_max_traceback(lenA, lenB, mat_max_alignement):
    list_max_alignement = []
    matrice_min = gm.zero_matrix(lenA + 1, lenB + 1)
    matrice_max = gm.zero_matrix(lenA + 1, lenB + 1)
    for x in range(0, lenB+1):
        for y in range(0, lenA+1):
            matrice_min[x][y] = mat_max_alignement[x][y][0]
            matrice_max[x][y] = mat_max_alignement[x][y][len(mat_max_alignement[x][y])-1]
    list_max_alignement.append(matrice_min)
    list_max_alignement.append(matrice_max)
    return list_max_alignement


# voir si on peut pas supp
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
    # voir si on  passe pas dans la fct d'en dessous peut etre plsu logique

    max_val = max(diag, down, right)
    matrice_score[i][j][0] = max_val  # retourne le max des 3 possibilité d'alignement
    if matrice_score[i][j][0] == right:
        mat_max[i][j].append('→')
        if right == diag:
            mat_max[i][j].append('↘')
        if right == down:
            mat_max[i][j].append('↓')
    elif matrice_score[i][j][0] == diag:
        mat_max[i][j].append('↘')
        if diag == down:
            mat_max[i][j].append('↓')
    else:
        mat_max[i][j].append('↓')
    # si on a creer un gap ou mets 1 comme deuxiemme indice de liste sinon 0
    if max(diag, down, right) == down or max(diag, down, right) == right:
        matrice_score[i][j][1] = 1
    else:
        matrice_score[i][j][1] = 0
    return matrice_score, mat_max

# modif ya pas longtemps a verifier avec d'autre test
def rempli_score(lenA, lenB, seqA, seqB, matrice_score, mat_max, liste_score, alignement):
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
            matrice_score, mat_max = rempli_symbole(i, j, diag, down, right, matrice_score, mat_max)
    return matrice_score, mat_max


# remonte la matrice de traceback pour creer un alignement selon needleman et wunsh
def nw_aligne(seqA, seqB, lenA, lenB, matrice_traceback):
    i = lenB
    j = lenA
    align1 = ""
    align2 = ""
    while i > 0 and j > 0:
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
            print("erreur alignement")
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
    return align1, align2


def sw_aligne(seqA, seqB, lenA, lenB, matrice_traceback, matrice_score):
    i = lenB
    j = lenA
    l = lenA
    c = lenB
    align1 = ""
    align2 = ""
    col = 0
    line = 0
    max_line = l
    max_col = c
    position_max = [0, 0]
    while i > 0:
        if matrice_score[i][j][0] > line:
            line = matrice_score[i][j][0]
            l = i
            max_line = l
        i -= 1
    i = lenB
    j = lenA
    while j > 0:
        if matrice_score[i][j][0] > col:
            col = matrice_score[i][j][0]
            c = j
            max_col = c
        j -= 1
    if col >= line:
        position_max[0] = max_line
        position_max[1] = max_col
        x = 0
        while x > max_col:
            align1 += seqA[x - 1]
            align2 += '-'
            matrice_traceback[max_line][x] = '↓'
            x -= 1
    else:
        position_max[0] = max_line
        position_max[1] = max_col
        x = 0
        while x > max_line:
            align1 += '-'
            align2 += seqB[x - 1]
            matrice_traceback[x][max_col] = '→'
            x -= 1
    score_final = matrice_score[max_line][max_col][0]
    i = lenB
    j = lenA
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


# creation d'alignement selon l'algorithme de needleman et wunsh
def matrix(seqA, seqB, liste_score, liste_symbole, type_alignement, algo_type):
    lenA = len(seqA)
    lenB = len(seqB)
    nb_alignement = 0
    dico_x_aligne, traceback_liste_mat = creation_alignement(lenA, lenB, nb_alignement, liste_score)
    score = dico_x_aligne[str(nb_alignement)]["matrice score"]
    if type_alignement is True:
        score, mat_max_traceback = rempli_score(lenA, lenB, seqA, seqB,score, traceback_liste_mat, liste_score, True)
    else:
        score, mat_max_traceback = rempli_score(lenA, lenB, seqA, seqB,score, traceback_liste_mat, liste_score, False)
    dico_x_aligne[str(nb_alignement)]["matrice score"] = score
    list_traceback = sep_mat_max_traceback(lenA, lenB, mat_max_traceback)
    for traceback2 in list_traceback:
        dico_x_aligne[str(nb_alignement)] = dico_x_aligne['0'].copy()
        dico_x_aligne[str(nb_alignement)]["matrice traceback"] = traceback2.copy()
        if algo_type is True:
            seqA_aligne, seqB_align = nw_aligne(seqA, seqB, lenA, lenB, traceback2)
            dico_x_aligne[str(nb_alignement)]["score final"] = \
                dico_x_aligne[str(nb_alignement)]["matrice score"][lenB][lenA][0]
        else:
            seqA_aligne, seqB_align, dico_x_aligne[str(nb_alignement)]["score final"] = \
                sw_aligne(seqA, seqB, lenA, lenB, traceback2, score)
        dico_x_aligne[str(nb_alignement)]["seqA aligne"] = seqA_aligne
        dico_x_aligne[str(nb_alignement)]["seqB aligne"] = seqB_align
        seq_symbole = symbole_alignement_fct(seqA_aligne, seqB_align, liste_symbole)
        dico_x_aligne[str(nb_alignement)]["seq symbole"] = seq_symbole
        dico_x_aligne[str(nb_alignement)]["score"] = calcul_score_match_gap_mismatch(liste_symbole, seq_symbole)
        nb_alignement += 1
    return dico_x_aligne, mat_max_traceback
