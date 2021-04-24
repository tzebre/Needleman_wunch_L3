# Import des différents fichiers de fonction
import MATHIEU_Theo_module.MATHIEU_Theo_matrix_prot_dic as mpd  # Fichier ou sont stocké les matrices de score
import MATHIEU_Theo_module.MATHIEU_Theo_gestion_matrice as gm  # Fichier qui gere la creation des matrices
from sys import setrecursionlimit  # Permet de gérer la limite de récursivité
setrecursionlimit(2000)


# Initialisation des premieres lignes et colonnes des matrices de trace et de score selon le type d'algorithme
def matrice_initialise(lenA, lenB, liste_score, type_algorithme):
    score_mat = gm.matrix_zero_list(lenA + 1, lenB + 1)
    traceback_mat = gm.matrix_str(lenA + 1, lenB + 1)
    traceback_mat[0][0] = ' '
    #  Premieres cases correspondantes au debut des sequences en ligne et en colonne
    if type_algorithme is True:  # Needleman
        # La matrice de score prend la valeur d'ouverture de gap et celle de trace la flèche correspondante
        score_mat[0][1] = liste_score[4]
        traceback_mat[0][1] = '⤏'  # Flèche en pointillé pour signifier qu'il s'agit d'une ouverture de gap
        score_mat[1][0] = liste_score[4]
        traceback_mat[1][0] = '⇣'
    else:  # Smith
        max_val = max(liste_score[4], 0)  # La matrice de score prend la valeur maximale entre ouverture de gap et 0
        score_mat[0][1] = max_val
        score_mat[1][0] = max_val
        if max_val == liste_score[4]:  # Si le max est l'ouverture de gap, la matrice de trace reçoit la flèche associée
            traceback_mat[0][1] = '→'
            traceback_mat[1][0] = '↓'
        else:  # Si le max est 0, la matrice de trace est remplie avec un caractère vide
            traceback_mat[0][1] = ' '
            traceback_mat[1][0] = ' '
    for i in range(2, lenA + 1):  # Reste de la première ligne
        if type_algorithme is True:
            score_mat[0][i] = score_mat[0][i - 1] + liste_score[5]
            traceback_mat[0][i] = '→'
        else:
            max_val = max((score_mat[0][i - 1] + liste_score[5]), 0)
            if max_val == (score_mat[0][i - 1] + liste_score[5]):
                traceback_mat[0][i] = '→'
            else:
                traceback_mat[0][i] = ' '
    for j in range(2, lenB + 1):  # Reste de la première colonne
        if type_algorithme is True:
            score_mat[j][0] = score_mat[j - 1][0] + liste_score[5]
            traceback_mat[j][0] = '↓'
        else:
            max_val = max((score_mat[j - 1][0] + liste_score[5]), 0)
            if max_val == (score_mat[j - 1][0] + liste_score[5]):
                traceback_mat[0][j] = '↓'
            else:
                traceback_mat[j][0] = ' '
    return score_mat, traceback_mat


# Creation d'un indice dans un dictionnaire qui regroupe les dictionnaire d'alignements
def creation_alignement(lenA, lenB, nb_alignement, liste_score, type_algorithme):
    dico_x_aligne = {}
    dico_aligne = {"seqA aligne": "", "seqB aligne": "", "score final": 0, "score": list, "seq symbole": "",
                   "matrice score": [], "liste traceback": []}
    dico_x_aligne[str(nb_alignement)] = dico_aligne
    dico_x_aligne[str(nb_alignement)]["matrice score"], traceback_mat = matrice_initialise(lenA, lenB,
                                                                                           liste_score, type_algorithme)
    return dico_x_aligne, traceback_mat


# Remonte la matrice de trace et retourne une liste avec les chemins possibles
def traverse_recursive(matrice_traces, col, row, liste_des_traces, trace):
    # Si la case ne contient pas de fleche (arrivé au bout), on ajoute la trace a la liste de trace possible
    if matrice_traces[row][col][0] not in "↘↓→⤏⇣":
        liste_des_traces.append(trace)
    # Remonte la matrice de trace en rappelant la fonction apres chaque déplacement
    else:
        for symbole in matrice_traces[row][col]:
            if len(trace) >= 1:
                #  Le and permet de ne pas quitter un gap sans passer par une case 'ouverture de gap'
                if symbole == '↘' and trace[-1] in '⤏⇣↘':
                    trace += symbole
                    liste_des_traces = traverse_recursive(matrice_traces, col - 1, row - 1, liste_des_traces, trace)
                    trace = trace[:-1]
                if symbole in '→⤏' and trace[-1] in '⤏⇣→↘':
                    trace += symbole
                    liste_des_traces = traverse_recursive(matrice_traces, col - 1, row, liste_des_traces, trace)
                    trace = trace[:-1]
                if symbole in '↓⇣' and trace[-1] in '⤏⇣↓↘':
                    trace += symbole
                    liste_des_traces = traverse_recursive(matrice_traces, col, row - 1, liste_des_traces, trace)
                    trace = trace[:-1]
            else:
                if symbole == '↘':
                    trace += symbole
                    liste_des_traces = traverse_recursive(matrice_traces, col - 1, row - 1, liste_des_traces, trace)
                    trace = trace[:-1]
                if symbole == '↓' or symbole == '⇣':
                    trace += symbole
                    liste_des_traces = traverse_recursive(matrice_traces, col, row - 1, liste_des_traces, trace)
                    trace = trace[:-1]
                if symbole == '→' or symbole == '⤏':
                    trace += symbole
                    liste_des_traces = traverse_recursive(matrice_traces, col - 1, row, liste_des_traces, trace)
                    trace = trace[:-1]
    return liste_des_traces


# Compare deux nucléotide a aligner et retourne le score associé
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


# Rempli la matrice de trace avec toute les directions possibles
def rempli_symbole(row, col, diag, down, right, matrice_score, mat_max, nw, d_o, r_o):
    # d_o et r_o sont True si le déplacement induit une ouverture de gap
    if d_o is True:
        fd = '⇣'
    else:
        fd = '↓'
    if r_o is True:
        fr = '⤏'
    else:
        fr = '→'
    if matrice_score[row][col] == right:  # Si le score maximal est possible depuis la case de gauche
        mat_max[row][col] += fr
        if right == diag:  # Si le score maximal est possible depuis la case en diagonale et gauche
            mat_max[row][col] += '↘'
        if right == down:
            mat_max[row][col] += fd  # Si le score maximal est possible depuis la case en haut en diagonale et a gauche
    elif matrice_score[row][col] == diag:  # Si le score maximal est possible depuis la case en diagonale
        mat_max[row][col] += '↘'
        if diag == down:  # Si le score maximal est possible depuis la case en haut en diagonale
            mat_max[row][col] += fd
    elif matrice_score[row][col] == down:  # Si le score maximal est possible depuis la case en haut
        mat_max[row][col] += fd
    else:
        if nw is False:  # Si aucune des direction correspond alors on ne mets aucune flèches
            mat_max[row][col] += ' '
    return matrice_score, mat_max


# Rempli la matrice de score et appelle la fonction de remplissage de la matrice de trace
def rempli_score(lenA, lenB, seqA, seqB, matrice_score, traceback_mat, liste_score, alignement, nw):
    if alignement is True:  # Alignement génomique on récupère la matrice de score correspondnate
        dico_score, liste_char = mpd.custom_dic_genomique(liste_score)
    else:  # Alignement protéique
        dico_score, liste_char = mpd.blosum62()
    for row in range(1, lenB + 1):  # Parcour de la matrice pour remplir les score
        for col in range(1, lenA + 1):
            AA = mpd.sort_char(seqA[col - 1], seqB[row - 1], liste_char)
            down = int
            right = int
            diag = matrice_score[row - 1][col - 1] + int(dico_score[AA[0]][AA[1]])
            # Case supérieur est un gap
            # Savoir si il y a un changement de brin pour le gap
            if '⇣' in traceback_mat[row - 1][col] or '↓' in traceback_mat[row - 1][col]:
                d_o = False
                down = matrice_score[row - 1][col] + liste_score[5]  # Décalage en dessous, extension gap
            else:  # Case supérieure n'est pas un gap
                d_o = True
                down = matrice_score[row - 1][col] + liste_score[4]  # Décalage en dessous, ouverture de gap
            # Case gauche est un gap
            # Savoir si il y a un changement de brin pour le gap
            if '⤏' in traceback_mat[row][col - 1] or '→' in traceback_mat[row][col - 1]:
                r_o = False
                right = matrice_score[row][col - 1] + liste_score[5]  # Décalage a droite extension gap
            else:
                r_o = True
                right = matrice_score[row][col - 1] + liste_score[4]  # Décalage a droite ouverture gap
            if nw is True:  # Algorithme de Needleman max de score entre les 3 cases possibles
                max_val = max([diag, down, right])
            else:  # Algorithme de Smith max de score entre les 3 cases possibles et 0
                max_val = max([diag, down, right, 0])
            matrice_score[row][col] = max_val
            matrice_score, traceback_mat = rempli_symbole(row, col, diag, down, right,
                                                          matrice_score, traceback_mat, nw, d_o, r_o)
    return matrice_score, traceback_mat


# Alignement selon l'algorithme de Needleman et Wunch grace a une chaine de caractère de flèche pour un chemin possible
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
        elif s_current == '↓' or s_current == '⇣':
            align1 += '-'
            align2 += seqB[i - 1]
            i -= 1
        elif s_current == '→' or s_current == '⤏':
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


# Alignement selon l'algorithme de Smith et Waterman grace a une chaine de caractère de flèche pour un chemin possible
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
        elif s_current == '↓' or s_current == '⇣':
            align1 += '-'
            align2 += seqB[i - 1]
            i -= 1
        elif s_current == '→' or s_current == '⤏':
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


# Retourne une chaine de caractère composé des symboles pour l'alignement
def symbole_alignement_fct(align1, align2, liste_symbole):
    symbole = ""
    for i in range(0, len(align1)):
        comparaison = symbole_compare(align1[i], align2[i], [0, 1, 2, 3, 4])
        symbole += liste_symbole[comparaison]
    return symbole


# Calcule des scores de l'alignement
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


# Retourne la liste des cases de score maximal pour le depart de l'algorithme de Smith et Waterman
def liste_depart_max(mat_score, seqA, seqB):
    co_score_max = []
    score_maximal = 0
    row = 0
    col = 0
    while row < len(seqB) + 1:
        while col < len(seqA) + 1:
            if mat_score[row][col] == score_maximal:
                co_score_max.append([])
                x = len(co_score_max)
                co_score_max[x - 1].append(row)
                co_score_max[x - 1].append(col)
            if int(mat_score[row][col]) > score_maximal:
                co_score_max = []
                score_maximal = int(mat_score[row][col])
                co_score_max.append([row, col])
            col += 1
        col = 0
        row += 1
    return co_score_max, score_maximal


#  Fonction qui gère tout l'alignement en appellant les autre fonction de ce fichier
def matrix(seqA, seqB, liste_score, liste_symbole, type_alignement, type_algorithme):
    lenA = len(seqA)
    lenB = len(seqB)
    nb_alignement = 0
    score_max = 0
    liste_des_traces = []
    liste_traceback = []
    depart_max = []
    liste_dico = []
    trace = ""
    seqA_align = ""
    seqB_align = ""
    dico_x_aligne, traceback_mat = creation_alignement(lenA, lenB, nb_alignement, liste_score, type_algorithme)
    score_mat = dico_x_aligne[str(nb_alignement)]["matrice score"]
    score_mat, traceback_mat = rempli_score(lenA, lenB, seqA, seqB, score_mat,
                                            traceback_mat, liste_score, type_alignement, type_algorithme)
    if type_algorithme is True:  # Algorithme de Needleman
        liste_traceback = traverse_recursive(traceback_mat, lenA, lenB, liste_des_traces, trace)
        for traceback in liste_traceback:  # Pour chaque trace possible on ajoute un alignement dans le dictionnaire
            dico_x_aligne[str(nb_alignement)] = dico_x_aligne['0'].copy()
            dico_x_aligne[str(nb_alignement)]["liste traceback"] = traceback
            seqA_align, seqB_align = nw_aligne(seqA, seqB, traceback)
            dico_x_aligne[str(nb_alignement)]["score final"] = \
                dico_x_aligne[str(nb_alignement)]["matrice score"][lenB][lenA]
            dico_x_aligne[str(nb_alignement)]["seqA aligne"] = seqA_align
            dico_x_aligne[str(nb_alignement)]["seqB aligne"] = seqB_align
            seq_symbole = symbole_alignement_fct(seqA_align, seqB_align, liste_symbole)
            dico_x_aligne[str(nb_alignement)]["seq symbole"] = seq_symbole
            dico_x_aligne[str(nb_alignement)]["score"] = calcul_score_match_gap_mismatch(liste_symbole, seq_symbole)
            nb_alignement += 1
        liste_dico.append(dico_x_aligne)
    else:  # Algorithme de Smith et Waterman
        depart_max, score_max = liste_depart_max(score_mat, seqA, seqB)
        for i, dep in enumerate(depart_max):  # On ajoute a une liste un dictionnaire par depart possible
            liste_des_traces = []
            liste_traceback = traverse_recursive(traceback_mat, dep[1], dep[0], liste_des_traces, trace)
            nb_alignement = 0
            for traceback in liste_traceback:  # Pour chaque trace possible on ajoute un alignement dans le dictionnaire
                dico_x_aligne[str(nb_alignement)] = dico_x_aligne['0'].copy()
                dico_x_aligne[str(nb_alignement)]["liste traceback"] = traceback
                dico_x_aligne[str(nb_alignement)]["score final"] = \
                    dico_x_aligne[str(nb_alignement)]["matrice score"][dep[0]][dep[1]]
                seqA_align, seqB_align = sw_aligne(seqA, seqB, dep[0], dep[1], traceback)
                dico_x_aligne[str(nb_alignement)]["seqA aligne"] = seqA_align
                dico_x_aligne[str(nb_alignement)]["seqB aligne"] = seqB_align
                seq_symbole = symbole_alignement_fct(seqA_align, seqB_align, liste_symbole)
                dico_x_aligne[str(nb_alignement)]["seq symbole"] = seq_symbole
                dico_x_aligne[str(nb_alignement)]["score"] = calcul_score_match_gap_mismatch(liste_symbole, seq_symbole)
                nb_alignement += 1
            tampon = dico_x_aligne.copy()
            liste_dico.append(tampon)
    return liste_dico, traceback_mat
