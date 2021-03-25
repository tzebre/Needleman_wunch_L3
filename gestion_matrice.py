def matrix_zero_list(col, ligne):
    matrice = []
    for x in range(ligne):
        matrice.append([])
        for y in range(col):
            matrice[x].append([0 for _ in range(2)])
    return matrice


def matrix_list_list(col, ligne):
    matrice = []
    for x in range(ligne):
        matrice.append([])
        for y in range(col):
            matrice[x].append([])
    return matrice


def zero_matrix(col, ligne):
    traceback_matrice = []
    for x in range(ligne):
        traceback_matrice.append([])
        for y in range(col):
            traceback_matrice[x].append(" ")
    return traceback_matrice
