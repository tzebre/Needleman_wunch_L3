def matrix_zero_list(col, row):
    matrice = []
    for x in range(row):
        matrice.append([])
        for y in range(col):
            matrice[x].append([0 for _ in range(2)])
    return matrice


def matrix_str(col, row):
    matrice = []
    for x in range(row):
        matrice.append([])
        for y in range(col):
            matrice[x].append('')
    return matrice


def zero_matrix(col, row):
    traceback_matrice = []
    for x in range(row):
        traceback_matrice.append([])
        for y in range(col):
            traceback_matrice[x].append(" ")
    return traceback_matrice
