# fichier de gestion de creation de matrice
# creation d'une matrice de liste de 2 indice = 0 pour la matrice de score
def matrix_zero_list(col, row):
    matrice = []
    for x in range(row):
        matrice.append([])
        for y in range(col):
            matrice[x].append([0 for _ in range(2)])
    return matrice


# creation d'une matrice de str vide
def matrix_str(col, row):
    matrice = []
    for x in range(row):
        matrice.append([])
        for y in range(col):
            matrice[x].append('')
    return matrice
