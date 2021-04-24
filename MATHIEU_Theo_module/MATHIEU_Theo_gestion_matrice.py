# Fichier pour la creation de matrices
# Creation d'une matrice avec des valeurs = 0 pour la matrice de score
def matrix_zero_list(col, row):
    matrice = []
    for x in range(row):
        matrice.append([])
        for y in range(col):
            matrice[x].append(0)
    return matrice


# Creation d'une matrice de chaine de caractère vide pour la matrice de traceback
def matrix_str(col, row):
    matrice = []
    for x in range(row):
        matrice.append([])
        for y in range(col):
            matrice[x].append('')
    return matrice
