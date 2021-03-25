# gestion des print

# print des different score de l'alignement
def print_score(score):
    print("Nombre de match :", score[0])
    print("Nombre de mismatch purine :", score[1])
    print("Nombre de mismatch pyrimidine :", score[2])
    print("Nombre de mismatch autre :", score[3])
    print("Nombre de mismatch total :", score[1] + score[2] + score[3])
    print("Nombre de gap ouvert :", score[4])
    print("Nombre de gap étendu :", score[5])
    print("Nombre de gap total :", score[4] + score[5])


# affichage d'une matrice (liste de liste)
def print_propre(mat, seqA, seqB):
    i = 0
    print('\t', end='')
    while i < len(seqA):
        print(seqA[i], '|', end=' ')
        i += 1
    print('')
    for x in range(len(mat)):
        if x >= 1:
            print(seqB[x-1], '|', end=' ')
        else:
            print(' ', '|', end=' ')
        for y in range(len(mat[x])):
            print(mat[x][y][0], '|', end=' ')
        print('')


def print_propre_list_test(mat, seqA, seqB):
    i = 0
    print('\t', end='')
    while i < len(seqA):
        print(seqA[i], '|', end=' ')
        i += 1
    print('')
    for x in range(len(mat)):
        if x >= 1:
            print(seqB[x-1], '|', end=' ')
        else:
            print(' ', '|', end=' ')
        for y in range(len(mat[x])):
            print(mat[x][y], '|', end=' ')
        print('')


# affichage d'une matrice de liste (liste de liste de liste)
def print_propre_list(mat, seqA, seqB):
    i = 0
    print('\t', end='')
    while i < len(seqA):
        print(seqA[i], '|', end=' ')
        i += 1
    print('')
    for x in range(len(mat)):
        if x >= 1:
            print(seqB[x-1], '|', end=' ')
        else:
            print(' ', '|', end=' ')
        for y in range(len(mat[x])):
            print(mat[x][y][0], '|', end=' ')
        print('')


# print les different alignement
def print_final(dico_x_aligne, seqA, seqB):
    print("==============================================================")
    print("Matrice de score : ")
    print_propre_list(dico_x_aligne['1']["matrice score"], seqA, seqB)
    print("==============================================================")
    deja_vue = ""
    for dico in dico_x_aligne:
        if deja_vue != dico_x_aligne[dico]["seq symbole"]:
            print("alignement possible n°", int(dico) + 1, sep='')
            print("Score total de l'alignement : ", dico_x_aligne[dico]["score final"], sep='')
            print(dico_x_aligne[dico]["seqA aligne"])
            print(dico_x_aligne[dico]["seq symbole"])
            print(dico_x_aligne[dico]["seqB aligne"])
            print("matrice de traceback : ")
            print_propre(dico_x_aligne[dico]["matrice traceback"], seqA, seqB)
            print_score(dico_x_aligne[dico]["score"])
            print("-------------------------------------------------------")
            deja_vue = dico_x_aligne[dico]["seq symbole"]
