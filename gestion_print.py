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


def print_propre(mat, seqA, seqB, type_matrice):
    i = 0
    print('\t', end='\t')
    while i < len(seqA):
        print(seqA[i].center(3), '|', sep='', end='')
        i += 1
    print()
    for row in range(len(mat)):
        if row >= 1:
            print(seqB[row - 1].center(3), '|', sep='', end='')
        else:
            print(''.center(3), '|', sep='', end='')
        for col in range(len(mat[row])):
            if type_matrice is True:  # score
                print(str(mat[row][col][0]).center(3), '|', sep='', end='')
            else:
                print(mat[row][col].center(3), '|', sep='', end='')
        print('')


# print les different alignement
def print_final(dico_x_aligne, seqA, seqB, mat_traceback):
    print("==============================================================")
    print("Matrice de score : ")
    print_propre(dico_x_aligne['1']["matrice score"], seqA, seqB, True)
    print('---------------------------------------------------------------')
    print("Matrice de traceback :")
    print_propre(mat_traceback, seqA, seqB, False)
    print("==============================================================")
    for dico in dico_x_aligne:
        print("alignement possible n°", int(dico) + 1, sep='')
        print("Score total de l'alignement : ", dico_x_aligne[dico]["score final"], sep='')
        print("liste de traceback : ")
        print(dico_x_aligne[dico]["liste traceback"][::-1])
        print(dico_x_aligne[dico]["seqA aligne"])
        print(dico_x_aligne[dico]["seq symbole"])
        print(dico_x_aligne[dico]["seqB aligne"])
        print_score(dico_x_aligne[dico]["score"])
        print("-------------------------------------------------------")
