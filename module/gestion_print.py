# Gestion des prints

# Print des différents score de l'alignement
def print_score(score, type_alignement):
    print("Nombre de match :", score[0])
    if type_alignement is True:
        print("Nombre de mismatch purine :", score[1])
        print("Nombre de mismatch pyrimidine :", score[2])
        print("Nombre de mismatch autre :", score[3])
        print("Nombre de mismatch total :", score[1] + score[2] + score[3])
    else:
        print("Nombre de mismatch :", score[1] + score[2] + score[3])
    print("Nombre de gap ouvert :", score[4])
    print("Nombre de gap étendu :", score[5])
    print("Nombre de gap total :", score[4] + score[5])


# Print des matrices de score et de traceback. Type matrice True = score False = traceback
def print_matrice(mat, seqA, seqB, type_matrice):
    i = 0
    print('\t', end='')
    while i < len(seqA):
        print(seqA[i].center(3), '|', sep='', end='')  # Permet de toujours avoir au minimum un écart de 3 entre deux |
        i += 1
    print()
    for row in range(len(mat)):
        if row >= 1:
            print(seqB[row - 1].center(3), '|', sep='', end='')
        else:
            print(''.center(3), '|', sep='', end='')
        for col in range(len(mat[row])):
            if type_matrice is True:  # Matrice de score
                print(str(mat[row][col][0]).center(3), '|', sep='', end='')
            else:  # Matrice de trace
                print(mat[row][col].center(3), '|', sep='', end='')
        print('')


# Print les différents alignement
def print_final(liste_dico, seqA, seqB, mat_traceback, type_alignement, type_algorithme):
    for d, dico_x_aligne in enumerate(liste_dico):
        if d < 1:
            print("==============================================================")
            print("Matrice de score : ")
            print_matrice(dico_x_aligne['0']["matrice score"], seqA, seqB, True)
            print('---------------------------------------------------------------')
            print("Matrice de traceback :")
            print_matrice(mat_traceback, seqA, seqB, False)
            print("==============================================================")
        for nb_alignement in range(0,len(dico_x_aligne)):
            dico = str(nb_alignement)
            if d < 1:
                if type_algorithme is False:
                    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
                    print("Depart possible n°", d + 1)
                print("alignement possible n°", int(dico) + 1, sep='')
                print("Score total de l'alignement : ", dico_x_aligne[dico]["score final"], sep='')
                print("Chemin de traceback : ")
                print(dico_x_aligne[dico]["liste traceback"][::-1])
                print(dico_x_aligne[dico]["seqA aligne"])
                print(dico_x_aligne[dico]["seq symbole"])
                print(dico_x_aligne[dico]["seqB aligne"])
                print_score(dico_x_aligne[dico]["score"], type_alignement)
            else:
                print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
                print("Depart possible n°", d + 1)
                print("liste de traceback : ")
                print(dico_x_aligne[dico]["liste traceback"][::-1])
                print(dico_x_aligne[dico]["seqA aligne"])
                print(dico_x_aligne[dico]["seq symbole"])
                print(dico_x_aligne[dico]["seqB aligne"])
                print_score(dico_x_aligne[dico]["score"], type_alignement)
        print("-------------------------------------------------------")
