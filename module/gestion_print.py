# Gestion des prints
import sys
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


# Save les different alignement dans un txt
def print_final_fichier(liste_dico, seqA, seqB, mat_traceback, type_alignement, file_name, liste_score):
    sys.stdout = open(file_name, 'w')
    for d, dico_x_aligne in enumerate(liste_dico):  # Chaque dictionnaire dans cette liste est une case de depart
        if d < 1:  # Si il s'agit du premier dictionnaire
            print('Paramètres : ')
            print('Score')
            if type_alignement is True:
                print('match    purine    pyrimidine    autre   gap ouvert    gap étendu')
                print('{0}\t{1}\t{2}\t{3}\t{4}\t{5}'.format(str(liste_score[0]), str(liste_score[1]),
                                                            str(liste_score[2]), str(liste_score[3]),
                                                            str(liste_score[4]), str(liste_score[5])))
            else:
                print('matrice    gap ouvert    gap étendu')
                print('BLOSUM62\t{0}\t{1}'.format(str(liste_score[4]), str(liste_score[5])))
            print("###############################################################")
            print("Matrice de score : ")
            print_matrice(dico_x_aligne['0']["matrice score"], seqA, seqB, True)
            print('"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""')

            print("Matrice de traceback :")
            print_matrice(mat_traceback, seqA, seqB, False)
            print("==============================================================")
        # Dico_x_aligne contient autan d'indice que de trace possible
        for nb_alignement in range(0, len(dico_x_aligne)):
            dico = str(nb_alignement)
            if int(dico) < 1:  # Si il s'agit du premier alignement
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
                print("alignement possible n°", int(dico) + 1, sep='')
                print("liste de traceback : ")
                print(dico_x_aligne[dico]["liste traceback"][::-1])
                print(dico_x_aligne[dico]["seqA aligne"])
                print(dico_x_aligne[dico]["seq symbole"])
                print(dico_x_aligne[dico]["seqB aligne"])
                print_score(dico_x_aligne[dico]["score"], type_alignement)
        print("-------------------------------------------------------------------------")
    sys.stdout.close()


# Print les différents alignement dans la console
def print_final(liste_dico, seqA, seqB, mat_traceback, type_alignement, type_algorithme, liste_score):
    for d, dico_x_aligne in enumerate(liste_dico):  # Chaque dictionnaire dans cette liste est une case de depart
        if d < 1:  # Si il s'agit du premier dictionnaire
            print('Paramètres : ')
            print('Score')
            if type_alignement is True:
                print('match    purine    pyrimidine    autre   gap ouvert    gap étendu')
                print('{0}\t{1}\t{2}\t{3}\t{4}\t{5}'.format(str(liste_score[0]), str(liste_score[1]),
                                                            str(liste_score[2]), str(liste_score[3]),
                                                            str(liste_score[4]), str(liste_score[5])))
            else:
                print('matrice    gap ouvert    gap étendu')
                print('BLOSUM62\t{0}\t{1}'.format(str(liste_score[4]), str(liste_score[5])))
            print("###############################################################")
            print("###############################################################")
            print("Matrice de score : ")
            print_matrice(dico_x_aligne['0']["matrice score"], seqA, seqB, True)
            print('"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""')
            print("Matrice de traceback :")
            print_matrice(mat_traceback, seqA, seqB, False)
            print("==============================================================")
        # Dico_x_aligne contient autan d'indice que de trace possible
        for nb_alignement in range(0, len(dico_x_aligne)):
            dico = str(nb_alignement)
            if int(dico) < 1:  # Si il s'agit du premier alignement
                if type_algorithme is False:
                    if d % 1 == 0 and d >= 1:
                        print(" ")
                        suite_input = input("appuyez sur entrée pour continuer (depart suivant), quitter avec ctrl+C")
                        suite_input = "  "
                        print("  ")
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
                if nb_alignement % 30 == 0 and nb_alignement >= 1:
                    print(" ")
                    suite_input = \
                        input("appuyez sur entrée pour continuer (30 alignements suivant), quitter avec ctrl+C")
                    suite_input = "  "
                    print("  ")
                print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
                print("alignement possible n°", int(dico) + 1, sep='')
                print("liste de traceback : ")
                print(dico_x_aligne[dico]["liste traceback"][::-1])
                print(dico_x_aligne[dico]["seqA aligne"])
                print(dico_x_aligne[dico]["seq symbole"])
                print(dico_x_aligne[dico]["seqB aligne"])
                print_score(dico_x_aligne[dico]["score"], type_alignement)
        print("-------------------------------------------------------------------------")