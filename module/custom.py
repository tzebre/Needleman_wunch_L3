# Fichier qui gère les inputs et la customisation du programme

# Retourne True, True/False pour un input y ou n et 'erreur', False si l'input n'est pas y ou n
def true_false(rep):
    if rep == 'y':
        return True, True
    elif rep == 'n':
        return False, True
    else:
        print("la réponse doit être y ou n")
        return 'erreur', False


# Retourne une séquence génomique ou protéique entrée a la main
def custom_seq_fct():
    seq: str = input("Sequence : ")
    return seq


# Retourne une séquence lue depuis un fichier Fasta
def lecture_fasta(seq):
    fasta = open(seq, 'r')
    seq = ""
    for line in fasta:
        if line[0] != '>':
            seq += line.replace('\n', '')
    return seq


# Retourne la séquence protéique ou génomique en majuscule et True si elle n'a pas d'erreurs
def netoyage_seq(seq, type_alignement):
    seq = seq.upper()
    if type_alignement is True:  # cas alignement genomique
        ok = "AaCcGgTtUuNn"
    else:  # cas alignement protéique
        ok = "AaRrNnDdCcQqEeGgHhIiLlKkMmFfPpSsTtWwYyVv"
    for nt in seq:
        if nt not in ok:
            print("erreur dans la sequence ", seq)
            return seq, False
    return seq, True


# Retourne la liste score_list avec les inputs des différents scores possibles de l'alignement
def custom_score_fct(score_list, type_alignement):
    diff_mis = bool
    penalite_input = bool
    if type_alignement is True:  # Pour un alignement génomique
        score_list[0] = int(input("score match : "))
        input_ok = False
        while input_ok is False:
            mis = input("Importance du type de mismatch ? ? y/n : ").lower().strip()
            diff_mis, input_ok = true_false(mis)
        if diff_mis is True:
            score_list[1] = int(input("score mismatch purine/purine : "))
            score_list[2] = int(input("score mismatch pyrimidine/pyrimidine : "))
            score_list[3] = int(input("score autre mismatch : "))
        else:
            score_list[1] = int(input("score de mismatch : "))
            score_list[2] = score_list[1]
            score_list[3] = score_list[1]
    else:  # Pour un alignement protéique
        print("choix de la matrice : que BLOSUM62 pour l'instant")
    input_ok = False
    while input_ok is False:
        penalite = input("Penalité de gap ouvert ? y/n : ").lower().strip()
        penalite_input , input_ok = true_false(penalite)
    if penalite_input is True:
        score_list[4] = int(input("score ouverture de gap : "))
        score_list[5] = int(input("score extension de gap : "))
    else:
        score_list[4] = int(input("score de gap : "))
        score_list[5] = score_list[4]
    return score_list


# Retourne la liste liste_symbole avec les inputs des différents symbole d'alignement possible
def custom_symbol(liste_symbole, custom_type):
    liste_symbole[0] = input("symbole de match : ")
    liste_symbole[4] = input("symbole gap : ")
    if custom_type is True:  # Seulement pour les d'alignements génomiques
        liste_symbole[1] = input("symbole mismatch purine/purine : ")
        liste_symbole[2] = input("symbole mismatch pyrimidine/pyrimidine : ")
    liste_symbole[3] = input("symbole autre mismatch : ")
    return liste_symbole


# Retourne les séquences, les listes de score/symbole et le type d'alignement
def custom_input(seqA, seqB, liste_score, liste_symbole):
    type_alignement, type_algorithme = True, True
    custom_seq_fasta = bool
    input_ok = False
    while input_ok is False:  # Temps que l'input n'est pas bon
        custom_type_input = input("alignement génomique (y) ou protéique (n) ? : ").lower().strip()
        type_alignement, input_ok = true_false(custom_type_input)  # True = génomique False = protéique
    input_ok = False
    while input_ok is False:
        custom_algo_input = input("type d'algorithme Needleman-Wunsch (y) ou Smith-Waterman (n) ? : ").lower().strip()
        type_algorithme, input_ok = true_false(custom_algo_input)  # True = Needleman False = Smith
    input_ok = False
    while input_ok is False:  # Temps que l'input n'est pas bon
        custom_seq_input = input("Custom séquences ? y/n : ").lower().strip()
        custom_seq, input_ok = true_false(custom_seq_input)
    input_ok = False
    if custom_seq is True:  # Si customisation des séquence True on demande si on importe depuis un fichier Fasta
        while input_ok is False:  # Temps que l'input n'est pas bon
            custom_seq_fasta_input = input("Depuis un fichier Fasta ? y/n : ").lower().strip()
            custom_seq_fasta, input_ok = true_false(custom_seq_fasta_input)
    if custom_seq is True:  # Customisation des séquences
        seqA = ""
        seqB = ""
        impA = False
        impB = False
        while impA is False:  # Recommence tant que la séquence est tronqué
            seqA = custom_seq_fct()
            if custom_seq_fasta is True:  # Lecture Fasta
                seqA = lecture_fasta(seqA)
            seqA, impA = netoyage_seq(seqA, type_alignement)
        while impB is False:  # Recommence tant que la séquence est tronquée
            seqB = custom_seq_fct()
            if custom_seq_fasta is True:  # Lecture Fasta
                seqB = lecture_fasta(seqB)
            seqB, impB = netoyage_seq(seqB, type_alignement)
    else:
        seqA, impA = netoyage_seq(seqA, type_alignement)
        seqB, impB = netoyage_seq(seqB, type_alignement)
    input_ok = False
    while input_ok is False:  # Temps que l'input n'est pas bon
        custom_score_input = input("Custom scores ? y/n : ").lower().strip()
        custom_score, input_ok = true_false(custom_score_input)
        if custom_score is True:  # Si le choix de customisation des score est True
            liste_score = custom_score_fct(liste_score, type_alignement)
    input_ok = False
    while input_ok is False:  # Temps que l'input n'est pas bon
        custom_symb_input = input("Custom symboles ? y/n : ").lower().strip()
        custom_symb, input_ok = true_false(custom_symb_input)
        if custom_symb is True:  # Si choix de customisation des symboles
            liste_symbole = custom_symbol(liste_symbole, type_alignement)
            if type_alignement is False:
                liste_symbole[1] = liste_symbole[3]
                liste_symbole[2] = liste_symbole[3]
    return seqA, seqB, liste_score, liste_symbole, type_alignement, type_algorithme
