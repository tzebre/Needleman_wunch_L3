# fichier qui gère les inputs et la customisation du programme

# return True et True ou False pour un input y ou n et 'erreur' et False si l'input n'est pas y ou n
def true_false(rep):
    if rep == 'y':
        return True, True
    elif rep == 'n':
        return False, True
    else:
        print("la reponsse doit etre True ou False")
        return 'erreur', False


# retourne une sequence genomique ou proteique entrée a la main
def custom_seq_fct():
    seq: str = input("Sequence : ")
    return seq


# retourne une sequence lu depuis un fichier fasta
def lecture_fasta(seq):
    fasta = open(seq, 'r')
    seq = ""
    for line in fasta:
        if line[0] != '>':
            seq += line.replace('\n', '')
    return seq


# retourne la sequence proteique ou genomique en majuscule seulement avec True si elle n'a pas d'erreur
def netoyage_seq(seq, type_alignement):
    seq = seq.upper()
    if type_alignement is True:  # cas alignement genomique
        ok = "AaCcGgTtUu"
    else:  # cas alignement proteique
        ok = "AaRrNnDdCcQqEeGgHhIiLlKkMmFfPpSsTtWwYyVv"
    for nt in seq:
        if nt not in ok:
            print("erreur dans la sequence ", seq)
            return seq, False
    return seq, True


# retourne la liste score_list avec les input des different score possible dans un alignement genomique ou proteique
def custom_score_fct(score_list, type_alignement):
    if type_alignement is True:  # seulement pour un aligenemnt génomique
        score_list[0] = int(input("score match : "))
        score_list[1] = int(input("score mismatch purine/purine : "))
        score_list[2] = int(input("score mismatch pyrimidine/pyrimidine : "))
        score_list[3] = int(input("score autre mismatch : "))
    else:
        print("choix de la matrice : que blosum62 pour l'instant")
    score_list[4] = int(input("score ouverture de gap : "))
    score_list[5] = int(input("score extension de gap : "))
    return score_list


# retourne la liste liste_symbole avec les input des different symbole d'alignement possible
def custom_symbol(liste_symbole, custom_type):
    liste_symbole[0] = input("symbole de match : ")
    liste_symbole[4] = input("symbole gap : ")
    if custom_type is True:  # seuelement en cas d'alignement genomique
        liste_symbole[1] = input("symbole mismatch purine/purine : ")
        liste_symbole[2] = input("symbole mismatch pyrimidine/pirymidine : ")
    liste_symbole[3] = input("symbole autre mismatch : ")
    return liste_symbole


# retourne les sequences, les liste de score/symbole et le type d'alignement
def custom_input(seqA, seqB, liste_score, liste_symbole):
    type_alignement, type_algorithme = True, True
    custom_seq_fasta = bool
    input_ok = False
    while input_ok is False:  # temps que l'input n'est pas bon
        custom_type_input = input("alignement genomique ou proteique ? y/n : ").lower().strip()
        type_alignement, input_ok = true_false(custom_type_input)  # True = genomique False = proteique

    input_ok = False
    while input_ok is False:
        custom_algo_input = input("type d'algorithme Needleman-Wunsch/Smith-Waterman ? y/n : ").lower().strip()
        type_algorithme, input_ok = true_false(custom_algo_input) # True = Needleman Famse = Smith

    input_ok = False
    while input_ok is False:
        custom_seq_input = input("custom sequences ? y/n : ").lower().strip()
        custom_seq, input_ok = true_false(custom_seq_input)
        if custom_seq is True:  # si customisation des sequence est True on demande si on import depuis un fasta
            custom_seq_fasta_input = input("from fasta ? y/n : ").lower().strip()
            custom_seq_fasta, input_ok = true_false(custom_seq_fasta_input)
        if custom_seq is True:  # customisation des sequences
            seqA = ""
            seqB = ""
            impA = False
            impB = False
            while impA is False:  # recommence tant que la sequence est tronqué
                seqA = custom_seq_fct()
                if custom_seq_fasta is True:  # lecture fasta
                    seqA = lecture_fasta(seqA)
                seqA, impA = netoyage_seq(seqA, type_alignement)
            while impB is False:  # recommence tant que la sequence est tronqué
                seqB = custom_seq_fct()
                if custom_seq_fasta is True:  # lecture fasta
                    seqB = lecture_fasta(seqB)
                seqB, impB = netoyage_seq(seqB, type_alignement)
        else:
            seqA, impA = netoyage_seq(seqA, type_alignement)
            seqB, impB = netoyage_seq(seqB, type_alignement)

    input_ok = False
    while input_ok is False:
        custom_score_input = input("custom scores ? y/n : ").lower().strip()
        custom_score, input_ok = true_false(custom_score_input)
        if custom_score is True:  # si le choix de customisation des score est True
            liste_score = custom_score_fct(liste_score, type_alignement)

    input_ok = False
    while input_ok is False:
        custom_symb_input = input("custom symbole ? y/n : ").lower().strip()
        custom_symb, input_ok = true_false(custom_symb_input)
        if custom_symb is True:  # si choix de customisation des symboles
            liste_symbole = custom_symbol(liste_symbole, type_alignement)

        return seqA, seqB, liste_score, liste_symbole, type_alignement, type_algorithme
