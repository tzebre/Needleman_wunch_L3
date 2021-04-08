from tkinter.filedialog import *
import module.custom as cst
import module.NW_fct as fct
import module.gestion_print as gpt

dico = {'liste_score':[], 'liste_symbole':[],'score_prot': "","seqA" :"", "seqB":"",'algo': bool, 'alignement':bool }
def total():
    fenetre = Tk()
    global value_algo
    global value_ali
    global fenetre_score
    global fenetre_symbole
    global fenetre_sequence
    global fenetre_validation
    global fenetre_algo_ali
    global valide_fin
    global valide
    global score
    global symb
    global prms_score
    global prms_symb
    global quiter
    fenetre_bas = Frame(fenetre)
    fenetre_bas.pack(side= "bottom", fill = BOTH)
    fenetre_algo_ali = Frame(fenetre, borderwidth=2, relief=GROOVE)
    fenetre_algo_ali.pack(side = "left", fill = BOTH)
    fenetre_score = Frame(fenetre, borderwidth=2, relief=GROOVE)
    fenetre_score.pack(side ="left", fill = BOTH)
    fenetre_symbole = Frame(fenetre, borderwidth=2, relief=GROOVE)
    fenetre_symbole.pack(side = "left", fill = BOTH)
    fenetre_sequence = Frame(fenetre_bas,borderwidth=2, relief=GROOVE)
    fenetre_sequence.pack(side = "bottom", fill = BOTH)
    fenetre_validation = Frame(fenetre_bas,borderwidth=2, relief=GROOVE)
    fenetre_validation.pack(fill = BOTH)
    algo_aliLF = LabelFrame(fenetre_algo_ali, text="Customisation des algorithme/alignement", labelanchor='n')
    algo_aliLF.pack(fill = BOTH)
    parametres = LabelFrame(fenetre_algo_ali, text="Parametres", labelanchor='n')
    parametres.pack(fill = BOTH)
    prms_score = Label(parametres, text="")
    prms_symb = Label(parametres, text="")
    value_algo = BooleanVar()
    value_ali = BooleanVar()
    Label(algo_aliLF, text="Type d'alignement").pack(fill = BOTH)
    Radiobutton(algo_aliLF, text="Genomique", variable=value_ali, indicatoron=0, value=True).pack()
    Radiobutton(algo_aliLF, text="Proteique", variable=value_ali, indicatoron=0, value=False).pack()
    Label(algo_aliLF, text="Choix de l'algorithme").pack(fill = BOTH)
    Radiobutton(algo_aliLF, text="Needleman et Wunch", variable=value_algo, indicatoron=0, value=True).pack()
    Radiobutton(algo_aliLF, text="Smith et Waterman", variable=value_algo, indicatoron=0, value=False).pack()
    valide = Button(algo_aliLF, text="APPLIQUER", command=choix_val, default = 'disable')
    valide.pack()
    valide_fin = Button(fenetre_validation, text="ALLIGNER", command=valid_final, default='disable')
    valide_fin.pack()
    Label(fenetre_bas, text = "Copyright © 2021 MATHIEU Theo - Tous droits réservés").pack(side = "bottom")
    fenetre.mainloop()

def valid_final():
    valide_fin.config(state ='disable')
    global haut
    global bas
    global resultat
    type_alignement = dico['alignement']
    type_algorithme = dico['algo']
    seqA = cst.lecture_fasta(dico['seqA'])
    seqA , bool= cst.netoyage_seq(seqA, type_alignement)
    seqB = cst.lecture_fasta(dico['seqB'])
    seqB, bool = cst.netoyage_seq(seqB, type_alignement)
    dico_x_aligne, mat_max_traceback = fct.matrix(seqA, seqB, dico['liste_score'], dico['liste_symbole'],
                                                  type_alignement, type_algorithme)
    fenetre = Tk()
    fenetre_haut = Frame(fenetre)
    fenetre_haut.pack(fill =BOTH)
    fenetre_bas = Frame(fenetre)
    fenetre_bas.pack(side="bottom",fill = BOTH)
    Label(fenetre_bas, text="Copyright © 2021 MATHIEU Theo - Tous droits réservés").pack()
    matrice = Frame(fenetre_haut)
    matrice.pack(side= 'left')
    alignement = Frame(fenetre_haut)
    alignement.pack()
    fenetre_score = Frame(matrice, borderwidth=2, relief=GROOVE)
    fenetre_score.pack()
    fenetre_trace = Frame(matrice, borderwidth=2, relief=GROOVE)
    fenetre_trace.pack()
    fenetre_scoreLF = LabelFrame(fenetre_score, text="matrice de score ", labelanchor='n')
    fenetre_scoreLF.pack()
    fenetre_traceLF = LabelFrame(fenetre_trace, text="matrice de trace ", labelanchor='n')
    fenetre_traceLF.pack()
    i=0
    mat_score = []
    mat_trace = []
    for row in range(len(seqB)+1):
        mat_score.append([])
        mat_trace.append([])
        for col in range(len(seqA)+1):
            mat_score[row] += (str(dico_x_aligne['0']["matrice score"][row][col][0]).center(3)+'|')
            mat_trace[row] += (str(mat_max_traceback[row][col]).center(3)+'|')
    while i < len(seqB)+1:
        Label(fenetre_scoreLF, text=mat_score[i]).pack()
        Label(fenetre_traceLF, text = mat_trace[i]).pack()
        i+=1
    ali = Frame(alignement)
    ali.pack()
    alibis = Frame(alignement)
    alibis.pack()
    aliter = Frame(alignement)
    aliter.pack()
    for dic in dico_x_aligne:
        no = ('alignement n°'+ str(int(dic)+1))
        if 3 < int(dic) < 6:
            no = LabelFrame(alibis, text=no, labelanchor='n')
        elif int(dic) > 6:
            no = LabelFrame(aliter, text=no, labelanchor='n')
        else :
            no = LabelFrame(ali, text=no, labelanchor='n')
        no.pack(side = 'left')
        total = ('score total : ' + str(dico_x_aligne['0']["score final"]))
        Label(no, text = str(total)).pack()
        Label(no, text = dico_x_aligne[dic]["liste traceback"][::-1]).pack()
        Label(no, text = dico_x_aligne[dic]["seqA aligne"]).pack()
        Label(no, text = dico_x_aligne[dic]["seq symbole"]).pack()
        Label(no, text = dico_x_aligne[dic]["seqB aligne"]).pack()
        score = dico_x_aligne[dic]['score']
        Label(no, text = ("Nombre de match : "+ str(score[0]))).pack()
        if type_alignement is True:
            Label(no, text = ("Nombre de mismatch purine : " + str(score[1]))).pack()
            Label(no, text = ("Nombre de mismatch pyrimidine : " + str(score[2]))).pack()
        Label(no, text = ("Nombre de mismatch autre : "+ str(score[3])))
        Label(no, text = ("Nombre de gap ouvert : "+str(score[4]))).pack()
        Label(no, text = ("Nombre de gap étendu : " + str(score[5]))).pack()
        Label(no, text = ("Nombre de gap total : "+ str(score[4] + score[5]))).pack()
    gpt.print_final(dico_x_aligne, seqA, seqB, mat_max_traceback)
    matrice.mainloop()

def geno():
    global valide_score_geno
    global valide_symb_geno
    global sc_match
    global sc_pu
    global sc_py
    global sc_mis
    global sc_ouv
    global sc_ext
    scoreLF = LabelFrame(fenetre_score, text="Score", labelanchor='n')
    scoreLF.pack()
    sc_match = IntVar()
    sc_pu = IntVar()
    sc_py = IntVar()
    sc_mis = IntVar()
    sc_ouv = IntVar()
    sc_ext = IntVar()
    Label(scoreLF, text="match").pack()
    sc_match = Spinbox(scoreLF,from_ = -50, to=50,  width=5)
    sc_match.pack()
    Label(scoreLF, text="pu/pu").pack()
    sc_pu = Spinbox(scoreLF, from_=-50, to=50, width=5)
    sc_pu.pack()
    Label(scoreLF, text="py/py").pack()
    sc_py = Spinbox(scoreLF, from_= -50, to = 50, width=5)
    sc_py.pack()
    Label(scoreLF, text="mismatch").pack()
    sc_mis = Spinbox(scoreLF, from_=-50, to=50, width=5)
    sc_mis.pack()
    Label(scoreLF, text="ouverture gap").pack()
    sc_ouv = Spinbox(scoreLF, from_= -50, to = 50, width=5)
    sc_ouv.pack()
    Label(scoreLF, text="extention gap").pack()
    sc_ext = Spinbox(scoreLF, from_=-50, to=50, width=5)
    sc_ext.pack()
    valide_score_geno = Button(scoreLF, text="APPLIQUER", command=recup_score_geno, default = 'disable')
    valide_score_geno.pack()
    symbLF = LabelFrame(fenetre_symbole, text="Symbole", labelanchor='n')
    symbLF.pack()
    global sy_match
    global sy_pu
    global sy_py
    global sy_mis
    global sy_gap
    sy_match = StringVar()
    sy_pu = StringVar()
    sy_py = StringVar()
    sy_mis = StringVar()
    sy_gap = StringVar()
    Label(symbLF, text="match").pack()
    Entry(symbLF, textvariable= sy_match, width=5, justify = 'center').pack()
    Label(symbLF, text="pu/pu").pack()
    Entry(symbLF, textvariable= sy_pu, width=5, justify = 'center').pack()
    Label(symbLF, text="py/py").pack()
    Entry(symbLF, textvariable= sy_py, width=5, justify = 'center').pack()
    Label(symbLF, text="mismatch").pack()
    Entry(symbLF, textvariable= sy_mis, width=5, justify = 'center').pack()
    Label(symbLF, text="gap").pack()
    Entry(symbLF, textvariable= sy_gap, width=5, justify = 'center').pack()
    valide_symb_geno = Button(symbLF, text="APPLIQUER", command=recup_symb_geno, default = 'disable')
    valide_symb_geno.pack()
    seq()


def recup_score_geno():
    score = []
    score.append(int(sc_match.get()))
    score.append(int(sc_pu.get()))
    score.append(int(sc_py.get()))
    score.append(int(sc_mis.get()))
    score.append(int(sc_ouv.get()))
    score.append(int(sc_ext.get()))
    dico['liste_score'] = score
    prms_score['text'] = score
    prms_score.pack()


def recup_symb_geno():
    symb = []
    symb.append(sy_match.get())
    symb.append(sy_pu.get())
    symb.append(sy_py.get())
    symb.append(sy_mis.get())
    symb.append(sy_gap.get())
    dico['liste_symbole'] = symb
    prms_symb['text'] = symb
    prms_symb.pack()

def recup_score_prot():
    score =[0,0,0,0]
    score.append(int(sc_ouv.get()))
    score.append(int(sc_ext.get()))
    dico['liste_score'] = score
    dico['score_prot'] = mat_score.get()
    prms_score['text'] = (str(score[4])+ '  '+str(score[5]) +'  '+dico['score_prot'])
    prms_score.pack()

def recup_symb_prot():
    symb = []
    symb.append(sy_match.get())
    symb.append(' ')
    symb.append(' ')
    symb.append(sy_mis.get())
    symb.append(sy_gap.get())
    dico['liste_symbole'] = symb
    prms_symb['text'] = (symb[0]+ '  '+ symb[3]+'  '+symb[4])
    prms_symb.pack()

def prot():
    global valide_score_prot
    global valide_symb_prot
    global sc_ouv
    global sc_ext
    global mat_score
    scoreLF = LabelFrame(fenetre_score, text="Score", labelanchor='n')
    scoreLF.pack()
    mat_score = StringVar()
    Label(scoreLF, text="matrice").pack()
    Radiobutton(scoreLF, text="BLOSUM62", variable=mat_score, indicatoron=0, value='blosum62').pack()
    Label(scoreLF, text="ouverture gap").pack()
    sc_ouv = Spinbox(scoreLF, from_= -50, to = 50, width=5)
    sc_ouv.pack()
    Label(scoreLF, text="extention gap").pack()
    sc_ext = Spinbox(scoreLF, from_=-50, to=50, width=5)
    sc_ext.pack()
    valide_score_prot = Button(scoreLF, text="APPLIQUER", command=recup_score_prot, default = 'disable')
    valide_score_prot.pack()
    symbLF = LabelFrame(fenetre_symbole, text="Symbole", labelanchor='n')
    symbLF.pack()
    global sy_match
    global sy_mis
    global sy_gap
    sy_match = StringVar()
    sy_mis = StringVar()
    sy_gap = StringVar()
    Label(symbLF, text="match").pack()
    Entry(symbLF, textvariable= sy_match, width=5, justify = 'center').pack()
    Label(symbLF, text="mismatch").pack()
    Entry(symbLF, textvariable= sy_mis, width=5, justify = 'center').pack()
    Label(symbLF, text="gap").pack()
    Entry(symbLF, textvariable= sy_gap, width=5, justify = 'center').pack()
    valide_symb_prot = Button(symbLF, text="APPLIQUER", command=recup_symb_prot, default = 'disable')
    valide_symb_prot.pack()
    seq()

def seq():
    filepath_seqA = askopenfilename(title="Ouvrir un fichier fasta pour la sequence A",filetypes=[('all files', '.*')])
    Label(fenetre_algo_ali, text = 'Path seqA :').pack()
    Label(fenetre_algo_ali, text = filepath_seqA).pack()
    filepath_seqB = askopenfilename(title="Ouvrir un fichier fasta pour la sequence B", filetypes=[('all files', '.*')])
    Label(fenetre_algo_ali, text='Path seqB :').pack()
    Label(fenetre_algo_ali, text = filepath_seqB).pack()
    dico['seqA']=filepath_seqA
    dico['seqB']=filepath_seqB


def choix_val():
    valide.config(state ='disable')
    type_algo = value_algo.get()
    type_ali = value_ali.get()
    dico['algo'] = type_algo
    dico['alignement'] = type_ali
    if type_ali == True:
        geno()
    else :
        print("prot")
        prot()

total()




