from tkinter import *
from tkinter.filedialog import *
import custom as cst
import NW_fct as fct
import gestion_print as gpt

dico = {'liste_score':[], 'liste_symbole':[],'score_prot': "","seqA" :"", "seqB":"",'algo': bool, 'alignement':bool, }
def total():
    fenetre = Tk()
    global value_algo
    global value_ali
    global fenetre_score
    global fenetre_symbole
    global fenetre_sequence
    global fenetre_validation
    global valide_fin
    global valide
    global score
    global symb
    fenetre_bas = Frame(fenetre)
    fenetre_bas.pack(side= "bottom")
    fenetre_algo_ali = Frame(fenetre, borderwidth=2, relief=GROOVE)
    fenetre_algo_ali.pack(side = "left")
    fenetre_score = Frame(fenetre, borderwidth=2, relief=GROOVE)
    fenetre_score.pack(side ="left")
    fenetre_symbole = Frame(fenetre, borderwidth=2, relief=GROOVE)
    fenetre_symbole.pack(side = "left")
    fenetre_sequence = Frame(fenetre_bas,borderwidth=2, relief=GROOVE)
    fenetre_sequence.pack(side = "bottom")
    fenetre_validation = Frame(fenetre_bas,borderwidth=2, relief=GROOVE)
    fenetre_validation.pack()
    algo_aliLF = LabelFrame(fenetre_algo_ali, text="customisation des algorithme/alignement", labelanchor='n')
    algo_aliLF.pack()
    value_algo = BooleanVar()
    value_ali = BooleanVar()
    Label(algo_aliLF, text="type d'alignement").pack()
    Radiobutton(algo_aliLF, text="genomique", variable=value_algo, indicatoron=0, value=True).pack()
    Radiobutton(algo_aliLF, text="proteique", variable=value_algo, indicatoron=0, value=False).pack()
    Label(algo_aliLF, text="choix de l'algorithme").pack()
    Radiobutton(algo_aliLF, text="Needleman et Wunch", variable=value_ali, indicatoron=0, value=True).pack()
    Radiobutton(algo_aliLF, text="Smith et Waterman", variable=value_ali, indicatoron=0, value=False).pack()
    valide = Button(algo_aliLF, text="APPLIQUER", command=choix_val, default = 'disable')
    valide.pack()
    valide_fin = Button(fenetre_validation, text="ALLIGNER", command=valid_final, default='disable')
    valide_fin.pack()
    fenetre.mainloop()

def valid_final():
    global haut
    global bas
    global resultat
    type_alignement = dico['alignement']
    type_algorithme = dico['algo']
    seqA = cst.lecture_fasta(dico['seqA'])
    seqB = cst.lecture_fasta(dico['seqB'])
    dico_x_aligne, mat_max_traceback = fct.matrix(seqA, seqB, dico['liste_score'], dico['liste_symbole'],
                                                  type_alignement, type_algorithme)
    fenetre = Tk()
    fenetre_bas = Frame(fenetre)
    fenetre_bas.pack(side="bottom")
    fenetre_score = Frame(fenetre, borderwidth=2, relief=GROOVE)
    fenetre_score.pack()
    fenetre_trace = Frame(fenetre, borderwidth=2, relief=GROOVE)
    fenetre_trace.pack()
    fenetre_scoreLF = LabelFrame(fenetre_score, text="matrice de score ", labelanchor='n')
    fenetre_scoreLF.pack()
    fenetre_traceLF = LabelFrame(fenetre_trace, text="matrice de trace ", labelanchor='n')
    fenetre_traceLF.pack()
    i=0
    while i < len(seqB)+1:
        Label(fenetre_scoreLF, text=dico_x_aligne['0']["matrice score"][i]).pack()
        Label(fenetre_traceLF, text = mat_max_traceback[i]).pack()
        i+=1

    fenetre.mainloop()
    gpt.print_final(dico_x_aligne, seqA, seqB, mat_max_traceback)

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
    sc_match = Spinbox(scoreLF, from_= -50, to = 50, width=5)
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
    valide_score_geno.config(state ='disable')
    score = []
    score.append(int(sc_match.get()))
    score.append(int(sc_pu.get()))
    score.append(int(sc_py.get()))
    score.append(int(sc_mis.get()))
    score.append(int(sc_ouv.get()))
    score.append(int(sc_ext.get()))
    dico['liste_score'] = score


def recup_symb_geno():
    valide_symb_geno.config(state ='disable')
    symb = []
    symb.append(sy_match.get())
    symb.append(sy_pu.get())
    symb.append(sy_py.get())
    symb.append(sy_mis.get())
    symb.append(sy_gap.get())
    dico['liste_symbole'] = symb

def recup_score_prot():
    valide_score_prot.config(state ='disable')
    score =[0,0,0,0]
    score.append(int(sc_ouv.get()))
    score.append(int(sc_ext.get()))
    dico['liste_score'] = score
    dico['score_prot'] = mat_score.get()

def recup_symb_prot():
    valide_symb_prot.config(state ='disable')
    symb = []
    symb.append(sy_match.get())
    symb.append(' ')
    symb.append(' ')
    symb.append(sy_mis.get())
    symb.append(sy_gap.get())
    dico['liste_symbole'] = symb

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
    seqLF = LabelFrame(fenetre_sequence, text="Sequence", labelanchor='n')
    seqLF.pack()
    filepath_seqA = askopenfilename(title="Ouvrir un fichier fasta pour la sequence A",filetypes=[('all files', '.*')])
    Label(seqLF, text = 'Path seqA :').pack()
    Label(seqLF, text = filepath_seqA).pack()
    filepath_seqB = askopenfilename(title="Ouvrir un fichier fasta pour la sequence B", filetypes=[('all files', '.*')])
    Label(seqLF, text='Path seqB :').pack()
    Label(seqLF, text = filepath_seqB).pack()
    dico['seqA']=filepath_seqA
    dico['seqB']=filepath_seqB



def choix_val():
    valide.config(state ='disable')
    type_algo = value_algo.get()
    type_ali = value_ali.get()
    dico['algo'] = type_algo
    dico['alignement'] = type_ali
    if type_algo == True:
        geno()
    else :
        prot()

total()




