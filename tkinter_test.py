from tkinter import *
from tkinter.filedialog import *

def total():
    fenetre = Tk()
    global value_algo
    global value_ali
    global score
    global fenetre_score
    global fenetre_symbole
    global fenetre_sequence
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
    algo_aliLF = LabelFrame(fenetre_algo_ali, text="customisation des algorithme/alignement", labelanchor='n')
    algo_aliLF.pack()
    value_algo = StringVar()
    value_ali = StringVar()
    Label(algo_aliLF, text="type d'alignement").pack()
    Radiobutton(algo_aliLF, text="genomique", variable=value_algo, indicatoron=0, value='y').pack()
    Radiobutton(algo_aliLF, text="proteique", variable=value_algo, indicatoron=0, value='n').pack()
    Label(algo_aliLF, text="choix de l'algorithme").pack()
    Radiobutton(algo_aliLF, text="Needleman et Wunch", variable=value_ali, indicatoron=0, value='y').pack()
    Radiobutton(algo_aliLF, text="Smith et Waterman", variable=value_ali, indicatoron=0, value='n').pack()
    Button(algo_aliLF, text="APPLIQUER", command=choix_val).pack()
    fenetre.mainloop()


def geno():
    scoreLF = LabelFrame(fenetre_score, text="Score", labelanchor='n')
    scoreLF.pack()
    Label(scoreLF, text="match").pack()
    Spinbox(scoreLF, from_= -20, to = 20, width=5).pack()
    Label(scoreLF, text="pu/pu").pack()
    Spinbox(scoreLF, from_=-20, to=20, width=5).pack()
    Label(scoreLF, text="py/py").pack()
    Spinbox(scoreLF, from_= -20, to = 20, width=5).pack()
    Label(scoreLF, text="mismatch").pack()
    Spinbox(scoreLF, from_=-20, to=20, width=5).pack()
    Label(scoreLF, text="ouverture gap").pack()
    Spinbox(scoreLF, from_= -20, to = 20, width=5).pack()
    Label(scoreLF, text="extention gap").pack()
    Spinbox(scoreLF, from_=-20, to=20, width=5).pack()

    symbLF = LabelFrame(fenetre_symbole, text="Symbole", labelanchor='n')
    symbLF.pack()
    match = StringVar()
    mis_pu = StringVar()
    mis_py = StringVar()
    mis = StringVar()
    ouv = StringVar()
    ext = StringVar()
    Label(symbLF, text="match").pack()
    Entry(symbLF, textvariable= match, width=5, justify = 'center').pack()
    Label(symbLF, text="pu/pu").pack()
    Entry(symbLF, textvariable= mis_pu, width=5, justify = 'center').pack()
    Label(symbLF, text="py/py").pack()
    Entry(symbLF, textvariable= mis_py, width=5, justify = 'center').pack()
    Label(symbLF, text="mismatch").pack()
    Entry(symbLF, textvariable= mis, width=5, justify = 'center').pack()
    Label(symbLF, text="ouverture gap").pack()
    Entry(symbLF, textvariable= ouv, width=5, justify = 'center').pack()
    Label(symbLF, text="extention gap").pack()
    Entry(symbLF, textvariable= ext, width=5, justify = 'center').pack()
    seq()

def seq():
    seqLF = LabelFrame(fenetre_sequence, text="Sequence", labelanchor='n')
    seqLF.pack()
    seqA_path = StringVar()
    filepath_seqA = askopenfilename(title="Ouvrir un fichier fasta pour la sequence A", filetypes=[('all files', '.*')])
    Label(seqLF, text = ('path seqA :' , filepath_seqA)).pack()
    seqB_path = StringVar()
    filepath_seqB = askopenfilename(title="Ouvrir un fichier fasta pour la sequence B", filetypes=[('all files', '.*')])
    Label(seqLF, text=('path seqB :', filepath_seqB)).pack()



def choix_val():
    a = value_algo.get()
    b = value_ali.get()
    print(a, b)
    if a == 'y':
        geno()


total()

