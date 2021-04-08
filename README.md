# Needleman_Wunch/Smith_Waterman
**Auteur : MATHIEU Theo**  
Programme réalisé en Python dans le cadre de la L3 B.I.S.M. pour l'U.E. bio-informatique.

## Table des matières
1. [Fonctionalité](#fonc)
2. [Principe d'alignement](#principe)
    1. [Needleman et Wunch](#NE)
    2. [Smith et Waterman](#SW)
3. [Notice d'utilisation](#notice)
4. [Fonctionement](#fonctionement)
    1. [Needleman et Wunch](#fNE)
    2. [Smith et Waterman](#fSW)
5. [Résultats](#resultats)
   1. [Fenêtre graphique](#graph)
    2. [Ligne de comande](#cmd)
6. [Améliorations](#amelioration)

## Fonctionalité  <a id="fonc"></a>
 - [x] Alignement selon l'algorithme de Needleman-Wunch ou Smith-Waterman 
    - De séquences protéiques ou génomiques 
    - Depuis un fichier Fasta ou une entrée manuelle de la séquence 
    - Choix des scores et des symboles d'alignement possible
    
 - [x] Utilisation en fenêtre graphique grâce à tkinter (peut lisible pour de tres longue sequences) ou en 
   commande console 

## Principe d'alignement <a id="principe"></a>
### Needleman et wunch <a id="NE"></a>
>L'algorithme Needleman – Wunsch est un algorithme utilisé en bio-informatique pour aligner des séquences protéiques ou 
nucléotidiques. C'était l'une des premières applications de la programmation dynamique pour comparer des séquences 
biologiques. L'algorithme a été développé par Saul B. Needleman et Christian D. Wunsch et publié en 1970.  Il est
également parfois appelé algorithme de correspondance optimal et technique d'alignement global.L'algorithme attribue un score à chaque alignement
possible, et le but de l'algorithme est de trouver tous les alignements possibles ayant le score le plus élevé.
[Source: Wikipedia](https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm)
### Smith et Waterman <a id="SW"></a>
>L'algorithme de Smith – Waterman effectue l'alignement des séquences locales; c'est-à-dire pour déterminer des régions 
similaires entre deux chaînes de séquences d'acide nucléique ou de séquences protéiques. 
Au lieu de regarder la séquence entière, l'algorithme de Smith – Waterman compare des segments de toutes les 
longueurs possibles et optimise la mesure de similarité. L'algorithme a été proposé pour la première fois par 
Temple F. Smith et Michael S. Waterman en 1981.Comme l'algorithme Needleman – Wunsch, dont il est une variante, 
Smith – Waterman est un algorithme de programmation dynamique.
La procédure de traçage commence à la cellule de la matrice ayant le score le plus élevé et se poursuit jusqu'à ce 
qu'une cellule avec un score zéro soit rencontrée, ce qui donne l'alignement local le plus élevé. En raison de sa 
complexité quadratique dans le temps et dans l'espace, il ne peut souvent pas être appliqué en pratique à des 
problèmes à grande échelle et est remplacé par des alternatives moins générales.
[Source: Wikipedia](https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm)
## Notice d'utilisation <a id="notice"></a>
```
# cloner le repository 
$ git clone https://github.com/tzebre/Needleman_wunch_L3.git

# Se placer dans le repository 
$ cd Needleman_wunch_l3
``` 
- Utilisation du programme en ligne de commande 
```
# Lancer le programme d'alignement avec 
$ python3 NW_main.py
```
- Utilisation du programme avec des fenêtres graphiques
```
# Instaler la librarie tkinter pour python3
$ sudo apt-get install python3-tk

# Lancer le programme d'alignement avec 
$ python3 NW_graphique.py
```
**Note** : Si vous décidez de prendre des séquences stockées dans des fichiers Fasta, les fichiers doivent se trouver dans le
même dossier que le programme. Ou donner leur PATH.

## Fonctionnement <a id="fonctionement"></a>
### Needleman et Wunch <a id="fNE"></a>
L'algorithme prend en entrée : 
  - 2 séquences génomiques ou protéiques.
    - Entrée à la main (seulement en utilisation ligne de commande) ou depuis un fichier Fasta
  - Alignement génomique   
      - Des scores (match, mismatch purine/putine, mismatch pyrimidine/pyrimidine, autre mismatch, 
        ouverture de gap, extention de gap), par défaut (2, 1, 1, -1, -10, -1)
      - Des symboles associés, par défaut ('|',':',':','!','')   
        
  - Alignement protéique 
      - Une matrice de score (Blosum62)
      - Des symboles associés, par défaut (match : '|', mismatch : '!', gap : ' ')
    
Pour l'exemple suivant, on utilise les séquences ci-dessous :
> ATGGCGT  
> ATGAGT  

Et les scores :
> Match : 2  
> Mismatch purine : 1  
> Mismatch pyrimidine : 1  
> Mismatch : -1  
> Ouverture de gap : -10  
> Extention de gap : -1
1) Création des matrices de score et de traceback sous cette forme :  

|   | ㅤ | A | T | G | G | C | G | T |
|---|---|---|---|---|---|---|---|---|  
| ㅤ |   |   |   |   |   |   |   |   |
| A |   |   |   |   |   |   |   |   |
| T |   |   |   |   |   |   |   |   |
| G |   |   |   |   |   |   |   |   |
| A |   |   |   |   |   |   |   |   |  
| G |   |   |   |   |   |   |   |   |
| T |   |   |   |   |   |   |   |   |  

**Note** : Les matrices en Python sont codées comme des listes de listes de liste. 
C'est à dire que chaque ligne est une liste remplie avec une liste par colonne. Chaque colonne 
est remplie par une liste qui représente la case.  
```py 
# tableau de 6 ligne et 5 colone
ligne = 6
col = 5
tableau = []
for x in range(ligne):
    tableau.append([])
    for y in range(col):
        tableau[x].append([])
```
|   |   | A | T |
|---|---|---|---|
|   | 0 |-10|-11|
| A |-10| X | Y |
| T |-11| Z |   |  

Score de X depuis :  
- Dessus : -10 + extention de gap(=-10) = - 11
- Gauche : 1 + extention de gap(=-1) = - 11
- Diagonale : 0 + match(=2) = 2  
  
Ici, le score maximum est 2. Il est obtenu lors d'un déplacement depuis la case en diagonale. 
La case X dans le tableau de score sera donc égale à 2. 
Cette même case dans le tableau de traceback sera remplie avec une flêche `↘` car le "chemin" depuis la case en 
diagonale apporte le meilleur score.  

Score de Y depuis :
- Dessus : -11 + extention de gap(=-1) =  - 12
- Gauche : X(=2) + ouverture de gap(=-10) = - 8 
  - Lors d'un déplacement depuis la case de gauche, on crée 
    un gap. La case de gauche n'étant pas un gap, ici on ouvre donc un nouveau gap. 
- Diagonale : -10 + autre mismatch(=-1) = - 11  

Ici le score maximum est - 8. Il est obtenu lors d'un déplacement depuis la case de gauche. 
La case Y dans le tableau de score sera donc égale a - 8. 
Cette même case dans le tableau de traceback sera remplie avec une flêche `→` car le "chemin" depuis la case 
à gauche apporte le meilleur score. 

**Note** : Dans le cas où deux chemins sont égaux, les deux flêches sont ajoutées. 
(avec `.append()` dans la liste qui correspond à la case) 

Pour trouver le ou les alignements optimaux selon l'algorithme de Needleman et Wunch,  
On remonte le tableau de trace depuis la case en bas à droite en suivant le sens des flêches.

**Note** : Afin de trouver tout les alignements optimaux, on remonte les cases du tableau grâce à de 
la programmation récursive.

### Smith et Waterman <a id="fSW"></a>
L'algorithme fonctionne de la même manière que le précédent. Cependant la case, prend le score maximal entre 
les trois directions possibles et 0. Dans le cas où la case à un score de 0, aucune flêche est inserée.  

Pour trouver le ou les alignements optimaux, on remonte les cases du tableau à partir de ou des cases avec le score maximal.  
On arrête l'alignement quand on trouve un score de 0 (il n'y aura donc pas de flêche).

## Resultats <a id="resultats"></a>
### Utilisation en mode graphique <a id="graph"></a>
Une fenêtre :
- Sur la partie gauche les matrices de score et de trace.
- Sur la partie droite les alignements ex-aequo  

**Note** : Même en utilisation graphique les résultats seront affichés dans la console (Pour être copiable)
### Utilisation en ligne de commande <a id="cmd"></a>
- Matrices de score et de trace 
- Alignement ex-aequo


## Amélioration <a id="amelioration"></a>
- Amélioration de l'utilisation graphique (web?)  
- Alignement en double par moment   
- Enregistrer les résultats dans un fichier   
- Ajout de matrice de substitution protéique
- Meilleur gestion des N comme nucleotide
- Bug non résolu : 
    - Lors de l'utilisation de Smith et Waterman le programme ne retourne que un des alignements 
      par point de depart avec pas la bonne trace d'affiché 


