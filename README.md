# Needleman_Wunch/Smith_Waterman
**Auteur : MATHIEU Theo**  
Programe realisé en python dans le cadre de la L3 BISM pour l'UE bioinformatique.

## Fonctionalité  
 - [x] Alignement selon l'algotihme de Needleman-wunch ou Smith-Waterman 
    - De sequence proteique ou genomique 
    - Depuis un fichier fasta ou un entrée manuelle de la séquence 
    - Choix des score et des symbole d'alignement possible  
    
## Principe d'alignement 
### Needleman et wunch 
>L'algorithme Needleman – Wunsch est un algorithme utilisé en bioinformatique pour aligner des séquences protéiques ou 
nucléotidiques. C'était l'une des premières applications de la programmation dynamique pour comparer des séquences 
biologiques. L'algorithme a été développé par Saul B. Needleman et Christian D. Wunsch et publié en 1970.  Il est
également parfois appelé algorithme de correspondance optimal et technique d'alignement global.L'algorithme attribue un score à chaque alignement
possible, et le but de l'algorithme est de trouver tous les alignements possibles ayant le score le plus élevé.
[Source: Wikipedia](https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm)
### Smith et Waterman 
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
## Notice d'utilisation
```
# cloner le repository 
$ git clone https://github.com/tzebre/Needleman_wunch_L3.git

# Se placer dans le repository 
$ cd Needleman_wunch_l3

# Lancer le programme d'alignement avec 
$ python3 NW_main.py
```
**Note** : Si vous decidez de prendre des sequence stockées dans des fichiers fasta. Les fichiers doivent se trouver dans le
meme dossier que le programme.

## Fonctionement 
### Needleman et wunch 
L'algorithme prend en entrée: 
  - 2 sequence genomique ou proteique.
    - Entrée a la main ou depuis un fichier fasta
  - Alignement genomique   
      - Des score (match, mismatch purine/putine, mismatch pyrimidine/pyrimidine, autre mismatch, 
        ouverture de gap, extention de gap), par defaut (2, 1, 1, -1, -10, -1)
      - Ainsi que les symboles associés, par defaut ('|',':',':','!','')   
        
  - Alignement proteique 
      - Une matrice de score (Blosum62)
      - Ainsi que des symbole associé, par defaut (match : '|', mismatch : '!', gap : ' ')
    
Pour l'exemple suivant on utilise les sequence suivante :
> ATGGCGT  
> ATGAGT  

Et le score par defaut :
1) Creation des matrice de score et de traceback sous cette forme 

|   | ㅤ | A | T | G | G | C | G | T |
|---|---|---|---|---|---|---|---|---|  
|ㅤ  |   |   |   |   |   |   |   |   |
| A |   |   |   |   |   |   |   |   |
| T |   |   |   |   |   |   |   |   |
| G |   |   |   |   |   |   |   |   |
| A |   |   |   |   |   |   |   |   |  
| G |   |   |   |   |   |   |   |   |
| T |   |   |   |   |   |   |   |   |  

**Note** : les matrices en python son codée comme des listes de listes de liste. 
C'est a dire que chaque ligne est une liste replie avec une liste par colone, et chaque colone 
est remplie par une liste qui represente la case.  
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
- dessu : -10 + extention de gap(=-10) = - 11
- gauche : 1 + extention de gap(=-1) = - 11
- diagonale : 0 + match(=2) = 2  
  
Ici le score maximum est 2 et il est obtenu lors d'un deplacement depuis la case en diagonale. 
la case X dans le tableau de score sera donc égale a 2. 
Cette même case dans le tableau de traceback sera rempli avec une fleche `↘`, car le "chemin" depuis la case en 
diagonale apporte le meilleur score.  

Score de Y depuis :
- dessu : -11 + extention de gap(=-1) =  - 12
- gauche : X(=2) + ouverture de gap(=-10) = - 8 
  - Lors d'un deplacement depuis la case de gauche on cree 
    un gap. La case gauche n'etant pas un gap, ici on ouvre donc un nouveau gap. 
- diagonale : -10 + autre mismatch(=-1) = - 11  

Ici le score maximum est - 8 et il est obtenu lors d'un deplacement depuis la case de gauche. 
la case y dans le tableau de score sera donc égale a - 8. 
Cette même case dans le tableau de traceback sera rempli avec une fleche `→`, car le "chemin" depuis la case 
a gauche apporte le meilleur score. 

**Note** : Dans le cas ou deux chemin sont egaux les deux fleche sont ajouté. 
(avec `.append()` dans la liste qui correspond a la case) 

### Smith et Waterman 







    

    


