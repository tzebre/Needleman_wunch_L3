'''mat_max = [[[' '], ['→'], ['→'], ['→'], ['→'], ['→'], ['→'], ['→']],
           [['↓'], ['↘'], ['→'], ['→'], ['→'], ['→', '↘'], ['→'], ['→']],
           [['↓'], ['↓'], ['↘'], ['↘'], ['→'], ['→'], ['→'], ['→']],
           [['↓'], ['↓'], ['↘', '↓'], ['↓'], ['↘'], ['→'], ['→'], ['→', '↘']],
           [['↓'], ['↓'], ['↘', '↓'], ['↓'], ['↘', '↓'], ['↘'], ['→', '↘'], ['↘']],
           [['↓'], ['↓'], ['↘', '↓'], ['↘'], ['↓'], ['↘', '↓'], ['↘'], ['↓']],
           [['↓'], ['↓'], ['↘'], ['↓'], ['↓'], ['↘', '↓'], ['↘'], ['→']],
           [['↓'], ['↓'], ['↓'], ['↘'], ['→', '↓'], ['↘', '↓'], ['↓'], ['↘']]]

col = 7
line = 7

def traverse_recursive(matrice_traces, col, row, liste_des_traces, trace):
    if matrice_traces[row][col][0] not in "↘↓→":
        if len(liste_des_traces) >= 1:
            liste_des_traces.append(trace[1:])
        else:
            liste_des_traces.append(trace)
        trace =''
    for symbole in matrice_traces[row][col]:
        if symbole == '↘':
            trace += '↘'
            liste_des_traces = traverse_recursive(matrice_traces,col-1,row-1,liste_des_traces,trace)
        if symbole == '↓':
            trace += '↓'
            liste_des_traces = traverse_recursive(matrice_traces,col,row-1,liste_des_traces,trace)
        if symbole == '→' :
            trace+='→'
            liste_des_traces = traverse_recursive(matrice_traces,col-1,row,liste_des_traces,trace)
    return liste_des_traces

liste_des_traces = []
trace = ""
test = traverse_recursive(mat_max,col,line,liste_des_traces, trace)
print(test)

def nw_aligne(seqA, seqB, trace):
    lenA = len(seqA)
    lenB = len(seqB)
    i = lenB
    j = lenA
    t = 0
    align1 = ""
    align2 = ""
    while i > 0 and j > 0:
        s_current = trace[t]
        if s_current == '↘':
            align1 += seqA[j - 1]
            align2 += seqB[i - 1]
            i -= 1
            j -= 1
        elif s_current == '↓':
            align1 += '-'
            align2 += seqB[i - 1]
            i -= 1
        elif s_current == '→':
            align1 += seqA[j - 1]
            align2 += '-'
            j -= 1
        else:
            print("erreur alignement")
            break
        t += 1
    while j > 0:
        align1 += seqA[j - 1]
        align2 += '-'
        j -= 1
        t+=1
    while i > 0:
        align1 += '-'
        align2 += seqB[i - 1]
        i -= 1
        t+=1
    align1 = align1[::-1]
    align2 = align2[::-1]
    return align1, align2
seqA = "gcatgcu"
seqB = "gattaca"
for i in test:
    print(i[::-1])
    a,b = nw_aligne(seqA,seqB,i)
    print(a)
    print(b)

print(len(seqA))'''

def liste_depart_max(mat_score,seqA,seqB):
    depart_max = []
    max_score = 0
    i = 0
    j = 0
    while i < len(seqB)+1:
        while j < len(seqA)+1:
            if mat_score[i][j][0] == max_score:
                depart_max.append([])
                x = len(depart_max)
                depart_max[x-1].append(i)
                depart_max[x-1].append(j)
            if mat_score[i][j][0] > max_score:
                x = len(depart_max)
                max_score = mat_score[i][j]
                depart_max[x-1].append(i)
                depart_max[x-1].append(j)
            j += 1
        i += 1
    return depart_max
seqA = "gcatgcu"
seqB = "gattaca"
score_mat = [[[0, 0], [-1, 1], [-2, 1], [-3, 1], [-4, 1], [-5, 1], [-6, 1], [-7, 1]], [[-1, 1], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0]], [[-2, 1], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0]], [[-3, 1], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0]], [[-4, 1], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0]], [[-5, 1], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0]], [[-6, 1], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0]], [[-7, 1], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0]]]
test = liste_depart_max(score_mat,seqA,seqB)
print(test)