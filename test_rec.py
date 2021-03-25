mat_max = [[[' '], ['→'], ['→'], ['→'], ['→'], ['→'], ['→'], ['→']],
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
    print(i)
    a,b = nw_aligne(seqA,seqB,i)
    print(a)
    print(b)
