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
    print(matrice_traces[row][col], row, col)
    if matrice_traces[row][col][0] not in "↘↓→":
        liste_des_traces.append(trace)
        trace =''
    for symbole in matrice_traces[row][col]:
        if symbole == '↘':
            trace += '↘'
            liste_des_traces = traverse_recursive(matrice_traces,col-1,row-1,liste_des_traces,trace)
            print(symbole)
        if symbole == '↓':
            trace += '↓'
            liste_des_traces = traverse_recursive(matrice_traces,col,row-1,liste_des_traces,trace)
            print(symbole)
        if symbole == '→' :
            trace+='→'
            liste_des_traces = traverse_recursive(matrice_traces,col-1,row,liste_des_traces,trace)
            print(symbole)
    return liste_des_traces

liste_des_traces = []
trace = ""
test = traverse_recursive(mat_max,col,line,liste_des_traces, trace)
print(test)