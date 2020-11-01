import numpy as np
def fprint(*args, **kwargs):
    print(*args, **kwargs)
    with open('ResultsAssignment4.txt', 'a') as file:
        print(*args, **kwargs, file=file)

def buildDCY(lines,n):
    Y = np.zeros((n, n))
    for key in lines:
        i = key[0]
        j = key[1]
        if i != j:
            Y[i][j] = 1 / lines[i, j]
            Y[j][i] = 1 / lines[i, j]

    for i in range(n):
        Y[i][i] = -Y[i].sum()

    return Y

def powerFlows(Y, angles,P, print=False):
    n=Y.shape
    PF=np.zeros(n)
    n=max(n)
    for i in range(n):
        for j in range(n):
            if i==j:
                PF[i,j]=P[i]
            else:
                PF[i, j] = abs(Y[i, j]) * (angles[i] - angles[j])
                if print and i>j and PF[i, j]!=0:
                    fprint('P',i+1,j+1,' = ',PF[i, j],sep="")
    return PF

def buildPTDF(Z,b,ik):
    nrOfbuses=len(Z)
    nrOfLines=len(ik)
    PTDF = np.zeros((nrOfLines,nrOfbuses))
    for i in range(nrOfLines):
        for k in range(nrOfLines):
            row = ik(i, k)
            if row!=0:
                for n in range(nrOfbuses):
                    PTDF[row, n] = b(i, k) * (Z[i, n] - Z[k, n])
    return PTDF

def printPTDF(PTDF,ik):
    lines=[str(line[0])+'-'+str(line[1]) for line in ik]
    fprint('PTDF MATRIX:')
    for line in range(PTDF.shape[1]):
        if line == 0:
            fprint('         ',end="")
        else:
            fprint('node', line,'    ',end="")

    fprint()
    for line in range(PTDF.shape[0]):
        fprint('line', lines[line], end="")
        for node in PTDF.shape[1]:
            fprint(PTDF[line, node])
        fprint()