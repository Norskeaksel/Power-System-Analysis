import numpy as np
def fprint(*args, **kwargs):
    print(*args, **kwargs)
    with open('ResultsAssignment4.txt', 'a') as file:
        print(*args, **kwargs, file=file)

def build_ik(keys,n):
    ik = -np.ones((n, n), dtype=int)
    for nr, (i, k) in enumerate(keys):
        ik[i, k] = nr
    return ik

def buildDCY(lines, n):
    Z=np.zeros((n, n))
    Y = np.zeros((n, n))
    for key in lines:
        i = key[0]
        j = key[1]
        if i != j:
            #Z[i,j] = lines[i, j]
            #Z[j,i] = lines[i, j]
            Y[i,j] = 1 / lines[i, j]
            Y[j,i] = 1 / lines[i, j]

    for i in range(n):
        #Z[i, i] = -Z[i].sum()
        Y[i,i] = -Y[i].sum()

    return Y#,Z

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

def build_b(Z,keys,n):
    b=np.zeros((n,n))
    for row,(i,k) in enumerate(keys):
        b[row,i]=(1/Z[i,k])

    b=np.maximum(b, b.transpose()) #make b symetrical
    return b

def buildPTDF(Z,B,ik,n):
    """
    The PTDF value from node n for the line between nodes i and k can then be calculated with PTDF_ik,n= B_ik(Zbus_in-Zbus_kn)
    """
    nrOfbuses=n
    nrOfLines=len(ik)
    PTDF = np.zeros((nrOfLines,nrOfbuses))
    for i in range(nrOfLines):
        for k in range(nrOfLines):
            row = ik[i, k]
            if row!=-1:
                for n in range(nrOfbuses):
                        PTDF[row, n] = B[i, k] * (Z[i, n] - Z[k, n])
    return PTDF

def printPTDF(PTDF,keys):
    lines=[]
    for i,k in keys:
        lines.append(f"line {i}-{k} ")
    fprint('PTDF MATRIX:')
    first=1
    for nodeNr in range(PTDF.shape[1]):
        if first:
            fprint('         ',end="")
            first=0
        fprint('node', nodeNr,end='      ')

    fprint()
    for lineNr in range(len(lines)):
        fprint(lines[lineNr], end="")
        for node in range(PTDF.shape[1]):
            val=PTDF[lineNr, node]
            space = "    "
            if val>=0:
                space = "     "
            fprint("%4.4f" %(PTDF[lineNr, node]),space, end="")
        fprint()

def newPowerFlows(PF,PTDF,newP):
