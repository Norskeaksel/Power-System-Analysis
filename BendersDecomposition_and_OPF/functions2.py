import numpy as np
def fprint(*args, **kwargs):
    print(*args, **kwargs)
    with open('ResultsAssignment5.txt', 'a') as file:
        print(*args, **kwargs, file=file)

def buildDCY(lines, n):
    Z=np.zeros((n, n))
    Y = np.zeros((n, n))
    for key in lines:
        i = key[0]
        j = key[1]
        if i != j:
            Y[i,j] = 1 / lines[i, j]
            Y[j,i] = 1 / lines[i, j]

    for i in range(n):
        Y[i,i] = -Y[i].sum()

    return Y

def IMML_angles(H,D0,i,k,change):
    """
    :param i: bus nr of line start
    :param k: bus nr of line end
    :param change: factor describing the new admittance at the line
    :return the new voltage angles of the system
    """
    Hinv=np.linalg.inv(H)
    M=np.zeros(shape=(len(H),1))
    M[i,0]=1
    M[k,0]=-1
    delta_h = H[i, k] - H[i, k] * change
    delta_h_inv=1/delta_h
    M_tran=np.transpose(M)
    z=M_tran @ Hinv @ M
    c=1/(delta_h_inv+z)
    #D0 = Hinv @ P
    deltaD=-Hinv @ M * c @ M_tran @ D0
    D=D0+deltaD

    return D

def powerFlows(Y, angles):
    n=Y.shape
    PF=np.zeros(n)
    PF_str=[]
    n=max(n)
    for i in range(n):
        for j in range(n):
            if i != j:
                PF[i, j] = abs(Y[i, j]) * (angles[i] - angles[j])
                if print and i<j and PF[i, j]!=0:
                    PF_str.append(f"P{i}{j} = {PF[i, j]}")

        PF[i,i]=sum(PF[i,:])
    return PF,PF_str

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
                    if k>i:
                        PTDF[row, n]*=-1

    return PTDF

def printPTDF(PTDF,lines):
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

def build_ik(keys,n):
    """
    :return: unsymetrical connection matrix containing the line numbers.
    """
    ik = -np.ones((n, n), dtype=int)
    for nr, (i, k) in enumerate(keys):
        ik[i, k] = nr
    return ik

def printFlowChange(loadChange,lines_str,flowChange):
    fprint('We get the following flow change when the load change is:',loadChange,':')
    for nr,line in enumerate(lines_str):
        fprint(line,':',flowChange[nr])

def updatePF(PF,lines_str,flowChange,n):
    for nr, line in enumerate(lines_str):
        ik=[int(s) for s in line if s.isdigit()]
        i=ik[0]
        k=ik[1]
        PF[i,k]+=flowChange[nr]
        PF[k, i] -= flowChange[nr]
    for i in range(n):
        PF[i,i]=0
        PF[i,i]=PF[i].sum()

    return PF

def pyomoResults(model,F):
    geneation = ["P0:", "P1:", "P2:", "P3:"]
    P = []
    fprint("\nResults:\nObjective value:\n", model.obj())
    fprint("Power generation:")
    for i in model.N:
        P.append(model.P[i].value)
        fprint(geneation[i], P[i])

    L = ["line 1-2: ", "line 2-3: ", "line 0-1: ", "line 0-2: "]
    fprint("Power flow:")
    for i in model.N:
        fprint(L[i], end='')
        a = 0
        for c in model.N:
            a += F[i][c] * model.delta[c].value
        fprint(a)

    return P