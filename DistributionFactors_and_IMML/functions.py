import numpy as np
def fprint(*args, **kwargs):
    print(*args, **kwargs)
    with open('ResultsAssignment4.txt', 'a') as file:
        print(*args, **kwargs, file=file)

def build_ik(keys,n):
    """
    :return: unsymetrical connection matrix containing the line numbers.
    """
    ik = -np.ones((n, n), dtype=int)
    for nr, (i, k) in enumerate(keys):
        ik[i, k] = nr
    return ik

def buildY(lines, n):
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

def powerFlows(Y, angles,P):
    n=Y.shape
    PF=np.zeros(n)
    PF_str=[]
    n=max(n)
    for i in range(n):
        for j in range(n):
            if i==j:
                PF[i,j]=P[i]
            else:
                PF[i, j] = abs(Y[i, j]) * (angles[i] - angles[j])
                if print and i<j and PF[i, j]!=0:
                    PF_str.append(f"P{i}{j} = {PF[i, j]}")
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

def newPowerFlows(PTDF,loadChange):
    flowChange=PTDF*loadChange
    print(flowChange)

def printFlowChange(loadChange,lines_str,flowChange):
    fprint('We get the following flow change when the load change is:',loadChange,':')
    for nr,line in enumerate(lines_str):
        fprint(line,':',flowChange[nr])


def IMML_angles(H,P,i,k,change):
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
    D0 = Hinv @ P
    deltaD=-Hinv @ M * c @ M_tran @ D0
    D=D0+deltaD

    return D