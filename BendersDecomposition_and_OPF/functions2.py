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