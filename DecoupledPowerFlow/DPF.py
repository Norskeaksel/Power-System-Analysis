import importlib
from copy import deepcopy
from DecoupledPowerFlow import newtonRapson3
importlib.reload(newtonRapson3)
from DecoupledPowerFlow.newtonRapson3 import *

def DPF(lines, X, PQsch, P, Q, V, D, Pnr, Qnr, slackbus, allowedMissmatch,taskNr,maxIteratins):
    buses = buildBuses(P, Q, V, D)
    if taskNr==2:
        fprint("Task 2 initialize system:\n")
    elif taskNr==3:
        fprint("Task 3 initialize system:\n")
    elif taskNr==4:
        fprint("Task 4 initialize system:\n")
    DPS = PowerSystem(lines, buses, slackbus, X, Pnr, Qnr, PQsch)
    DPS.buildJacobian()
    if taskNr==2:
        fprint("Task 2a., primal FDPF solution:\n")
    elif taskNr==3:
        fprint("Task 3, primal FDPF solution:\n")
    elif taskNr==4:
        fprint("Task 4 primal FDPF solution:\n")
    DPS2,c2=primalDPF(deepcopy(DPS), allowedMissmatch,maxIteratins)
    if taskNr==2:
        fprint("Task 2b., dual FDPF solution:\n")
    elif taskNr==3:
        fprint("Task 3, dual FDPF solution:\n")
    elif taskNr==4:
        fprint("Task 4, dual FDPF solution:\n")
    DPS3,c3=dualDPF(deepcopy(DPS),allowedMissmatch,maxIteratins)
    c4="NaN"
    DPS4=None
    if taskNr==2:
        fprint("Task 2c., standard FDPF solution:\n")
        DPS4,c4=standardDPF(deepcopy(DPS),allowedMissmatch,maxIteratins)

    return [DPS2,DPS3,DPS4],[c2,c3,c4]

def primalDPF(DPS, allowedMissmatch,maxIteratins):
    c=0
    oldDeviation = 1e5
    while 1 and c<maxIteratins:
        c+=1
        PQsch = DPS.PQsch
        PQk = DPS.PQk
        deltaPQ = PQsch - PQk
        deltaP=deltaPQ[:len(deltaPQ)//2]

        Hinv=np.linalg.inv(DPS.H)
        deltaD=np.matmul(Hinv,deltaP)
        newD=DPS.getD()+deltaD
        DPS.setD(newD)
        LeqInv=np.linalg.inv(DPS.Leq)

        #Update equations to find new delta Q
        DPS.PFequations()
        PQk = DPS.PQk
        deltaPQ = PQsch - PQk
        deltaQ = deltaPQ[len(deltaPQ)//2:]
        deltaV=np.matmul(LeqInv,deltaQ)
        newV=DPS.getV()+deltaV
        DPS.setV(newV)
        DPS.PFequations()
        maxDeviation = max(abs(DPS.PQsch - DPS.PQk))

        DPS.DVk=np.append(deltaD,deltaV)
        DPS.print(c)
        if maxDeviation < allowedMissmatch:
            break
        elif oldDeviation < maxDeviation:
            fprint("Algorithm not converging with following load:\n", DPS.PQsch)
            c = -1
            break
        oldDeviation = maxDeviation+0.1
    return DPS,c

def dualDPF(DPS, allowedMissmatch,maxIteratins):
    c=0
    oldDeviation=1e5
    while 1 and c<maxIteratins:
        c+=1
        PQsch = DPS.PQsch
        PQk = DPS.PQk
        deltaPQ = PQsch - PQk
        deltaQ = deltaPQ[len(deltaPQ) // 2:]

        Linv=np.linalg.inv(DPS.L)
        deltaV=np.matmul(Linv,deltaQ)
        newV=DPS.getV()+deltaV
        DPS.setV(newV)
        HeqInv = np.linalg.inv(DPS.Heq)

        # Update equations to find new delta P
        DPS.PFequations()
        PQk = DPS.PQk
        deltaPQ = PQsch - PQk
        deltaP=deltaPQ[:len(deltaPQ)//2]
        deltaD = np.matmul(HeqInv, deltaP)
        newD = DPS.getD() + deltaD
        DPS.setD(newD)
        DPS.PFequations()
        maxDeviation = max(abs(DPS.PQsch - DPS.PQk))

        DPS.DVk = np.append(deltaD, deltaV)
        DPS.print(c)
        if maxDeviation < allowedMissmatch:
            break
        elif oldDeviation < maxDeviation:
            fprint("Algorithm not converging with following load:\n", DPS.PQsch)
            c = -1
            break
        oldDeviation = maxDeviation+0.1

    return DPS,c

def standardDPF(DPS, allowedMissmatch,maxIteratins):
    c = 0
    oldDeviation=1e5
    while 1 and c < maxIteratins:
        c += 1
        PQsch = DPS.PQsch
        PQk = DPS.PQk
        deltaPQ = PQsch - PQk
        deltaP = deltaPQ[:len(deltaPQ) // 2]
        deltaQ = deltaPQ[len(deltaPQ) // 2:]

        Hinv = np.linalg.inv(DPS.H)
        Linv = np.linalg.inv(DPS.L)
        deltaD = np.matmul(Hinv, deltaP)
        newD = DPS.getD() + deltaD
        DPS.setD(newD)

        Linv=np.linalg.inv(DPS.L)
        deltaV=np.matmul(Linv,deltaQ)
        newV=DPS.getV()+deltaV
        DPS.setV(newV)

        DPS.PFequations()
        maxDeviation = max(abs(DPS.PQsch - DPS.PQk))
        DPS.DVk = np.append(deltaD, deltaV)
        DPS.print(c)
        if maxDeviation < allowedMissmatch:
            break
        elif oldDeviation<maxDeviation:
            fprint("Algorithm not converging with following load:\n",DPS.PQsch)
            c=-1
            break
        oldDeviation = maxDeviation+0.1
    return DPS,c
