import importlib
from copy import deepcopy
from DecoupledPowerFlow import newtonRapson3
importlib.reload(newtonRapson3)
from DecoupledPowerFlow.newtonRapson3 import *

def cout(buses):
    for i in buses:
        fprintResults("bus", i, buses[i])
    fprintResults()

def DPF(lines, X, PQsch, P, Q, V, D, Pnr, Qnr, slackbus, allowedMissmatch):
    buses = buildBuses(P, Q, V, D)
    fprint("Task 2 initialize system:\n")
    DPS = PowerSystem(lines, buses, slackbus, X, Pnr, Qnr, PQsch)
    DPS.buildJacobian()
    fprint("Task 2a., primal FDPF solution:\n")
    DPS2=primalDPF(deepcopy(DPS), allowedMissmatch)

    fprint("Task 2b., dual FDPF solution:\n")
    DPS3=dualDPF(deepcopy(DPS),allowedMissmatch)

    fprint("Task 2c., standard FDPF solution:\n")
    DPS4=standardDPF(deepcopy(DPS),allowedMissmatch)
    return DPS2,DPS3

def primalDPF(DPS, allowedMissmatch):
    c=0
    while 1 and c<10:
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

    return DPS

def dualDPF(DPS, allowedMissmatch):
    c=0
    while 1 and c<10:
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

    return DPS

def standardDPF(DPS, allowedMissmatch):
    c = 0
    while 1 and c < 10:
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