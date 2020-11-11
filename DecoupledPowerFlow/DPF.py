import importlib
from copy import deepcopy
from DecoupledPowerFlow import newtonRapson3

importlib.reload(newtonRapson3)
from DecoupledPowerFlow.newtonRapson3 import *


def DPF(lines, X, PQsch, buses, Pnr, Qnr, slackbus, allowedMissmatch, taskNr, maxIteratins):
    fprint(f"Task {taskNr} initialize system:\n")
    DPS = PowerSystem(lines, buses, slackbus, X, Pnr, Qnr, PQsch)
    DPS.buildJacobian()
    if taskNr == 2:
        fprint("Task 2a., primal FDPF solution:\n")
    else:
        fprint(f"Task {taskNr}, primal FDPF solution:\n")

    DPS2, c2 = primalDPF(deepcopy(DPS), allowedMissmatch, maxIteratins)
    if taskNr == 2:
        fprint("Task 2b., dual FDPF solution:\n")
    else:
        fprint(f"Task {taskNr}, dual FDPF solution:\n")

    DPS3, c3 = dualDPF(deepcopy(DPS), allowedMissmatch, maxIteratins)
    c4 = "NaN"
    DPS4 = None
    if taskNr == 2:
        fprint("Task 2c., standard FDPF solution:\n")
        DPS4, c4 = standardDPF(deepcopy(DPS), allowedMissmatch, maxIteratins)

    return [DPS2, DPS3, DPS4], [c2, c3, c4]


def primalDPF(DPS, allowedMissmatch, maxIteratins):
    c = 0
    oldDeviation = 1e5
    while 1 and c < maxIteratins:
        c += 1
        PQsch = DPS.PQsch
        PQk = DPS.PQk
        deltaPQ = PQsch - PQk
        deltaP = deltaPQ[:len(deltaPQ) // 2]

        #Find new angles
        Hinv = np.linalg.inv(DPS.H)
        deltaD = np.matmul(Hinv, deltaP)
        newD = DPS.getD() + deltaD
        DPS.setD(newD)
        LeqInv = np.linalg.inv(DPS.Leq)

        # Update equations to find new delta Q
        DPS.PFequations()
        PQk = DPS.PQk
        deltaPQ = PQsch - PQk
        deltaQ = deltaPQ[len(deltaPQ) // 2:]

        #Find new voltages
        deltaV = np.matmul(LeqInv, deltaQ)
        newV = DPS.getV() + deltaV
        DPS.setV(newV)

        #Calculate power injections
        DPS.PFequations()
        maxDeviation = max(abs(DPS.PQsch - DPS.PQk))

        #record correction vector
        DPS.DVk = np.append(deltaD, deltaV)
        DPS.print(c)
        if maxDeviation < allowedMissmatch:
            break
        elif oldDeviation < maxDeviation:
            fprint("Algorithm not converging with following load:\n", DPS.PQsch)
            c = -1
            break
        oldDeviation = maxDeviation + 0.1 # stop algorithm if improvement is less than 0.1
    return DPS, c


def dualDPF(DPS, allowedMissmatch, maxIteratins):
    c = 0
    oldDeviation = 1e5
    while 1 and c < maxIteratins:
        c += 1
        PQsch = DPS.PQsch
        PQk = DPS.PQk
        deltaPQ = PQsch - PQk
        deltaQ = deltaPQ[len(deltaPQ) // 2:]

        # Find new voltages
        Linv = np.linalg.inv(DPS.L)
        deltaV = np.matmul(Linv, deltaQ)
        newV = DPS.getV() + deltaV
        DPS.setV(newV)
        HeqInv = np.linalg.inv(DPS.Heq)

        # Update equations to find new delta P
        DPS.PFequations()
        PQk = DPS.PQk
        deltaPQ = PQsch - PQk
        deltaP = deltaPQ[:len(deltaPQ) // 2]

        # Find new angles
        deltaD = np.matmul(HeqInv, deltaP)
        newD = DPS.getD() + deltaD
        DPS.setD(newD)

        DPS.PFequations()
        maxDeviation = max(abs(DPS.PQsch - DPS.PQk))

        # record correction vector
        DPS.DVk = np.append(deltaD, deltaV)
        DPS.print(c)
        if maxDeviation < allowedMissmatch:
            break
        elif oldDeviation < maxDeviation:
            fprint("Algorithm not converging with following load:\n", DPS.PQsch)
            c = -1
            break
        oldDeviation = maxDeviation + 0.1 # stop algorithm if improvement is less than 0.1

    return DPS, c


def standardDPF(DPS, allowedMissmatch, maxIteratins):
    c = 0
    oldDeviation = 1e5
    while 1 and c < maxIteratins:
        c += 1
        PQsch = DPS.PQsch
        PQk = DPS.PQk
        deltaPQ = PQsch - PQk
        deltaP = deltaPQ[:len(deltaPQ) // 2]
        deltaQ = deltaPQ[len(deltaPQ) // 2:]

        #Find new angles
        Hinv = np.linalg.inv(DPS.H)
        Linv = np.linalg.inv(DPS.L)
        deltaD = np.matmul(Hinv, deltaP)
        newD = DPS.getD() + deltaD
        DPS.setD(newD)

        #Find new voltages without using the newfound angles
        Linv = np.linalg.inv(DPS.L)
        deltaV = np.matmul(Linv, deltaQ)
        newV = DPS.getV() + deltaV
        DPS.setV(newV)

        # Update equations to find new delta P
        DPS.PFequations()
        maxDeviation = max(abs(DPS.PQsch - DPS.PQk))

        # record correction vector
        DPS.DVk = np.append(deltaD, deltaV)
        DPS.print(c)
        if maxDeviation < allowedMissmatch:
            break
        elif oldDeviation < maxDeviation:
            fprint("Algorithm not converging with following load:\n", DPS.PQsch)
            c = -1
            break
        oldDeviation = maxDeviation + 0.1 # stop algorithm if improvement is less than 0.1
    return DPS, c
