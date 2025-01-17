import numpy as np
from math import *


def fprint(*args, **kwargs):
    print(*args, **kwargs)
    with open('ResultsAssignment3.txt', 'a') as file:
        print(*args, **kwargs, file=file)

class Bus:
    def __init__(self, p=0.0, q=0.0, v=1, d=0):
        self.p = p
        self.q = q
        self.v = v
        self.d = d

    def __str__(self):
        return str(self.__dict__)


class PowerSystem:
    def __init__(self, lines: dict, buses: dict, slackbus: int, X: list, Pnr: list, Qnr: list, PQsch: np.array):
        self.buses = buses
        self.lines = lines
        self.n = len(buses)
        self.slackbus = slackbus
        self.X = X
        self.Pnr = Pnr
        self.Qnr = Qnr
        self.PQsch = PQsch
        self.deltaPQ = PQsch
        self.buildYbus()
        self.buildB()
        self.PFequations()
        self.S=[]
        self.H=[]
        self.N=[]
        self.M=[]
        self.L=[]
        self.Heq=[]
        self.Leq=[]

        fprint("ITERATION NR: 0")
        fprint("Ybus magnitude:")
        fprint(self.Ymag, '\n')
        fprint("Ybus angles:")
        fprint(self.Ytheta, '\n')
        fprint('Net injections:')
        for i in self.buses:
            fprint("bus", i, self.buses[i])

        fprint()
        self.printMissmatchVector()
        fprint()

    def buildYbus(self):
        n = self.n
        lines = self.lines
        Ybus = np.zeros((n, n), dtype=np.complex_)
        for key in lines:
            i = key[0]
            j = key[1]
            if i != j:
                Ybus[i][j] = -1 / lines[i, j]
                Ybus[j][i] = -1 / lines[i, j]

        for i in range(n):
            Ybus[i][i] = -Ybus[i].sum()
        self.Ybus = Ybus
        self.Ymag = abs(Ybus)
        self.Ytheta = np.angle(Ybus)
        return Ybus

    def buildB(self):
        n = self.n
        lines = self.lines
        B = np.zeros((n, n))
        for key in lines:
            i = key[0]
            j = key[1]
            if i != j:
                B[i][j] = -1 / lines[i, j].imag
                B[j][i] = -1 / lines[i, j].imag

        for i in range(n):
            B[i][i] = -B[i].sum()
        self.B = B
        return B

    def PFequations(self):
        buses = self.buses
        Y = self.Ymag
        O = self.Ytheta
        n = self.n
        Pnr = self.Pnr
        Qnr = self.Qnr
        s = self.slackbus
        V = [buses[i].v for i in range(n)]
        D = [buses[i].d for i in range(n)]
        Peq = np.full(n, -1.)

        for i in Pnr:
            Peq[i] = sum(V[i] * V[j] * Y[i][j] * cos(O[i][j] - D[i] + D[j]) for j in range(n))
            buses[i].p = Peq[i]

        buses[s].p = sum(V[s] * V[j] * Y[s][j] * cos(O[s][j] - D[s] + D[j]) for j in range(n))

        PQk = [i for i in Peq if i != -1]
        Qeq = np.full(n, -1.)
        for i in Qnr:
            Qeq[i] = -sum(V[i] * V[j] * Y[i][j] * sin(O[i][j] - D[i] + D[j]) for j in range(n))
            buses[i].q = Qeq[i]

        buses[s].q = -sum(V[s] * V[j] * Y[s][j] * sin(O[s][j] - D[s] + D[j]) for j in range(n))

        PQk.extend([i for i in Qeq if i != -1])
        self.PQk = PQk

    def buildS(self):
        Snr=list(set(self.Pnr).intersection(self.Qnr))
        Seq = np.full(self.n, 0.j)
        PQk=self.PQk
        Psize=len(self.Pnr)
        for i in Snr:
            Seq[i]=PQk[i]+1j*PQk[i+Psize]

        S = [i for i in Seq if i != 0]
        self.S=S


    def dPi_dDk(self, V, Y, O, D, i, j, k):
        """Equation to devivate on Dk:
        V[i] * V[j] * Y[i][j] * cos(O[i][j] - D[i] + D[j])"""

        if i == j:
            return 0
        if k == i:
            return V[i] * V[j] * Y[i][j] * sin(O[i][j] - D[i] + D[j])
        if k == j:
            return -V[i] * V[j] * Y[i][j] * sin(O[i][j] - D[i] + D[j])
        else:
            return 0

    def dPi_dVk(self, V, Y, O, D, i, j, k):
        """Equation to devivate on Vk:
        V[i] * V[j] * Y[i][j] * cos(O[i][j] - D[i] + D[j])"""

        if i == j and k == i:
            return 2 * V[i] * Y[i][j] * cos(O[i][j] - D[i] + D[j])
        if k == i:
            return V[j] * Y[i][j] * cos(O[i][j] - D[i] + D[j])
        if k == j:
            return V[i] * Y[i][j] * cos(O[i][j] - D[i] + D[j])
        else:
            return 0

    def dQi_dDk(self, V, Y, O, D, i, j, k):
        # Equation to devivate on Dk:
        # V[i] * V[j] * Y[i][j] * sin(O[i][j] - D[i] + D[j])

        if i == j:
            return 0
        if k == i:
            return -V[i] * V[j] * Y[i][j] * cos(O[i][j] - D[i] + D[j])
        if k == j:
            return V[i] * V[j] * Y[i][j] * cos(O[i][j] - D[i] + D[j])
        else:
            return 0

    def dQi_dVk(self, V, Y, O, D, i, j, k):
        """Equation to devivate on Vk:
        V[i] * V[j] * Y[i][j] * sin(O[i][j] - D[i] + D[j])"""

        if i == j and k == i:
            return 2 * V[i] * Y[i][i] * sin(O[i][i])
        if k == i:
            return V[j] * Y[i][j] * sin(O[i][j] - D[i] + D[j])
        if k == j:
            return V[i] * Y[i][j] * sin(O[i][j] - D[i] + D[j])
        else:
            return 0

    def buildJacobian(self):
        buses = self.buses
        Y = self.Ymag
        O = self.Ytheta
        n = self.n
        Pnr = self.Pnr
        Qnr = self.Qnr
        Dnr = Pnr
        Vnr = Qnr
        V = [buses[i].v for i in range(n)]
        D = [buses[i].d for i in range(n)]
        jacobian = []
        #jacobian elements
        Pdim=len(Pnr)
        Qdim=len(Qnr)
        H=np.zeros((Pdim,Pdim))
        N=np.zeros((Pdim,Qdim))
        M=np.zeros((Qdim,Pdim))
        L=np.zeros((Qdim,Qdim))
        Heq=np.zeros((Pdim,Pdim))
        Leq=np.zeros((Qdim,Qdim))
        B=self.B
        for i in Pnr:
            dProws = []
            for k in Dnr:
                dPi_dDi = sum(self.dPi_dDk(V, Y, O, D, i, j, k) for j in range(n))
                dProws.append(dPi_dDi)
                H[i,k]=dPi_dDi
                Heq[i,k]=sum(self.dPi_dDk(V, abs(B), O, D, i, j, k) for j in range(n))

            for k in Qnr:
                dPi_dVi = sum(self.dPi_dVk(V, Y, O, D, i, j, k) for j in range(n))
                dProws.append(dPi_dVi)
                N[i,k]=dPi_dVi

            jacobian.append(dProws)

        for i in Qnr:
            dQrows = []
            for k in Pnr:
                dQi_dDi = -sum(self.dQi_dDk(V, Y, O, D, i, j, k) for j in range(n))
                dQrows.append(dQi_dDi)
                M[i,k]=dQi_dDi

            for k in Qnr:
                dQi_dVi = -sum(self.dQi_dVk(V, Y, O, D, i, j, k) for j in range(n))
                dQrows.append(dQi_dVi)
                L[i,k]=dQi_dVi
                Leq[i, k]=-sum(self.dQi_dVk(V, abs(B), O, D, i, j, k) for j in range(n))

            jacobian.append(dQrows)

        self.jacobian = np.array(jacobian)
        self.H=H
        self.N=N
        self.M=M
        self.L=L
        self.Heq=Heq
        self.Leq=Leq
        Pdim=H.shape
        Qdim=L.shape
        Pzeros=np.zeros(Pdim)
        Qzeros=np.zeros(Qdim)
        Ppart=np.hstack([H, Pzeros])
        Qpart=np.hstack([Qzeros,L])
        self.DPFjacobian=np.vstack([Ppart,Qpart])

    def printMissmatchVector(self):
        fprint('Missmatch Vector:')
        c = 0
        deltaPQ = self.PQsch - self.PQk
        for i in self.Pnr:
            fprint('Delta P', i, " = ", deltaPQ[c], sep="")
            c += 1

        for i in self.Qnr:
            fprint('Delta Q', i, " = ", deltaPQ[c], sep="")
            c += 1

        fprint()

    def printCorrectionVector(self):
        fprint('Correction Vector:')
        DVk = self.DVk
        c = 0
        for i in self.Pnr:
            fprint('Delta D', i, " = ", DVk[c], sep="")
            c += 1

        for i in self.Pnr:
            fprint('Delta V', i, " = ", DVk[c], sep="")
            c += 1

        fprint()

    def print(self, itNr):
        # fprint(len(self.lines),"lines",len(self.buses),"buses. Admittance matrix:")
        # fprint(np.around(self.Ybus,2))
        fprint("ITERATION NR:", itNr)
        fprint('Net injections:')
        for i in self.buses:
            fprint("bus", i, self.buses[i])

        fprint()
        self.printMissmatchVector()
        self.printCorrectionVector()
        fprint()

    def iteration(self, itNr):
        self.buildJacobian()

        PQsch = self.PQsch
        PQk = self.PQk
        buses = self.buses
        n = self.n
        s = self.slackbus
        deltaPQ = PQsch - PQk
        DVk = np.linalg.solve(self.jacobian, deltaPQ)
        c = 0
        for i in self.Pnr:
            # buses[i].p = PQk[c]
            buses[i].d += DVk[c]
            c += 1

        for i in self.Qnr:
            # buses[i].q = PQk[c]
            buses[i].v += DVk[c]
            c += 1

        V = [buses[i].v for i in range(n)]
        D = [buses[i].d for i in range(n)]

        self.DVk = DVk
        self.buses = buses
        self.PFequations()

    def buildPredictionVector(self):
        jacobian=self.jacobian
        zero=np.zeros(len(jacobian))
        zero[-1]=1
        self.predictionVector=np.linalg.solve(jacobian, zero)
        return self.predictionVector

    def takePredictionStep(self,ba,step):
        predictionVector=self.predictionVector
        self.PQsch-=step*ba
        halfWay=len(predictionVector)//2
        addDeltas=step*predictionVector[0:halfWay]
        addVoltages=step*predictionVector[halfWay:]
        for i in range(halfWay):
            self.buses[i].d+=addDeltas[i]
            self.buses[i].v+=addVoltages[i]

        self.PFequations() #Recalculate PF equations with new V and D

    def getV(self):
        return [self.buses[i].v for i in range(self.n) if i!=self.slackbus]

    def getD(self):
        return [self.buses[i].d for i in range(self.n) if i!=self.slackbus]

    def setV(self,V):
        for idx,val in enumerate(V):
            self.buses[idx].v=val

    def setD(self,D):
        for idx,val in enumerate(D):
            self.buses[idx].d=val