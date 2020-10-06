import numpy as np
from math import *


# def round_expr(expr, num_digits):
#    return expr.xreplace({n : round(n, num_digits) for n in expr.atoms(sy.Number)})

class Bus:
    def __init__(self, p=0.0, q=0.0, v=1, d=0):
        self.p = p
        self.q = q
        self.v = v
        self.d = d

    def __str__(self):
        return str(self.__dict__)


class Line:
    def __init__(self, i: int, j: int, r: float, x: float):
        self.line = {(i, j): complex(r, x)}  # line from i to j has impedance

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
        self.dnr = [i for i in Pnr]
        self.vnr = [i for i in Qnr]
        self.buildYbus()
        print("ITERATION 0:")
        self.PFequations()

    def buildYbus(self):
        n = self.n
        lines = self.lines
        Ybus = np.zeros((n, n), dtype=np.complex_)
        for key in lines:
            i = key[0]
            j = key[1]
            if i != j:
                Ybus[i][j] = 1 / lines[i, j]
                Ybus[j][i] = 1 / lines[i, j]

        for i in range(n):
            Ybus[i][i] = -Ybus[i].sum()
        self.Ybus = Ybus
        self.Ymag = abs(Ybus)
        self.Ytheta = np.angle(Ybus)
        return Ybus

    def PFequations(self):
        buses = self.buses
        Y = self.Ymag
        O = self.Ytheta
        n = self.n
        Pnr = self.Pnr
        Qnr = self.Qnr
        dnr = self.dnr
        vnr = self.vnr
        js = len(Pnr) + len(Qnr)

        V = [buses[i].v for i in range(n)]
        D = [buses[i].d for i in range(n)]

        jacobian = []
        Peq = np.full(n, -1)
        for i in Pnr:
            dProws = []
            Peq[i] = sum(V[i] * V[j] * Y[i][j] * cos(O[i][j] - D[i] + D[j]) for j in range(n))
            for j in dnr:
                dProws.append(sum(V[i] * V[j] * Y[i][j] * sin(O[i][j] - D[i] + D[j]) for j in range(n)))

            for j in vnr:
                dProws.append(sum(V[j] * Y[i][j] * cos(O[i][j] - D[i] + D[j]) for j in range(n)))

            jacobian.append(dProws)

        PQk = [i for i in Peq if i != -1]
        Qeq = np.full(n, -1)
        for i in Qnr:
            dQrows = []
            Qeq[i] = -sum(V[i] * V[j] * Y[i][j] * sin(O[i][j] - D[i] + D[j]) for j in range(n))
            for j in dnr:
                dQrows.append(sum(V[i] * V[j] * Y[i][j] * cos(O[i][j] - D[i] + D[j]) for j in range(n)))

            for j in vnr:
                dQrows.append(-sum(V[j] * Y[i][j] * sin(O[i][j] - D[i] + D[j]) for j in range(n)))

            jacobian.append(dQrows)

        PQk.extend([i for i in Qeq if i != -1])
        self.PQk = np.array(PQk)
        self.Peq = Peq
        self.Qeq = Qeq
        self.jacobian = jacobian

    """
        # Calculate all power flow equations, not just the ones needed
        dPi_dDi=[]
        dPi_dVi=[]
        dQi_dDi=[]
        dQi_dVi=[]
        for i in range(n):
            Pi=sum(V[i]*V[j]*Y[i][j]*cos(O[i][j]-D[i]+D[j]) for j in range(n))
            Qi=-sum(V[i]*V[j]*Y[i][j]*sin(O[i][j]-D[i]+D[j]) for j in range(n))
            P.append(Pi)
            Q.append(Qi)

            dPi_dDi.append(sum(V[i]*V[j]*Y[i][j]*sin(O[i][j]-D[i]+D[j]) for j in range(n)))
            dPi_dVi.append(sum(V[j]*Y[i][j]*cos(O[i][j]-D[i]+D[j]) for j in range(n)))
            dQi_dDi.append(sum(V[i]*V[j]*Y[i][j]*cos(O[i][j]-D[i]+D[j]) for j in range(n)))
            dQi_dVi.append(-sum(V[j]*Y[i][j]*sin(O[i][j]-D[i]+D[j]) for j in range(n)))"""

    def iteration(self):
        PQsch = self.PQsch
        PQk = self.PQk
        buses=self.buses
        deltaPQ = PQsch - PQk
        #DVk = np.linalg.solve(self.jacobian, deltaPQ)
        c = 0
        for i in self.Pnr:
            buses[i].p=PQk[c]
            #buses[i].d = DVk[c]
            c += 1

        for i in self.Qnr:
            buses[i].q=PQk[c]
            #buses[i].v = DVk[c]
            c += 1

        self.buses=buses
        self.print()

    def addLines(self, extra):
        self.lines.update(extra)

    def addBus(self, extra):
        self.buses.update(extra)
        self.n += len(extra)
        self.buildYbus()

    def print(self):
        # print(len(self.lines),"lines",len(self.buses),"buses. Admittance matrix:")
        # print(np.around(self.Ybus,2))
        print("Jacobian:")
        print(np.around(self.jacobian, 2))
        print("Net injections:")
        c = 0
        PQk = self.PQk
        for i in self.Pnr:
            print('P', i, " = ", PQk[c], sep="")
            c += 1

        for i in self.Pnr:
            print('Q', i, " = ", PQk[c], sep="")
            c += 1

        for i in self.buses:
            print("bus",i,self.buses[i])

        """
        print("angles:")
        for idx, val in enumerate(self.buses):
            print("Delta",idx," = ",val.d,sep="")

        print("Voltages:")
        for idx, val in enumerate(self.buses):
            print("V",idx," = ",val.v,sep="")"""
        # self.printPower(self.P,'P')
        # self.printPower(self.Q,'Q')
