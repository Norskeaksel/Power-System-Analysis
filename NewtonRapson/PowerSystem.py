#Test commit from pycharm
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
        self.buildYbus()
        print("ITERATION 0:")
        self.PFequations()
        self.buildJacobian()
        self.numeric()

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

    def dPi_dDk(self,V,Y,O,D,i,j,k):
        """Equation to devivate on Dk:
        V[i] * V[j] * Y[i][j] * cos(O[i][j] - D[i] + D[j])"""

        if i==j:
            return 0
        if k==i:
            return V[i] * V[j] * Y[i][j] * sin(O[i][j] - D[i] + D[j])
        if k==j:
            return -V[i] * V[j] * Y[i][j] * sin(O[i][j] - D[i] + D[j])
        else:
            return 0

    def dPi_dVk(self,V,Y,O,D,i,j,k):
        """Equation to devivate on Vk:
        V[i] * V[j] * Y[i][j] * cos(O[i][j] - D[i] + D[j])"""

        if i==j and k==i:
            return 2*V[i]* Y[i][j] * cos(O[i][j] - D[i] + D[j])
        if k==i:
            return V[j] * Y[i][j] * cos(O[i][j] - D[i] + D[j])
        if k==j:
            return V[i] * Y[i][j] * cos(O[i][j] - D[i] + D[j])
        else:
            return 0

    def dQi_dDk(self,V,Y,O,D,i,j,k):
        # Equation to devivate on Dk:
        # V[i] * V[j] * Y[i][j] * sin(O[i][j] - D[i] + D[j])

        if i==j:
            return 0
        if k==i:
            return -V[i] * V[j] * Y[i][j] * cos(O[i][j] - D[i] + D[j])
        if k==j:
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

    def PFequations(self):
        buses = self.buses
        Y = self.Ymag
        O = self.Ytheta
        n = self.n
        Pnr = self.Pnr
        Qnr = self.Qnr
        V = [buses[i].v for i in range(n)]
        D = [buses[i].d for i in range(n)]
        Peq = np.full(n, -1.)

        for i in Pnr:
            Peq[i] = sum(V[i] * V[j] * Y[i][j] * cos(O[i][j] - D[i] + D[j]) for j in range(n))

        PQk = [i for i in Peq if i != -1]
        Qeq = np.full(n, -1.)
        for i in Qnr:
            Qeq[i] = -sum(V[i] * V[j] * Y[i][j] * sin(O[i][j] - D[i] + D[j]) for j in range(n))

        PQk.extend([i for i in Qeq if i != -1])

        self.PQk=PQk

    def buildJacobian(self):
        buses = self.buses
        Y = self.Ymag
        O = self.Ytheta
        n = self.n
        Pnr = self.Pnr
        Qnr = self.Qnr
        Dnr=Pnr
        Vnr=Qnr
        V = [buses[i].v for i in range(n)]
        D = [buses[i].d for i in range(n)]

        jacobian = []
        for i in Pnr:
            dProws = []
            for k in Dnr:
                dPi_dDi=sum(self.dPi_dDk(V, Y, O, D, i, j, k) for j in range(n))
                dProws.append(dPi_dDi)

            for k in Qnr:
                dPi_dVi=sum(self.dPi_dVk(V, Y, O, D, i, j, k) for j in range(n))
                dProws.append(dPi_dVi)

            jacobian.append(dProws)

        for i in Qnr:
            dQrows = []
            for k in Pnr:
                dQi_dDi=-sum(self.dQi_dDk(V, Y, O, D, i, j, k) for j in range(n))
                dQrows.append(dQi_dDi)

            for k in Qnr:
                dQi_dVi=-sum(self.dQi_dVk(V, Y, O, D, i, j, k) for j in range(n))
                dQrows.append(dQi_dVi)

            jacobian.append(dQrows)

        self.jacobian=np.array(jacobian)

    def numeric(self):
        from copy import deepcopy
        def f(i, eps, d_or_v):
            backup = deepcopy(self.buses)
            if d_or_v == 'v':
                self.buses[i].v += eps
            else:
                self.buses[i].d += eps
            self.PFequations()
            self.buses = deepcopy(backup)
            return np.array(self.PQk)

        eps = 1e-5
        jacobian = []
        for i in self.Pnr:
            jacobian.append((f(i, eps / 2, 'd') - f(i, -eps / 2, 'd')) / eps)
        for i in self.Qnr:
            jacobian.append((f(i, eps / 2, 'v') - f(i, -eps / 2, 'b')) / eps)

        self.numericJacobian = np.array(jacobian).transpose()
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
        DVk = np.linalg.solve(self.jacobian, deltaPQ)
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
        print("Numeric Jacobian:")
        print(np.around(self.numericJacobian, 2))
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


