import numpy as np
from math import *
import importlib

from NewtonRapson import Settings1
importlib.reload(Settings1)
from NewtonRapson.Settings1 import *

global fprinting
fprinting = 1


def stop_fprinting():
    global fprinting
    fprinting = 0


def fprint(*args, **kwargs):
    if Task1 and fprinting:
        print(*args, **kwargs)
        with open('ResultsAssignment1.txt', 'a') as file:
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
    def __init__(self, lines: dict, buses: dict, slackbus: int, X: list, Pnr: list, Qnr: list, PQsch: np.array,
                 useNumeric=False):
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
        self.PFequations()

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
        for i in Pnr:
            dProws = []
            for k in Dnr:
                dPi_dDi = sum(self.dPi_dDk(V, Y, O, D, i, j, k) for j in range(n))
                dProws.append(dPi_dDi)

            for k in Qnr:
                dPi_dVi = sum(self.dPi_dVk(V, Y, O, D, i, j, k) for j in range(n))
                dProws.append(dPi_dVi)

            jacobian.append(dProws)

        for i in Qnr:
            dQrows = []
            for k in Pnr:
                dQi_dDi = -sum(self.dQi_dDk(V, Y, O, D, i, j, k) for j in range(n))
                dQrows.append(dQi_dDi)

            for k in Qnr:
                dQi_dVi = -sum(self.dQi_dVk(V, Y, O, D, i, j, k) for j in range(n))
                dQrows.append(dQi_dVi)

            jacobian.append(dQrows)

        self.jacobian = np.array(jacobian)

    def numeric(self, useNumeric=False):
        if useNumeric == False:
            return
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
            jacobian.append((f(i, eps / 2, 'v') - f(i, -eps / 2, 'v')) / eps)

        self.numericJacobian = np.array(jacobian).transpose()

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

    def print(self, itNr, showNumeric=False):
        fprint("ITERATION NR:", itNr)
        fprint("Jacobian:")
        fprint(self.jacobian, '\n')
        if showNumeric:
            fprint("Numeric Jacobian:")
            fprint(self.numericJacobian, '\n')

        fprint('Net injections:')
        for i in self.buses:
            fprint("bus", i, self.buses[i])

        fprint()
        self.printMissmatchVector()
        self.printCorrectionVector()
        fprint()

    def iteration(self, itNr):
        self.buildJacobian()
        # self.numeric() Calculate the Jacobian numericly to check correctness

        PQsch = self.PQsch
        PQk = self.PQk
        buses = self.buses
        n = self.n
        s = self.slackbus
        deltaPQ = PQsch - PQk
        DVk = np.linalg.solve(self.jacobian, deltaPQ)
        c = 0
        # update angles and voltages based on J^-1*missmatch
        for i in self.Pnr:
            buses[i].d += DVk[c]
            c += 1

        for i in self.Qnr:
            buses[i].v += DVk[c]
            c += 1

        # update system before next iteration
        self.DVk = DVk
        self.buses = buses
        self.PFequations()
