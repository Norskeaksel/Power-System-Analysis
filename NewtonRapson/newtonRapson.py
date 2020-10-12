import importlib
import shutil
import os
import matplotlib as plt

from NewtonRapson import PowerSystem

importlib.reload(PowerSystem)
from NewtonRapson.PowerSystem import *

open('Results.txt', 'w').close()


def buildBuses(P, Q, V, D):
    buses = {}
    for i in range(len(P)):
        buses[i] = Bus(P[i], Q[i], V[i], D[i])

    return buses


def newtonRapson(lines, X, PQsch, P, Q, V, D, Pnr, Qnr, slackbus, allowedMissmatch):

    buses = buildBuses(P, Q, V, D)

    PS = PowerSystem(lines, buses, slackbus, X, Pnr, Qnr, PQsch)
    lastEffectDeviation=1e5
    i = 0
    while 1:
        i += 1
        PS.iteration(i)
        PS.print(i)
        maxActiveDeviation = max(abs(P[j] - PS.buses[j].p) for j in Pnr)
        maxReactiveDeviation = max(abs(Q[j] - PS.buses[j].q) for j in Qnr)
        maxEffectDeviation = max(maxActiveDeviation, maxReactiveDeviation)
        if lastEffectDeviation < maxEffectDeviation:
            fprint("Algorithm not convergin with the following load:")
            fprint(buses)
            break
        if maxEffectDeviation < allowedMissmatch:
            break


    filename = os.path.basename('Results.txt')
    dest = os.path.join('NewtonRapson', filename)
    shutil.move('Results.txt', dest)
