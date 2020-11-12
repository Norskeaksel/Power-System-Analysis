import importlib
import os, sys

try:
    open('ResultsAssignment2.txt', 'w').close()
except:
    pass
sys.path.append('..')

from ContinuationPowerFlow import newtonRapson2
importlib.reload(newtonRapson2)
from ContinuationPowerFlow.newtonRapson2 import *


def z(r, x):
    return complex(r, x)


def printPredictionVector(predictionVector):
    fprint("\nPrediction vector:", predictionVector, "meaning:")
    c = 0
    for i in Pnr:
        fprint("Delta D", i, " = ", predictionVector[c], sep="")
        c += 1

    for i in Qnr:
        fprint("Delta V", i, " = ", predictionVector[c], sep="")
        c += 1

    fprint("Delta S", " = ", predictionVector[c], '\n', sep="")


# defining the Ybus elements
r01 = 0.1
x01 = 0.2
r02 = 0.05
x02 = 0.25
r12 = 0.05
x12 = 0.15
lines = {}
lines[0, 1] = z(r01, x01)
lines[0, 2] = z(r02, x02)
lines[1, 2] = z(r12, x12)

# defining the known states
d2 = 0
v2 = 1

# defining the power injections
P0sch = -0.8
Q0sch = -0.5
P1sch = -0.4
Q1sch = -0.5
P2sch = 0
Q2sch = 0
PQsch = np.array([P0sch, P1sch, Q0sch, Q1sch])

# defining the inital values of unknowns
d0 = 0
d1 = 0
v0 = 1
v1 = 1

X = [d0, d1, v0, v1]

# defining which corresponding powerflow equations to use
Pnr = [0, 1]
Qnr = [0, 1]
Snr = [0, 1]
# defining the buses in the system

P = np.array([P0sch, P1sch, P2sch])
Q = np.array([Q0sch, Q1sch, Q2sch])
V = np.array([v0, v1, v2])
D = np.array([d0, d1, d2])

slackbus = 2
allowedMissmatch = 1e-5

# define increase
ba = np.array([0.3, 0.7, 0, 0])
fprint("Task 1:")
PS = newtonRapson2(lines, X, PQsch, P, Q, V, D, Pnr, Qnr, slackbus, allowedMissmatch)
buses = PS.buses
fprint("Base case conditions assuming flat start:\n")
for i in buses:
    fprint("bus", i, buses[i])

# Task 2
fprint('\nTask 2:')
oneCol = len(Pnr + Qnr)  # which column the extended jacobian should contain a 1
jacobian=PS.extendJacobian(ba, oneCol)
fprint('Extended Jacobian:')
fprint(jacobian)
predictionVector = PS.buildPredictionVector()
printPredictionVector(predictionVector)


# Task 3
fprint('\nTask3:')
step = 0.3
PS.takePredictionStep(ba, step)
P = PS.PQsch[:len(PS.PQsch) // 2]
Q = PS.PQsch[len(PS.PQsch) // 2:]
# iterate until solution is found for new load:
i = 0
PS.print(i)
while 1:
    i += 1
    PS.CPFiteration(ba, oneCol)
    PS.print(i)
    maxActiveDeviation = max(abs(P[j] - PS.buses[j].p) for j in Pnr)
    maxReactiveDeviation = max(abs(Q[j] - PS.buses[j].q) for j in Qnr)
    maxEffectDeviation = max(maxActiveDeviation, maxReactiveDeviation)
    if maxEffectDeviation < allowedMissmatch:
        break

    lastEffectDeviation = maxEffectDeviation

fprint('Final bus values:')
buses = PS.buses
for i in buses:
    fprint("bus", i, buses[i])

# Task 4
fprint('\nTask 4:')
predictionVector = PS.buildPredictionVector()
printPredictionVector(predictionVector)

# find step that corresponds with load with 30% increase
step = -sum(P) * 0.3
newP = P - step * ba[:len(Pnr)]
# Update loads
for idx, val in enumerate(newP):
    PS.buses[idx].p = val

# Task 5
fprint('Task 5:')
PS.takePredictionStep(ba, step)

# find bus with larges rate of change
maxV = 0
maxVIdx = 0
for idx, val in enumerate(predictionVector[:-1]):
    if abs(val) > maxV:
        maxV = abs(val)
        maxVIdx = idx

oneCol = maxVIdx
P = PS.PQsch[:len(PS.PQsch) // 2]
Q = PS.PQsch[len(PS.PQsch) // 2:]

i = 0
# iterate until solution is found for new load:
oldDeviation=1e5
while 1:
    i += 1
    PS.CPFiteration(ba, oneCol)
    PS.print(i)
    maxActiveDeviation = max(abs(P[j] - PS.buses[j].p) for j in Pnr)
    maxReactiveDeviation = max(abs(Q[j] - PS.buses[j].q) for j in Qnr)
    maxEffectDeviation = max(maxActiveDeviation, maxReactiveDeviation)
    if maxEffectDeviation >= oldDeviation or i == 10:
        break
    oldDeviation=maxEffectDeviation

fprint('Task 5, Final bus values:')
buses = PS.buses
for i in buses:
    fprint("bus", i, buses[i])

print("All done. Results saved to ResultAssignment2.txt")
