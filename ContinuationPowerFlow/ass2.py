import importlib
import os
try:
    open('ResultsAssignment2.txt', 'w').close()
except:
    pass
os.chdir("..")  # Necesseary due to a Pycharm bug
from ContinuationPowerFlow import newtonRapson2

importlib.reload(newtonRapson2)
from ContinuationPowerFlow.newtonRapson2 import *
from matplotlib import pyplot as plt

def z(r, x):
    return complex(r, x)

def printPredictionVector(predictionVector,taskNr):
    fprintResults("\nTask", taskNr,"Prediction vector:", predictionVector, "meaning:")
    c = 0
    for i in Pnr:
        fprintResults("Delta D", i, " = ", predictionVector[c], sep="")
        c += 1

    for i in Qnr:
        fprintResults("Delta V", i, " = ", predictionVector[c], sep="")
        c += 1

    fprintResults("Delta S", i, " = ", predictionVector[c], '\n', sep="")

# defining the Ybus elements
r12 = 0.1
x12 = 0.2
r13 = 0.05
x13 = 0.25
r23 = 0.05
x23 = 0.15
lines = {}
lines[0, 1] = z(r12, x12)
lines[0, 2] = z(r13, x13)
lines[1, 2] = z(r23, x23)

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

PS = newtonRapson2(lines, X, PQsch, P, Q, V, D, Pnr, Qnr, slackbus, allowedMissmatch, ba)
buses=PS.buses
fprintResults("Task 1, Base case conditions assuming flat start:\n")
for i in buses:
    fprintResults("bus",i,buses[i])



#Task 2
oneCol=len(Pnr+Qnr) # which column the extended jacobian should contain a 1
PS.extendJacobian(ba,oneCol)
predictionVector=PS.buildPredictionVector()
printPredictionVector(predictionVector,2)

#Task 3
step=0.3
PS.takePredictionStep(ba,step)
P=PS.PQsch[:len(PS.PQsch)//2]
Q=PS.PQsch[len(PS.PQsch)//2:]
# iterate until solution is found for new load:
i=0
while 1:
    i += 1
    PS.CPFiteration(ba,oneCol)
    maxActiveDeviation = max(abs(P[j] - PS.buses[j].p) for j in Pnr)
    maxReactiveDeviation = max(abs(Q[j] - PS.buses[j].q) for j in Qnr)
    maxEffectDeviation = max(maxActiveDeviation, maxReactiveDeviation)
    if maxEffectDeviation < allowedMissmatch:
        break

    lastEffectDeviation = maxEffectDeviation


fprintResults('Task 3, New bus values:')
buses=PS.buses
for i in buses:
    fprintResults("bus",i,buses[i])
    
#Task 4
predictionVector=PS.buildPredictionVector()
printPredictionVector(predictionVector,len(X))

#find step that corresponds with load with 30% increase
step=-sum(P)*0.3
newP=P-step*ba[:len(Pnr)]
for idx,val in enumerate(newP):
    PS.buses[idx].p=val

#Task 5
PS.takePredictionStep(ba,step)

#find bus with larges rate of change
maxV=0
maxVIdx=0
for idx,val in enumerate(predictionVector[:-1]):
    if abs(val) > maxV:
        maxV = abs(val)
        maxVIdx=idx

oneCol=maxVIdx
P=PS.PQsch[:len(PS.PQsch)//2]
Q=PS.PQsch[len(PS.PQsch)//2:]

i=0
while 1:
    i += 1
    PS.CPFiteration(ba,oneCol)
    maxActiveDeviation = max(abs(P[j] - PS.buses[j].p) for j in Pnr)
    maxReactiveDeviation = max(abs(Q[j] - PS.buses[j].q) for j in Qnr)
    maxEffectDeviation = max(maxActiveDeviation, maxReactiveDeviation)
    if maxEffectDeviation < allowedMissmatch or i==10:
        break

fprintResults('Task 5, Final bus values:')
buses=PS.buses
for i in buses:
    fprintResults("bus",i,buses[i])


filename = os.path.basename('ResultsAssignment2.txt')
dest = os.path.join(assignmentName, filename)
shutil.move('Results.txt', dest)

print("All done")
