import importlib
import os, sys
import shutil
from copy import deepcopy

sys.path.append('..')

from DistributionFactors_and_IMML import functions

importlib.reload(functions)
from DistributionFactors_and_IMML.functions import *

try:
    open('ResultsAssignment4.txt', 'w').close()
except:
    pass

fprint('Task 1:\nThe AC powerflow equations are as follows:')
fprint('Pi= sum(Vi * Vj * Yij * cos(Oij - Di + Dj) for j in range(n))')
fprint('Qi= -sum(Vi * Vj * Yij * sin(Oij - Di + Dj) for j in range(n))')
fprint('Where V are voltage magnitudes, Y is admittance magnitude and O is admittance angle\n')
fprint('In DC power flow we make 3 critical assumptions')
fprint('1. Line resistances are negligible compared to line reactances (RL<<XL, for all lines).')
fprint('This means that we set B=-Y(all except slack row and column), that we disregard grid losses and that Oij becomes pi/2')
fprint('2. All bus voltages are constant at 1pu')
fprint('3. Voltage angle differences between buses are smal')
fprint('This assumption means that we can linearize the sine and cosine terms in the AC power flow equations')
fprint('With these assumtions, our new DC power flow equation become:')
fprint('Pi=sum(Bij*(Di-Dj) for j in range(n)). The term Bij*(Di-Dj) = the power flow Pij.\n')

# define nr of buses
n = 4
# defining the Ybus elements
x01 = 0.2
x02 = 0.1
x12 = 0.25
x23 = 0.25
# defining the loads
P0 = -1.25
P1 = -0.4
P2 = -0.6
P = np.array([P0, P1, P2])
lines = {}
lines[0, 1] = x01
lines[0, 2] = x02
lines[1, 2] = x12
lines[2, 3] = x23
slackbus = 3

keys = lines.keys()
# define line numbers
ik = build_ik(keys, n) #a connection matrix with line numbers 0-n or -1 if no connection
lines_str = []
for i, k in keys:
    lines_str.append(f"line {i}-{k} ")

fprint('With the buildDCY function and by removing the slack row and column we get that the B matrix in our system is:')
Y = buildY(lines, n)
B = np.delete(-Y, slackbus, 0)
B = np.delete(B, slackbus, 1)
fprint(B)

# Task 2
angles = np.linalg.solve(B, P) # angles=B^-1*P
angles = np.append(angles, 0)
P = np.append(P, -sum(P)) # assume all generation comes from the last node
PF, PF_str = powerFlows(Y, angles, P)
fprint('\nTask 2:\nThe voltage angles in the load flow solutions are:')
fprint(angles)
fprint('The power flows are given by the following matrix:')
fprint(PF)
fprint('From this matrix we see power generation/demand in each bus as well as the following power flows:')
for i in PF_str:
    fprint(i)

# Task 3
fprint('\nTask 3:')
fprint(
    'Power Transfer Distribution Factors (PTDF) represents the change in real power transfer that occurs on transmission lines when their corresponding buses change their power injections.')
fprint('Each row in the PTDF corresponds to a line and each column corresponds to a bus')
# The  Zbus matrix  is constructed  by  adding  "+1"  to  the  diagonal  element  corresponding  to  the  slack-node  in  Y,  followed  by  an  inverse  operation.
Z = deepcopy(Y)
Z[slackbus, slackbus] += 1
Z = np.linalg.inv(Z)

PTDF = buildPTDF(Z, Y, ik, n)
printPTDF(PTDF, lines_str)

# Task 4
loadChange = [-0.5, 0, 0, 0]
flowChange = (PTDF * loadChange).sum(axis=1)
fprint('\nTask 4:')
printFlowChange(loadChange, lines_str, flowChange)

# Task 5
fprint('\nTask 5:')
fprint('The original power flows were:')
fprint(PF_str)
loadChange = [-0.5, 0.3, 0, 0]
flowChange = (PTDF * loadChange).sum(axis=1)
printFlowChange(loadChange, lines_str, flowChange)

P += loadChange
P = P[:-1]
angles = np.linalg.solve(B, P)
angles = np.append(angles, 0)
P = np.append(P, -sum(P))
fprint('When we solve the load flow with this new load we get power flows')
newPF, newPF_str = powerFlows(Y, angles, P)
fprint(newPF_str)
fprint(
    'We see by adding the changes to the original flows that we get almost the exact same solution as when we resolve with the new load')

# Part 2

# Task 1

fprint('\nPart 2:')
fprint('\nTask 1:')
fprint(
    'The Inverse Matrix Modification Lemma is a technique for solving the new load flow cases that occur when the connections in a solved power system is changed, without the need for resolving evererything')

P = np.array([P0, P1, P2])
i=0
k=1
change=0
angles = IMML_angles(B, P, i, k,change)

fprint('\nTask 2:\nThe voltage angles in the load flow solutions using IMML are:')
fprint(angles)

lines[0, 1] = np.inf
Y = buildY(lines, n)
B = np.delete(-Y, slackbus, 0)
B = np.delete(B, slackbus, 1)
anglesDC = np.linalg.solve(B, P)
anglesDC = np.append(anglesDC, 0)
P = np.append(P, -sum(P))
PF, PF_str = powerFlows(Y, anglesDC, P)

fprint('The voltage angles using DC power flow is')
fprint(anglesDC)
fprint(
    'These angles are identical to the ones found with the IMML approach, which means that the load flow solution is:')
fprint(PF)
fprint('From this matrix we see power generation/demand in each bus as well as the following power flows:')
fprint('P01 =', PF[0, 1])
for i in PF_str:
    fprint(i)

# Task 3
P = np.array([P0, P1, P2])
lines[0, 1] = x01
Y = buildY(lines, n)
B = np.delete(-Y, slackbus, 0)
B = np.delete(B, slackbus, 1)
i=0
k=2
change=0.5
angles = IMML_angles(B, P, i, k, change)

fprint('\nTask 3:\nThe voltage angles in the load flow solutions using IMML are:')
fprint(angles)

lines[0, 2] *= 2
Y = buildY(lines, n)
B = np.delete(-Y, slackbus, 0)
B = np.delete(B, slackbus, 1)
anglesDC = np.linalg.solve(B, P)
anglesDC = np.append(anglesDC, 0)
P = np.append(P, -sum(P))
PF, PF_str = powerFlows(Y, anglesDC, P)

fprint('The voltage angles using DC power flow is')
fprint(anglesDC)
fprint(
    'These angles are identical to the ones found with the IMML approach, which means that the load flow solution is:')
fprint(PF)
fprint('From this matrix we see power generation/demand in each bus as well as the following power flows:')
for i in PF_str:
    fprint(i)
