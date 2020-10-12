import importlib
from NewtonRapson import PowerSystem

importlib.reload(PowerSystem)
from NewtonRapson.PowerSystem import *

def z(r, x):
    return complex(r, x)


def buildBuses(P, Q, V, D):
    buses = {}
    for i in range(len(P)):
        buses[i] = Bus(P[i], Q[i], V[i], D[i])

    return buses

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
d2=0
v2=1

# defining the power injections
P0sch = -0.8
Q0sch = -0.5
P1sch = -0.4
Q1sch = -0.5
P2sch = 0
Q2sch = 0

# defining the inital values of unknowns
d0=0
d1=0
v0=1
v1=1

X=[d0,d1,v0,v1]

# defining which corresponding powerflow equations to use
Pnr=[0,1]
Qnr=[0,1]

#defining the buses in the system

P = np.array([P0sch, P1sch, P2sch])
Q = np.array([Q0sch, Q1sch, Q2sch])
V = np.array([v0, v1, v2])
D = np.array([d0, d1, d2])

buses = buildBuses(P, Q, V, D)

PQsch=np.array([P0sch,P1sch,Q0sch,Q1sch])
slackbus=2
PS = PowerSystem(lines, buses,slackbus,X,Pnr,Qnr,PQsch)
final=4
for i in range(1,final+1):
    PS.iteration(i)
    PS.print(i)

