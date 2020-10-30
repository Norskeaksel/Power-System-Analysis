import importlib
import os
try:
    open('ResultsAssignment3.txt', 'w').close()
except:
    pass
os.chdir("..")  # Necesseary due to a Pycharm bug
from DecoupledPowerFlow import newtonRapson3
importlib.reload(newtonRapson3)
from DecoupledPowerFlow.newtonRapson3 import *
from DecoupledPowerFlow import DPF
importlib.reload(DPF)
from DecoupledPowerFlow.DPF import *

def z(r, x):
    return complex(r, x)

# defining the Ybus elements
r12 = 0.1 #0.05
x12 = 0.2
r13 = 0.05
x13 = 0.2 #0.1
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
P0sch = -1
Q0sch = -0.5
P1sch = -0.5
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
allowedMissmatch = 1e-3

# define increase
ba = np.array([0.3, 0.7, 0, 0])

#fprint("Task 1, Base case conditions assuming flat start:\n")
#PS = newtonRapson3(lines, X, PQsch, P, Q, V, D, Pnr, Qnr, slackbus, allowedMissmatch)

# Task 2

primalDPS,dualDPS=DPF(lines, X, PQsch, P, Q, V, D, Pnr, Qnr, slackbus, allowedMissmatch)

# Task 3

filename = os.path.basename('ResultsAssignment3.txt')
dest = os.path.join("DecoupledPowerFlow", filename)
shutil.move('Results.txt', dest)