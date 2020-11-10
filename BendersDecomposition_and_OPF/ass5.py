import importlib
import os
import shutil
import pandas as pd
import numpy as np
from copy import deepcopy

df = pd.read_excel('systemData.xlsx',converters={'Line':str})
print(df)
os.chdir("..")  # Necesseary due to a Pycharm bug
from BendersDecomposition_and_OPF import functions2
importlib.reload(functions2)
from BendersDecomposition_and_OPF.functions2 import *
from BendersDecomposition_and_OPF import OPF_Solver
importlib.reload(OPF_Solver)
from BendersDecomposition_and_OPF.OPF_Solver import *
from BendersDecomposition_and_OPF import subproblem
importlib.reload(subproblem)
from BendersDecomposition_and_OPF.subproblem import *
try:
    open('ResultsAssignment5.txt', 'w').close()
except:
    pass


n=4 # nr of buses
slackbus=3
# Parameters
costs = df["GenCosts"].tolist()
loads = df["Loads"].tolist()
transCap = df["Capacity"].tolist()
reactance = df["Reactance"].tolist()

lines = {}
for idx,ik in enumerate(df["Line"].tolist()):
    i=int(ik[0])
    k=int(ik[1])
    lines[i,k]=reactance[idx]

B=buildDCY(lines, n)

# create connection matrix F which is the load flow equations described in matrix form with a zero diagonal
#P_ij = -B_ij(delta_i-delta_j)
#==>
#P_12 = 0 * delta_0 -B_12 *delta_1  + B_12 *delta_2 + 0
#P_23 =0 + 0 + delta_0 * B_01 + B_01 *delta_1  + 0 *delta_2
#P_01 = delta_0 * B_01 + B_01 *delta_1  + 0 *delta_2 + 0
#P_02 =  delta_0 * B_02 + 0 *delta_1  + B_02 *delta_2 + 0



F = np.array([[0,          -B[1,2],     B[1,2],     0], #P12
              [0,          0,          -B[2,3],     B[2,3]], #P23
              [-B[0,1],    B[0,1],      0,          0], #P01
              [-B[0,2],    0,           B[0,2],     0], #P02
              ])
fprint('Task 1:')
fprint('System matrix:')
fprint(B)

model=solve(n,B,F,costs,loads,transCap)
P=pyomoResults(model,F)

#Task 2
fprint('\nTask 2:')
model.dual.display()
fprint("Disclamer: I don't know how to change the key names so we must count the constraints set in the code to see what the keys correspond to. Doing that, we find that")
fprint('constraints 2-5 show the dual value of the P+L=B*Delta constraints and represents cost of increasing the load (marginal cost) in bus 0-3 respectivly.')
dual=[4.0,3.5,3.0,2.0]
fprint('Dual values:',dual)
fprint('We see that the dual values equal the operating costs in their buses except for bus 1 where the operating cost is 5, but the dual value is 3.5.')
fprint('This means that we can increase the load in bus 1 without using it\'s generator.')

# Part 2
fprint('\nPart 2:')
# Task 1
B=np.delete(B, slackbus, 0)
B= np.delete(B, slackbus, 1)
D0=[-model.delta[i].value for i in model.N if i!=slackbus]
#Predict angles
angles=IMML_angles(B,D0,1,2,0.5)
angles=np.append(angles,0)
#Update lines
lines[1,2]*=2
Y=buildDCY(lines,n)
# Calculate flows in new system with predicted angles
PF,PF_str=powerFlows(Y,angles)
fprint('\nTask 2:\nThe voltage angles in the load flow solutions using IMML are:')
fprint(angles)
fprint('With these angles we get the following load flow solution:')
fprint(PF)
fprint('From this matrix we see the net power generation/demand in each bus as well as the following power flows:')
for i in PF_str:
    fprint(i)
fprint('We see that the IMML approach gives an infeasible solution as line 0-2 is overloaded')

fprint('\nTask 2 see file subproblem.py')

# Task 3
B=-Y
B=np.delete(B, slackbus, 0)
B= np.delete(B, slackbus, 1)


# calculate PTDF matrix
Z=deepcopy(Y)
Z[slackbus,slackbus]+=1
Z=np.linalg.inv(Z)

keys=lines.keys()
# define line numbers
ik=build_ik(keys,n)
lines_str = []
for i, k in keys:
    lines_str.append(f"line {i}-{k} ")

PTDF=buildPTDF=buildPTDF(Z,Y,ik,n)
#printPTDF(PTDF,lines_str)

sub=subSolve(transCap[1],PF,PTDF,1)
# New generation
deltaP=np.zeros(n)
deltaP+= np.array([sub.Pup[i].value for i in sub.N])
deltaP-= np.array([sub.Pdown[i].value for i in sub.N])
P+=deltaP
flowChange=(PTDF*deltaP).sum(axis=1)
printFlowChange(deltaP,lines_str,flowChange)
PF = updatePF(PF,lines_str,flowChange,n)
fprint('\nTask 3:')
fprint('Final solution:')
fprint('Loads:')
fprint(loads)
fprint('Generation:')
fprint(P)
fprint('Flows:')
fprint(PF)
for i in range(n):
    for k in range(n):
        if k>i and PF[i,k]!=0:
            fprint('Line ',i,'-',k,' = ',PF[i,k],sep="")

missmatch=np.zeros(n)
for i in range(n):
    missmatch[i]=loads[i]+P[i]-PF[i,i]

#print(missmatch)

#Task 4
fprint('\nTask 4:')
fprint('Add constraint sum(model.P[i]-P[i])*dk_dp[i] for i in model.N) <=0')
Ks=sub.obj()
rc=np.array([i[-1] for i in sub.rc.items()])[len(sub.rc)//2:]
subCost=np.ones(n)
dk_dp=subCost-rc
dk_dp[1]=-1
model.constraints.add(sum((model.P[i]-P[i])*dk_dp[i] for i in model.N) <=0)

#Task 5

opt = SolverFactory("gurobi")
opt.solve(model, load_solutions=True)
fprint('\nTask 5, base case:')
pyomoResults(model,F)

fprint('\nTask 5, New case:')
B=buildDCY(lines, n)
F = np.array([[0,          -B[1,2],     B[1,2],     0], #P12
              [0,          0,          -B[2,3],     B[2,3]], #P23
              [-B[0,1],    B[0,1],      0,          0], #P01
              [-B[0,2],    0,           B[0,2],     0], #P02
              ])
model=solve(n,B,F,costs,loads,transCap,P,dk_dp)
pyomoResults(model,F)

filename = os.path.basename('ResultsAssignment5.txt')
dest = os.path.join("BendersDecomposition_and_OPF", filename)
shutil.move('ResultsAssignment5.txt', dest)
