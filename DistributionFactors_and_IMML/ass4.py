import importlib
import os
import shutil
from copy import deepcopy
os.chdir("..")  # Necesseary due to a Pycharm bug
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
fprint('This means that we set Y=-Y, that we disregard grid losses and that Oij becomes pi/2')
fprint('2. All bus voltages are constant at 1pu')
fprint('3. Voltage angle differences between buses are smal')
fprint('This assumption means that we can linearize the sine and cosine terms in the AC power flow equations' )
fprint('With these assumtions, our new DC power flow equation become:')
fprint('Pi=sum(Bij*(Di-Dj) for j in range(n)). The term Bij*(Di-Dj) = the power flow Pij\n')

#define nr of buses
n=4
# defining the Ybus elements
x12 = 0.2
x13 = 0.1
x23 = 0.25
x34 = 0.25
# defining the loads
P0=-1.25
P1=-0.4
P2=-0.6
P=np.array([P0,P1,P2])
lines = {}
lines[0, 1] = x12
lines[0, 2] = x13
lines[1, 2] = x23
lines[2, 3] = x34
slackbus=3

keys=lines.keys()
# define line numbers
ik=build_ik(keys,n)

fprint('With the buildDCY function we get that the Y matrix in our system is:')
Y=buildDCY(lines, n)
B=np.delete(-Y, slackbus, 0)
B= np.delete(B, slackbus, 1)
fprint(B)

#Task 2
angles=np.linalg.solve(B, P)
angles=np.append(angles,0)
P=np.append(P,-sum(P))
PF=powerFlows(Y,angles,P)
fprint('\nTask 2:\nThe voltage angles in the load flow solutions are:')
fprint(angles)
fprint('The power flows are given by the following matrix:')
fprint(PF)
fprint('From this matrix we see power generation/demand in each bus as well as the following power flows:')
PF=powerFlows(Y,angles,P,True)

#Task 3
fprint('\nTask 3:')
fprint('Power Transfer Distribution Factors (PTDF) represents the change in real power transfer that occurs on transmission lines when their corresponding buses change their power injections.')
fprint('Each row in the PTDF corresponds to a line and each column corresponds to a bus')
# the  Zbus matrix  is constructed  by  adding  "+1"  to  the  diagonal  element  corresponding  to  the  slack-node  in  Y,  followed  by  an  inverse  operation.
Z=deepcopy(Y)
Z[slackbus,slackbus]+=1
Z=np.linalg.inv(Z)


#b=build_b(Z,keys,n)
#print(b)
PTDF=buildPTDF=buildPTDF(Z,Y,ik,n)
printPTDF(PTDF,keys)
filename = os.path.basename('ResultsAssignment4.txt')
dest = os.path.join("DistributionFactors_and_IMML", filename)
shutil.move('ResultsAssignment4.txt', dest)