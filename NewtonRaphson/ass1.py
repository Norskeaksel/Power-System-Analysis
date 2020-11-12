import importlib
from matplotlib import pyplot as plt

from NewtonRaphson import newtonRaphson
importlib.reload(newtonRaphson)
from NewtonRaphson.newtonRaphson import *



try:
    open('ResultsAssignment1.txt', 'w').close()
except:
    pass


def z(r, x):
    return complex(r, x)


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

# defining the buses in the system
P = np.array([P0sch, P1sch, P2sch])
Q = np.array([Q0sch, Q1sch, Q2sch])
V = np.array([v0, v1, v2])
D = np.array([d0, d1, d2])

slackbus = 2
allowedMissmatch = 1e-5
buses = newtonRapson(lines, X, PQsch, P, Q, V, D, Pnr, Qnr, slackbus, allowedMissmatch)
stop_fprinting()

# I assume we are to plot the voltage as a function of active power demand as that's whats changing
if Task2:
    buses = newtonRapson(lines, X, PQsch, P, Q, V, D, Pnr, Qnr, slackbus,
                         allowedMissmatch)  # only necessary if skipping task 1 in Settings1.py
    deltaP = 0.2
    bus0Share = 0.3
    bus1Share = 0.7
    calcFlatStart = True
    calcOldStart = True
    oldStartBuses = 0
    flatStartBuses = 0
    newV = [buses[i].v for i in buses]
    newD = [buses[i].d for i in buses]

    load = [sum(i for i in P)]
    flatBus1V = [buses[0].v]
    flatBus2V = [buses[1].v]
    oldBus1V = [buses[0].v]
    oldBus2V = [buses[1].v]
    while 1: # repeat until the algorithm fails to converge for both previous load and flat start
        # increase the bus loads according to specifications
        P[0] -= bus0Share * deltaP
        P[1] -= bus1Share * deltaP
        PQsch[0] = P[0]
        PQsch[1] = P[1]
        load.append(sum(i for i in P))
        # solve system for increased load
        if calcOldStart:
            oldStartBuses = newtonRapson(lines, X, PQsch, P, Q, newV, newD, Pnr, Qnr, slackbus, allowedMissmatch)

        if calcFlatStart:
            flatStartBuses = newtonRapson(lines, X, PQsch, P, Q, V, D, Pnr, Qnr, slackbus, allowedMissmatch)

        try:
            oldBus1V.append(oldStartBuses[0].v)
            oldBus2V.append(oldStartBuses[1].v)
            newV = [oldStartBuses[i].v for i in buses]
            newD = [oldStartBuses[i].d for i in buses]
        except:  # Algorithm not converging
            oldLoad = [-load[i] for i in range(len(oldBus1V))]
            plt.plot(oldLoad, oldBus1V, label="Bus 1")
            plt.plot(oldLoad, oldBus2V, label="Bus 2")
            plt.ylabel('Voltage [pu]')
            plt.xlabel('Active flatLoad [pu]')
            plt.title('Previous start Voltage given active system load')
            plt.legend()
            filename = os.path.basename('VprevActive.png')
            plt.savefig(filename)
            plt.show()
            calcOldStart = False
            oldStartBuses = Bus()

        try:
            flatBus1V.append(flatStartBuses[0].v)
            flatBus2V.append(flatStartBuses[1].v)
        except:  # Algorithm not converging
            flatLoad = [-load[i] for i in range(len(flatBus1V))]
            plt.plot(flatLoad, flatBus1V, label="Bus 1")
            plt.plot(flatLoad, flatBus2V, label="Bus 2")
            plt.ylabel('Voltage [pu]')
            plt.xlabel('Active flatLoad [pu]')
            plt.title('Flat start Voltage given active system load')
            plt.legend()
            filename = os.path.basename('VflatActive.png')
            plt.savefig(filename)
            plt.show()
            calcFlatStart = False
            flatStartBuses = Bus()

        if calcOldStart == 0 and calcFlatStart == 0:
            break

#same as Task 2 but for reactive load
if Task3:
    # reset
    PQsch = np.array([P0sch, P1sch, Q0sch, Q1sch])
    P = np.array([P0sch, P1sch, P2sch])
    Q = np.array([Q0sch, Q1sch, Q2sch])
    V = np.array([v0, v1, v2])
    D = np.array([d0, d1, d2])
    buses = newtonRapson(lines, X, PQsch, P, Q, V, D, Pnr, Qnr, slackbus, allowedMissmatch)

    deltaQ = 0.2
    calcFlatStart = True
    calcOldStart = True
    oldStartBuses = 0
    flatStartBuses = 0
    newV = [buses[i].v for i in buses]
    newD = [buses[i].d for i in buses]

    load = [sum(i for i in Q)]
    flatBus1V = [buses[0].v]
    flatBus2V = [buses[1].v]
    oldBus1V = [buses[0].v]
    oldBus2V = [buses[1].v]
    while 1:
        Q[0] -= 0.3 * deltaQ
        Q[1] -= 0.7 * deltaQ
        PQsch[2] = Q[0]
        PQsch[3] = Q[1]
        load.append(sum(i for i in Q))
        if calcOldStart:
            oldStartBuses = newtonRapson(lines, X, PQsch, P, Q, newV, newD, Pnr, Qnr, slackbus, allowedMissmatch)

        if calcFlatStart:
            flatStartBuses = newtonRapson(lines, X, PQsch, P, Q, V, D, Pnr, Qnr, slackbus, allowedMissmatch)

        try:
            oldBus1V.append(oldStartBuses[0].v)
            oldBus2V.append(oldStartBuses[1].v)
            newV = [oldStartBuses[i].v for i in buses]
            newD = [oldStartBuses[i].d for i in buses]
        except:  # Algorithm not converging
            oldLoad = [-load[i] for i in range(len(oldBus1V))]
            plt.plot(oldLoad, oldBus1V, label="Bus 1")
            plt.plot(oldLoad, oldBus2V, label="Bus 2")
            plt.ylabel('Voltage [pu]')
            plt.xlabel('Reactive flatLoad [pu]')
            plt.title('Previous start Voltage given reactive system load')
            plt.legend()
            filename = os.path.basename('VprevReactive.png')
            # dest = os.path.join(assignmentName, filename)
            plt.savefig(filename)  # Not working
            plt.show()
            calcOldStart = False
            oldStartBuses = Bus()

        try:
            flatBus1V.append(flatStartBuses[0].v)
            flatBus2V.append(flatStartBuses[1].v)
        except:  # Algorithm not converging
            flatLoad = [-load[i] for i in range(len(flatBus1V))]
            plt.plot(flatLoad, flatBus1V, label="Bus 1")
            plt.plot(flatLoad, flatBus2V, label="Bus 2")
            plt.ylabel('Voltage [pu]')
            plt.xlabel('Reactive flatLoad [pu]')
            plt.title('Flat start Voltage given reactive system load')
            plt.legend()
            filename = os.path.basename('VflatReactive.png')
            # dest = os.path.join(assignmentName, filename)
            plt.savefig(filename)
            plt.show()
            calcFlatStart = False
            flatStartBuses = Bus()

        if calcOldStart == 0 and calcFlatStart == 0:
            break

print("All done. Results saved to ResultAssignment1.txt")
