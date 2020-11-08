import pyomo.environ as pyo
from pyomo.opt import SolverFactory
def solve(B,F,costs,loads,transCap):

    # initialize model
    model = pyo.ConcreteModel()
    model.N = pyo.Set(ordered=True, initialize=[0, 1, 2, 3])  # Set of the 4 nodes in our system
    model.P = pyo.Var(model.N,domain=pyo.NonNegativeReals)  # Make 4 P-variables that represent power generated at each node
    model.delta = pyo.Var(model.N)  # Make 4 variables that represent the voltage angles at each node
    model.constraints = pyo.ConstraintList()  # Define the list of all constraints

    # Define the objective function to be the sum of the costs of the generators times their respective power production
    def ObjectiveFunction(model):
        return sum(model.P[i] * costs[i] for i in model.N)


    model.obj = pyo.Objective(rule=ObjectiveFunction, sense=pyo.minimize)

    # Set slack bus angle to zero
    model.constraints.add(model.delta[3] == 0)

    # Define constraint P+L=B*Delta
    for i in model.N:
        lhs = model.P[i] + loads[i]
        rhs = sum(B[i][c] * model.delta[c] for c in model.N)
        model.constraints.add(lhs == rhs)

    # Define constraint F*Delta<=FLmax
    for i in model.N:
        lhs=sum(F[i][c] * model.delta[c] for c in model.N)
        model.constraints.add(lhs <= transCap[2-i])
        model.constraints.add(lhs >= -transCap[2-i])

    # Solve model
    opt = SolverFactory("gurobi")
    model.dual = pyo.Suffix(direction=pyo.Suffix.IMPORT)
    results = opt.solve(model, load_solutions=True)
    return model