import pyomo.environ as pyo
from pyomo.opt import SolverFactory
def subSolve(cap, PF,PTDF,line):

    # initialize sub
    sub = pyo.ConcreteModel()
    sub.N = pyo.Set(ordered=True, initialize=[0, 1, 2, 3])  # Set of the 4 nodes in our system
    sub.Pup = pyo.Var(sub.N,domain=pyo.NonNegativeReals)  # Make 4 P-variables that represent increased power generation at each node
    sub.Pdown = pyo.Var(sub.N, domain=pyo.NonNegativeReals)  # Make 4 P-variables that represent decreased power generation at each node
    sub.constraints = pyo.ConstraintList()  # Define the list of all constraints

    # Define the objective function to be the sum of change in power production
    def ObjectiveFunction(sub):
        return sum(sub.Pup[i] + sub.Pdown[i] for i in sub.N)

    sub.obj = pyo.Objective(rule=ObjectiveFunction, sense=pyo.minimize)
    # Load balance constraint
    sub.constraints.add(sum(sub.Pup[i]-sub.Pdown[i] for i in sub.N) == 0)

    # Constraint for not violation of line 0-2
    lhs=-cap - PF[0, 2]
    mid= sum((sub.Pup[i] - sub.Pdown[i])*PTDF[line,i] for i in sub.N)
    rhs= cap - PF[0, 2]
    sub.constraints.add(pyo.inequality(lhs,mid,rhs))

    # Solve model
    opt = SolverFactory("gurobi")
    sub.dual = pyo.Suffix(direction=pyo.Suffix.IMPORT)
    sub.rc = pyo.Suffix(direction=pyo.Suffix.IMPORT)
    results = opt.solve(sub, load_solutions=True)
    return sub

