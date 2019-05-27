from table_func_approximator_Lagrange import approximator
from fenics import *
import sympy as sym
import mshr

f = open('data/table_func.txt', 'r')

f_code = approximator(f)

def q(u):
    "Return nonlinear coefficient"
    return 1 + u**2

# Use SymPy to compute f from the manufactured solution u
x, y = sym.symbols('x[0], x[1]')
u = 1 + x + 2*y
f = - sym.diff(q(u)*sym.diff(u, x), x)*sym.diff(q(u)*sym.diff(u, y), y)
u_code = sym.printing.ccode(u)

# Create mesh and define function space
domain =  mshr.Rectangle(Point(0., 0.), Point(8., 8.))
mesh = mshr.generate_mesh(domain, 32)

V = FunctionSpace(mesh, 'P', 1)

# Define boundary condition
u_D = Expression(u_code, degree=2)

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, u_D, boundary)

# Define variational problem
u = Function(V)  # Note: not TrialFunction!
v = TestFunction(V)
f = Expression(f_code, degree=2)
F = q(u)*dot(grad(u), grad(v))*dx - f*v*dx

# Compute solution
solve(F == 0, u, bc)


#--------------------------Ploting pylab------------------------
#get array componets and triangulation :
v = u.compute_vertex_values(mesh)
x = mesh.coordinates()[:,0]
y = mesh.coordinates()[:,1]
t = mesh.cells()

import pylab as plb

ax = plb.axes()
cm = plb.get_cmap('viridis')

ax.tricontourf(x, y, t, v, 10, cmap = cm)
ax.triplot(x, y, t, '-', color='k', lw=0.2, alpha=0.4)

# Output in the file
# print("Lenear_Poisson.pdf")
# plb.savefig("results/lenear_Poisson.%s" % "pdf", bbox_inches="tight")
print("Lenear_Poisson.png")
plb.savefig("results/nonlenear_Poisson.%s" % "png", bbox_inches="tight")
