from __future__ import print_function

# Warning: from fenics import * will import both `sym` and
# `q` from FEniCS. We therefore import FEniCS first and then
# overwrite these objects.
from fenics import *

def q(u):
    "Return nonlinear coefficient"
    return 1 + u**2

# Use SymPy to compute f from the manufactured solution u
import sympy as sym
x, y = sym.symbols('x[0], x[1]')
u = 1 + x + 2*y
f = - sym.diff(q(u)*sym.diff(u, x), x) - sym.diff(q(u)*sym.diff(u, y), y)
f = sym.simplify(f)
u_code = sym.printing.ccode(u)
f_code = sym.printing.ccode(f)
print("f = ", f_code)
print("u = ", u_code)

# Create mesh and define function space
mesh = UnitSquareMesh(8, 8)
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
# get array componets and triangulation :
v = u.compute_vertex_values(mesh)
x = mesh.coordinates()[:,0]
y = mesh.coordinates()[:,1]
t = mesh.cells()

from pylab import *
from mpl_toolkits.axes_grid1 import make_axes_locatable

ax = plt.axes()

cm = get_cmap('viridis')
c = ax.tricontourf(x, y, t, v, 10, cmap = cm)
p = ax.triplot(x, y, t, '-', color='k', lw=0.2, alpha=0.4)

# Output in the file
for ext in ["png", "pdf"]:
    print("Nonlinear_Poisson.%s" % (ext,))
    plt.savefig("Nonlinear_Poisson.%s" % (ext,), bbox_inches="tight")


# #--------------------------Ploting PVD------------------------
# # Plot solution and mesh
# plot(u)
# plot(mesh)
#
# # Save solution to file in VTK format
# vtkfile = File('poisson/solution.pvd')
# vtkfile << u
#
# # Compute maximum error at vertices
# vertex_values_u_D = u_D.compute_vertex_values(mesh)
# vertex_values_u = u.compute_vertex_values(mesh)
