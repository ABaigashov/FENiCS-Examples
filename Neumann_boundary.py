from __future__ import print_function
from fenics import *

# Create mesh and define function space
mesh = UnitSquareMesh(8, 8)
V = FunctionSpace(mesh, 'P', 1)
  
# Define Dirichlet boundary condition
u_D = Expression('1 + x[0]*x[0] + x[1]*x[1]', degree = 2)

tol = -1e16

def boundary_D(x, on_boundary):
    if on_boundary:
        if near(x[0], 0, tol) or near(x[1], 1, tol):
            return True
        else:
            return False
    else:
        return False

# Define Neuman boundary condition
g = Expression('0', degree = 1)

bc = DirichletBC(V, u_D, boundary_D)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(-6.0)
a = dot(grad(u), grad(v))*dx
L = f*v*dx - g*v*ds

# Compute solution
u = Function(V)
solve(a == L, u, bc)

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
    print("Neumann_boundary.%s" % (ext,))
    plt.savefig("Neumann_boundary.%s" % (ext,), bbox_inches="tight")
