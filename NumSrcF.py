from fenics import *
import numpy as np
from mshr import *

# solve -laplace(u) = f in [0, 1]^2 with homog. Dirichlet boundary conditions
# and funny 'numerical' rhs

#mesh = RectangleMesh(-1,-1,1,1,10,10,'right')

domain=Rectangle(Point(-1,-1), Point(10, 10))

mesh=generate_mesh(domain, 32)


# source term calculation
F = FunctionSpace(mesh, 'CG', 1)  # so that dofs are only in mesh vertices
n = F.dim()

print(n)



d = mesh.geometry().dim()

print(d)


f = Function(F)

print(f)




# dofmap approach
F_dof_coordinates = F.tabulate_dof_coordinates().reshape(n,d)

print(F_dof_coordinates[10])

F_dof_coordinates.resize((n, d))
F_dof_x = F_dof_coordinates[:, 0]
F_dof_y = F_dof_coordinates[:, 1]

f_values = np.zeros(n)
f_values = abs(F_dof_x + F_dof_y)
f.vector()[:] = f_values
plot(f, title='Source term')


print(F_dof_coordinates[1000])
print(f_values[1000])

# # vertex to dofmap approach
# vertex_values = np.zeros(mesh.num_vertices())
# for vertex in vertices(mesh):
#   x = vertex.x(0)
#   y = vertex.x(1)
#   vertex_values[vertex.index()] = abs(x + y)
#
# f.vector()[:] = vertex_values[dof_to_vertex_map(F)]
# plot(f, title='Source term')

# variational problem definition
V = FunctionSpace(mesh, 'CG', 2)

u = TrialFunction(V)
v = TestFunction(V)

a = inner(grad(u), grad(v))*dx
L = f*v*dx
bc = DirichletBC(V, Constant(0.), DomainBoundary())

u = Function(V)
solve(a == L, u, bc)

plot(u, title='Solution')
interactive()
