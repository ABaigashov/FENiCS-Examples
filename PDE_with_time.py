from fenics import *
import numpy as np

T = 20.0 # final time
num_steps = 100 # number of time steps
dt = T/num_steps # time step size
alpha = 30 # parameter alpha
beta = 12 # parameter beta

# Create mesh and define function space
nx = ny = 64
mesh = UnitSquareMesh(nx, ny)
V = FunctionSpace(mesh, 'P', 1)

# Define boundary condition
u_D = Expression('1 + x[0]*x[0] + alpha*x[1]*x[1] + beta*t',
                 degree=2, alpha=alpha, beta=beta, t=0)

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, u_D, boundary)

# Define initial value
u_n = interpolate(u_D, V)
#u_n = project(u_D, V)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(beta - 2 - 2*alpha)

F = u*v*dx + dt*dot(grad(u), grad(v))*dx - (u_n + dt*f)*v*dx
a, L = lhs(F), rhs(F)

# Time-stepping
u = Function(V)
t = 0
for n in range(num_steps):
    # Update current time
    t += dt
    u_D.t = t
    # Compute solution
    solve(a == L, u, bc)
    # Plot solution

    #--------------------------Ploting pylab------------------------
    #get array componets and triangulation :
    v = u.compute_vertex_values(mesh)
    x = mesh.coordinates()[:,0]
    y = mesh.coordinates()[:,1]
    tt = mesh.cells()

    import pylab as plb

    ax = plb.axes()
    cm = plb.get_cmap('viridis')
    #
    # ax.tricontourf(x, y, tt, v, 10, cmap = cm)
    # ax.triplot(x, y, tt, '-', color='k', lw=0.2, alpha=0.4)

    ax.tricontourf(x, y, tt, v, 10, cmap = cm)
    #ax.triplot(x, y, tt, '-', color='k', lw=0.2, alpha=0.4)

    # Output in the file
    # print("Lenear_Poisson.pdf")
    # plb.savefig("results/lenear_Poisson.%s" % "pdf", bbox_inches="tight")
    print("heat_eauation{}.png".format(n))
    plb.savefig("results/heat_eauation{}.%s".format(n) % "png", bbox_inches="tight")

    # # Compute error at vertices
    # u_e = interpolate(u_D, V)
    # error = np.abs(u_e.vector().array() - u.vector().array()).max()
    # print('t = %.2f: error = %.3g' % (t, error))
    # # Update previous solution
    u_n.assign(u)