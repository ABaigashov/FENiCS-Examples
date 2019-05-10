#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 10 13:44:58 2019

@author: astrobfu
"""

from fenics import *
from mshr import *
from math import sin, cos, pi

R=10

b=1

a=0.001



domain = Sphere(Point(0, 0,0), R)

cylinder = Sphere(Point(0, 0,0), b) - Sphere(Point(0, 0,0), a)



mesh = generate_mesh(domain, 32)

#plot(mesh)

x = mesh.coordinates()[:,0]
y = mesh.coordinates()[:,1]
z=mesh.coordinates()[:,2]
t = mesh.cells()

import pylab as plb

ax = plb.axes()
cm = plb.get_cmap('viridis')

ax.triplot(x, y, z, t, '-', color='k', lw=0.2, alpha=0.4)

# Output in the file
# print("Lenear_Poisson.pdf")
# plb.savefig("results/lenear_Poisson.%s" % "pdf", bbox_inches="tight")
print("Lenear_Poisson.png")
plb.savefig("results/Horror.%s" % "png", bbox_inches="tight")
