import numpy as np
import sympy as smp

f0=0.45

R0=5

x=smp.Symbol('x')

N=100

f=f0*(1-x*x/(R0*R0))


density = []
pressure = []
enthalpy = []

f = open('EOS.txt', 'r')

lines = f.read().split()

f.close()

for i in range(0, len(lines)-2, 3):
    enthalpy.append(float(lines[i]))
    density.append(float(lines[i+1]))
    pressure.append(float(lines[i+2]))

print(enthalpy[151], density[151], pressure[151])
