from scipy.interpolate import lagrange
import numpy as np
import sympy as sym

f = open('table_func.txt', 'r')
lines = f.read().split()
f.close()

x = sym.symbols('x[0]')

abscissa = []
ordinate = []
for i in range(0, len(lines)-1, 2):
    abscissa.append(float(lines[i]))
    ordinate.append(float(lines[i+1]))

poly = lagrange(abscissa, ordinate)



f_code = sym.printing.ccode(poly)

print(f_code)
