#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 10 13:44:58 2019

@author: astrobfu
"""

import fenics as fen


R=10

b=1

a=0.001



domain =fen.Circle(Point(0, 0), R)

cylinder = fen.Circle(Point(0, 0), b) - fen.Circle(Point(0, 0), a)

domain.set_subdomain(1, cylinder)

mesh = fen.generate_mesh(domain, 32)

plot(mesh)