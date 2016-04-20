#! /usr/bin/env python

# Flexural wavelength calculator

# 2D

h = 60000 # m
rho_f = 917. # ice
rho_m = 3300.
drho = rho_m - rho_f
E=65E9
nu=0.25
g = 9.8

D = (E * h**3) / (1 - nu**2)

alpha1D = (4*D/(drho * g))**.25
alpha2D =   (D/(drho * g))**.25

lambda1D = alpha1D * 2*3.14159
lambda2D = alpha2D * 2*3.14159

print "1D:", lambda1D/1000, 'km'
print "2D:", lambda2D/1000, 'km'
