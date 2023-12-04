#! /usr/bin/env python

# Flexural wavelength calculator

# 2D

Te = 30000  # m
rho_f = 1000.0  # water
rho_m = 3300.0  # mantle
drho = rho_m - rho_f
E = 65e9  # Youong's modulus
nu = 0.25  # Poisson's ratio
g = 9.8

D = (E * Te**3) / (12 * (1 - nu**2))

alpha1D = (4 * D / (drho * g)) ** 0.25
alpha2D = (D / (drho * g)) ** 0.25

lambda1D = alpha1D * 2 * 3.14159
lambda2D = alpha2D * 2 * 3.14159

print("")

print("1D:")
print("Flexural wavelength:", lambda1D / 1000, "km")
print("Distance to first zero-crossing:", 0.375 * lambda1D / 1000, "km")
print("Flexural parameter:", alpha1D / 1000, "km")

print("")

print("2D:")
print("Flexural wavelength:", lambda2D / 1000, "km")
print("Distance to first zero-crossing:", 0.375 * lambda2D / 1000, "km")
print("Flexural parameter:", alpha2D / 1000, "km")

print("")
