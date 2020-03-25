#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:../lib

import pyctmm
import sys

stack = pyctmm.create_stack(2, 633e-9, 0)
print('\n', stack)

pyctmm.set_ind(stack, 0, 1, 0)
pyctmm.set_ind(stack, 1, 1.515, 0)
pyctmm.set_d(stack, 0, 0)
pyctmm.set_d(stack, 1, 0)

print(pyctmm.get_ind(stack, 0))
print(pyctmm.get_ind(stack, 1))
print(pyctmm.get_d(stack, 0))
print(pyctmm.get_d(stack, 1))

pyctmm.evaluate(stack)

matrix = pyctmm.get_matrix(stack)
print(matrix)
fresnel_coefs = pyctmm.get_amplitude(stack)
print(fresnel_coefs)
power_coefs = pyctmm.get_power(stack)
print(power_coefs)
power_phase_coefs = pyctmm.get_power_phase(stack)
print(power_phase_coefs)