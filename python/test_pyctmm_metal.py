#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:../lib

import pyctmm
import sys

stack = pyctmm.create_stack(3, 633e-9, 0)

pyctmm.set_ind(stack, 0, 1, 0)
pyctmm.set_ind(stack, 1, 3, -2.999)
pyctmm.set_ind(stack, 2, 1, 0)

pyctmm.set_d(stack, 0, 0)
pyctmm.set_d(stack, 1, 50e-9)
pyctmm.set_d(stack, 2, 0)

pyctmm.evaluate(stack)

power_coefs = pyctmm.get_power(stack)
print('\n', power_coefs)
