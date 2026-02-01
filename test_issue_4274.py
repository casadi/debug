"""Test for issue #4274: Windows binaries crash on thread map of a Function with IO"""
from casadi import *

x = MX.sym('x')

c = Function('x', [x], [x.printme(4)**2], ["x"], ["out"])
local_inputs = {"x": 5}

F = c.map(2, "thread")

wide_inputs = {}
for k, v in local_inputs.items():
    wide_inputs[k] = horzcat(local_inputs[k], local_inputs[k] * 1.2)

print("wide_inputs:", wide_inputs)
result = F(**wide_inputs)
print("result:", result)
print("SUCCESS: No crash occurred!")
