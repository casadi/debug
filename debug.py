from casadi import *

x = SX.sym('x', Sparsity.diag(4))

print(SX(-2))

print(SX(1,1)**(-2))

x = MX.sym('x', Sparsity.diag(4))

print(MX(-2))
print(x**(-2))
