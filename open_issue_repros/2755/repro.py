import casadi as ca
try:
    x=ca.SX.sym('x',2)
    s=ca.nlpsol('s','bonmin',{'x':x,'f':ca.sumsqr(x-1.3)},{'discrete':[True,False],'print_time':False})
    r=s(x0=0); print(r['x'])
    ok = abs(float(r['x'][0])-1)<1e-6
except Exception as e:
    print("EXC", e); ok=False
print("BUG ABSENT" if ok else "BUG PRESENT")
