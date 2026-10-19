import casadi as ca
x=ca.MX.sym('x',2)
try:
    s=ca.nlpsol('s','fatrop',{'x':x,'f':ca.sumsqr(x-1)},{'structure_detection':'none'})
    r=s(); print(r['x'])
    print("BUG ABSENT")
except Exception as e:
    print(e); print("BUG PRESENT")
