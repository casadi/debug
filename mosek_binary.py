import os
os.sched_setaffinity(0, sorted(os.sched_getaffinity(0))[:4])
import casadi as ca
import numpy as np
import subprocess, os
ca.DM.rng(1)
n,m,ms=10,8,5
x=ca.MX.sym('x',n)
M=ca.DM.rand(n,n); H=0.5*(M@M.T)
f=100*ca.DM.rand(n,1); A=ca.DM.rand(m,n)
bupper=20*ca.DM.rand(m,1); blower=-20*ca.DM.rand(m,1)
s=ca.qpsol('solver','mosek',{'f':0.5*(x.T@H@x)+f.T@x,'x':x,'g':A@x},{'discrete':[1]*ms+[0]*(n-ms)})
args=dict(lbx=[0]*ms+[-10]*(n-ms),ubx=[1]*ms+[10]*(n-ms),lbg=blower,ubg=bupper)
s.generate('solver.c')
r=os.environ['MOSEKDIR']
subprocess.run(['gcc','-O3','-fPIC','-shared','solver.c','-I'+r+'/h','-L'+r+'/bin','-Wl,-rpath,'+r+'/bin','-lmosek64','-o','solver.so'],check=True)
g=ca.external('solver','./solver.so')
s(**args)
print('CASADI',ca.__file__,flush=True)
for i in range(4):
 a=s(**args)['x'].full().ravel(); b=g(**args)['x'].full().ravel()
 print(i,'native', a,flush=True)
 print(i,'maxdiff',max(abs(a-b)), 'generated',b,flush=True)
