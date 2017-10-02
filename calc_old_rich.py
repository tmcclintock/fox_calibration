import numpy as np

h = 0.670435
#Msun/h
lMs = np.array([13.0, 13.1, 13.2, 13.4, 13.6, 13.8, 14.0, 14.2, 14.5])#, 15.0, 16.0])
Ms = 10**lMs/h #Msun

def rich(M):
    Mp = 10**14.371
    F = 1.12
    return 30. * (M/Mp)**(1./F)

lams = rich(Ms)
for lM, lam in zip(lMs, lams)[3:]:
    print "log10(M)=%.2f   lambda=%.0f"%(lM - np.log10(h), lam)
