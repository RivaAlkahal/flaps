import numpy as np
from constants import *

###############################################################################

def compute_elemental_pressure(nel,p,iconP,solve_stokes):

    pc=np.zeros(nel,dtype=np.float64)  

    if solve_stokes:
       for iel in range(0,nel):
           pc[iel]=(p[iconP[0,iel]]+p[iconP[1,iel]]+p[iconP[2,iel]]+p[iconP[3,iel]])/4

       print(spacing+" -> pc (m,M) %e %e " %(np.min(pc),np.max(pc)))

    return pc

###############################################################################
