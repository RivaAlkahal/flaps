###############################################################################
import numpy as np
from constants import *
from basis_functions import *

###############################################################################

def project_pressure_on_Q2(NV,nel,mV,mP,xV,zV,p,iconP,iconV,debug,istep,solve_stokes):

    q=np.zeros(NV,dtype=np.float64)

    if solve_stokes:
 
       rVnodes=[-1, 1,1,-1, 0,1,0,-1,0]
       sVnodes=[-1,-1,1, 1,-1,0,1, 0,0]

       for iel in range(0,nel):
           for i in range(0,mV):
               inode=iconV[i,iel]
               rq=rVnodes[i]
               sq=sVnodes[i]
               NNNP=NNN(rq,sq,'Q1')
               q[inode]=np.dot(p[iconP[0:mP,iel]],NNNP[0:mP])
           #end for

           #q[iconV[6,iel]]=0.5*(p[iconP[2,iel]]+p[iconP[3,iel]])

       #end for

       print(spacing+" -> q (m,M) %e %e " %(np.min(q),np.max(q)))

       if debug:
          np.savetxt('q_'+str(istep)+'.ascii',np.array([xV,zV,q]).T,fmt='%1.4e')

    return q

###############################################################################
