###############################################################################
import numpy as np
from constants import *
import os
###############################################################################

def compute_dynamic_topography(output_dir,nel,NV,viscosity_nodal,viscosity_elemental,errc,err1,err2,err3,q,pc,solve_stokes,\
                               surfaceV,cmbV,inner_element,outer_element,innerQ2,outerQ2,istep,theta,thetac,\
                               uppermantle_rho,lowermantle_rho,rho_core,g0,bottom_p_avrg):

    dyn_topo_elemental=np.zeros(nel,dtype=np.float64)
    dyn_topo_nodal1=np.zeros(NV,dtype=np.float64)
    dyn_topo_nodal2=np.zeros(NV,dtype=np.float64)
    dyn_topo_nodal3=np.zeros(NV,dtype=np.float64)

    if solve_stokes:

       for i in range(0,NV):
           if surfaceV[i]:
              dyn_topo_nodal1[i]= -(2*viscosity_nodal[i]*err1[i]-q[i])/(uppermantle_rho*g0) 
              dyn_topo_nodal2[i]= -(2*viscosity_nodal[i]*err2[i]-q[i])/(uppermantle_rho*g0) 
              dyn_topo_nodal3[i]= -(2*viscosity_nodal[i]*err3[i]-q[i])/(uppermantle_rho*g0) 
           if cmbV[i]:
              dyn_topo_nodal1[i]= -(2*viscosity_nodal[i]*err1[i]-(q[i]-bottom_p_avrg))/((lowermantle_rho-rho_core)*g0) 
              dyn_topo_nodal2[i]= -(2*viscosity_nodal[i]*err2[i]-(q[i]-bottom_p_avrg))/((lowermantle_rho-rho_core)*g0) 
              dyn_topo_nodal3[i]= -(2*viscosity_nodal[i]*err3[i]-(q[i]-bottom_p_avrg))/((lowermantle_rho-rho_core)*g0) 

       for iel in range(0,nel):
           if outer_element[iel]:
              dyn_topo_elemental[iel]= -(2*viscosity_elemental[iel]*errc[iel]-pc[iel])/(uppermantle_rho*g0) 
           if inner_element[iel]:
              dyn_topo_elemental[iel]= -(2*viscosity_elemental[iel]*errc[iel]-(pc[iel]-bottom_p_avrg))/((lowermantle_rho-rho_core)*g0) 

       np.savetxt(os.path.join(output_dir,'DT1_R1_'+str(istep)+'.ascii'),np.array([theta[innerQ2],dyn_topo_nodal1[innerQ2]]).T,fmt='%1.16e')
       np.savetxt(os.path.join(output_dir,'DT1_R2_'+str(istep)+'.ascii'),np.array([theta[outerQ2],dyn_topo_nodal1[outerQ2]]).T,fmt='%1.16e')

       np.savetxt(os.path.join(output_dir,'DT2_R1_'+str(istep)+'.ascii'),np.array([theta[innerQ2],dyn_topo_nodal2[innerQ2]]).T,fmt='%1.16e')
       np.savetxt(os.path.join(output_dir,'DT2_R2_'+str(istep)+'.ascii'),np.array([theta[outerQ2],dyn_topo_nodal2[outerQ2]]).T,fmt='%1.16e')

       np.savetxt(os.path.join(output_dir,'DT3_R1_'+str(istep)+'.ascii'),np.array([theta[innerQ2],dyn_topo_nodal3[innerQ2]]).T,fmt='%1.16e')
       np.savetxt(os.path.join(output_dir,'DT3_R2_'+str(istep)+'.ascii'),np.array([theta[outerQ2],dyn_topo_nodal3[outerQ2]]).T,fmt='%1.16e')

       np.savetxt(os.path.join(output_dir,'DTc_R1_'+str(istep)+'.ascii'),np.array([thetac[inner_element],dyn_topo_elemental[inner_element]]).T,fmt='%1.16e')
       np.savetxt(os.path.join(output_dir,'DTc_R2_'+str(istep)+'.ascii'),np.array([thetac[outer_element],dyn_topo_elemental[outer_element]]).T,fmt='%1.16e')

       print(spacing+" -> dyn_topo_nodal1 surface (m,M) %e %e | nel= %d" %(np.min(dyn_topo_nodal1[surfaceV]),np.max(dyn_topo_nodal1[surfaceV]),nel))
       print(spacing+" -> dyn_topo_nodal2 surface (m,M) %e %e | nel= %d" %(np.min(dyn_topo_nodal2[surfaceV]),np.max(dyn_topo_nodal2[surfaceV]),nel))
       print(spacing+" -> dyn_topo_nodal3 surface (m,M) %e %e | nel= %d" %(np.min(dyn_topo_nodal3[surfaceV]),np.max(dyn_topo_nodal3[surfaceV]),nel))
       print(spacing+" -> dyn_topo_nodal1 cmb (m,M) %e %e | nel= %d" %(np.min(dyn_topo_nodal1[cmbV]),np.max(dyn_topo_nodal1[cmbV]),nel))
       print(spacing+" -> dyn_topo_nodal2 cmb (m,M) %e %e | nel= %d" %(np.min(dyn_topo_nodal2[cmbV]),np.max(dyn_topo_nodal2[cmbV]),nel))
       print(spacing+" -> dyn_topo_nodal3 cmb (m,M) %e %e | nel= %d" %(np.min(dyn_topo_nodal3[cmbV]),np.max(dyn_topo_nodal3[cmbV]),nel))

    return dyn_topo_elemental,dyn_topo_nodal1,dyn_topo_nodal2,dyn_topo_nodal3

###############################################################################
