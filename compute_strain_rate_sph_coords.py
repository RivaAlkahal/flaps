###############################################################################
import numpy as np
from constants import *
import os
###############################################################################
# Note that there is a pb when we run the model in plane strain. In that 
# case I set the err,ert,ett to zero for x<0 (t -> theta). See manual.

def compute_strain_rate_sph_coords(output_dir,solve_stokes,nel,NV,theta,thetac,istep,\
                                   inner_element,outer_element,innerQ2,outerQ2,\
                                   exx1,exx2,exx3,exxc,\
                                   exz1,exz2,exz3,exzc,\
                                   ezz1,ezz2,ezz3,ezzc):

    errc=np.zeros(nel,dtype=np.float64)  
    err1=np.zeros(NV,dtype=np.float64)  
    err2=np.zeros(NV,dtype=np.float64)  
    err3=np.zeros(NV,dtype=np.float64)  

    ettc=np.zeros(nel,dtype=np.float64)  
    ett1=np.zeros(NV,dtype=np.float64)  
    ett2=np.zeros(NV,dtype=np.float64)  
    ett3=np.zeros(NV,dtype=np.float64)  

    ertc=np.zeros(nel,dtype=np.float64)  
    ert1=np.zeros(NV,dtype=np.float64)  
    ert2=np.zeros(NV,dtype=np.float64)  
    ert3=np.zeros(NV,dtype=np.float64)  

    if solve_stokes:

       for i in range(0,NV):

           err1[i]=exx1[i]*np.sin(theta[i])**2+2*exz1[i]*np.sin(theta[i])*np.cos(theta[i])+ezz1[i]*np.cos(theta[i])**2
           err2[i]=exx2[i]*np.sin(theta[i])**2+2*exz2[i]*np.sin(theta[i])*np.cos(theta[i])+ezz2[i]*np.cos(theta[i])**2
           err3[i]=exx3[i]*np.sin(theta[i])**2+2*exz3[i]*np.sin(theta[i])*np.cos(theta[i])+ezz3[i]*np.cos(theta[i])**2

           ett1[i]=exx1[i]*np.cos(theta[i])**2-2*exz1[i]*np.sin(theta[i])*np.cos(theta[i])+ezz1[i]*np.sin(theta[i])**2
           ett2[i]=exx2[i]*np.cos(theta[i])**2-2*exz2[i]*np.sin(theta[i])*np.cos(theta[i])+ezz2[i]*np.sin(theta[i])**2
           ett3[i]=exx3[i]*np.cos(theta[i])**2-2*exz3[i]*np.sin(theta[i])*np.cos(theta[i])+ezz3[i]*np.sin(theta[i])**2

           ert1[i]=(exx1[i]-ezz1[i])*np.sin(theta[i])*np.cos(theta[i])+exz1[i]*(-np.sin(theta[i])**2+np.cos(theta[i])**2)
           ert2[i]=(exx2[i]-ezz2[i])*np.sin(theta[i])*np.cos(theta[i])+exz2[i]*(-np.sin(theta[i])**2+np.cos(theta[i])**2)
           ert3[i]=(exx3[i]-ezz3[i])*np.sin(theta[i])*np.cos(theta[i])+exz3[i]*(-np.sin(theta[i])**2+np.cos(theta[i])**2)

       for iel in range(0,nel):
           errc[iel]=exxc[iel]*np.sin(thetac[iel])**2+2*exzc[iel]*np.sin(thetac[iel])*np.cos(thetac[iel])+ezzc[iel]*np.cos(thetac[iel])**2
           ettc[iel]=exxc[iel]*np.cos(thetac[iel])**2-2*exzc[iel]*np.sin(thetac[iel])*np.cos(thetac[iel])+ezzc[iel]*np.sin(thetac[iel])**2
           ertc[iel]=(exxc[iel]-ezzc[iel])*np.sin(thetac[iel])*np.cos(thetac[iel])+exzc[iel]*(-np.sin(thetac[iel])**2+np.cos(thetac[iel])**2)

       np.savetxt(os.path.join(output_dir,'errc_R1_'+str(istep)+'.ascii'),np.array([thetac[inner_element],errc[inner_element]]).T,fmt='%1.4e')
       np.savetxt(os.path.join(output_dir,'errc_R2_'+str(istep)+'.ascii'),np.array([thetac[outer_element],errc[outer_element]]).T,fmt='%1.4e')
       np.savetxt(os.path.join(output_dir,'err1_R1_'+str(istep)+'.ascii'),np.array([theta[innerQ2],err1[innerQ2]]).T,fmt='%1.4e')
       np.savetxt(os.path.join(output_dir,'err1_R2_'+str(istep)+'.ascii'),np.array([theta[outerQ2],err1[outerQ2]]).T,fmt='%1.4e')
       np.savetxt(os.path.join(output_dir,'err2_R1_'+str(istep)+'.ascii'),np.array([theta[innerQ2],err2[innerQ2]]).T,fmt='%1.4e')
       np.savetxt(os.path.join(output_dir,'err2_R2_'+str(istep)+'.ascii'),np.array([theta[outerQ2],err2[outerQ2]]).T,fmt='%1.4e')
       np.savetxt(os.path.join(output_dir,'err3_R1_'+str(istep)+'.ascii'),np.array([theta[innerQ2],err3[innerQ2]]).T,fmt='%1.4e')
       np.savetxt(os.path.join(output_dir,'err3_R2_'+str(istep)+'.ascii'),np.array([theta[outerQ2],err3[outerQ2]]).T,fmt='%1.4e')

       print(spacing+" -> errc (m,M) %e %e | nel= %d" %(np.min(errc),np.max(errc),nel))
       print(spacing+" -> err1 (m,M) %e %e | nel= %d" %(np.min(err1),np.max(err1),nel))
       print(spacing+" -> err2 (m,M) %e %e | nel= %d" %(np.min(err2),np.max(err2),nel))
       print(spacing+" -> err3 (m,M) %e %e | nel= %d" %(np.min(err3),np.max(err3),nel))

       print(spacing+" -> ettc (m,M) %e %e | nel= %d" %(np.min(ettc),np.max(ettc),nel))
       print(spacing+" -> ett1 (m,M) %e %e | nel= %d" %(np.min(ett1),np.max(ett1),nel))
       print(spacing+" -> ett2 (m,M) %e %e | nel= %d" %(np.min(ett2),np.max(ett2),nel))
       print(spacing+" -> ett3 (m,M) %e %e | nel= %d" %(np.min(ett3),np.max(ett3),nel))

       print(spacing+" -> ertc (m,M) %e %e | nel= %d" %(np.min(ertc),np.max(ertc),nel))
       print(spacing+" -> ert1 (m,M) %e %e | nel= %d" %(np.min(ert1),np.max(ert1),nel))
       print(spacing+" -> ert2 (m,M) %e %e | nel= %d" %(np.min(ert2),np.max(ert2),nel))
       print(spacing+" -> ert3 (m,M) %e %e | nel= %d" %(np.min(ert3),np.max(ert3),nel))


    return err1,err2,err3,errc,ert1,ert2,ert3,ertc,ett1,ett2,ett3,ettc

###############################################################################
