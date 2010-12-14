c
c    Library VRpack
c
c    A collection of functions for the saft-vr equation 
c    designed to be used in the saft_dft surface tension 
c    programs.
c
c    Original 24 January 2001 
c              GJ Gloor  Imperial College London
c
c************************************************************************
      real*8 function  eta_eff_sw(lambda,eta)
C************************************************************************
C
C  returns eta effective for a given lambda and eta
c
c
      implicit none
c
      real*8 lambda,eta,neff
      real*8 A(1:3,1:3)
      integer i,j
C
c     
      a(1,1)=2.25855d0
      a(1,2)=-1.50349d0
      a(1,3)=0.249434d0
      a(2,1)=-0.669270d0
      a(2,2)=1.40049d0
      a(2,3)=-0.827739d0
      a(3,1)=10.1576d0
      a(3,2)=-15.0427d0
      a(3,3)=5.30827d0
c
      neff = 0
      do 10 i=1,3
       do 20 j=1,3
        neff= neff+a(i,j)*lambda**(j-1)*eta**i
20     continue 
10    continue
      eta_eff_sw=neff
      return
      end
c*************************************************************************
      real*8 function deta_eff_dn(lambda,eta)
c**************************************************************************
      implicit none
      real*8 lambda,eta,neff
      real*8 A(1:3,1:3)
      integer i,j
C
c
      a(1,1)=2.25855
      a(1,2)=-1.50349
      a(1,3)=0.249434
      a(2,1)=-0.669270
      a(2,2)=1.40049
      a(2,3)=-0.827739
      a(3,1)=10.1576
      a(3,2)=-15.0427
      a(3,3)=5.30827
c
      neff = 0
      do 10 i=1,3
       do 20 j=1,3
        neff=neff+dfloat(i)*a(i,j)*lambda**(j-1)*eta**(i-1)
20     continue
10    continue
      deta_eff_dn=neff
      return
      end
c******************************************************************************
      real*8 function deta_eff_dlam(lambda,eta)
c*****************************************************************************
      implicit none
      real*8 lambda,eta,neff
      real*8 A(1:3,1:3)
      integer i,j
C
c
      a(1,1)=2.25855
      a(1,2)=-1.50349
      a(1,3)=0.249434
      a(2,1)=-0.669270
      a(2,2)=1.40049
      a(2,3)=-0.827739
      a(3,1)=10.1576
      a(3,2)=-15.0427
      a(3,3)=5.30827

      neff = 0
       do 20 i= 1,3
        neff= neff+(a(i,2)+2.d0*a(i,3)*lambda)*eta**i
20     continue
      deta_eff_dlam=neff
      return
      end
C******************************************************************************
      real*8 function d2eta_eff_dn2(lambda, eta)
c******************************************************************************
c
c     Calculates the second derivative with respect to eta of the
c     first term in the Barker-Henderson perturbation expanson
c     of the attractive energy for a square well fluid.
c
      implicit none 
      dimension a(1:3,1:3)
      real*8 a
      real*8 eta,lambda,epsilon
      real*8 cte1,cte2,cte3,cte4
      real*8 term1,term2,PI
      real*8 eta2,eta3,lambda2,lambda3
      Integer*2 j
      COMMON/pidata/pi
c
      a(1,1)=2.25855
      a(1,2)=-1.50349
      a(1,3)=0.249434
      a(2,1)=-0.669270
      a(2,2)=1.40049
      a(2,3)=-0.827739
      a(3,1)=10.1576
      a(3,2)=-15.0427
      a(3,3)=5.30827
c
      eta2=eta*eta
      eta3=eta2*eta
      lambda2=lambda*lambda
      lambda3=lambda2*lambda
c
      term1=0.0d0
      term2=0.0d0
      do 10 j=1,3
        term1=term1+a(2,j)*lambda**(j-1)
        term2=term2+a(3,j)*lambda**(j-1)
10    continue
c   
      d2eta_eff_dn2=(term1*2.d0+6.d0*term2*eta)
c
      return
      end
c******************************************************************************
      real*8 function d2eta_eff_dlamdn(lambda,eta)
c*****************************************************************************
      implicit none
      real*8 lambda,eta,neff
      real*8 A(1:3,1:3)
      integer i,j
C
c
      a(1,1)=2.25855
      a(1,2)=-1.50349
      a(1,3)=0.249434
      a(2,1)=-0.669270
      a(2,2)=1.40049
      a(2,3)=-0.827739
      a(3,1)=10.1576
      a(3,2)=-15.0427
      a(3,3)=5.30827

      neff = 0
       do 20 i= 1,3
        neff=neff+dfloat(i)*(a(i,2)+2.d0*a(i,3)*lambda)*eta**(i-1)
20     continue
      d2eta_eff_dlamdn=neff
      return
      end
C
c*************************************************************************
      real*8 function g_hs(eta)
c*************************************************************************
      implicit none
      real*8 eta
C
C     returns gHS(1;neff)
c     j. chem. phys. 106 4168
C
      g_hs=(1.0d0-0.5d0*eta)/(1.0d0-eta)**3.0d0
      return
      end 
c******************************************************************************
      real*8 function dg_hs_dn(eta)
c******************************************************************************
C     calculates the first derivative wrt eta of the 
c     hard sphere g
      implicit none 
      real*8 eta
c
      dg_hs_dn=(2.5d0-eta)/(1.0d0-eta)**4.0d0
      return
      end 
c******************************************************************************
      real*8 function d2g_hs_dn2(eta)
c******************************************************************************
      implicit none
      real*8 eta
      d2g_hs_dn2=(9.d0-3.d0*eta)/(1.d0-eta)**5.d0
      return
      end    
c******************************************************************************
      real*8 function g_sw(eta,lambda,epsilon)
******************************************************************************
c
c     calculates the square well g at contact using the saft vr model
      implicit none
      real*8 eta,lambda,da1_deta,da1_dlam
      real*8 term1,beta,epsilon,g_hs
      real*8 rho_l,rho_v,temp,pressure,cp_bulk
      common/thermo/rho_l,rho_v,temp,pressure,cp_bulk
      beta=1d0/temp
c
      term1=da1_deta(eta,lambda,epsilon)-lambda/3.0d0/eta*
     &       da1_dlam(eta,lambda,epsilon)
c
      g_sw=g_hs(eta)+0.25d0*beta*term1
c
      return
      end
c******************************************************************************
      real*8 function dg_sw_dn(eta,lambda,epsilon)
c******************************************************************************
      implicit none
      real*8 eta,lambda,epsilon
      real*8 rho_l,rho_v,temp,pressure,cp_bulk
      real*8 beta,term1,term2,term3
      real*8 da1_dlam,d2a1_dlamdn,d2a1_dn2,dg_hs_dn
c
      common/thermo/rho_l,rho_v,temp,pressure,cp_bulk
      beta=1.d0/temp
c
      term1=d2a1_dn2(eta,lambda,epsilon)
      term2=lambda/3.d0/eta/eta*da1_dlam(eta,lambda,epsilon)
      term3=lambda/3.d0/eta*d2a1_dlamdn(eta,lambda,epsilon)
c
      dg_sw_dn=dg_hs_dn(eta)+0.25d0*beta*(term1+term2-term3)
c
      return 
      end
c**************************************************************************
      real*8 function K_hs(eta)
c**************************************************************************
c
c  returns the compressibility of a hard sphere fluid using
c  the Percus- Yevek model      
c
      implicit none
      real*8  eta
      k_hs=(1.d0-eta)**4.0d0/(1.0d0+4.0d0*(eta+eta**2))
      return
      end
c**************************************************************************
      real*8 function dkhs_deta(eta)
c**************************************************************************
      implicit none
      real*8 eta,numerator,denominator,term1,term2
c 
      term1=(1.0d0-eta)**3.d0
      term2=(8.0d0+20.0d0*eta+8.0d0*eta**2.d0)
      numerator=term1*term2
      denominator=(1.0d0+4.0d0*eta+4.0d0*eta**2.d0)**2.d0
c
      dkhs_deta=-numerator/denominator
c
      end
c*******************************************************
      real*8 function a1_sw(eta,lambda,epsilon)
c********************************************************
      implicit none
      real*8 eta,lambda,epsilon,cte1
      real*8 g_hs,neff,eta_eff_sw
c
      neff=eta_eff_sw(lambda,eta)
      cte1=-4d0*eta*epsilon*(lambda**3d0-1d0)
c
      a1_sw=cte1*(g_hs(neff))
c
      return
      end
c***********************************************************
      real*8 function da1_deta(eta,lambda,epsilon)
c***********************************************************
      implicit none

      real*8 eta,lambda,epsilon, dg_hs_dn, deta_eff_dn 
      real*8 dg_deta,cte1,neff,g_hs,eta_eff_sw
      cte1=-4.0d0*epsilon*(lambda**3.d0-1.d0)
      neff=eta_eff_sw(lambda,eta)
c
      da1_deta=cte1*(g_hs(neff)+eta*dg_hs_dn(neff)*
     &         deta_eff_dn(lambda,eta))
c
      return
      end
c***********************************************************
      real*8 function d2a1_dn2(eta,lambda,epsilon)
c***********************************************************
      implicit none

      real*8 eta,lambda,epsilon
      real*8 dg_deta,cte1,neff, d2g_hs_dn2 
      real*8 dg_hs_dn,eta_eff_sw, deta_eff_dn, d2eta_eff_dn2
      cte1=-4.0d0*epsilon*(lambda**3.d0-1.d0)
      neff=eta_eff_sw(lambda,eta)
c
      d2a1_dn2=cte1*(2.0*dg_hs_dn(neff)
     &        *deta_eff_dn(lambda,eta)
     &      +d2g_hs_dn2(neff)
     &        *deta_eff_dn(lambda,eta)**2.0d0 *eta
     &      +eta*dg_hs_dn(neff)*d2eta_eff_dn2(lambda,eta))
c
      return
      end
c******************************************************************************
      real*8 function da1_dlam(eta,lambda,epsilon)
c******************************************************************************
      implicit none

      real*8 eta,lambda,epsilon
      real*8 g_hs,dg_hs_dn,deta_eff_dlam
      real*8 cte1,neff,eta_eff_sw
      real*8 term1,term2
      
c
      neff=eta_eff_sw(lambda,eta)
      cte1=-4.d0*epsilon*eta
      term1=3.d0*lambda**2.d0*g_hs(neff)
      term2=(lambda**3.d0-1.d0)*dg_hs_dn(neff)*deta_eff_dlam(lambda,eta)
c
      da1_dlam=cte1*(term1+term2)
c
      return 
      end
c****************************************************************************
      real*8 function d2a1_dlamdn(eta,lambda,epsilon)
c****************************************************************************
c
c     Calculates the second derivative with respect to lambda X eta of the 
c     first term of the Barker-Henderson perturbation expansion of 
c     the attractive energy for a square well fluid
c
      implicit none
      real*8 eta,lambda,epsilon
      real*8 term1,term2,term3,term4,term5
      real*8 cte1,cte2,neff,cte3
      real*8 g_hs,dg_hs_dn,dneff_dlam,dneff_deta
      real*8 d2eta_dn2,d2eta_eff_dlamdn,eta_eff_sw,deta_eff_dn
      real*8 deta_eff_dlam,d2g_hs_dn2 
c
      cte1=-4.0d0*epsilon
      cte2=lambda**3.d0-1.d0
      cte3=3.d0*lambda**2.d0
      neff=eta_eff_sw(lambda,eta)
c
      term1= g_hs(neff)+eta*dg_hs_dn(neff)*deta_eff_dn(lambda,eta)
      term2=dg_hs_dn(neff)*deta_eff_dlam(lambda,eta)
      term3=d2g_hs_dn2(neff)*deta_eff_dn(lambda,eta)
     &     *deta_eff_dlam(lambda,eta)*eta
      term4=eta*dg_hs_dn(neff)*d2eta_eff_dlamdn(lambda,eta)
c
      d2a1_dlamdn=cte1*(cte3*term1+cte2*(term2+term3+term4))
c
      return
      end
c**********************************************************
      real*8 function a2_sw(eta,lambda,epsilon)
c**********************************************************
      implicit none
      real*8 eta,lambda,epsilon
      real*8 da1_deta,k_hs
c
      a2_sw=0.5d0*epsilon*k_hs(eta)*eta
     &       *da1_deta(eta,lambda,epsilon)
c
      return
      end
c***********************************************************
      real*8 function dfmonoa1(eta,lambda,epsilon)
c***********************************************************
      implicit none
      real*8 eta,lambda,epsilon
      real*8 g_hs,dg_hs_dn,rho,dg_deta 
      real*8 cte1,term1,pi,neff,eta_eff_sw,deta_eff_dn 
      common/pidata/pi
c 
      neff=eta_eff_sw(lambda,eta)
      rho=6.0d0/pi*eta
      cte1= -4.d0*epsilon*(lambda**3.d0-1.d0)
c
      term1=eta*(2.d0*g_hs(neff)+ eta*dg_hs_dn(neff)
     &      *deta_eff_dn(lambda,eta)) 
c
      dfmonoa1=cte1*term1
c
      end
c************************************************************
      real*8 function dfmonoa2(eta,lambda,epsilon)
c***********************************************************
      implicit none
      real*8 eta,lambda,epsilon
      real*8 term1,term2,term3,eta_eff_sw,neff
      real*8 k_hs,da1_deta,d2a1_dn2,dkhs_deta
      real*8 deta_eff_dn
c
      neff=eta_eff_sw(lambda,eta)
      term1=K_hs(eta)*da1_deta(eta,lambda,epsilon)
      term2=dkhs_deta(eta)*da1_deta(eta,lambda,epsilon)*eta
c    & *deta_eff_dn(lambda,eta)
      term3=k_hs(eta)*d2a1_dn2(eta,lambda,epsilon)*eta
      dfmonoa2=epsilon*eta*(2.0d0*term1+term2+term3)/2.0d0
      end

