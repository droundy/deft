cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real*8 function  f_assoc(epsilon_hb, k_hb,eta, radius, temp)
c
c       Calculate the contribution due to the association term in the
c       Euler-Lagrange equation for the equilibrium density profile and
c       the excess Helmholtz free energy due to the association for
c       evaluating the surface tension
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       implicit none
c
       integer ndmax,ndim
       parameter(ndmax=41000,ndim=ndmax-1)
c
       dimension fraction(0:ndmax+1)
c
       real*8 rho, radius
       real*8 pi
       real*8 sigma,epsilon_mf,epsilon_hb,k_hb,epsilon,lambda
       real*8 temp
       real*8 eta,eta2,eta3,cte2,cte3
       real*8 numerator,denominator
       real*8 f_hb
       real*8 ghs,del,x_frac,x_assoc,dg_hs
       real*8 term1,term2
       real*8 fraction
       real*8 g_sw,dg_sw_dn 
c
       integer ndata,nzinv
       integer i
       integer type,n_sites,att 
c
c       Define an useful constant
c
       pi= acos(-1.0d0)          
c     Work out the actual density in ordinary units
       rho = eta/(4d0*pi*radius**3d0/3d0)
c
c       Calculate the fraction of non-bonded water molecules associated
c       at a given site
c
       x_frac=x_assoc(epsilon_hb, k_hb,eta, radius, temp)
c
c       Calculate the excess Helmholtz free energy due to the association
c
       f_assoc=4d0*(dlog(x_frac)-0.5d0*x_frac + 0.5d0)
c       f_assoc=rho*f_assoc

       return
       end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real*8 function  x_assoc(epsilon_hb, k_hb,eta, radius, temp)
c
c       Calculate the contribution due to the association term in the
c       Euler-Lagrange equation for the equilibrium density profile and
c       the excess Helmholtz free energy due to the association for
c       evaluating the surface tension
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       implicit none
c
       integer ndmax,ndim
       parameter(ndmax=41000,ndim=ndmax-1)
c
       dimension fraction(0:ndmax+1)
c
       real*8 rho, radius
       real*8 pi
       real*8 sigma,epsilon_mf,epsilon_hb,k_hb,epsilon,lambda
       real*8 temp
       real*8 eta,eta2,eta3,cte2,cte3
       real*8 numerator,denominator
       real*8 f_hb
       real*8 ghs,del,x_frac,dg_hs
       real*8 term1,term2
       real*8 fraction
       real*8 g_sw,dg_sw_dn 
c
       integer ndata,nzinv
       integer i
       integer type,n_sites,att 
c
c       Define an useful constant
c
       pi= acos(-1.0d0)          
c
c       Calculate f_hb factor
       f_hb=exp(epsilon_hb/temp)-1.0d0
c
c       Calculate the pair radial distribution function of hard-spheres
c       at the contact length
c
       eta2=eta*eta
       eta3=eta2*eta
c
c       Calculate the delta integral for association
c
       del=k_hb*f_hb*g_sw(eta,lambda,epsilon)
c     Work out the actual density in ordinary units
       rho = eta/(4d0*pi*radius**3d0/3d0)
c
c       Calculate the fraction of non-bonded water molecules associated
c       at a given site
c
       x_frac=-1.0d0+dsqrt(1.0d0+8.0d0*rho*del)
       x_frac=x_frac/(4.0d0*rho*del)

       x_assoc = x_frac
       return
       end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real*8 function  delta(epsilon_hb, k_hb,eta, temp)
c
c       Calculate the contribution due to the association term in the
c       Euler-Lagrange equation for the equilibrium density profile and
c       the excess Helmholtz free energy due to the association for
c       evaluating the surface tension
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       implicit none
c
       real*8 f_hb
       real*8 sigma,epsilon_mf,m,epsilon_hb,k_hb,epsilon,lambda
       real*8 rho_l,rho_v,temp,pressure,cp_bulk
       real*8 eta,eta2,eta3,cte2,cte3
       real*8 numerator,denominator
       real*8 term1,term2
       real*8 fraction
       real*8 g_sw,dg_sw_dn 
c
       integer ndata,nzinv
       integer i
       integer type,n_sites,att 
c
       f_hb=exp(epsilon_hb/temp)-1.0d0
c
c       Calculate the delta integral for association
c
       delta=k_hb*f_hb*g_sw(eta,temp,lambda,epsilon)
       return
       end
