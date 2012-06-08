c This is a program I (David Roundy) wrote by very slightly modifying
c the code sent me by George Jackson, just so it would dump interesting
c and relevant data which I can compare to my other code, which tries to
c implement a functional based on his work.

       program eta_effective_vs_eta
 
c XXXXXXXXXXXXXXXXXXXXXXXXXXX
c
c       Program to calculate the interfacial properties of water-like
c       molecules. The program is written to solve the local
c       density functional theory proposed by Blas, Martin del Rio,
c       de Miguel and Jackson, based on the Wertheim's theory.
c
c       Water molecules are modeled as spherical hard-spheres with
c       four associating symmetric sites. Attractive interactions
c       are treated by means of the attractive part of the Yukawa
c       potential. In bulk, the attractive free energs is accounted for
c       a van der Waals mean-field term, where the constant 'a' is
c       given the bulk-integral of the attractive part of the Yukawa
c       potential.
c
c       This program is written in reduced units: the attractive
c       mean-field energy (epsilon_mf) and the molecular
c       volume (b=pi*sigma/6).
c
       implicit none
       real*8 sigma,epsilon_mf,m,epsilon_hb,k_hb,epsilon,lambda
       real*8 rho_l,rho_v,temp,pressure,cp_bulk
       real*8 zmin,zmax,deltaz
       real*8 pi
       real*8 rho_old,rho_new
       real*8 xcsfit
       real*8 tolerance,tol
       real*8 lambinv,lambinv2
       real*8 zi
       real*8 delta
       real*8 fraction
c
       real*8 v1,v2,v3,v4,v5
       real*8 eta
       real*8 rho_min,rho_max
       real*8 wmin,wmax,width
       integer i

c Here I define the types of some functions...
       real*8 a1_sw, eta_eff_sw
c
c
c       COMMONS
c
       common/parameters1/sigma,epsilon_mf,m
       common/parameters2/epsilon_hb,k_hb,epsilon,lambda
       common/thermo/rho_l,rho_v,temp,pressure,cp_bulk
       common/prof_param/zmin,zmax,deltaz
       common/pidata/pi
       common/profiles/rho_old,rho_new
       common/csfit/xcsfit
       common/bonding/fraction
c
c
       character*10 prof
       character*30 name
c
       open(unit=10,file='fortran/saft_dft.in',status='old')
c
       read(10,*)sigma,epsilon_mf,m
       read(10,*)epsilon_hb,k_hb,lambda
       read(10,*)rho_l,rho_v
       read(10,*)temp,pressure,cp_bulk

       do i=0,16
          eta = i/32d0
          write(*,"(2F15.10)") eta, eta_eff_sw(lambda, eta)
       end do
       stop
       end

c
c Here we include the little library
c
      include '../fortran/vrpack.f'
