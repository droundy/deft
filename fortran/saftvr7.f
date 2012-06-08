cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       program saftvr_dft
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
c       NOTE: the liquid phase is supossed to be at the LEFT side of
c       the interface in the coordinate system used in the program.
c
c       Last modif.: April 12th, 2000, Huelva (SPAIN), F. J. Blas
c
c       Last modif.: April 24th, 2000, Huelva (SPAIN), F. J. Blas
c
c       Last modif.: April 25th, 2000, Huelva (SPAIN), F. J. Blas
c            (first succesful compilation)
c
c       Last modif.: May 17th, 2000, Huelva (SPAIN), F. J. Blas
c            (density profiles are consistent)
c
c       Last modif.: May 22th, 2000, Huevla (SPAIN), F. J. Blas
c            (the program calculates the surface tension as well)
c
c       Last modif.: May 23th, 2000, Huelva (SPAIN), F. J. Blas
c            (the program can run continuously at different
c            temperatures along the LV coexistence curve)
c
c       Last modif.: May 25th, 2000, Huelva (SPAIN), F. J. Blas
c            (the program has been modified to calculate the
c             10-90 interfacial thickness)
c
c       Last modif.: July 28th, 2000, London (U.K.), F. J. Blas
c            (postdoctoral research at Imperial College of Science
c             technology and medicine under the supervision of
c             Dr. George Jackson)
c
c            (the program has been modified to account, in a local
c             way, for the chain contribution)
c
c       Last Modf.:  November 9th, 2000, London(U.K.) G.J. Gloor
c              ( addition of square well potential, changed DELTAZ to 
c                NZINV in input file to correct integer conversion
c                error on UNIX systems)
c
c
c      Last Modif.:   December 11, 2000, London(U.K.) G.J. Gloor
c              ( added corrections to surface tension for 
c                square well potential._)
c      
c      Last Modif:    February 19, 2001 London (U.K.) G.J. Gloor
c                ( added the energy and derivative terms for 
c                  the saft-vr model in file VRPACK.F, 
c                 changed name of program, added vr terms to 
c                 NEW_PROFILE and SURFACE_TENSION routines.
c
c
c      Last Modif:  March 1, 2001 London (U.K.)  G.J. Gloor
c                  ( added the vr formulation to the attractive term)
c                  
c
c      Last Modif: July 10, 2001 London (U.K.)  G.J. Gloor
c                  ( added newton ralphson routine in place of 
c                    bisection (XROOTN) in order to speed up convergence
c                    of rho_new(i) solution.
c
c      Last Modif: September 4, 2001 London (U.K.) G.J. Gloor
c                  ( added 2 more association models,
c                    2:1 model for linear primary alcohols and amines
c                    and 
c                    3:1 model for ammonia
c                    additionally added a routine to 'comb' the output
c                    of the profile files to reduce the disk space requirement)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       implicit none
c
       integer ndmax,ndim
       parameter(ndmax=41000,ndim=ndmax-1)
c
       dimension rho_old(0:ndmax+1),rho_new(0:ndmax+1)
       dimension xcsfit(4)
       dimension fraction(0:ndmax+1)
c
       dimension v1(100),v2(100),v3(100),v4(100),v5(100)
       dimension name(100)
c
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
       real*8 gamma
c
       real*8 v1,v2,v3,v4,v5
       real*8 rho_min,rho_max
       real*8 wmin,wmax,width
c
       integer ndata,nzinv
       integer menu,type,n_sites,att
       integer niter,niter_max
       integer i,j
c
       integer nruns
c
c       COMMONS
c
       common/parameters1/sigma,epsilon_mf,m
       common/parameters2/epsilon_hb,k_hb,epsilon,lambda
       common/thermo/rho_l,rho_v,temp,pressure,cp_bulk
       common/prof_param/zmin,zmax,deltaz
       common/pidata/pi
       common/integer1/ndata,nzinv
       common/profiles/rho_old,rho_new
       common/csfit/xcsfit
       common/bonding/fraction
       common/surf_tension/gamma
       common/syst/type,n_sites,att
c
c
       character*10 prof
       character*30 name
c
       open(unit=10,file='saft_dft.in',status='old')
c
       read(10,*)sigma,epsilon_mf,m
       read(10,*)epsilon_hb,k_hb,lambda
       read(10,*)rho_l,rho_v
       read(10,*)temp,pressure,cp_bulk
       read(10,*)zmin,zmax,nzinv
       read(10,*)niter_max,tolerance
       read(10,*)menu
       read(10,*)type
       read(10,*)att
c
       if((type.ne.1).and.(type.ne.2).and.(type.ne.3).
     &      and.(type.ne.4).and.(type.ne.5).and.
     &       (type.ne.6)) then
          write(*,*)'Problems with the selection'
          stop
       end if
c
c       Calculate the number of sites
c
       if(type.eq.1) then
          n_sites=0
       else if(type.eq.2) then
          n_sites=1
       else if(type.eq.3) then
          n_sites=2
       else if(type.eq.4) then
          n_sites=4
       else if(type.eq.5) then
          n_sites=3
       else if(type.eq.6) then
          n_sites=4
       end if
c
c       Check calculation type
c
       if((menu.eq.0).or.(menu.eq.1)) then
          nruns=1
          v1(1)=temp
          v2(1)=rho_l
          v3(1)=rho_v
          v4(1)=pressure
          v5(1)=cp_bulk
       end if
c
c       Read the file that contains the thermodynamic conditions
c       needed for solving the Euler-Lagrange equation (DFT)
c
       if(menu.eq.3) then
          open(unit=15,file='saft_dft.run',status='old')
          read(15,*)nruns
          do j=1,nruns
             read(15,100)name(j)
             read(15,*)v1(j),v2(j),v3(j),v4(j),v5(j)
          end do
       end if
c
       pi=acos(-1.0d0)
c
c       Calculation of the dispersive energy of the Yukawa potential
c
        If(att.eq.1) then
C        *** Yukawa potential   **
          lambinv=1.0d0/lambda
          lambinv2=lambinv*lambinv
          epsilon=epsilon_mf/(12.0d0*(lambinv+lambinv2))
        else
C        *** Square well **
          epsilon=epsilon_mf/(4.0d0*(lambda**3.d0-1.d0))
c           epsilon=epsilon_mf
        endif
c
        epsilon_hb=epsilon_hb*epsilon
c
c       Number of grid points
       deltaz=1.0d0/((dfloat(nzinv)))
c
       ndata=int((zmax-zmin)*nzinv) 
       write(6,*) ndata
c
       if(ndata.gt.ndmax)then
          write(*,*)'Problems with NDATA: too much grid points!!!'
          stop
       end if
c
c       nzinv=1+int(1.0d0/deltaz)
C        nzinv=20
c
c       Open files
c
      if (menu.ne.3)then
        open(unit=12,file='saft_dft.out',status='unknown')
        open(unit=13,file='saft_dft.bond',status='unknown')
      endif
        open(unit=14,file='saft_dft.surf',status='unknown')
        open(unit=15,file='saft_dft.iter',status='unknown')
c
       write(14,*)'#pwdpwdTemperature / Thickness / Surface Tension'
c
c       Loop over all the 'nruns' different thermodynamic conditions
c
       do j=1,nruns
c
          temp=v1(j)
          rho_l=v2(j)
          rho_v=v3(j)
          pressure=v4(j)
          cp_bulk=v5(j)
c
c       Boundary conditions over the density profile
c
          rho_old(0)=rho_l
          rho_old(ndata+1)=rho_v
c
c       Positions at which the cubic splines are calculated
c
          xcsfit(1)=-deltaz
          xcsfit(2)=0.d0
          xcsfit(3)=deltaz
          xcsfit(4)=2.0d0*deltaz
c
c       Initial profile from a tanh function
c
          if((menu.eq.0).or.(menu.eq.3))then
             call guesses
          end if
c
c       Read a previous equilibrium density profile
c
          if(menu.eq.1) then
             open(unit=11,file='saft_dft.out',status='old')
             do i=1,ndata
                read(11,*)zi,rho_old(i)
             end do
          end if
c
c       Calculate the bulk integrals
c
           if (att.eq.1) then
              call bulks_yk
           else
              call bulks_sw
           endif

c
c       Calculate the constants from the attractive integral
c
           if (att.eq.1) then
             call cte_phi_attr_yk
           else
             call cte_phi_attr_sw
           endif

c
c       Write the density profiles at different thermodynamic condtions
c
          if(menu.eq.3) then
             open(unit=19+j,file=name(j),status='unknown')
          end if
c
c       Iteration loop
c
          do niter=1,niter_max
c
c       HS contribution
c
             call f_hs
c
c       Chain contribution
c
            call f_chain
c
c       Association contribution (SAFT)
c
             call f_assoc
c
c       Attrative integral
c
             call phi_attractive
c
c       Calculate the new profile using an iterative procedure
c
             call new_profile
c
c       Check the convergence of the solution
c
             tol=tolerance
c
             do i=1,ndata
                delta=abs(rho_new(i)-rho_old(i))
                if(delta.gt.tol) then
                   tol=delta
                end if
             end do
c
             write(15,*)'Run:',j,' Iter.:',niter,' Max. dif.:',tol
             write(6,*)'Run:',j,' Iter.:',niter,' Max. dif.:',tol
c
             if(tol.le.tolerance) goto 9
c
c       Save the new configuration for next iteration
c
             do i=1,ndata
                rho_old(i)=rho_new(i)
             end do
c
          end do
c
c       No solution
c
          if (menu.ne.3) then
            write(12,*)'Number maximum of iterations exceeded!!!'
          endif
c
9         if (menu.ne.3) then 
            write(12,*)'Number of iterations:',niter
          endif
c
c       Write the final equilibrium density profile
c
          do i=1,ndata
             zi=zmin+float(i-1)*deltaz
             if(menu.eq.3) then
                if(mod(i,10).eq.0) then
                  write(19+j,*) zi,rho_new(i)
                endif
             else
                 write(12,*)zi,rho_new(i)
             end if
          end do
c
c       Write the fraction of non-bonded molecules at a given site
c       and the total fraction of non-bonded molecules
c
          if (menu.ne.3) then
          do i=1,ndata
            if(mod(i,10).eq.0) then
              zi=zmin+float(i-1)*deltaz
              if((type.eq.1).or.(type.eq.2))then
                 write(13,*)zi,fraction(i),fraction(i)
              else if(type.eq.3) then
                 write(13,*)zi,fraction(i),fraction(i)*2
              else if(type.eq.4) then
                 write(13,*)zi,fraction(i),fraction(i)*4
              end if
            endif
          end do
          endif
c
c       Calculate the interface width from the density profile
c       Actually, we use the definition of the 10-90 interfacial thickness
c       of the LV free interface
c
          rho_min=rho_v+0.10*(rho_l-rho_v)
          rho_max=rho_v+0.90*(rho_l-rho_v)
c
c       Set the counter 'i' equal to one
c
          i=1
c
c       Calculation of the upper border
c
10          zi=zmin+float(i-1)*deltaz
          if(rho_new(i).eq.rho_max) then
             wmax=zi
             go to 20
          end if
c
          if(rho_new(i).lt.rho_max) then
             wmax=zi-0.5d0*deltaz
             go to 20
          end if
c
c       Increase the counter
c
          i=i+1
          go to 10
c
c       Calculation of the lower border
c
20          zi=zmin+float(i-1)*deltaz
          if(rho_new(i).eq.rho_min) then
             wmin=zi
             goto 30
          end if
c
          if(rho_new(i).lt.rho_min) then
             wmin=zi-0.5*deltaz
             go to 30
          end if
c
c       Increase the counter
c
          i=i+1
          goto 20
c
30          width=wmin-wmax
c
c       Calculate the surface tension
c
          call surface_tension
c
c       Write the interfacial thickness and the surface tension
c
          write(14,*)temp,width,gamma
c
       end do
c
100       format(a30)
c
       stop
       end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       subroutine guesses
c
c       This subroutine generates the initial density profile
c       from the density of liquid and gas phases in bulk.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       implicit none
c
       integer ndmax,ndim
       parameter(ndmax=41000,ndim=ndmax-1)
c
       dimension rho_old(0:ndmax+1),rho_new(0:ndmax+1)
c
       real*8 rho_l,rho_v,temp,pressure,cp_bulk
       real*8 zmin,zmax,deltaz
       real*8 rho_old,rho_new
       real*8 cte1,cte2,zi
c
       integer ndata,nzinv
       integer i
c
c       COMMONS
c
       common/thermo/rho_l,rho_v,temp,pressure,cp_bulk
       common/integer1/ndata,nzinv
       common/prof_param/zmin,zmax,deltaz
       common/profiles/rho_old,rho_new
c
       cte1=0.5d0*(rho_l+rho_v)
       cte2=-0.5d0*(rho_l-rho_v)
c
c       Loop over the grid points
c
       do i=1,ndata
          zi=zmin+float(i-1)*deltaz
          rho_old(i)=cte1+cte2*tanh(zi)
          rho_new(i)=rho_old(i)
       end do
c
       return
       end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        subroutine bulks_yk
c
c       This subroutine calculates the attractive integrals performed
c       over the two bulk regions, from -oo (infinity) to zmin, and from
c       zmax to +oo (infinity).
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        implicit none
c
        integer ndmax,ndim
        parameter(ndmax=41000,ndim=ndmax-1)
c
        dimension bulk(ndmax)
c
       real*8 sigma,epsilon_mf,m,epsilon_hb,k_hb,epsilon,lambda
       real*8 zmin,zmax,deltaz
       real*8 bulk
       real*8 pi
       real*8 cte1,cte2,zi

c
        integer i
        integer ndata,nzinv
c
c       COMMONS
c
       common/parameters1/sigma,epsilon_mf,m
       common/parameters2/epsilon_hb,k_hb,epsilon,lambda
       common/integer1/ndata,nzinv
       common/prof_param/zmin,zmax,deltaz
       common/bulk_integral/bulk
       common/pidata/pi
c
c
c       Define some useful constants to calculate the integrals
c       associated to the bulk regions (liquid and vapor)
c
        cte1=(2.0d0*pi*epsilon)/(lambda*epsilon_mf)
        cte2=cte1/lambda
c
c       Loop over the grid points
c
        do i=1,ndata
           zi=dfloat(i-1)*deltaz
c
c       Condition over zi: if zi>1 uses phi_attr>;  if zi<1 uses phi_attr<
c
           if(zi.gt.1.0d0) then
              bulk(i)=-cte2*dexp(-lambda*(zi-1.0d0))
           else
              bulk(i)=-cte1*(1.0d0-zi)-cte2
           end if
        end do
c
        return
        end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        subroutine bulks_sw
c
c       This subroutine calculates the attractive integrals performed
c       over the two bulk regions, from -oo (infinity) to zmin, and from
c       zmax to +oo (infinity).
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        implicit none
c
        integer ndmax,ndim
        parameter(ndmax=41000,ndim=ndmax-1)
c
        dimension bulk(ndmax)
c
       real*8 sigma,epsilon_mf,m,epsilon_hb,k_hb,epsilon,lambda
       real*8 zmin,zmax,deltaz
       real*8 bulk
       real*8 pi
       real*8 cte1,cte2,zi
c
c
        integer i
        integer ndata,nzinv
c
c       COMMONS
c
c
       common/parameters1/sigma,epsilon_mf,m
       common/parameters2/epsilon_hb,k_hb,epsilon,lambda
       common/integer1/ndata,nzinv
       common/prof_param/zmin,zmax,deltaz
       common/bulk_integral/bulk
       common/pidata/pi
c
c       Define some useful constants to calculate the integrals
c       associated to the bulk regions (liquid and vapor)
c
        cte1=-pi*epsilon/epsilon_mf

c
c       Loop over the grid points
c
        do i=1,ndata
           zi=dfloat(i-1)*deltaz
c
c       Condition over zi: if zi>1 uses phi_attr>;  if zi<1 uses phi_attr<
c
           if(zi.lt.1.0d0) then
               bulk(i)=cte1*(-(lambda**2.0d0-1.0d0)*zi+2.0d0/3.0d0*
     &              lambda**3.0-2.0/3.0)

           else if((zi.ge.1.0d0).and.(zi.le.lambda)) then 

              bulk(i)=cte1*(2.0d0/3.0d0*lambda**3.0+zi**3.0/3.0-
     &            lambda**2.0d0*zi)
           else if (zi.gt.lambda) then

              bulk(i)=0.0d0
           end if

        end do
c
        return
        end
c

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        subroutine cte_phi_attr_yk
c
c       This subroutine calculates the constants multiplying
c       the coeffients coef(j) obtained from the cubic spline
c       interpolation of the density profile in each grid point
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        implicit none
c
       real*8 zmin,zmax,deltaz
       real*8 delt2,delt3,delt4
       real*8 sigma,epsilon_mf,m,epsilon_hb,k_hb,epsilon,lambda
       real*8 pi
       real*8 c1_intb,c2_intb,c3_intb,c4_intb
       real*8 cte1
c
c
        integer i,j,kalpha
        integer ndata,nzinv
c
c       COMMONS
c
       common/integer1/ndata,nzinv
       common/prof_param/zmin,zmax,deltaz
       common/deltas/delt2,delt3,delt4
       common/parameters1/sigma,epsilon_mf,m
       common/parameters2/epsilon_hb,k_hb,epsilon,lambda
       common/pidata/pi
       common/cintb/c1_intb,c2_intb,c3_intb,c4_intb
c
c
c       Definition of some useful constants
c
        delt2=deltaz*deltaz
        delt3=delt2*deltaz
        delt4=delt3*deltaz
c
        cte1=(2.0d0*pi*epsilon)/(lambda*epsilon_mf)
c
c       Constants that appear in integrals type-b (those that are
c       expressed in terms of phi_attr<)
c
c       NOTE: In our particular case, the constants (after the integrals
c       are performed) have no dependency with the indexes 'i' and 'j',
c       and due to that are numbers (and not vector matrix)
c
        c1_intb=-cte1*deltaz
        c2_intb=-cte1*delt2/2.0d0
        c3_intb=-cte1*delt3/3.0d0
        c4_intb=-cte1*delt4/4.0d0
c
c       The interval [zmin,zmax] is divided in three subintervals,
c       according to impose the conditions on the 'j' index
c
c       Values of zi close to the lower border (zmin)
c
        do i=1,nzinv+1
c
c       Loop over the 'j' index: points with j<(i+nzinv) are closer than
c       1 from zmin, and constants from integral type
c       b should be used (these constants are calculated above);
c       points with j>(i+nzinv-1) are in the opposite situation, and
c       since j>i, constants are calculated from integral type
c       a2
c
           do j=1,ndata-1
              kalpha=i-j-1
              if(j.ge.(i+nzinv)) then
                 call coef_inta2(kalpha)
              end if
           end do
        end do
c
c       Intermediate values of zi (not close to zmin neither zmax)
c
        do i=nzinv+2,ndata-nzinv-1

c
c       Loop over the 'j' index (similar comments that previous one)
c
           do j=1,ndata-1
              kalpha=i-j-1

              if(j.le.(i-nzinv-1)) then
                 call coef_inta1(kalpha)
              else if(j.ge.(i+nzinv)) then
                 call coef_inta2(kalpha)
              end if
           end do
        end do
c
c       Values of zi close to the upper border (zmax)
c
        do i=ndata-nzinv,ndata
           do j=1,ndata-1
              kalpha=i-j-1
              if(j.le.(i-nzinv-1)) then
                 call coef_inta1(kalpha)
              end if
           end do
        end do
c
        return
        end
c

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        subroutine cte_phi_attr_sw
c
c       This subroutine calculates the constants multiplying
c       the coeffients coef(j) obtained from the cubic spline
c       interpolation of the density profile in each grid point
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        implicit none
c
       real*8 zmin,zmax,deltaz
       real*8 delt2,delt3,delt4
       real*8 sigma,epsilon_mf,m,epsilon_hb,k_hb,epsilon,lambda
       real*8 pi
       real*8 c1_intb,c2_intb,c3_intb,c4_intb
       real*8 cte1

c
        integer i,j,kalpha
        integer ndata,nzinv,nlaminv
c
c       COMMONS
c
       common/integer1/ndata,nzinv
       common/prof_param/zmin,zmax,deltaz
       common/deltas/delt2,delt3,delt4
        common/parameters1/sigma,epsilon_mf,m
        common/parameters2/epsilon_hb,k_hb,epsilon,lambda
       common/pidata/pi
       common/cintb/c1_intb,c2_intb,c3_intb,c4_intb
c
c
c       Definition of some useful constants
c
        delt2=deltaz*deltaz
        delt3=delt2*deltaz
        delt4=delt3*deltaz
  
c
        cte1=(pi*epsilon)*((lambda**2d0-1.0d0)/epsilon_mf)
c
c       Constants that appear in integrals type-b (those that are
c       expressed in terms of phi_attr<)
c
c       NOTE: In our particular case, the constants (after the integrals
c       are performed) have no dependency with the indexes 'i' and 'j',
c       and due to that are numbers (and not vector matrix)
c
        c1_intb=-cte1*deltaz
        c2_intb=-cte1*delt2/2.0d0
        c3_intb=-cte1*delt3/3.0d0
        c4_intb=-cte1*delt4/4.0d0
c
c       The interval [zmin,zmax] is divided in three subintervals,
c       according to impose the conditions on the 'j' index
c
c       Values of zi close to the lower border (zmin)
        nlaminv=int(lambda*nzinv)
        write(6,*) nlaminv
c
        do i=1,nzinv+1
c
c       Loop over the 'j' index: points with j<(i+nzinv) are closer than
c       1 from zmin, and constants from integral type
c       b should be used (these constants are calculated above);
c       points with j>(i+nzinv-1) are in the opposite situation, and
c       since j>i, constants are calculated from integral type
c       a2
c
           do j=1,ndata-1
              kalpha=i-j-1
              if((j.ge.(i+nzinv)).and.(j.le.(i+nlaminv))) then
               call coef_inta2_sw(kalpha)
              end if
           end do
        end do
c
c       Intermediate values of zi (not close to zmin neither zmax)
c
        do i=nzinv+2,ndata-nzinv-1

c
c       Loop over the 'j' index (similar comments that previous one)
c
           do j=1,ndata-1
              kalpha=i-j-1

              if((j.le.(i-nzinv-1)).and.(j.ge.(i-nlaminv-1))) then
               call coef_inta1_sw(kalpha)
              else if((j.ge.(i+nzinv)).and.(j.le.(i+nlaminv))) then
               call coef_inta2_sw(kalpha)
              end if
           end do
        end do
c
c       Values of zi close to the upper border (zmax)
c
        do i=ndata-nzinv,ndata
           do j=1,ndata-1
              kalpha=i-j-1
              if((j.le.(i-nzinv-1)).and.(j.ge.(i-nlaminv-1))) then
              call coef_inta1_sw(kalpha)
              end if
           end do
        end do
c
        return
        end
c


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       subroutine coef_inta1(kalpha)
c
c       This subrotine calculates the constants multiplying the
c       coeficients of the cubic spline interpolation of the
c       density profile
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       implicit none
c
       integer ndmax,ndim
       parameter(ndmax=41000,ndim=ndmax-1)
c
       dimension c1_inta1(-ndim:ndim),c2_inta1(-ndim:ndim)
       dimension c3_inta1(-ndim:ndim),c4_inta1(-ndim:ndim)
c
       real*8 zmin,zmax,deltaz
       real*8 delt2,delt3,delt4
       real*8 sigma,epsilon_mf,m,epsilon_hb,k_hb,epsilon,lambda
       real*8 lambinv,lambinv2,lambinv3
       real*8 expolambda,expolambda1,expolambda2
       real*8 c1_inta1,c2_inta1,c3_inta1,c4_inta1
       real*8 pi
       real*8 cte1,alp,exp_alp
       real*8 c1inta1,c2inta1,c3inta1,c4inta1
c
       integer kalpha
c
c       COMMONS
c
       common/prof_param/zmin,zmax,deltaz
       common/deltas/delt2,delt3,delt4
        common/parameters1/sigma,epsilon_mf,m
        common/parameters2/epsilon_hb,k_hb,epsilon,lambda
       common/lambdas1/lambinv,lambinv2,lambinv3
       common/lambdas2/expolambda,expolambda1,expolambda2
       common/cinta1/c1_inta1,c2_inta1,c3_inta1,c4_inta1
       common/pidata/pi
c
c       Define some useful constants
c
       lambinv=1.0d0/lambda
       lambinv2=lambinv*lambinv
       lambinv3=lambinv2*lambinv
c
       expolambda=exp(lambda)
       expolambda1=exp(lambda*deltaz)
       expolambda2=1.0/expolambda1
c
       cte1=(2.0d0*pi*epsilon*lambinv2*expolambda)/epsilon_mf
c
c       alp is the variable 'a' defined in the theory as a=(i-j-1)*deltaz
c
       alp=float(kalpha+1)*deltaz
c
       exp_alp=exp(-lambda*alp)
c
c       Calculate the constants cj_inta1
c
       c1inta1=expolambda1-1.0d0
       c2inta1=expolambda1*(deltaz+lambinv*(-1.0d0+expolambda2))
       c3inta1=expolambda1*(delt2-2.0d0*deltaz*lambinv+
     &           2.0d0*lambinv2*(1.0d0-expolambda2))
       c4inta1=expolambda1*(delt3-3.0d0*delt2*lambinv+
     &           6.0d0*deltaz*lambinv2+6.0d0*lambinv3*(-1.0d0+
     &           expolambda2))
c
c       Multiply by the adequate normalization constants
c
       c1_inta1(kalpha)=-cte1*exp_alp*c1inta1
       c2_inta1(kalpha)=-cte1*exp_alp*c2inta1
       c3_inta1(kalpha)=-cte1*exp_alp*c3inta1
       c4_inta1(kalpha)=-cte1*exp_alp*c4inta1
c
       return
       end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       subroutine coef_inta2(kalpha)
c
c       This subrotine calculates the constants multiplying the
c       coeficients of the cubic spline interpolation of the
c       density profile
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       implicit none
c 
       integer ndmax,ndim
       parameter(ndmax=41000,ndim=ndmax-1)
c
       dimension c1_inta2(-ndim:ndim),c2_inta2(-ndim:ndim)
       dimension c3_inta2(-ndim:ndim),c4_inta2(-ndim:ndim)
c
       real*8 zmin,zmax,deltaz
       real*8 delt2,delt3,delt4
       real*8 sigma,epsilon_mf,m,epsilon_hb,k_hb,epsilon,lambda
       real*8 lambinv,lambinv2,lambinv3
       real*8 expolambda,expolambda1,expolambda2
       real*8 c1_inta2,c2_inta2,c3_inta2,c4_inta2
       real*8 pi
       real*8 cte1,alp,exp_alp
       real*8 c1inta2,c2inta2,c3inta2,c4inta2
c
       integer kalpha
c
c       COMMONS
c
       common/prof_param/zmin,zmax,deltaz
       common/deltas/delt2,delt3,delt4
        common/parameters1/sigma,epsilon_mf,m
        common/parameters2/epsilon_hb,k_hb,epsilon,lambda
       common/lambdas1/lambinv,lambinv2,lambinv3
       common/lambdas2/expolambda,expolambda1,expolambda2
       common/cinta2/c1_inta2,c2_inta2,c3_inta2,c4_inta2
       common/pidata/pi
c
c       Define some useful constants
c
        lambinv=1.0d0/lambda
        lambinv2=lambinv*lambinv
        lambinv3=lambinv2*lambinv
c
        expolambda=exp(lambda)
        expolambda1=exp(lambda*deltaz)
        expolambda2=1.d0/expolambda1
c
       cte1=(2.0d0*pi*epsilon*lambinv2*expolambda)/epsilon_mf
c
c       alp is the variable 'a' defined in the theory as a=(i-j-1)*deltaz
c
       alp=float(kalpha+1)*deltaz
c
       exp_alp=exp(lambda*alp)
c
c       Calculate the constants
c
       c1inta2=expolambda2-1.d0
       c2inta2=expolambda2*(deltaz+lambinv*(1.0d0-expolambda1))
       c3inta2=expolambda2*(delt2+2.0d0*deltaz*lambinv+
     &          2.0d0*lambinv2*(1.0d0-expolambda1))
       c4inta2=expolambda2*(delt3+3.0d0*delt2*lambinv+
     &           6.0d0*deltaz*lambinv2+6.0d0*lambinv3*(1.0d0-
     &           expolambda1))
c
c       Multiply by the adequate normalization constants
c
       c1_inta2(kalpha)=cte1*exp_alp*c1inta2
       c2_inta2(kalpha)=cte1*exp_alp*c2inta2
       c3_inta2(kalpha)=cte1*exp_alp*c3inta2
       c4_inta2(kalpha)=cte1*exp_alp*c4inta2
c 
       return
       end
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        subroutine coef_inta1_sw(kalpha)
c
c       This subrotine calculates the constants multiplying the
c       coeficients of the cubic spline interpolation of the
c       density profile
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        implicit none
c
        integer ndmax,ndim
        parameter(ndmax=41000,ndim=ndmax-1)
c
        dimension c1_inta1(-ndim:ndim),c2_inta1(-ndim:ndim)
        dimension c3_inta1(-ndim:ndim),c4_inta1(-ndim:ndim)
c
       real*8 zmin,zmax,deltaz
       real*8 delt2,delt3,delt4
       real*8 sigma,epsilon_mf,m,epsilon_hb,k_hb,epsilon,lambda
       real*8 lambinv,lambinv2,lambinv3
       real*8 expolambda,expolambda1,expolambda2
       real*8 c1_inta1,c2_inta1,c3_inta1,c4_inta1
       real*8 pi
       real*8 cte1,alp,exp_alp
       real*8 c1inta1,c2inta1,c3inta1,c4inta1
       real*8 alp2lam2 
c
c
        integer kalpha
c
c       COMMONS
c

       common/prof_param/zmin,zmax,deltaz
       common/deltas/delt2,delt3,delt4
        common/parameters1/sigma,epsilon_mf,m
        common/parameters2/epsilon_hb,k_hb,epsilon,lambda
       common/lambdas1/lambinv,lambinv2,lambinv3
       common/lambdas2/expolambda,expolambda1,expolambda2
       common/cinta1/c1_inta1,c2_inta1,c3_inta1,c4_inta1
       common/pidata/pi

c
c       Define some useful constants
c
        cte1=pi*epsilon/epsilon_mf
        delt2 = deltaz*deltaz
        delt3 = deltaz*delt2
        delt4 = deltaz*delt3


c       alp is the variable 'a' defined in the theory as a=(i-j-1)*deltaz
c
        alp=(dfloat(kalpha+1)*deltaz)
c
        alp2lam2=alp*alp-lambda*lambda
c
c       Calculate the constants cj_inta1
C
C        if (abs(alp).le.(lambda)) then
c       Multiply by the adequate normalization constants
c
        c1_inta1(kalpha)=cte1*deltaz*(alp2lam2-alp*
     &      deltaz+1.0d0/3.0d0*delt2)
        c2_inta1(kalpha)=cte1*delt2*(1.0d0/2.0d0*alp2lam2-2.0d0/3.0d0*
     +     alp*deltaz+1.0d0/4.0d0*delt2)
        c3_inta1(kalpha)=cte1*delt3*(1.0d0/3.0d0*alp2lam2-1.0d0/2.0d0*
     +     alp*deltaz+1.0d0/5.0d0*delt2)
        c4_inta1(kalpha)=cte1*delt4*(1.0d0/4.0d0*alp2lam2-2.0d0/5.0d0*
     +     alp*deltaz+1.0d0/6.0d0*delt2)

c        end if
c
        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        subroutine coef_inta2_sw(kalpha)
c
c       This subrotine calculates the constants multiplying the
c       coeficients of the cubic spline interpolation of the
c       density profile
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        implicit none
c
        integer ndmax,ndim
        parameter(ndmax=41000,ndim=ndmax-1)
c
        dimension c1_inta2(-ndim:ndim),c2_inta2(-ndim:ndim)
        dimension c3_inta2(-ndim:ndim),c4_inta2(-ndim:ndim)
c
       real*8 zmin,zmax,deltaz
       real*8 delt2,delt3,delt4
       real*8 sigma,epsilon_mf,m,epsilon_hb,k_hb,epsilon,lambda
       real*8 lambinv,lambinv2,lambinv3
       real*8 expolambda,expolambda1,expolambda2
       real*8 c1_inta2,c2_inta2,c3_inta2,c4_inta2
       real*8 pi
       real*8 cte1,alp,exp_alp
       real*8 c1inta2,c2inta2,c3inta2,c4inta2
       real*8 alp2lam2
c

c
        integer kalpha
c
c       COMMONS
c
       common/prof_param/zmin,zmax,deltaz
       common/deltas/delt2,delt3,delt4
        common/parameters1/sigma,epsilon_mf,m
        common/parameters2/epsilon_hb,k_hb,epsilon,lambda
       common/lambdas1/lambinv,lambinv2,lambinv3
       common/lambdas2/expolambda,expolambda1,expolambda2
       common/cinta2/c1_inta2,c2_inta2,c3_inta2,c4_inta2
       common/pidata/pi
c

c
c       Define some useful constants
c
        cte1=pi*epsilon/epsilon_mf
        delt2 = deltaz*deltaz
        delt3 = deltaz*delt2
        delt4 = deltaz*delt3


c       alp is the variable 'a' defined in the theory as a=(i-j-1)*deltaz
c
        alp=(float(kalpha+1)*deltaz)
c
        alp2lam2=alp*alp-lambda*lambda
c
c       Calculate the constants cj_inta1
c
c        if(abs(alp).le.(lambda)) then
c       Multiply by the adequate normalization constants
c
        c1_inta2(kalpha)=cte1*deltaz*(alp2lam2-alp*deltaz+1.0d0/3.0d0
     &      *delt2)
        c2_inta2(kalpha)=cte1*delt2*(1.0d0/2.0d0*alp2lam2-2.0d0/3.0d0*
     +    alp*deltaz+1.0d0/4.0d0*delt2)
        c3_inta2(kalpha)=cte1*delt3*(1.0d0/3.0d0*alp2lam2-1.0d0/2.0d0*
     +     alp*deltaz+1.0d0/5.0d0*delt2)
        c4_inta2(kalpha)=cte1*delt4*(1.0d0/4.0d0*alp2lam2-2.0d0/5.0d0*
     +    alp*deltaz+1.0d0/6.0d0*delt2)
         
c        end if
c
        return
        end
c

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       subroutine phi_attractive
c
c       This subroutine calculates the integral associated with
c       the attractive part of the intermolecular interaction
c       potential
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       implicit none
c
       integer ndmax,ndim
       parameter(ndmax=41000,ndim=ndmax-1)
c
       dimension phi_attr(ndmax)
       dimension a(1:25000)
       dimension rhoj(-1000:1000)
       dimension rho_old(0:ndmax+1),rho_new(0:ndmax+1)
c
       real*8 sigma,epsilon_mf,m,epsilon_hb,k_hb,epsilon,lambda
       real*8 rho_intb,rho_inta1,rho_inta2
       real*8 phi_attr
       real*8 rho_integral,bulkt
       real*8 rho_old,rho_new
       real*8 rho_l,rho_v,temp,pressure,cp_bulk
       real*8 a
       real*8 rhoj
       real*8 etaav,g_hs,simpson,p_att,pi,etaj
       real*8 eetai,eetaj,eetaav
c
       integer ndata,nzinv
       integer i,j,lmin,lmax,k
       integer index_bulk1,index_bulk2
       integer kalpha
       real*8 zmin,zmax,deltaz,etaefi,etaefj
       real*8 etaef,eta_eff_sw,gswav,eta,etai
       real*8 dfmonoa2,dmi,cte1
       real*8 phi,dg_hs_dn,deta_eff_dn,d2eta_eff_dn2
       real*8 d2g_hs_dn2,dkhs_deta,da1_deta,k_hs,a1
c       real*8 d2eta_eff_dn2
c
c       COMMONS
c
       common/parameters1/sigma,epsilon_mf,m
       common/parameters2/epsilon_hb,k_hb,epsilon,lambda
       common/integer1/ndata,nzinv
       common/rhoint/rho_intb,rho_inta1,rho_inta2
       common/phi/phi_attr
       common/thermo/rho_l,rho_v,temp,pressure,cp_bulk 
       common/prof_param/zmin,zmax,deltaz 
       common/profiles/rho_old,rho_new 
       pi= acos(-1.0d0)   
c
c       Set the attractive integral equal to zero
c
       do i=1,ndata
          phi_attr(i)=0.d0
       end do
       cte1=1.d0/(4.d0*(lambda*lambda*lambda-1))
c
c
       do i=1,ndata
         lmin=-int(lambda*float(nzinv))-1
	 lmax=+int(lambda*float(nzinv))+1
	 do j=lmin,lmax
	    k=j-lmin+1
	   if((i+j).le.1) then
	     rhoj(j)=rho_l
	   elseif ((i+j).ge.(ndata-1)) then
	      rhoj(j)=rho_v
	   elseif (((i+j).gt.1).and.((i+j).lt.(ndata-1))) then    
	      rhoj(j)=rho_old(i+j)
	   endif 
	   etaj=pi/6.d0*rhoj(j)
	   etai=pi/6.0d0*rho_old(i)  
	   etaav=(rhoj(j)+rho_old(i))/2.0d0*pi/6.0d0
	   eetai=eta_eff_sw(lambda,etai)
	   eetaj=eta_eff_sw(lambda,etaj)
	   eetaav=eta_eff_sw(lambda,etaav)
	   gswav=(g_hs(eetai)+g_hs(eetaj))/2.d0
c	   dmi=dfmonoa2(etai,lambda,epsilon)
	   
c	   if ((j).ne.0) then
	   phi=
     &   p_att(0.0d0,(dfloat(j)*deltaz),lambda,1.0d0,epsilon)
             a(k)=(g_hs(eetaav)*etaj+
     &        g_hs(eetaav)*etaj
     &        +dg_hs_dn(eetaav)*
     &        deta_eff_dn(lambda,etaav)*etaj*etai
     &        )
     &        *phi/2.d0
     
     
c         variation for a2sw:

cc
c     
c            a(k)=a1+0.25d0*(dkhs_deta(etai)*etai*a1*2.0d0*temp/cte1
c     &       +k_hs(etai)*a1*2.0d0*temp/cte1
c     &       +k_hs(etai)*etai*(2.d0*g_hs(eetaav)+
c     &       2.d0*dg_hs_dn(eetaav)*deta_eff_dn(lambda,etaav)*
c     &       (etaj+etai)+
c     &       (d2g_hs_dn2(eetaav)*deta_eff_dn(lambda,etaav)**2.d0+
c c    &        dg_hs_dn(eetaav)*d2eta_eff_dn2(lambda,etaav))*
c c    &       etai*etaj)*phi*temp/cte1+
c     &       dkhs_deta(etaj)*etaj*a1*2.0d0*temp/cte1
c     &       +k_hs(etaj)*a1*2.0d0*temp/cte1
c     &       +k_hs(etaj)*etaj*(2.d0*g_hs(eetaav)+
c     &       2.d0*dg_hs_dn(eetaav)*deta_eff_dn(lambda,etaav)*
c     &       (etaj+etai)+
c     &       (d2g_hs_dn2(eetaav)*deta_eff_dn(lambda,etaav)**2.d0+
c     &        dg_hs_dn(eetaav)*d2eta_eff_dn2(lambda,etaav))*
c     &       etai*etaj)*phi*temp/cte1)
c     &         /temp*cte1/2.d0
c             else 
c             a(k)=0.d0
c             endif
c     &   - (dfmonoa2(etai,lambda,epsilon)+
c     &     -dfmonoa2(etai,lambda,epsilon)
c     &       dmi)/2.0d0
c     &      /temp*phi*cte1/2.d0
c     &   * 0.5d0
	 enddo
	 phi_attr(i)=6.d0/pi*simpson((lmax-lmin),deltaz,a)

	enddo
	return
 	end   

c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       function bulkt(index1,index2)
c
       implicit none
c
       integer ndmax,ndim
       parameter(ndmax=41000,ndim=ndmax-1)
c
       dimension bulk(ndmax)
c 
       real*8 bulk,cte1,pi
       real*8 rho_l,rho_v,temp,pressure,cp_bulk
       real*8 bulkt
c
       integer index1,index2
c
c       COMMONS
c
       common/bulk_integral/bulk
       common/thermo/rho_l,rho_v,temp,pressure,cp_bulk
       common/pdata/pi 
c
       cte1=pi/6.0d0
       bulkt=rho_l*bulk(index1)+rho_v*bulk(index2)
    
c
       return
       end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       subroutine rho_integralb(i,j)
c
c       Calculate the intregral type b
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       implicit none
c
       integer ndmax,ndim
       parameter(ndmax=41000,ndim=ndmax-1)
c
       dimension xcsfit(4)
       dimension rho_old(0:ndmax+1),rho_new(0:ndmax+1)
       dimension y(4),coef(4)
c
       real*8 xcsfit
       real*8 c1_intb,c2_intb,c3_intb,c4_intb
       real*8 rho_old,rho_new
       real*8 rho_intb,rho_inta1,rho_inta2
       real*8 y,coef
c
       integer i,j
c
c       COMMONS
c
       common/csfit/xcsfit
       common/cintb/c1_intb,c2_intb,c3_intb,c4_intb
       common/profiles/rho_old,rho_new
       common/rhoint/rho_intb,rho_inta1,rho_inta2
c
c       Update the vector y(4) before calling the fitting subroutine
c
       y(1)=rho_old(j-1)
       y(2)=rho_old(j)
       y(3)=rho_old(j+1)
       y(4)=rho_old(j+2)
c
       call polcoef(xcsfit,y,4,coef)
c
c       Calculate the integral type b
c
       rho_intb=coef(1)*c1_intb+coef(2)*c2_intb+
     &           coef(3)*c3_intb+coef(4)*c4_intb
c
       return
       end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       subroutine rho_integrala1(i,j)
c 
c       Calculate the intregral type a1
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       implicit none
c
       integer ndmax,ndim
       parameter(ndmax=41000,ndim=ndmax-1)
c
       dimension xcsfit(4)
       dimension c1_inta1(-ndim:ndim),c2_inta1(-ndim:ndim)
       dimension c3_inta1(-ndim:ndim),c4_inta1(-ndim:ndim)
       dimension rho_old(0:ndmax+1),rho_new(0:ndmax+1)
       dimension y(4),coef(4)
c
       real*8 xcsfit
       real*8 c1_inta1,c2_inta1,c3_inta1,c4_inta1
       real*8 rho_old,rho_new
       real*8 rho_intb,rho_inta1,rho_inta2
       real*8 y,coef
c
       integer i,j
c
c       COMMONS
c
       common/csfit/xcsfit
       common/cinta1/c1_inta1,c2_inta1,c3_inta1,c4_inta1
       common/profiles/rho_old,rho_new
       common/rhoint/rho_intb,rho_inta1,rho_inta2
c
c       Update the vector y(4) before calling the fitting subroutine
c
       y(1)=rho_old(j-1)
       y(2)=rho_old(j)
       y(3)=rho_old(j+1)
       y(4)=rho_old(j+2)
c
       call polcoef(xcsfit,y,4,coef)
c
c       Calculate the integral type b
c
       rho_inta1=coef(1)*c1_inta1(i)+coef(2)*c2_inta1(i)+
     &            coef(3)*c3_inta1(i)+coef(4)*c4_inta1(i)
c
       return
       end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       subroutine rho_integrala2(i,j)
c
c       Calculate the intregral type a2
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       implicit none
c
       integer ndmax,ndim
       parameter(ndmax=41000,ndim=ndmax-1)
c
       dimension xcsfit(4)
       dimension c1_inta2(-ndim:ndim),c2_inta2(-ndim:ndim)
       dimension c3_inta2(-ndim:ndim),c4_inta2(-ndim:ndim)
       dimension rho_old(0:ndmax+1),rho_new(0:ndmax+1)
       dimension y(4),coef(4)
c
       real*8 xcsfit
       real*8 c1_inta2,c2_inta2,c3_inta2,c4_inta2
       real*8 rho_old,rho_new
       real*8 rho_intb,rho_inta1,rho_inta2
       real*8 y,coef
c
       integer i,j
c
c       COMMONS
c
       common/csfit/xcsfit
       common/cinta2/c1_inta2,c2_inta2,c3_inta2,c4_inta2
       common/profiles/rho_old,rho_new
       common/rhoint/rho_intb,rho_inta1,rho_inta2
c
c       Update the vector y(4) before calling the fitting subroutine
c
       y(1)=rho_old(j-1)
       y(2)=rho_old(j)
       y(3)=rho_old(j+1)
       y(4)=rho_old(j+2)
c
       call polcoef(xcsfit,y,4,coef)
c
c       Calculate the integral type a2
c
       rho_inta2=coef(1)*c1_inta2(i)+coef(2)*c2_inta2(i)+
     &          coef(3)*c3_inta2(i)+coef(4)*c4_inta2(i)
c
       return
       end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       subroutine f_hs
c
c       Calculate the contribution due to the hard-sphere term in the
c       Euler-Lagrange equation for the equilibrium density profile. It
c       also calculated the excess Helmholtz free energy of a HS
c       system in order to obtain the surface tension
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       implicit none
c
       integer ndmax,ndim
       parameter(ndmax=41000,ndim=ndmax-1)
c
       dimension rho_old(0:ndmax+1),rho_new(0:ndmax+1)
       dimension fhs(ndmax),df_hs(ndmax)
c
       real*8 rho_old,rho_new
       real*8 sigma,epsilon_mf,m,epsilon_hb,k_hb,epsilon,lambda
       real*8 fhs,df_hs
       real*8 pi
       real*8 cte1,eta,eta2,eta3
       real*8 numerator1,denominator1
       real*8 numerator2,denominator2
       real*8 eta_eff_sw,neff,a1,g_hs
c 
       integer ndata,nzinv
       integer i
c
c       COMMONS
c
       common/profiles/rho_old,rho_new
       common/parameters1/sigma,epsilon_mf,m
       common/parameters2/epsilon_hb,k_hb,epsilon,lambda
       common/hard_sphere/fhs,df_hs
       common/pidata/pi
       common/integer1/ndata,nzinv
c
c       Define useful constant
c
       cte1=pi/6.0d0
       do i=1,ndata
          eta=cte1*rho_old(i)
          eta2=eta*eta
          eta3=eta2*eta
c
          numerator1=8.0d0*eta-9.0d0*eta2+3.0d0*eta3
          denominator1=(1.0d0-eta)**3d0
c
          numerator2=4.0d0*eta-3.0d0*eta2
          denominator2=(1.0d0-eta)**2.d0
c
c       Calculate the excess Helmholtz free energy of HS system
c
          fhs(i)=rho_old(i)*(numerator2/denominator2)
c
c       Calculate the contribution to the Euler-Langrange equation
c
          df_hs(i)=(m*numerator1)/denominator1
c
       end do
c       print *,'rho=',rho_old(10),rho_old(11)
c       print *,'eta=',eta
c       print *,'temp*df_hs(10)=',0.09*df_hs(10)
c
       return
       end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
       subroutine f_chain
c
c       Calculate the contribution due to the chain formation in the
c       Euler-Lagrange equation for the equilibrium density profile
c       and the excess Helmholtz free energy due to the association
c       for evaluating the surface tension. The free energy functional
c       associated to this effect (connectivity of the chain) is
c       obtained using the local density approximation (LDA), as
c       in the HS term. This means that the functioanl is expressed
c       as the intregral, over the whole volume, of the excess
c       Helmholtz free energy of the uniform fluid with density
c       rho(r).
c
c
c       This routine has been modified for the saft-vr scheme
c       GG 25 jan 2001
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       implicit none
c
       integer ndmax,ndim
       parameter(ndmax=41000,ndim=ndmax-1)
c
       dimension rho_old(0:ndmax+1),rho_new(0:ndmax+1)
       dimension fchain(ndmax),df_chain(ndmax)
c
       real*8 rho_old,rho_new
       real*8 fchain,df_chain
       real*8 pi
       real*8 sigma,epsilon_mf,m,epsilon_hb,k_hb,epsilon,lambda
       real*8 cte1,eta,eta2,eta3
       real*8 numerator,denominator
       real*8 g_hs,dg_hs
       real*8 term1,term2
       real*8 g_sw,beta,d2a1_dn2,dn_drho
       real*8 rho_l,rho_v,temp,pressure,cp_bulk
       real*8 dg_sw_dn
c
       integer ndata,nzinv
       integer i
c
c       COMMONS
c
       common/profiles/rho_old,rho_new
       common/chain/fchain,df_chain
       common/pidata/pi
       common/parameters1/sigma,epsilon_mf,m
       common/parameters2/epsilon_hb,k_hb,epsilon,lambda
       common/integer1/ndata,nzinv
       common/thermo/rho_l,rho_v,temp,pressure,cp_bulk

c
c       Define useful constant
c
       cte1=pi/6.0d0
       beta=1d0/temp
c      beta=1d0
       dn_drho=cte1
c
       do i=1,ndata
          eta=cte1*rho_old(i)
          eta2=eta*eta
          eta3=eta2*eta
c
c
c       Calculate the radial distribution function of
c       HS at the contact length
c
c       to the chain formation at LDA level
c
         fchain(i)=(1.0d0-m)*(dlog(g_sw(eta,lambda,epsilon))-
     &    1.0d0/(4.0d0*temp*(lambda**3.0d0-1.0d0))) 
         fchain(i)=(rho_old(i)*fchain(i))/m
c
c       Calculate the contribution to the Euler-Lagrange equation
c
c
c       Derivative of the radial distribution function
c       of HS at the contact length with respect to
c       the packing fraction
c
c
c       Calculate the derivative of the Helmholtz free
c       energy due to the chain formation with respect
c       to the density profile
c
        df_chain(i)=(1.0d0-m)*(dlog(g_sw(eta,lambda,epsilon))+
     &    1.0d0/(g_sw(eta,lambda,epsilon))*
     &    dg_sw_dn(eta,lambda,epsilon)*eta-
     &    1.0d0*epsilon/temp)
c
       end do
c
c       print *,'df_chain(10)=',0.09*df_chain(10)
c
       return
       end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       subroutine f_assoc
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
       dimension rho_old(0:ndmax+1),rho_new(0:ndmax+1)
       dimension fassoc(ndmax),df_assoc(ndmax)
       dimension fraction(0:ndmax+1)
c
       real*8 rho_old,rho_new
       real*8 fassoc,df_assoc
       real*8 pi
       real*8 sigma,epsilon_mf,m,epsilon_hb,k_hb,epsilon,lambda
       real*8 rho_l,rho_v,temp,pressure,cp_bulk
       real*8 cte1,eta,eta2,eta3,cte_sw,cte2,cte3
       real*8 numerator,denominator
       real*8 f_hb
       real*8 ghs,delta,x_frac,dg_hs,ddelta
       real*8 term1,term2
       real*8 fraction
       real*8 g_sw,dg_sw_dn 
c
       integer ndata,nzinv
       integer i
       integer type,n_sites,att
c
c       COMMONS
 
 
 
c
       common/profiles/rho_old,rho_new
       common/association/fassoc,df_assoc
       common/pidata/pi
        common/parameters1/sigma,epsilon_mf,m
        common/parameters2/epsilon_hb,k_hb,epsilon,lambda
       common/integer1/ndata,nzinv
       common/thermo/rho_l,rho_v,temp,pressure,cp_bulk
       common/bonding/fraction
       common/syst/type,n_sites,att
c
c       Define an useful constant
c
       cte1=pi/6.0d0
c
c       Calculate f_hb factor
       cte_sw = 1.d0/(4.0d0*((lambda**3)-1.0d0))
       f_hb=exp(epsilon_hb/temp)-1.0d0
c
       do i=1,ndata
c
c       Calculate the pair radial distribution function of hard-spheres
c       at the contact length
c
          eta=cte1*rho_old(i)
          eta2=eta*eta
          eta3=eta2*eta
c
c       Calculate the delta integral for association
c
          delta=k_hb*f_hb*g_sw(eta,lambda,epsilon)
c
c       Calculate the fraction of non-bonded water molecules associated
c       at a given site
c
          if(type.eq.1) then
             x_frac=1.0d0
          else if((type.eq.2).or.(type.eq.3)) then
             x_frac=-1.0d0+dsqrt(1.0d0+4.0d0*(rho_old(i)/m)*delta)
             x_frac=x_frac/(2.0d0*(rho_old(i)/m)*delta)
          else if(type.eq.4) then
             x_frac=-1.0d0+dsqrt(1.0d0+8.0d0*(rho_old(i)/m)*delta)
             x_frac=x_frac/(4.0d0*(rho_old(i)/m)*delta)
          else if (type.eq.5) then
             cte2=1.0d0/rho_old(i)*m/delta
             cte3=cte2-1.0d0
             x_frac=(-cte3+dsqrt(cte3*cte3+8.0d0*cte2))/4.0d0
          else if(type.eq.6) then
             cte2=1.0d0/rho_old(i)*m/delta
             cte3=cte2-2.0d0
             x_frac=(-cte3+dsqrt(cte3*cte3+12.0d0*cte2))/6.0d0
          end if
c
c       Calculate the fraction of non-bonded molecules at a given
c       site
c
          fraction(i)=x_frac
c
c       Calculate the excess Helmholtz free energy due to the association
c
          if(type.le.4) then
            fassoc(i)=dfloat(n_sites)*(dlog(x_frac)-0.5d0*x_frac)+
     &  0.5d0* dfloat(n_sites)
            fassoc(i)=(rho_old(i)*fassoc(i)/m)
          else if(type.eq.5) then
            fassoc(i)=2.0d0*dlog(x_frac)+dlog(2.0d0*x_frac-1.0d0)-2.0d0*
     &  x_frac+2.0d0
            fassoc(i)=(rho_old(i)*fassoc(i)/m)
          else if(type.eq.6) then
            fassoc(i)=3.0d0*dlog(x_frac)+dlog(3.0d0*x_frac-2.0d0)-3.0d0*
     &  x_frac+3.0d0
            fassoc(i)=(rho_old(i)*fassoc(i)/m)
          end if
c
c       Calculate the derivative of the pair radial distribution function
c       with respect to the density profile
c
c
c       Calculate the derivative of the delta intregal with respect
c       to the density profile
c
         ddelta=k_hb*f_hb*dg_sw_dn(eta,lambda,epsilon)
c
c       Calculate the derivative of the excess Helmholtz free energy
c       due to the association with respect to the density profile
c
          term1=(m*fassoc(i))/rho_old(i)
        if (type.le.4) then
           term2=-0.5d0*dfloat(n_sites)*(1.0d0-x_frac)*
     &     (1.0d0+eta*ddelta/delta) 
        else if(type.eq.5) then
           term2=rho_old(i)/m*(2.0d0/x_frac+2.0d0
     &       /(2.0d0*x_frac-1.d0)-2.0d0)
     &     *(-x_frac*(2.d0*x_frac-1.d0)*(eta*ddelta+delta)/(1.d0+
     &       4.0d0*rho_old(i)/m*x_frac*delta-rho_old(i)/m*delta))
        else if(type.eq.6) then
           term2=rho_old(i)/m*(3.d0/x_frac+3.d0
     &       /(3.d0*x_frac-2.d0)-3.d0)
     &     *(-x_frac*(3.d0*x_frac-2.d0)*(eta*ddelta+delta)/(1.d0+
     &       6.d0*rho_old(i)/m*x_frac*delta-2.d0*rho_old(i)/m*delta))
        endif
c        Write(*,*) 'delta', delta,ddelta
        df_assoc(i)=(term1+term2)
         
c
       end do
c         write(*,*) fassoc(1),df_assoc(1)
       return
       end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       subroutine new_profile
c
c       Calculate the new density profile using, in a selfconsistent way,
c       the Euler-Lagrange equation for the equilibrium density profile. In
c       fact, this method is based on the Picard method for solving
c       algebraic equations by using an iterative method
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       implicit none
c
       integer ndmax,ndim
       parameter(ndmax=41000,ndim=ndmax-1)
c
       dimension rho_old(0:ndmax+1),rho_new(0:ndmax+1)
       dimension fhs(ndmax),fchain(ndmax),fassoc(ndmax)
       dimension df_hs(ndmax),df_chain(ndmax)
       dimension df_assoc(ndmax),phi_attr(ndmax)
c
       real*8 ghs,etaeff_sw
       real*8 rho_old,rho_new
       real*8 fhs,df_hs,fchain,df_chain,fassoc,df_assoc,phi_attr
       real*8 rho_l,rho_v,temp,pressure,cp_bulk
       real*8 sigma,epsilon_mf,m
       real*8 chpot,xrootn,neff,pi
       real*8 epsilon_hb,k_hb,epsilon,lambda
       real*8 eta,eta_eff_sw,g_hs,da1_deta,dg_hs_dn
       real*8 dfmonoa1,dfmonoa2,betat,cte1
       real*8 deta_eff_dn
c
       integer ndata,nzinv
       integer i
c
c
c       COMMONS
c
       common/profiles/rho_old,rho_new
       common/hard_sphere/fhs,df_hs
       common/chain/fchain,df_chain
       common/association/fassoc,df_assoc
       common/phi/phi_attr
       common/thermo/rho_l,rho_v,temp,pressure,cp_bulk
       common/integer1/ndata,nzinv
       common/parameters1/sigma,epsilon_mf,m
       common/parameters2/epsilon_hb,k_hb,epsilon,lambda
       common/pidata/pi

c
c       Loop over all the grid points
c
        betat=(4.0d0*(lambda**3.0d0-1.0d0))
c       do i=1,ndata
c       end do
        do i=1,ndata
c
c       The attractive contribution should by multiplied by 'm', the
c       chain length, because the equilibrium equations
c       have obtained by minimizing with respect to the density
c       profiles of chains, and all the terms in phi_attr are
c       given in terms of monomeric density
c
           eta=pi/6.0d0*rho_old(i)
          neff=eta_eff_sw(lambda,eta)
          chpot=cp_bulk
     &     -m*phi_attr(i)
c     &     -6.d0/pi*eta*eta*0.5d0*dg_hs_dn(neff)*
c     &        deta_eff_dn(lambda,eta)
     &    -temp*df_chain(i)
     &    -temp*df_assoc(i)
     &    - m*dfmonoa2(eta,lambda,epsilon)/temp
C           write(6,*)
c     &    - m*dfmonoa1(eta,lambda,epsilon)
C           pause
     
c     &    - 2.0d0*m*eta
c        if(mod(i,100).eq.0) then 
c         write((i+10),*) i,phi_attr(i),dfmonoa1(eta,lambda,epsilon)
c     &   ,dfmonoa2(eta,lambda,epsilon)/temp,g_hs(neff)
c        endif
c       write(*,*) dfmonoa1(eta,lambda,epsilon), dfmonoa2(eta,lambda,
c    &     epsilon)/temp,temp*df_hs(i), log(rho_old(i))
c       Find the solution of cp_hs(rho(i))=chpot
c
c         rho_new(i)=xroot(chpot,temp)
c          write(*,*) i
          rho_new(i)=xrootn(rho_old(i),chpot,temp)
c          rho_new(i)=rho_l*exp(excess_cp/temp-df_hs(i)-df_assoc(i)-
c     &      m* phi_attr(i)/temp)
c
c          writE(77,*)I,RHO_NEW(I)
c      write(*,*) temp*df_assoc(100), phi_attr(100)
c       pause
        if (I.eq.10) then
	 write(55,*) g_hs(neff)
	endif 

       end do
       open(unit=55,file='test3.dat',access='append')
       write(55,*) temp, phi_attr(10),phi_attr(ndata-10)
       eta=pi/6.0*rho_old(10)
       write(55,*) dfmonoa1(eta,lambda,epsilon),dfmonoa2(eta,
     & lambda,epsilon)
     
       
       close(55)
       return
       end
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       real*8 function XrootN(rho_old,chpot,temp)

       implicit none

       real*8 rho_old,chpot,temp
       real*8 rtn,dx,ctol,chp_hs,dchp_hs,xroot
       integer j,jmax


       jmax=50
       ctol=1.d-13

       rtn=rho_old
       do 80 j=1,jmax
         dx=(chp_hs(rtn,temp)-chpot)/dchp_hs(rtn,temp)
         rtn=rtn-dx 
         if (rtn.lt.1d-14)  rtn=1d-14
c       write(*,*) j,rtn,chpot,chp_hs(rtn,temp),dchp_hs(rtn,temp)
c       pause
         if (dabs(chp_hs(rtn,temp)-chpot).le.ctol) goto 90
80     continue
       write(*,*) 'jmax exceeded, try bisection'
       xrootn=xroot(chpot,temp)
90     xrootn=rtn
       return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       real*8 function dchp_hs(rho,temp)

       implicit none
       real*8 rho,eta,temp,pi
       real*8 sigma,epsilon_mf,m
       common/pidata/pi
       common/parameters1/sigma,epsilon_mf,m

       eta=pi/6.d0*rho
       dchp_hs=(m*2.d0*(4.d0-eta)/(1.d0-eta)**4
     & +1.d0/eta)*temp*pi/6.0d0
       return
       end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       function xroot(chpot,temp)
c
       implicit none
c
       real*8 xroot,chpot,temp
       real*8 x1,x2,xhalf,f1,f2,fhalf,prod,facc,chp_hs,ftol
c
       integer i,imax
c
c       We fix the the values that bracket all the solutions
       parameter(imax=200,ftol=1.d-10)
c
       xroot=0.0d0
       x1=1.0e-12
       x2=1.20
c
c
       f1=chp_hs(x1,temp)-chpot
       f2=chp_hs(x2,temp)-chpot
c
c       In any case check the brackets
c
       prod=f1*f2
       if(prod.gt.0) then
          print*,chpot
          write(6,*) f1,f2
          print *,'x1=',x1
          print *,'x2=',x2
          pause 'ERROR IN BISECTION'
       end if
c
       do i=1,imax
          xhalf=0.5d0*(x1+x2)
          fhalf=chp_hs(xhalf,temp)-chpot
c
c       check convergence
          if(dabs(fhalf).lt.ftol) then
             xroot=xhalf
             return
          end if
c
          prod=fhalf*f2
          if(prod.gt.0.d0) then
c
c       root in [x1,xh]
c
             x2=xhalf
             f2=fhalf
          else
c
c       root in [xh,x2]
c
             x1=xhalf
             f1=fhalf
          end if
c
       end do
c
       pause 'BISECTION EXCEEDING MAXIMUN ITERATIONS'
c
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       real*8 function chp_hs(rho,temp)
c
c       implicit none
c
       real*8 sigma,epsilon_mf,m,epsilon_hb,k_hb,epsilon,lambda
       real*8 pi
c
       real*8 rho,temp
       real*8 cte1
       real*8 eta,eta2,eta3
       real*8 numerator,denominator
c
c       COMMONS
c
       common/parameters1/sigma,epsilon_mf,m
       common/parameters2/epsilon_hb,k_hb,epsilon,lambda
       common/pidata/pi
c
       cte1=pi/6.0d0
c
        eta=cte1*rho
        eta2=eta*eta
        eta3=eta2*eta
c
        numerator=8.0d0*eta-9.0d0*eta2+3.0d0*eta3
        denominator=(1.0d0-eta)**3
c
c       Calculate the HS chemical potential
c
c       rho=0.65308743
c       eta=0.34195578
        chp_hs=temp*(dlog(eta/m)+(m*(numerator/denominator)))
c       print *,'rho=',rho
c       print *,'eta=',eta
c       print *,temp*log(rho/m)
c       print *,temp*(m*(numerator/denominator))
c
       return
       end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       subroutine surface_tension
c
c       This subroutine calculates the surface tension of the free
c       LV interface using the Simpson's rule
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       implicit none
c
       integer ndmax,ndim
       parameter(ndmax=41000,ndim=ndmax-1)
c
       dimension rho_old(0:ndmax+1),rho_new(0:ndmax+1)
       dimension fhs(ndmax),df_hs(ndmax)
       dimension fchain(ndmax),df_chain(ndmax)
       dimension fassoc(ndmax),df_assoc(ndmax)
       dimension phi_attr(ndmax),bulk(ndmax)
c
       dimension f(41000)
c
       real*8 rho_old,rho_new
       real*8 fhs,df_hs,fchain,df_chain,fassoc,df_assoc
       real*8 phi_attr,bulk
       real*8 rho_l,rho_v,temp,pressure,cp_bulk
       real*8 zmin,zmax,deltaz
       real*8 pi
       real*8 sigma,epsilon_mf,m,epsilon_hb,k_hb,epsilon,lambda
       real*8 lambinv,lambinv2,lambinv3
c
       real*8 f_ideal,f,simpson
       real*8 gamma,a1_sw,a2_sw,neff,eta_eff_sw,g_hs
       real*8 rho_l2,rho_v2,rho_lv,eta
       real*8 cte1,cte2,cte3,expo1
       real*8 correction1,correction2,correction3
       real*8 correction4,correction5
       real*8 a_one
c
       integer ndata,nzinv
       integer i,index
       integer type,n_sites,att
c
c       COMMONS
c
       common/profiles/rho_old,rho_new
       common/hard_sphere/fhs,df_hs
       common/chain/fchain,df_chain
       common/association/fassoc,df_assoc
       common/phi/phi_attr
       common/bulk_integral/bulk
       common/thermo/rho_l,rho_v,temp,pressure,cp_bulk
       common/prof_param/zmin,zmax,deltaz
       common/pidata/pi
       common/parameters1/sigma,epsilon_mf,m
       common/parameters2/epsilon_hb,k_hb,epsilon,lambda
       common/lambdas1/lambinv,lambinv2,lambinv3
       common/integer1/ndata,nzinv
       common/surf_tension/gamma
       common/syst/type,n_sites,att

c
       gamma=0.0d0
c
c       Transform the reduced pressure taken from the bulk program
c       (p*b/epsilon_mf) into reduced pressure with respect
c       the same molecular parameters used in
c       SAFT-DFT (p*sigma/epsilon_mf)
c
       pressure=6.0d0*pressure/pi 
c
c       Loop over the grid points
c
c       open(unit=99,file='test2.txt',status='unknown')
       cte1=4.00d0*(lambda**3-1.0d0)
       do i=1,ndata
c          write(*,*) i,ndata
          eta=pi/6.d0*rho_old(i)
          f_ideal= (rho_old(i)/m)*(log(eta/m)-1.0d0)
          neff=eta_eff_sw(lambda,eta)
          f(i)=temp*cte1*(f_ideal+fhs(i))
     &         +rho_old(i)*cte1*
     &         a2_sw(eta,lambda,epsilon)/temp
     &         + cte1*temp*fchain(i)
     &         +cte1*temp*fassoc(i) 
     &         -cte1*cp_bulk*(rho_old(i)/m)+cte1*pressure
     &         +0.50d0*rho_old(i)*cte1*a_one(i)
c     &         +cte1*rho_old(i)*eta
c
c       write(99,*) i,f(i)
        enddo
        open(unit=40,file='fassoc.txt',access='append')
        write(40,*) fassoc(1),fassoc(ndata-1)
        close(40)
c
c       Calculate the main contribution to the surface tension
c
       gamma=simpson(ndata,deltaz,f)
c
c       Calculate the corrections to the surface tension
c
c       Define some useful constants
c
c       cte1=pi*lambinv*epsilon*(0.5d0+lambinv+lambinv2)
       cte2=pi*lambinv3*epsilon
c
       rho_l2=rho_l*rho_l
       rho_v2=rho_v*rho_v
       rho_lv=rho_l*rho_v
c
       if(att.eq.1) then
c
c

ci        ***Yukawa Potential******
c
c

          expo1=exp(-lambda*(zmax-zmin-1.0d0))
c
c          cte1=pi*lambinv*epsilon*(0.5d0+lambinv+lambinv2)
          cte2=pi*lambinv3*epsilon

c         Correction 1
c 
          correction1=cte1*rho_l2
c
c         Correction 2
c
          do i=1,ndata
            f(i)=rho_old(i)*bulk(i)
          end do
c
          correction2=0.5d0*rho_l*simpson(ndata,deltaz,f)
c
c         Correction 3
c
          correction3=cte2*rho_lv*expo1
c
c         Correction 4
c
          do i=1,ndata
            index=ndata-i+1
            f(i)=rho_old(i)*bulk(index)
          end do
c
          correction4=0.5d0*rho_v*simpson(ndata,deltaz,f)
c
c          Correction 5
c
          correction5=cte1*rho_v2
c
          write(6,*) correction1,correction2,correction3,
     &               correction4, correction5,
     &               gamma

c
       else
          cte2=pi*(1.0/8.0*(lambda**4-1.0))
c
c         Correction 1
c
          correction1=rho_l2*cte2
c
c         Correction 2
c
          do i=1,ndata
            f(i)=rho_old(i)*bulk(i)
          end do
c
          correction2=cte1*0.5d0*rho_l*simpson(ndata,deltaz,f)
c
c          Correction 3
c
          cte3=zmin-zmax
          correction3=0.0d0
C ******** note: for the square well correction 3 must be 0
C          as long as zmax-zmin is greater than lambda which 
c           should always be the case.
c
c         Correction 4
c
          do i=1,ndata
            index=ndata-i+1
            f(i)=rho_old(i)*bulk(index)
          end do



          correction4=cte1*0.5d0*rho_v*simpson(ndata,deltaz,f)
c
c         Correction 5
c
          correction5=rho_v2*cte2
            cte1=pi/6.d0
c      correction1=a1_sw(cte1*rho_l,lambda,epsilon)/temp
c      correction2=a1_sw(cte1*rho_v,lambda,epsilon)/temp
c      correction4=a2_sw(cte1*rho_l,lambda,epsilon)/temp**2
c      correction5=a2_sw(cte1*rho_v,lambda,epsilon)/temp**2
          write(6,*) correction1,correction2,correction3,
     &               correction4, correction5,
     &               gamma
c
       endif

       gamma=gamma+correction1+correction2+2.0d0*correction3+
     &         correction4+correction5

    
       return
       end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       subroutine polcoef(x,y,n,cof)
c
c       Given arrays X and Y of length N containing a tabuled
c       function Yj = f(Xj), this routine returns an array of
c       coefficients COF, also of length N, such that:
c                         j-1
c       Yj = sum_j (COFj*Xj   )
c
c       Taken from Num. Rec. (1st. edition), pg.93
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       implicit none
c
c       Variables shared with the main program
c
       integer n
       real*8 x(n),y(n),cof(n)
c
       integer nmax,i,j,k
       parameter(nmax=10)
       real*8 s(nmax),phi,ff,b
c
       do i=1,n
          s(i)=0.0d0
          cof(i)=0.0d0
       end do
c
       s(n)=-x(1)
c
       do i=2,n
          do j=n+1-i,n-1
             s(j)=s(j)-x(i)*s(j+1)
          end do
          s(n)=s(n)-x(i)
       end do
c
       do j=1,n
          phi=float(n)
          do k=n-1,1,-1
             phi=float(k)*s(k+1)+x(j)*phi
          end do
c
          ff=y(j)/phi
          b=1.0d0
c
          do k=n,1,-1
             cof(k)=cof(k)+b*ff
             b=s(k)+x(j)*b
          end do
       end do
c
       return
       end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       real*8 function simpson(n,h,f)
c
c       Calculate de integral of f(n) using the alternative
c       extended Simpson's rule (Numerical Recipes, pg.108)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       implicit none
c
       dimension f(4100)
c
       real*8 h,f,sum
c
       integer n,i
c
c       Initialize the integral using the first and last four terms
c
       sum=(17.0d0*(f(1)+f(n))+59.0d0*(f(2)+f(n-1))+
     &      43.0d0*(f(3)+f(n-2))+49.0d0*(f(4)+f(n-3)))/48.0d0
c        sum=(5.d0*(f(1)+f(n))+13.d0*(f(2)+f(n-1)))/12.d0
c
c       Make the sum over the resting terms
c
       do i=5,n-4
          sum=sum+f(i)
       end do
c
c       Normalize by the increment
c
       simpson=h*sum
c
       return
       end
       
c---------------------------------------------------------------------------------
	real*8 function p_att(z1,z2,lambda,sigma,epsilon)
c 
c        Calculates the planar attractive potential at a
c        distance z1-z2	
c---------------------------------------------------------------------------------	       
	implicit none
	real*8 z1,z2,lambda,sigma,epsilon
	real*8 z12,pi
	pi=3.14159265358973d0
c	
	z12=dabs(z1-z2)
	if (z12.le.sigma) then
	  p_att=-pi*epsilon*(lambda*lambda-1.0d0)
	elseif ((z12.gt.sigma).and.(z12.le.lambda)) then
	  p_att=-pi*epsilon*(lambda*lambda-z12*z12)
	elseif (z12.gt.lambda) then
	  p_att=0.0d0
	endif
	return
	end
	
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       real*8 function a_one(i)
c
c       This subroutine calculates the integral associated with
c       the attractive part of the intermolecular interaction
c       potential
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       implicit none
c
       integer ndmax,ndim
       parameter(ndmax=41000,ndim=ndmax-1)
c
       dimension phi_attr(ndmax)
       dimension a(1:25000)
       dimension rhoj(-1000:1000)
       dimension rho_old(0:ndmax+1),rho_new(0:ndmax+1)
c
       real*8 sigma,epsilon_mf,m,epsilon_hb,k_hb,epsilon,lambda
       real*8 rho_intb,rho_inta1,rho_inta2
       real*8 phi_attr
       real*8 rho_integral,bulkt
       real*8 rho_old,rho_new
       real*8 rho_l,rho_v,temp,pressure,cp_bulk
       real*8 a
       real*8 rhoj
       real*8 etaav,g_hs,simpson,p_att,pi,etaj
       real*8 eetai,eetaj,eetaav
c
       integer ndata,nzinv
       integer i,j,lmin,lmax,k
       integer index_bulk1,index_bulk2
       integer kalpha
       real*8 zmin,zmax,deltaz,etaefi,etaefj
       real*8 etaef,eta_eff_sw,gswav,eta,etai
       real*8 dfmonoa2,dmi,cte1
       real*8 phi,dg_hs_dn,deta_eff_dn
c
c       COMMONS
c
       common/parameters1/sigma,epsilon_mf,m
       common/parameters2/epsilon_hb,k_hb,epsilon,lambda
       common/integer1/ndata,nzinv
       common/rhoint/rho_intb,rho_inta1,rho_inta2
       common/phi/phi_attr
       common/thermo/rho_l,rho_v,temp,pressure,cp_bulk 
       common/prof_param/zmin,zmax,deltaz 
       common/profiles/rho_old,rho_new 
       pi= acos(-1.0d0)   
c
c       Set the attractive integral equal to zero
c
      cte1=1.d0/(4.d0*(lambda*lambda*lambda-1))
c
           lmin=-int(lambda*float(nzinv))-1
	 lmax=+int(lambda*float(nzinv))+1

	 do j=lmin,lmax
	    k=j-lmin+1
	   if((i+j).le.1) then
	     rhoj(j)=rho_l
	   elseif ((i+j).ge.(ndata-1)) then
	      rhoj(j)=rho_v
	   elseif (((i+j).gt.1).and.((i+j).lt.(ndata-1))) then    
	      rhoj(j)=rho_old(i+j)
	   endif 
	   etaj=pi/6.d0*rhoj(j)
	   etai=pi/6.0d0*rho_old(i)  
	   etaav=(rhoj(j)+rho_old(i))/2.0d0*pi/6.0d0
	   eetai=eta_eff_sw(lambda,etai)
	   eetaj=eta_eff_sw(lambda,etaj)
	   eetaav=eta_eff_sw(lambda,etaav)
	   gswav=(g_hs(eetai)+g_hs(eetaj))/2.d0
c	   dmi=dfmonoa2(etai,lambda,epsilon)
	   
c	   if ((j).ne.0) then
	   phi=
     &   p_att(0.0d0,(dfloat(j)*deltaz),lambda,1.0d0,epsilon)
             a(k)=g_hs(eetaav)*etaj
     &        *phi
c     &        + g_hs(eetai)*phi-0.5d0*(dg_hs_dn(eetai)*
c     &        deta_eff_dn(lambda,etai))*etai*etai)/2.d0
c             else 
c             a(k)=0.d0
c             endif
cc     &   + (dfmonoa2(etai,lambda,epsilon)+
c     &     +pi/6.d0*(dfmonoa2(etaj,lambda,epsilon)
c     &       dmi)/2.0d0
c     &      /temp/(float(lmax-lmin)))
c     &   * 0.5d0
	 enddo
	 a_one=6.d0/pi*simpson((lmax-lmin),deltaz,a)
	return
 	end   

c	      

      include 'vrpack.f'

