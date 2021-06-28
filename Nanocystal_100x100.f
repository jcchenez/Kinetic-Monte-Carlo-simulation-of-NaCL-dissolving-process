      PROGRAM MAIN
      implicit none
      real*8  NAvoga,pi,delta
      real*8  Evib,Erot,Norm,k_Bolz,eV2Joule,Tn
      real*8  vel0,vel,dvel,alpha,integcoeff
      real*8  RTn_up,RTn_down,Ein,reactionProb
      real*8  R_EinTn,RTn,EnergyVel,H2mass,aLattice
      real*8  Rateconst,Pressure,Rmoleconst,en,TotRateconst
      real*8  temperature,timestep,s2fs
      real*8  xa,x,y,b,c,d,rx,symmetryfactor,prefactor,barrier
      integer v,j,i,m,n,k,atomtype,imax,unitcell,step,timeend,iatom
      integer totatomt
      character*8 chatomt0,chx 

c      parameter(symmetryfactor=4.00,imax=8,unitcell=4)
      
       parameter(symmetryfactor=4.00,imax=8,unitcell=100)




      parameter(NAvoga=6.02205E+23)
      parameter(Rmoleconst=8.31441)
      parameter(k_Bolz= 1.38066E-23)
      parameter(eV2Joule=1.60217E-19)
      parameter(H2mass=1.66054E-27*2.00)
      parameter(temperature=348.15)
      parameter(timestep=1.00)
      parameter(timeend=999999)
      parameter (delta=2.8500)
      parameter (s2fs=1.0d15)
      real*8  ratecoeff,delxsqrt2
      dimension atomtype(imax), prefactor(imax), barrier(imax),
     &          ratecoeff(imax) 
      real*8  xt0,yt0,zt0
      dimension xt0(unitcell),yt0(unitcell),zt0(unitcell),
     &              chatomt0(unitcell,unitcell,unitcell)
      character*8 charatomt(unitcell**3)
      character*8 charatomt2(unitcell**3)
      real*8  coordt(unitcell**3,3)
      real*8  coordt2(unitcell**3,3)
      real*8  dis,tottime
      real*8  time2propagation,totfreq 
      real*8   weight(imax)
      integer   typelabel(imax+1,unitcell**3)
      real*8  region(imax),setregion 
      real*8  randnumber1,randnumber2,randnumber3
      integer  bondnearest(unitcell**3) 
      integer  bondnearestnb(unitcell**3)
      integer  itypy2remov
      integer  delatom(timeend),idel
      integer  ncount(imax)
      integer  ncornerCl,nvacancyCl,nstepsiteCl,nsurfaceCl
      integer  ncornerNa,nvacancyNa,nstepsiteNa,nsurfaceNa
      integer  ninside,nothers,integ
      integer  nsorted, nkind

      open(11, file='traj.xyz',status='replace')

c      pi=acos(-1.00)
      call  Lecture(temperature,atomtype, barrier, ratecoeff)
      call  GenerateCrystal(unitcell,chatomt0,xt0,yt0,zt0,iatom)

C print rate constants
      DO i=1, imax
         print*, atomtype(i), barrier(i), ratecoeff(i)
      ENDDO


      step=0

      open(22, file='geoinput.xyz', status='old')
      read(22,*) integ
      read(22,*) chx     

      if(integ.ne.iatom) then
        print*,'error----------- reading'
c        #EXIT
      endif

      DO i=1, iatom
        read(22,*) charatomt(i),
     &             coordt(i,1),coordt(i,2),
     &             coordt(i,3)
c       write(33,*) iatom,charatomt(i),coordt(i,1),
c     &              coordt(i,2),coordt(i,3)
      ENDDO               

      totatomt=iatom
      tottime=0.0

      delxsqrt2=delta*sqrt(2.0)
c     4.03050865276332088908
c-------sorting the atom connections in the list
100   continue 
     
      ncornerCl=0
      nvacancyCl=0
      nstepsiteCl=0
      nsurfaceCl=0

      ncornerNa=0
      nvacancyNa=0
      nstepsiteNa=0
      nsurfaceNa=0


      ninside=0
      nothers=0
      
      DO i=1, totatomt
         bondnearest(i) =0
         bondnearestnb(i) =0
         nkind=99999

         DO j=1,totatomt

           if(abs(coordt(i,1)-coordt(j,1)).lt.4.50
     &   .and.abs(coordt(i,2)-coordt(j,2)).lt.4.50
     &   .and.abs(coordt(i,3)-coordt(j,3)).lt.4.50)
     &     then
            

           dis=sqrt(
     &      (coordt(i,1)-coordt(j,1))**2
     &     +(coordt(i,2)-coordt(j,2))**2 +
     &      (coordt(i,3)-coordt(j,3))**2 )
           else 
             goto 105
           endif  

           if( abs(dis-2.850).le.0.0010) then  
             bondnearest(i)=bondnearest(i)+1
           endif

           if(bondnearest(i).ge.6) then 
               goto 110 
           endif    

           if( abs(dis-delxsqrt2).le.0.0010 ) then
             bondnearestnb(i) =bondnearestnb(i) +1
           endif
105        continue           
         ENDDO
110      continue 

c         print*, i
c add periodic conditions for x,-x, y,-y, xy,-xy,-x-y,x-y

       if(bondnearest(i).eq.6
c    & .or.bondnearest(i).eq.5.and.bondnearestnb(i).gt.8
     &  ) then
           ninside=ninside+1
           typelabel(9,ninside)=i
           nkind=9
           goto 200
       endif

 
 
       if(bondnearest(i).eq.5
     & .and.bondnearestnb(i).eq.8.and.charatomt(i).eq.'Cl')
     & then
           nsurfaceCl=nsurfaceCl+1
           typelabel(1,nsurfaceCl)=i
           nkind=1 
           goto 200 
       endif

       if(bondnearest(i).eq.5
     & .and.bondnearestnb(i).eq.8.and.charatomt(i).eq.'Na')
     &  then
           nsurfaceNa=nsurfaceNa+1
           typelabel(2,nsurfaceNa)=i
           nkind=2
           goto 200
       endif

       if(bondnearest(i).eq.5
     & .and.bondnearestnb(i).le.7.and.charatomt(i).eq.'Cl')
     & then
           nvacancyCl=nvacancyCl+1
           typelabel(3,nvacancyCl)=i
           nkind=3
           goto 200
       endif

       if(bondnearest(i).eq.5
     & .and.bondnearestnb(i).le.7.and.charatomt(i).eq.'Na')
     &  then
            nvacancyNa=nvacancyNa+1
            typelabel(4,nvacancyNa)=i
            nkind=4
            goto 200
       endif

      if(bondnearest(i).eq.4.and.charatomt(i).eq.'Cl')
     & then
           nstepsiteCl=nstepsiteCl+1
           typelabel(5,nstepsiteCl)=i
           nkind=5
           goto 200
      endif
      if(bondnearest(i).eq.4.and.charatomt(i).eq.'Na')
     & then
           nstepsiteNa=nstepsiteNa+1
           typelabel(6,nstepsiteNa)=i
           nkind=6
           goto 200
      endif


       if(bondnearest(i).le.3.and.charatomt(i).eq.'Cl')
     &  then
           ncornerCl=ncornerCl+1
           typelabel(7,ncornerCl)=i
           nkind=7
           goto 200
       endif


       if(bondnearest(i).le.3.and.charatomt(i).eq.'Na')
     & then
           ncornerNa=ncornerNa+1
           typelabel(8,ncornerNa)=i
           nkind=8
           goto 200
       endif 
         
c      if (step.eq.0) then
c        print*,i,coordt(i,1),coordt(i,2),coordt(i,3),
c     &      bondnearest(i), bondnearestnb(i), step,nkind
c      endif
200   continue         
      ENDDO     
     
      nsorted=ncornerCl+ ncornerNa+nvacancyCl+nvacancyNa+ 
     &        nstepsiteCl+nstepsiteNa 
     &       +nsurfaceCl+nsurfaceNa+
     &        ninside
      print*,'sorting results :',nsorted
      print*, 'corner ',ncornerCl, ncornerNa 
      print*, 'vacancy',nvacancyCl,nvacancyNa
      print*, 'step   ',nstepsiteCl,nstepsiteNa
      print*, 'surface',nsurfaceCl,nsurfaceNa
      print*, 'inside and others',    ninside 

      ncount(1)=nsurfaceCl
      ncount(2)=nsurfaceNa
      ncount(3)=nvacancyCl
      ncount(4)=nvacancyNa
      ncount(5)=nstepsiteCl
      ncount(6)=nstepsiteNa
      ncount(7)=ncornerCl
      ncount(8)=ncornerNa


c=======================================
c          No.   barrier (eV)             rate constants (second -1)
c          1  0.780000000000000        21.7711036946430     
c          2  0.750000000000000        61.2623525601113     
c          3  0.260000000000000        779693348.276845     
c          4  0.100000000000000        186808185930.564     
c          5  0.250000000000000        1088130295.85611     
c          6  0.130000000000000        68726572289.6966     
c          7  0.170000000000000        14539611507.6294     
c          8  5.000000000000000E-002   944462327368.562     

c        sorting results :        1000 atoms
c   no.    corner     vacancy    bridge     surface     inside        others
c          8           0          96         384         512           0
c------------------------------------
c                   Cl          Na
c  corner            4           4
c  vacancy           0           0
c  step              48          48
c  surface           192         192
c  inside and others         512

c======================================= 
c----------start propagate the dynamics
c solving the master equation:
c    
      totfreq=0.0
      DO i=1,imax
         totfreq=totfreq+ratecoeff(i)*ncount(i)
      ENDDO
      time2propagation=1/totfreq
      tottime=tottime+time2propagation
      
 
      print*,'time2propagation&totaltime (-s)',time2propagation, 
     &          tottime, '  (-s)' 
      
      region=0.0
      setregion=0.0
      DO i=1,imax
         weight(i)=ratecoeff(i)*ncount(i)/totfreq
         setregion=setregion+weight(i) 
         region(i)=setregion
         print*,'weight of atom type:',i, weight(i), region(i)
      ENDDO  


      randnumber1=rand()
c      randnumber1=(rand()+rand())/2
c      randnumber1=(rand()+rand())/2

      print*, randnumber1 
c----------
c time2propagation  1.391379168589234E-013 (-s)
c weight of atom type:           1  5.816037150321096E-010
c weight of atom type:           2  1.636591894481438E-009
c weight of atom type:           3   0.00000000000000     
c weight of atom type:           4   0.00000000000000     
c weight of atom type:           5  7.267208766552160E-003
c weight of atom type:           6  0.458998660859645     
c weight of atom type:           7  8.092045028438357E-003
c weight of atom type:           8  0.525642083127169     
c----------
c------find which type of atom to be removed--------
      if(randnumber1.le.region(1)) then 
          itypy2remov=1
          goto 400
      endif

      if(randnumber1.gt.region(1).and.randnumber1.le.region(2))
     &then
        itypy2remov=2
        goto 400
      endif  

      if(randnumber1.gt.region(2).and.randnumber1.le.region(3))
     &then
        itypy2remov=3
        goto 400
      endif

      if(randnumber1.gt.region(3).and.randnumber1.le.region(4))
     &then
        itypy2remov=4
        goto 400
      endif  

      if(randnumber1.gt.region(4).and.randnumber1.le.region(5))
     & then
        itypy2remov=5
        goto 400
      endif  

      if(randnumber1.gt.region(5).and.randnumber1.le.region(6))
     & then
        itypy2remov=6
        goto 400
      endif  

      if(randnumber1.gt.region(6).and.randnumber1.le.region(7))
     & then
         itypy2remov=7
         goto 400
      endif   

      if(randnumber1.gt.region(7).and.randnumber1.le.region(8))
     &then
        itypy2remov=8
        goto 400
      endif 

400   continue      
       print*, 'Atom type to remove ',itypy2remov
        
     
c-----find which exact atom to be removed from the type found above
      randnumber2=rand()
      print*, randnumber2
      print*, 'no. in the type',ncount(itypy2remov)
      delatom(step)=int(randnumber2*ncount(itypy2remov))+1  
      print*, delatom(step)
      print*, 'atom to be deleted at step ', step, '  is: ',
     & typelabel(itypy2remov,delatom(step))
      idel=typelabel(itypy2remov,delatom(step))

c------------------------------------------------------
c generate the geometry for time step:    step+1 
c------------------------------------------------------       
c       DO i=1,imax+1
c       DO i=1, ncount(itypy2remov)
       print*, 'idel =', idel, totatomt    
       totatomt=totatomt-1
       step=step+1
c       DO i=1, idel-1
c          charatomt2(step,i)=charatomt(step-1,i)
c          DO  j=1,3
c           coordt2(step,i,j)=coordt(step-1,i,j)
c          ENDDO
c       ENDDO
       DO i=idel, totatomt
          charatomt2(i)=charatomt(i+1)
          DO  j=1,3
            coordt2(i,j)=coordt(i+1,j)
          ENDDO
       ENDDO
       DO i=idel, totatomt
           charatomt(i)=charatomt2(i)
           DO  j=1,3
             coordt(i,j)=coordt2(i,j)   
           ENDDO
       ENDDO    
             
       print*,'   '
       print*,'   '
       print*,'   '
       print*,'   '
       print*,'========================================'
       print*,'   '
       print*, totatomt
       print*, 'step = ',step, 'time = ',tottime*s2fs
     

        
       write(11,*) 'step = ',step, 'timestep=', time2propagation,
     &                     'tottime = ',tottime*s2fs
     

       if(mod(step,5000)==0) then       
       DO i=1,totatomt
          print*, charatomt(i),
     &     coordt(i,1),coordt(i,2),coordt(i,3)
       ENDDO
       endif

 
       if(step.le.timeend) goto 100







  




   
      
      end
c---------------------------------------------------------------
      SUBROUTINE GenerateCrystal(tot, chatom,x0,y0,z0,natom)
      implicit none

      integer i,j,k,l,m,n,layer,layerNaCl,layerwater,natom,tot
      real*8  x0,y0,z0,delta
      dimension x0(tot),y0(tot),z0(tot)
      parameter (delta=2.850, layerwater=0)
c      parameter( tot=8)
      character*8 charNa,CharCl,CharH,CharO,chatom
      dimension chatom(tot,tot,tot)
c      common /inputgeo/chatom,x0,y0,z0

      open(2, file='geoinput.xyz', status='replace')

      write(2,*)  '   ',tot*tot*tot
c     #          +      tot*tot*3*layerwater
      write(2,*)  'i =        0, time =        0.000, E =     0.0'

      natom=0
      DO layer = 1, tot
         z0(layer)=(layer-1)*delta
         DO i = 1, tot
            y0(i)=(i-1)*delta
            DO j = 1, tot
               x0(j)=(j-1)*delta 
               natom=natom+1
               if ( mod(layer-1,2)==0 .and.
     *              mod(i-1,2)==0  .and.
     *              mod(j-1,2)==0 ) then
                    chatom(j,i,layer)='Cl'
c                  print*, 'Cl  ', x(j), y(i), z(layer)
                  write(2,*) chatom(j,i,layer), x0(j), y0(i), z0(layer)
               elseif( mod(layer-1,2)==0 .and.
     *              mod(i-1,2)==1  .and.
     *              mod(j-1,2)==1 ) then
                    chatom(j,i,layer)='Cl'
c                  print*, 'Cl  ', x(j), y(i), z(layer)
                  write(2,*) chatom(j,i,layer), x0(j), y0(i), z0(layer)
               elseif( mod(layer-1,2)==1 .and.
     *              mod(i-1,2)==0  .and.
     *              mod(j-1,2)==1 ) then
                    chatom(j,i,layer)='Cl'
c                  print*, 'Cl  ', x(j), y(i), z(layer)
                  write(2,*) chatom(j,i,layer), x0(j), y0(i), z0(layer)
               elseif( mod(layer-1,2)==1 .and.
     *              mod(i-1,2)==1  .and. 
     *              mod(j-1,2)==0 ) then
                    chatom(j,i,layer)='Cl'
c                  print*, 'Cl  ', x(j), y(i), z(layer)
                  write(2,*) chatom(j,i,layer), x0(j), y0(i), z0(layer)
               else
                  chatom(j,i,layer)='Na'
c                  print*, 'Na  ', x(j), y(i), z(layer)   
                  write(2,*) chatom(j,i,layer), x0(j), y0(i), z0(layer)
               endif
            ENDDO
          ENDDO
      ENDDO

c      layer =tot
c      z0(layer)=(layer-1)*delta
c      DO i = 1, tot
c            y0(i)=(i-1)*delta
c            DO j = 1, int(tot/2)
c               x0(j)=(j-1)*delta
c               natom=natom+1
c               if ( mod(layer-1,2)==0 .and.
c     *              mod(i-1,2)==0  .and.
c     *              mod(j-1,2)==0 ) then
c                    chatom(j,i,layer)='Cl'
c                  print*, 'Cl  ', x(j), y(i), z(layer)
c                  write(2,*) chatom(j,i,layer), x0(j), y0(i), z0(layer)
c               elseif( mod(layer-1,2)==0 .and.
c     *              mod(i-1,2)==1  .and.
c     *              mod(j-1,2)==1 ) then
c                    chatom(j,i,layer)='Cl'
cc                  print*, 'Cl  ', x(j), y(i), z(layer)
c                  write(2,*) chatom(j,i,layer), x0(j), y0(i), z0(layer)
c               elseif( mod(layer-1,2)==1 .and.
c     *              mod(i-1,2)==0  .and.
c     *              mod(j-1,2)==1 ) then
c                    chatom(j,i,layer)='Cl'
cc                  print*, 'Cl  ', x(j), y(i), z(layer)
c                  write(2,*) chatom(j,i,layer), x0(j), y0(i), z0(layer)
c               elseif( mod(layer-1,2)==1 .and.
c     *              mod(i-1,2)==1  .and.
c     *              mod(j-1,2)==0 ) then
c                    chatom(j,i,layer)='Cl'
cc                  print*, 'Cl  ', x(j), y(i), z(layer)
c                  write(2,*) chatom(j,i,layer), x0(j), y0(i), z0(layer)
c               else
c                  chatom(j,i,layer)='Na'
cc                  print*, 'Na  ', x(j), y(i), z(layer)   
c                  write(2,*) chatom(j,i,layer), x0(j), y0(i), z0(layer)
c               endif
c        ENDDO
c      ENDDO

      
      close(2)
      end
c----------------------------------------------------------------------------


      SUBROUTINE Lecture(temp, atomtype, barrier, ratecoeff)
      implicit none
      include 'parameters.inc'
c      real*8  NAvoga,Rmoleconst,k_Bolz,eV2Joule,H2mass
c      implicit none
      integer i, j, k, imax, atomtype
      real*8   temp, prefactor, barrier,
     &                  ratecoeff, pi, viben
      parameter(imax=8)
      dimension atomtype(imax), prefactor(imax), barrier(imax),
     &          ratecoeff(imax)
      character*80 chara1
      character*8  chara2

c      common /setTEMP/temp
c      common /inputpara/ atomtype, barrier, ratecoeff
      pi=acos(-1.00)
      OPEN(10,FILE='input_improved_barrier.dat',STATUS='OLD')
      read(10,*) chara1

      DO i=1, imax
         read(10,*) atomtype(i), chara2, prefactor(i), barrier(i)
c         print*, atomtype(i), chara2, prefactor(i), barrier(i)
      ENDDO
      CLOSE(10)

      DO j=1, imax
          ratecoeff(j)=(1.0/(prefactor(j)*1.d-15))*
     &            exp((-barrier(j)*eV2Joule) / (k_Bolz*temp))

      Viben= (1.00/(prefactor(j)*1.d-15))*4.13558d-15
      print*,'Vib energy: ',Viben,' eV',Viben*8065.73,' cm-1'
      ENDDO

      END

c-------------------------------------------------------------------------

