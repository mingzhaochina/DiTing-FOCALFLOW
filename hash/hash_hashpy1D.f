c Sample main driver program for using the focal mechanism subroutines.  
c Uses P polarities and S/P amplitude ratios.

      include 'param.inc'
      include 'rot.inc'

c variables for storing earthquake input information  
      integer icusp                 ! event ID
      real qlat,qlon,qdep           ! location
      integer iyr,imon,idy,ihr,imn  ! origin time, year, month, day, hour, minute
      real qsec                     ! origin time, seconds
      real seh, sez                 ! location error, horizontal & vertical 
c
c variables for polarity information, input to HASH subroutines
      character*5 sname(npick0)                        ! station name
      character*5 iname                                ! station name
      character*1 pickpol                              ! polarity pick : U, u, +, D, d, or - ; onset, I or E
      integer p_pol(npick0),qp                         ! polarity pick (+/-1), and reversal (+/-1)
      real sp_ratio(npick0),spout(npick0)              ! log10(S/Pratio), S/P ratio
      real p_azi_mc(npick0,nmc0),p_the_mc(npick0,nmc0) ! azimuth and takeoff angle for each trial
      integer index(nmc0)                              ! index into velocity models, for each trial
      real qdep2(nmc0)                                 ! new depth, for each trail
      integer nmc                                      ! number of trails with different azimuth and take-off angles
      integer npol,nppl,nspr                           ! number of observations, P polarities, and S/P ratios
c
c variables for set of acceptable mechanisms, output of HASH subroutines
      integer nout2                                    ! number of acceptable mechanisms returned
      integer nmult                                    ! number of solutions (1 if no mulitples)
      real str_avg(5),dip_avg(5),rak_avg(5)            ! solution(s)
      real f1norm(3,nmax0),f2norm(3,nmax0)             ! normal vectors to the fault planes
      real strike2(nmax0),dip2(nmax0),rake2(nmax0)     ! strike, dip and rake
      real var_est(2,5),var_avg(5)                     ! variance of each plane, average
      real mfrac(5),stdr(5),mavg(5)                    ! fraction misfit polarities, station distribution       
      real prob(5)                                     ! probability true mechanism is "close" to preferred solution(s)
      character*1 qual(5)                              ! solution quality rating
      character*10 flag                                ! flag for best or perturbed take-off angle
c
c control parameters
      real del, delmax                                 ! maximum station distance
      real dang,dang2                                  ! grid angle to sample focal sphere
      integer maxout                                   ! max number of acceptable mechanisms output
      real badfrac                                     ! assumed rate of polarity error (fraction)
      real cangle                                      ! definition of "close" == 45 degrees
      real qbadfac                                     ! assumed noise in amplitude ratios, log10 (0.3 for factor of 2)
c
c file names
      character*100 outfile1,outfile2,outfile3,corfile,fpfile
      character*100 stfile,plfile,ampfile
      
      degrad=180./3.1415927
      rad=1./degrad
      
      print *,'Enter polarity and S/P ratio input file name from hashpy'       
      read (*,'(a)') fpfile
      open (11,file=fpfile, status='old')

      print *,'Enter station file name (NonLinLoc format)'       
      read (*,'(a)') stfile

      print *,'Enter output file name for preferred focal mechanisms'
      read (*,'(a)') outfile1
      open (13,file=outfile1)

      print *,'Enter output file name for other planes'
      read (*,'(a)') outfile2
      open (14,file=outfile2)

      print *,'Enter output file name for take-off angles'
      read (*,'(a)') outfile3
      open (15,file=outfile3)

      print *,'Enter grid angle for focal mech search, in degrees 
     &  (max ',dang0,')'
      read *,dang
      dang2=max(dang0,dang) ! don't do finer than dang0

      print *,'Enter number of trials (e.g., 30)'
      read *,nmc

      print *,'Enter maxout for focal mech. output (e.g., 500)'
      read *,maxout

      print *,'Enter maximum allowed source receiver distance'
      read *,delmax

      print *,'Enter number of polarities assumed bad'
      read *,nmismax

      print *,'Enter the assumed noise in amplitude ratios, log10  
     &  (e.g. 0.3 for a factor of 2)'
      read *,qbadfac

      print *,'Enter angle for computing mechanisms probability, 
     &         in degrees (e.g., 45)'
      read *,cangle

      print *,'Enter probability threshold for multiples (e.g., 0.1)'
      read *,prob_max

c make tables of takeoff angles for various velocity models
      ntab=10
      call MK_TABLE(ntab)

c read in earthquake location, etc. (hashpy format)
120   continue
      read (11,125,end=505) iyr,imon,idy,ihr,imn,qsec,qlat,
     &                qlon,qdep,seh,sez,icusp
125   format (i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,f5.2,1x,f10.6,1x,f10.6,1x,
     &        f9.6,1x,f5.2,1x,f5.2,1x,a46)

      aspect=cos(qlat/degrad)

c choose a new source location and velocity model for each trial 
      qdep2(1)=qdep
      index(1)=1
      do nm=2,nmc
        call RAN_NORM(val)
        qdep2(nm)=qdep+sez*val    ! randomly perturbed source depth
        index(nm)=mod(nm,ntab)+1  ! index used to choose velocity model
      end do

c read in polarities
      k=1
      nspr=0
130   continue
      read (11,*,end=140) iname,pickpol,qp,s2p
       if (iname.eq.'     ')  goto 140 ! end of data for this event
       sname(k)=iname
       spout(k)=s2p
c NonLinLoc station format
c        call GETSTAT_NLL(stfile,sname(k),
c     &               flat,flon,felv)
c Geena's station format
        call GETSTAT_GL(stfile,sname(k),
     &               flat,flon,felv)
        if (flat.eq.999.) go to 130
        dx=(flon-qlon)*111.2*aspect
        dy=(flat-qlat)*111.2
        range=sqrt(dx**2+dy**2)
        if (range.gt.delmax) go to 130
        qazi=90.-atan2(dy,dx)*degrad
        if (qazi.lt.0.) qazi=qazi+360.
        if (pickpol.eq.'U'.or.
     &                    pickpol.eq.'u'.or.
     &                    pickpol.eq.'+') then
          p_pol(k)=1
        else if (pickpol.eq.'D'.or.
     &                    pickpol.eq.'d'.or.
     &                    pickpol.eq.'-') then
          p_pol(k)=-1
        else
          goto 130
        end if
       if (s2p.ne.0) then
        sp_ratio(k)=log10(s2p)
        nspr=nspr+1
       end if
        do 105 nm=1,nmc  ! find azimuth and takeoff angle for each trial
          p_azi_mc(k,nm)=qazi
          call GET_TTS(index(nm),range,qdep2(nm),
     &                 p_the_mc(k,nm),iflag)
105     continue
        k=k+1
      goto 130
140   continue
      nppl=k-1

cc view polarity data
      do k=1,nppl
        print *,k,sname(k),p_azi_mc(k,1),p_the_mc(k,1),p_pol(k)
      end do

      if (nppl.lt.1) then
        print *,'*** warning - no p-wave polarity data for event',
     &            icusp
      end if
      if (nspr.lt.1) then
        print *,'*** warning - no s/p amplitude ratios for event',
     &            icusp
      end if
      
      nextra=0 !max(nint(nppl*badfrac*0.5),0)
      qmismax=nspr*qbadfac !max(nspr*qbadfac,2.0)                    
      qextra=max(nspr*qbadfac*0.5,2.0)

      call FOCALAMP_MC(p_azi_mc,p_the_mc,sp_ratio,p_pol,nppl,nmc,
     &    dang2,nmax0,nextra,nmismax,qextra,qmismax,nf2,strike2,dip2,
     &    rake2,f1norm,f2norm)
      nout2=min(nmax0,nf2)  ! number mechs returned from sub
      nout1=min(maxout,nf2)  ! number mechs to return
      
c find the probable mechanism from the set of acceptable solutions          
      call MECH_PROB(nout2,f1norm,f2norm,cangle,prob_max,nmult,
     &        str_avg,dip_avg,rak_avg,prob,var_est)           

      do 390 imult=1,nmult
      
      var_avg(imult)=(var_est(1,imult)+var_est(2,imult))/2.
      print *,'cid = ',icusp,imult,'  mech = ',
     &          str_avg(imult),dip_avg(imult),rak_avg(imult)

c find misfit for prefered solution
      call GET_MISF_AMP(npol,p_azi_mc,p_the_mc,sp_ratio,
     &      p_pol,str_avg(imult),dip_avg(imult),rak_avg(imult),
     &      mfrac(imult),mavg(imult),stdr(imult))
      
c solution quality rating, completely ad-hoc - make up your own!
      if ((prob(imult).gt.0.8).and.(var_avg(imult).le.25)) then
        qual(imult)='A'
      else if ((prob(imult).gt.0.6).and.(var_avg(imult).le.35)) then
        qual(imult)='B'
      else if ((prob(imult).gt.0.5).and.(var_avg(imult).le.45)) then
        qual(imult)='C'
      else
        qual(imult)='D'
      end if

390   continue

400   continue
       
      do i=1,nmult
      write (13,*) nint(str_avg(i)),nint(dip_avg(i)),nint(rak_avg(i)),
     & nint(var_avg(i)), qual(i)
c ,qual(i)
      end do
     
c output set of acceptable mechanisms
      do 500 ic=1,nout1
      write (14,550) strike2(ic),dip2(ic),rake2(ic),
     & f1norm(1,ic),f1norm(2,ic),f1norm(3,ic),
     & f2norm(1,ic),f2norm(2,ic),f2norm(3,ic)

550   format (f6.1,2x,f6.1,2x,f6.1,2x,
     & f8.5,2x,f8.5,2x,f8.5,2x,f8.5,2x,f8.5,2x,f8.5)

500   continue

c output takeoff angles
      do k=1,nppl
      do m=1,nmc
        if (m.eq.1) then
            flag='best_angle'
        else
            flag='pert_angle'
        end if
        write (15,560) sname(k),p_azi_mc(k,m),p_the_mc(k,m),p_pol(k),
     & spout(k),flag

560   format (a5,1x,f6.2,2x,f6.2,2x,i2,2x,f8.4,2x,a10)
      end do
      end do

      
505   continue
      close(11)
      close(12)
      close(13)
      close(14)
      stop
      end
