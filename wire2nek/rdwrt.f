      subroutine read_rea

      include 'SIZE'

      character*80  line
      character*32  key, args(50)
      character*80  s80


      call blank (line,      80)
      call blank (args,   50*32)
      call rzero (param,    200)

C Open .rea parameters file
      open(unit=9,file=reaname_in,status='old',err=58)

C Initialize default parameter description
      call init_cparam

C Read in Parameters
      read(9,*) 
      read(9,*,err=59) vnekton
      read(9,*,err=59) ndimrea
      read(9,*,err=59) nparam_in
      do ip=1,nparam_in
        read(9,*,err=59) param(ip)
      enddo

      rewind(9)
      read(9,*) 
      read(9,*) 
      read(9,*) 
      read(9,*) 

      do ip = 1, nparam_in
        call blank(s80,80)
        read(9,'(19x,a80)',err=59) s80 
        if (nindx1(cparam(ip),' ',1).eq.0) then  !add parameter description if one does not exist 
          call chcopy(cparam(ip),s80,80)
        endif
      enddo

      nparam = max(nparam,nparam_in)
      npscal = param(23)  

      read(9,*,err=59) npsdata                  ! skip passive scalar data
      do i=1,npsdata
        read(9,*)
      enddo

      read(9,*,err=59) nlogic 

      if (nlogic.gt.0) then
        do ilogic=1,nlogic

          read(9,'(a80)',err=59) line
          call parse(line,' ',args,nargs)
          do iarg=1,nargs 
            key = args(iarg) 
            if (key(1:length(key)).eq.'IFFLOW') then
              read(args(1),*) ifflow
            elseif (key(1:length(key)).eq.'IFHEAT') then
              read(args(1),*) ifheat
            elseif (key(1:length(key)).eq.'IFTRAN') then
              read(args(1),*) iftran           
c            elseif (key(1:length(key)).eq.'IFADVC') then   !  IFADV=T and IFTMSH=F by default
c              do iiarg=1,npscal+2                          !  in write_rea_top
c                read(args(iiarg),*) ifadv(iiarg)           !
c              enddo                                        !
c            elseif (key(1:length(key)).eq.'IFTMSH') then   !
c              do iiarg=1,npscal+3                          !
c                read(args(iiarg),*) iftmsh(iiarg-1)        !
c              enddo
            elseif (key(1:length(key)).eq.'IFAXIS') then
              read(args(1),*) ifaxis
              ifaxiso = .true.
            elseif (key(1:length(key)).eq.'IFSTRS') then
              read(args(1),*) ifstrs
              ifstrso = .true.
            elseif (key(1:length(key)).eq.'IFSPLIT') then
              read(args(1),*) ifsplit
              ifsplito = .true.
            elseif (key(1:length(key)).eq.'IFMGRID') then
              read(args(1),*) ifmgrid
              ifmgrido = .true.
            elseif (key(1:length(key)).eq.'IFMODEL') then
              read(args(1),*) ifmodel
              ifmodelo = .true.
            elseif (key(1:length(key)).eq.'IFKEPS') then
              read(args(1),*) ifkeps
              ifkepso = .true.
            elseif (key(1:length(key)).eq.'IFMVBD') then
              read(args(1),*) ifmvbd
              ifmvbdo = .true.
            elseif (key(1:length(key)).eq.'IFCHAR') then
              read(args(1),*) ifchar
              ifcharo = .true.
            elseif (key(1:length(key)).eq.'IFLOMACH') then
              read(args(1),*) iflomach
              iflomacho = .true.
            elseif (key(1:length(key)).eq.'IFLOWMACH') then
              read(args(1),*) iflomach
              iflomacho = .true.
            elseif (key(1:length(key)).eq.'IFLOMA') then
              read(args(1),*) iflomach
              iflomacho = .true.
            endif
          enddo
        enddo
      endif

      read(9,'(a80)',err=59) line
      call parse(line,' ',args,nargs)
      read(args(1),*) xfac
      read(args(2),*) yfac
      read(args(3),*) xzero 
      read(args(4),*) yzero

      close(9)

      return

58    write(6,*) ' Error opening file ', reaname_in
59    write(6,*) ' Error reading parameters from file.'
      stop

      end

c-----------------------------------------------------------------------

      subroutine write_rea_top

      include 'SIZE'

      character*3 cbc3

      real  pcond(11), prhocp(11)

      character*10 my_fmt1, my_fmt2

C Modify logical switches, npscal and nflds according to input.exodus file
      ifflow = .true.
      ifheat = .true.
      npscal = 0             ! we overwrite npscal value from rea file
      nflds = 1

C Setup some custom write formats
      write(my_fmt1,'(a, i0, a)')  '(',npscal+2,'l2,a46)'   ! consistent with prenek
      write(my_fmt2,'(a, i0, a)')  '(',npscal+3,'l2,a43)'
c      write(my_fmt1,'(a, i0, a)')  '(',nflds,'l2,a46)'
c      write(my_fmt2,'(a, i0, a)')  '(',nflds+1,'l2,a43)'

C Check B.C.'s to set logical switches
      ifmvbd =    .false.
      iftmsh(0) = .false.
      nsides = 2*num_dim
      do iel=1,num_elem
        do iside=1,nsides 
          cbc3 = cbc(iside,iel,1)
          if (cbc3(1:1).eq.'M'.or.cbc3(1:1).eq.'m') ifmvbd = .true.
c         Test for moving mesh in solid
          cbc3 = cbc(iside,iel,2)
          if (cbc3(1:1).eq.'M'.or.cbc3(1:1).eq.'m') iftmsh(0) = .true.
        enddo
      enddo
      
      if (iflomach) ifsplit = .true.

C Write out parameter stuff
      write(10,*) '****** PARAMETERS *****'
      write(10,*) vnekton, ' NEKTON VERSION '
      write(10,*) num_dim, ' DIMENSIONAL RUN'
      write(10,*) nparam,  ' PARAMETERS FOLLOW'
      
      do ip=1,nparam
        l=length(cparam(ip)) 
        write(10,1049) param(ip), ip, cparam(ip)(1:l) 
      enddo
 1049 format(g14.6,1x,'p',i3.3,' ',a)

      write(10,*)'       4  Lines of passive scalar data follows:',
     & ' 2 CONDUCT; 2 RHOCP'
      call rone(pcond,11)
      call rone(prhocp,11)
      write(10,'(5g14.6)')(PCOND (I),I=3,11)
      write(10,'(5g14.6)')(PRHOCP(I),I=3,11)

      write(10,*) nlogic, '  LOGICAL SWITCHES FOLLOW'
      write(10,*) ifflow, '     IFFLOW'
      write(10,*) ifheat, '     IFHEAT'
      write(10,*) iftran, '     IFTRAN'
C Set IFADV=T and IFTMSH=F by default for all fields (always npscal+2)
      iftmsh(0) = .false.
      do i=1,npscal+2
        ifadv(i)  = .true.
        iftmsh(i) = .false.
      enddo
      write(10,my_fmt1)  (ifadv(i),i=1,npscal+2), 
     & '    IFNAV & IFADVC (convection in P.S. fields)'
      write(10,my_fmt2)  (iftmsh(i),i=0,npscal+2),        
     & '  IFTMSH (IF mesh for this field is T mesh)'
      
      if (ifaxiso)   write(10,*) ifaxis,   '     IFAXIS'
      if (ifstrso)   write(10,*) ifstrs,   '     IFSTRS'
      if (ifsplito)  write(10,*) ifsplit,  '     IFSPLIT'
      if (ifmgrido)  write(10,*) ifmgrid,  '     IFMGRID'
      if (ifmodelo)  write(10,*) ifmodel,  '     IFMODEL'
      if (ifkepso)   write(10,*) ifkeps,   '     IFKEPS'
      if (ifmvbdo)   write(10,*) ifmvbd,   '     IFMVBD'
      if (ifcharo)   write(10,*) ifchar,   '     IFCHAR'
      if (iflomacho) write(10,*) iflomach, '     IFLOMACH'

      write(10,'(4G14.6,'' XFAC,YFAC,XZERO,YZERO'')')
     &  xfac, yfac, xzero, yzero  


      return
      end

c-----------------------------------------------------------------------

      subroutine write_rea_bot

      include 'SIZE'


      write(10,*) '            0 PRESOLVE/RESTART OPTIONS *****'
      write(10,*) '            0        INITIAL CONDITIONS *****'
      write(10,*) '***** DRIVE FORCE DATA ***** BODY FORCE, FLOW, Q'
      write(10,*) '            1                 ',
     &                                'Lines of Drive force data follow'
      write(10,*) 'C'
      write(10,*) '***** Variable Property Data *****',
     &                                ' Overrrides Parameter data.'
      write(10,*) '            1 Lines follow.'
      write(10,*) '            0 PACKETS OF DATA FOLLOW'
      write(10,*) '***** HISTORY AND INTEGRAL DATA *****' 
      write(10,*) '            0   POINTS.  Hcode, I,J,H,IEL'

      nspecs = npscal + 5

      write(10,*) '***** OUTPUT FIELD SPECIFICATION *****'
      write(10,'(12x,i2,3x,a)') nspecs,'SPECIFICATIONS FOLLOW'
      write(10,*) 'T      COORDINATES'
      write(10,*) 'T      VELOCITY'
      write(10,*) 'T      PRESSURE'
      write(10,*) 'T      TEMPERATURE'
      write(10,*) 'F      TEMPERATURE GRADIENT'
      write(10,'(i2,3x,a)') npscal,'   PASSIVE SCALARS'
      do i=1,npscal
        write(10,'(a,i2)') ' T PS ', i 
      enddo
      write(10,*) '***** OBJECT SPECIFICATION *****'
      write(10,*) '     0 Surface Objects'
      write(10,*) '     0 Volume  Objects'
      write(10,*) '     0 Edge    Objects'
      write(10,*) '     0 Point   Objects'

      return
      end

c-----------------------------------------------------------------------

      subroutine init_cparam

      include 'SIZE'

      call blank(cparam,80*200)

      cparam( 1)='DENSITY'
      cparam( 2)='VISCOS'
      cparam( 3)='BETAG'
      cparam( 4)='GTHETA'
      cparam( 5)='PGRADX'
      cparam( 6)='FLOWRATE'
      cparam( 7)='RHOCP'
      cparam( 8)='CONDUCT'
      cparam( 9)='QVOL'
      cparam(10)='FINTIME'
      cparam(11)='NSTEPS'
      cparam(12)='DT'
      cparam(13)='IOCOMM'
      cparam(14)='IOTIME'
      cparam(15)='IOSTEP'
      cparam(16)='PSSOLVER: 0=default'
      cparam(17)='AXIS'
      cparam(18)='GRID < 0 --> # cells on screen'
      cparam(19)='INTYPE'
      cparam(20)='NORDER'
      cparam(21)='DIVERGENCE'
      cparam(22)='HELMHOLTZ'
      cparam(23)='NPSCAL'

      cparam(24)='TOLREL'
      cparam(25)='TOLABS'
      cparam(26)='COURANT/NTAU'
      cparam(27)='TORDER'
      cparam(28)='TORDER: mesh velocity (0: p28=p27)'

      cparam(29)='= magnetic visc if > 0, = -1/Rm if < 0'
      cparam(30)='> 0 ==> properties set in uservp()'
      cparam(31)='NPERT: #perturbation modes'
      cparam(32)='#BCs in re2 file, if > 0'


      cparam(41)='1-->multiplicative SEMG'
      cparam(42)='0=gmres/1=pcg'
      cparam(43)='0=semg/1=schwarz'
      cparam(44)='0=E-based/1=A-based prec.'
      cparam(45)='Relaxation factor for DTFS'
      cparam(46)='reserved'
      cparam(47)='vnu: mesh matieral prop.'
      cparam(49)='reserved'
      cparam(52)='IOHIS'

      cparam(54)='|p54|=1,2,3-->fixed flow rate dir=x,y,z'
      cparam(54)='fixed flow rate dir: |p54|=1,2,3=x,y,z'
      cparam(55)='vol.flow rate (p54>0) or Ubar (p54<0)'

      cparam(59)='!=0 --> full Jac. eval. for each el.'
      cparam(60)='!=0 --> init. velocity to small nonzero'

      cparam(62)='>0 --> force byte_swap for output'
      cparam(63)='=8 --> force 8-byte output'
      cparam(64)='=1 --> perturbation restart'

      cparam(65)='#iofiles (eg, 0 or 64); <0 --> sep. dirs'
      cparam(66)='output : <0=ascii, else binary'
      cparam(67)='restart: <0=ascii, else binary'
      cparam(68)='iastep: freq for avg_all'
      cparam(68)='iastep: freq for avg_all (0=iostep)'

      cparam(74)='verbose Helmholtz'


      cparam(84)='!=0 --> sets initial timestep if p12>0'
      cparam(85)='dt ratio if p84 !=0, for timesteps>0'
      cparam(86)='reserved' !=0 --> use skew-symm form (convection)'

      cparam(89)='reserved'

      cparam(93)='Number of previous pressure solns saved'
      cparam(94)='start projecting velocity after p94 step'
      cparam(95)='start projecting pressure after p95 step'

      cparam(99)='dealiasing: <0--> off/3--> old/4--> new'
      cparam(100)='reserved'
      cparam(101)='Number of additional modes to filter'
      cparam(102)='Dump out divergence at each time step'
      cparam(103)='weight of stabilizing filter (.01)'

      cparam(107)='!=0 --> add to h2 array in hlmhotz eqn'

      cparam(116)='!=0: x elements for fast tensor product'
      cparam(117)='!=0: y elements for fast tensor product'
      cparam(118)='!=0: z elements for fast tensor product'
c
c
      nparam=118

      return
      end

c-----------------------------------------------------------------------

      integer function nindx1(s1,s2,l2)
C
C     Return index of first character Not equal to S2
C
      CHARACTER*80 S1,S2
C
      N1=80-L2+1
      NINDX1=0
      IF (N1.LT.1) return
C
      DO 100 I=1,N1
         I2=I+L2-1
         IF (S1(I:I2).NE.S2(1:L2)) THEN
            NINDX1=I
            return
         endif
  100 CONTINUE
C
      return
      end

c-----------------------------------------------------------------------

      INTEGER FUNCTION LENGTH(STRING)

C Returns length of string ignoring trailing blanks

      CHARACTER*(*) STRING

      DO I = LEN(STRING), 1, -1
        IF(STRING(I:I) .NE. ' ') GO TO 20
      ENDDO
20    LENGTH = I
      END

c-----------------------------------------------------------------------

      SUBROUTINE PARSE (LINE,DELIM,ARGMNT,COUNT)

      CHARACTER*80  LINE         ! text line containing arguments
      CHARACTER*32  ARGMNT(50)    ! list of command line arguments
      CHARACTER*1  DELIM         ! list of command line arguments

      INTEGER SIZE(50)         ! number of chars in each argument
      INTEGER COUNT            ! number of arguments passed
      INTEGER LEFT,RIGHT       ! endpoints of arguments
C
      LEFT = 1
      RIGHT = 0
      COUNT = 0
C
      DO 10 I=1,50
      ARGMNT(I) = ' '
   10 CONTINUE
C
      DO 11 I=1,80
      IF (LINE(I:I) .NE. ' ' .AND. LINE(I:I) .NE. DELIM) THEN
        RIGHT = I
      ELSE IF (I .NE. 1 .AND. LINE(I-1:I-1) .NE. ' ' .AND.
     &         LINE(I-1:I-1) .NE. DELIM) THEN
        COUNT = COUNT + 1
        ARGMNT(COUNT) = LINE(LEFT:RIGHT)
        SIZE(COUNT) = RIGHT - LEFT + 1
        LEFT = I+1
      ELSE
        LEFT = I+1
      ENDIF
   11 CONTINUE
C
      RETURN
      END

c-----------------------------------------------------------------------
