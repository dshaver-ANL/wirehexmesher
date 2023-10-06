c-----------------------------------------------------------------------
      program wires
      include 'SIZE'
      integer nreps
      real partial
      logical ifperiodic
      nreps   = 1    !number of full spans
      partial = 0.0  !total length = (nreps+partial)*spanlength
      ifperiodic = .true.
      reaname_in="base.rea"
      reaname_out="wire_out.rea"
C Read .rea file to copy parameters
      call read_rea
c Load elements and bc's
      call load_convert(ifperiodic)
      write(*,*) "done convert"
c Load elements and bc's for pin
c      call load_convert_pin
c Generate a rea file with new geometry elements, curved sides and bc's
      call replicate(nreps,partial,ifperiodic)
      call test_range
      call gen_rea
      
c Output generated mesh to fld file
c      call outpost(xm1,ym1,zm1,pr,t,'   ')
c     Exit
      write(6,'(a)') '-------------------------------------'
      write(6,'(a)') ' Wire-mesh conversion finished ! '
      write(6,'(a)') '-------------------------------------'
c      stop
      end 
c-----------------------------------------------------------
      subroutine replicate(nreps,partial,ifperiodic)
      include 'SIZE'
      integer nreps,partial_elem,num_elem,blk_lay,elm_blk
      real partial
      logical ifperiodic
      real spanlen,z1_min,z1_max

      open(194,file='data_v5f.out',form='formatted',status='old')
      write(6,*) "*** Loading input data from mesher ***" 
      read(194,*) rlayer, rblk_lay, relm_blk, rnel_ftf, r_in
      close(194)
c   12 format(5F)

c     Parameters num_elem MUST be smaller than the parameter you are compiled with

      z1_min= 1000000.0
      z1_max=-1000000.0

      do i=1,num_elem 
      do j=1,27
        if (z1_max.lt.zm1(j,1,1,i)) z1_max=zm1(j,1,1,i)
        if (z1_min.gt.zm1(j,1,1,i)) z1_min=zm1(j,1,1,i)
      enddo
      enddo

      spanlen=z1_max-z1_min
  
      write(6,*) "============================================="
      write(*,*) "length of 1 span = ",spanlen
      num_elem1=num_elem   

      write(*,*) 'here 1',partial,num_elf,num_elem1
      partial_elem = partial*num_elem
      num_elem2=num_elem*nreps+partial_elem
   
      if (num_elem2.gt.max_num_elem) then
        write(6,*)" TOO MANY REPS ... check SIZE"
        stop
      endif 
      num_elem=num_elem2
      num_elf = num_elem2

      layer_base=rlayer
      layer=layer_base*nreps
      blk_lay=rblk_lay
      elm_blk=relm_blk

      nel_layer = blk_lay*elm_blk

      write(6,*)"* Layers: ",layer
      write(6,*)"* Blocks per Layer: ",blk_lay
      write(6,*)"* Elements per block: ",elm_blk


      do k=1,(nreps-1)
        write(*,*) 'rep', k
        do i=1,num_elem1
          do j=1,27
            xm1(j,1,1,i+k*num_elem1)=xm1(j,1,1,i)
            ym1(j,1,1,i+k*num_elem1)=ym1(j,1,1,i)
            zm1(j,1,1,i+k*num_elem1)=zm1(j,1,1,i)+dble(k)*spanlen
          enddo
          do j=1,4
            cbc(j,i+k*num_elem1,1)=cbc(j,i,1)
            cbc(j,i+k*num_elem1,2)=cbc(j,i,2)
            bc(5,j,i+k*num_elem1,1)=bc(5,j,i,1)
          enddo
        enddo
      enddo


c     THis section duplicates a fraction of a full pitch in addition
      if (partial.gt.0.0) then
        write(*,*) 'DOING PARTIAL'
        if(ifperiodic) write(*,*) 
     &     "Partial span incompatible with periodic BCs"
        ifperiodic=.false.
        k=nreps
        do i=1,nint(num_elem1*partial)
          do j=1,27
            xm1(j,1,1,i+k*num_elem1)=xm1(j,1,1,i)
            ym1(j,1,1,i+k*num_elem1)=ym1(j,1,1,i)
            zm1(j,1,1,i+k*num_elem1)=zm1(j,1,1,i)+dble(k)*spanlen
          enddo
          do j=1,4
            cbc(j,i+k*num_elem1,1)=cbc(j,i,1)
            cbc(j,i+k*num_elem1,2)=cbc(j,i,2)
            bc(5,j,i+k*num_elem1,1)=bc(5,j,i,1)
          enddo
        enddo
      endif


c      do k=1,(nreps-1)   
c      do i=1,num_elem1
c      do j=1,27  
c      xm1(j,1,1,i+k*num_elem1)=xm1(j,1,1,i)
c      ym1(j,1,1,i+k*num_elem1)=ym1(j,1,1,i)
c      zm1(j,1,1,i+k*num_elem1)=zm1(j,1,1,i)+dble(k)*spanlen
c      enddo
c      do j=1,4
c      cbc(j,i+k*num_elem1,1)=cbc(j,i,1)
c      cbc(j,i+k*num_elem1,2)=cbc(j,i,2)
c      enddo
c      enddo
c      enddo     

      do i=1,num_elem
        cbc(5,i,1)='E  '
        cbc(5,i,2)='E  '
        cbc(6,i,1)='E  '
        cbc(6,i,2)='E  '
      enddo

      last_rep   = num_elem1*(nreps-1+partial)
      last_layer = (layer_base-1)*nel_layer 

      write(*,*) "here", last_rep, last_layer, nel_layer, num_elem1

      do i=1,nel_layer
        if(ifperiodic) then  !needs to be tested
          cbc(5,i,1)='P  '
          cbc(5,i,2)='P  '
          bc(1,5,i,1)=i+last_rep+last_layer
          bc(2,5,i,1)=6
          bc(1,5,i,2)=i+last_rep+last_layer
          bc(2,5,i,2)=6
          bc(5,5,i,1)=0
          cbc(6,i+last_rep+last_layer,1)='P  '
          cbc(6,i+last_rep+last_layer,2)='P  '
          bc(1,6,i+last_rep+last_layer,1)=i
          bc(2,6,i+last_rep+last_layer,1)=5
          bc(1,6,i+last_rep+last_layer,2)=i
          bc(2,6,i+last_rep+last_layer,2)=5
          bc(5,6,i+last_rep+last_layer,1)=0
        else
          cbc(5,i,1)='v  '
          cbc(5,i,2)='t  '
          bc(5,5,i,1)=3
          cbc(6,i+last_rep+last_layer,1)='O  '
          cbc(6,i+last_rep+last_layer,2)='I  '
          bc(5,6,i+last_rep+last_layer,1)=4
        endif
      enddo
 
      write(*,*) "num_elem", num_elem, "num_elem1", num_elem1,
     &    "num_elem2", num_elem2

      return
      end  
c-----------------------------------------------------------------------
      subroutine test_range
      include 'SIZE'
      real*8 z1_max, z1_min, y1_max, y1_min, x1_max, x1_min     
 
      x1_min= 1000000.0
      x1_max=-1000000.0
      y1_min= 1000000.0
      y1_max=-1000000.0
      z1_min= 1000000.0
      z1_max=-1000000.0

      do i=1,num_elem 
      do j=1,27
        if (x1_max.lt.xm1(j,1,1,i)) x1_max=xm1(j,1,1,i)
        if (x1_min.gt.xm1(j,1,1,i)) x1_min=xm1(j,1,1,i)
        if (y1_max.lt.ym1(j,1,1,i)) y1_max=ym1(j,1,1,i)
        if (y1_min.gt.ym1(j,1,1,i)) y1_min=ym1(j,1,1,i)
        if (z1_max.lt.zm1(j,1,1,i)) z1_max=zm1(j,1,1,i)
        if (z1_min.gt.zm1(j,1,1,i)) z1_min=zm1(j,1,1,i)
      enddo
      enddo
  
      write(6,*) "============================================="
      write(6,*) "X range: ", x1_min," - ",x1_max
      write(6,*) "Y range: ", y1_min," - ",y1_max
      write(6,*) "Z range: ", z1_min," - ",z1_max
      write(6,*) "============================================="

      return
      end
 
c-----------------------------------------------------------------------
c
      subroutine load_convert(ifperiodic)
c
c     This routine loads the wire mesh and assigns the coordinates to
c     a lx1^3 mesh (hex27). 
c     Boundary conditions are also loaded through files 
c     using an appropriate convention   
c----------------------------------------------------------------------
      include 'SIZE'

      integer ie,f,f1
      integer layer, nel_layer, nel_pin
      integer elm_blk
      integer blk_lay
      integer i_inout
      real*8 period
      real*8 xx(27)
      real*8 yy(27)
      real*8 zz(27)
      real*8 xbc(6)
      real*8 ybc(6)
      integer block(1000000)
      integer block_L, el_rm, iee,skip
      real*8 pi,d1,d2,d3,d4,d5
      logical ifperiodic

      open(194,file='data_v5f.out',form='formatted',status='old')
      write(6,*) "*** Loading input data from mesher ***" 
      read(194,*) rlayer, rblk_lay, relm_blk, rnel_ftf, r_in
      close(194)

c     Inlet/outlet? Put to 1 to change 
      i_inout=1.0 !r_in
      if(ifperiodic) i_inout=-1.0 
c     Parameters num_elem MUST be smaller than the parameter you are compiled with
      layer=rlayer
      blk_lay=rblk_lay
      elm_blk=relm_blk 
      nel_ftf=rnel_ftf 
      nel_layer=blk_lay*elm_blk
      num_elem=layer*(nel_layer+nel_ftf) 
      write(6,*)"* Layers: ",layer
      write(6,*)"* Blocks per Layer: ",blk_lay
      write(6,*)"* Elements per block: ",elm_blk
      write(6,*)"* Elements in additional layer (if any): ",nel_ftf

      el_rm = 0
      block = block*0

c      do bl = 1,blk_lay
c        if(bl.eq.43.or.bl.eq.44.or.bl.eq.14.or.bl.eq.15) then
c           do e=1,elm_blk
c             block((block_L-1)*nel_layer+(bl-1)*elm_blk+e)=1
c             el_rm = el_rm+1
c           enddo
c         endif
c      enddo 
      write(*,*) el_rm," elements removed"
      write(*,*) num_elem, "original element"
c      layer=120
c      blk_lay=372
c      elm_blk=36
c      nel_ftf   = 5*6*12*layer
c      num_elem=layer*blk_lay*elm_blk + nel_ftf
      num_dim=3 
c      nel_layer = blk_lay*elm_blk
c      nel_ftf   = 5*6*12 
c      nel_pin   = 61*6*elm_blk
         
      open(195,file='BC_v5f.out',form='formatted',status='old')
      open(196,file='BCcon_v5f.out',form='formatted',status='old')
      open(197,file='xmesh_v5f.out',form='formatted',status='old')
      open(198,file='ymesh_v5f.out',form='formatted',status='old')
      open(199,file='zmesh_v5f.out',form='formatted',status='old')


      iee = 1
      do ie=1,num_elem
        if (mod(ie,10000).eq.0)
     &                       write(6,*) " ** Loaded ",ie," elements ..."

        do kk=1,27
          xx(kk)=0.0
          yy(kk)=0.0
          zz(kk)=0.0
        enddo

        do kk=1,6
          xbc(kk)=0.0
          ybc(kk)=0.0
        enddo

        read(195,*) (xbc(kk),kk=1,6)
        read(197,*) (xx(kk),kk=1,27)
        read(198,*) (yy(kk),kk=1,27)
        read(199,*) (zz(kk),kk=1,27)

        P = 1.189 
        pi = 4.*atan(1.)
        d1 = 4*P
        d2 = d1*sin(pi/3)
        d3 = 1.6*P 
        d4 = d1-P
        d5 = d4*sin(pi/3)

        block_s = 74 
        block_f = 77 

        skip=0

        if (skip.eq.0) then
 
          do kk=1,27
            xm1(kk,1,1,iee)=xx(kk)
            ym1(kk,1,1,iee)=yy(kk)
            zm1(kk,1,1,iee)=zz(kk)
          enddo
     
          do kk=1,6
            cbc(kk,iee,1)='E  '
            cbc(kk,iee,2)='E  '
          enddo

c         do j=1,5
c           bc(j,kk,ie,1)=0.0
c           bc(j,kk,ie,2)=0.0
c         enddo

          if (xbc(4).gt.0.0) then
            cbc(4,iee,1)='W  '
            cbc(4,iee,2)='I  '
c           cbc(4,ie,2)='f  '
            bc(5,4,iee,1) = 1
          endif
          if (xbc(2).gt.0.0) then
            cbc(2,iee,1)='W  '
            cbc(2,iee,2)='I  '  !outer wall
            bc(5,2,iee,1) = 2
c           cbc(2,ie,2)='f  '
          endif

          period=elm_blk*blk_lay

          if ((layer.gt.1).and.(ifperiodic)) then  
            if (xbc(5).gt.0.0) then
              if (ie.le.nel_layer) then
                period=nel_layer*(layer-1)
              else
                period=nel_ftf*(layer-1)
              endif 
              cbc(5,iee,1)='P  '
              cbc(5,iee,2)='P  '
              bc(1,5,iee,1)=ie+period-el_rm
              bc(2,5,iee,1)=6.0
              bc(1,5,iee,2)=ie+period-el_rm
              bc(2,5,iee,2)=6.0
              bc(5,5,iee,1) = 0
            endif
            if (xbc(6).gt.0.0) then
              if (ie.le.(nel_layer*layer)) then
                period=nel_layer*(layer-1)
              else
                period=nel_ftf*(layer-1)
              endif
              cbc(6,iee,1)='P  '
              cbc(6,iee,2)='P  '
              bc(1,6,iee,1)=ie-period+el_rm
              bc(2,6,iee,1)=5.0
              bc(1,6,iee,2)=ie-period+el_rm
              bc(2,6,iee,2)=5.0
              bc(5,6,iee,1) = 0
            endif
          else
            if (xbc(5).gt.0.0) then
              cbc(5,iee,1)='v  '
              cbc(5,iee,2)='t  '
              bc(5,5,iee,1) = 3
            endif
            if (xbc(6).gt.0.0) then 
              cbc(6,iee,1)='O  '
              cbc(6,iee,2)='I  '
              bc(5,6,iee,1) = 4
            endif 
          endif  
        iee = iee+1
        else
          write(*,*) "element removed"
          el_rm = el_rm+1
c         read(195,11) (xbc(kk),kk=1,6)
c         read(197,10) (xx(kk),kk=1,27)
c         read(198,10) (yy(kk),kk=1,27)
c         read(199,10) (zz(kk),kk=1,27)

        endif ! remove elements

      enddo !Main loop

      num_elem = num_elem - el_rm
      num_elf=num_elem

      close(195)
      close(196)
      close(197)
      close(198)
      close(199)

      return 
      end
C--------------------------------------------------------------------
c
      subroutine load_convert_pin
c
c     This routine loads the solid mesh and assigns the coordinates to
c     a lx1^3 mesh (hex27). 
c     Boundary conditions are also loaded through files 
c     using an appropriate convention   
c----------------------------------------------------------------------
      include 'SIZE'

      integer ie,f,f1
      integer layer
      integer elm_blk, nel_layer, ne_f
      integer blk_lay,period

      real*8 rlayer, rmlay
      real*8 xx(27)
      real*8 yy(27)
      real*8 zz(27)

      ne_f=num_elem
      num_elf=ne_f
      open(194,file='pdata_v5f.out',form='formatted',status='old')
      write(6,*) "*** Loading solid input data from mesher ***"
      read(194,*) rlayer, rmlay
      close(194)

      num_elem  = rmlay !rlayer
      nel_layer = rlayer !rmlay 
      num_dim=3


      write(6,*)"*** Also loading",num_elem,"solid  elements ***"

      if (nid.eq.0) then
      open(197,file='pxmesh_8-14.out',form='formatted',status='old')
      open(198,file='pymesh_8-14.out',form='formatted',status='old')
      open(199,file='pzmesh_8-14.out',form='formatted',status='old')
      endif

      do ie=1,num_elem
        if (mod(ie,10000).eq.0)
     &  write(6,*) "Loaded ",ie," elements ..."

      do kk=1,27
      xx(kk)=0.0
      yy(kk)=0.0
      zz(kk)=0.0
      enddo

      read(197,*) (xx(kk),kk=1,27)
      read(198,*) (yy(kk),kk=1,27)
      read(199,*) (zz(kk),kk=1,27)

c   10 format(27F)
c   11 format(2F)

      do kk=1,27
      xm1(kk,1,1,ne_f+ie)=xx(kk)
      ym1(kk,1,1,ne_f+ie)=yy(kk)
      zm1(kk,1,1,ne_f+ie)=zz(kk)
      enddo

      do kk=1,6
      cbc(kk,ne_f+ie,2)='E  '
      enddo

      cbc(4,ne_f+ie,2)='f  '
      cbc(2,ne_f+ie,2)='E  '

      period=num_elem-nel_layer

      if (ie.le.nel_layer) then
       cbc(5,ne_f+ie,2)='I  '
      endif

      if (ie.gt.period) then
       cbc(6,ne_f+ie,2)='I  '
      endif

      enddo

      if (nid.eq.0) then
      close(197)
      close(198)
      close(199)
      endif

      num_elem=num_elf+num_elem

      return
      end
C--------------------------------------------------------------------


C--------------------------------------------------------------------

      subroutine gen_rea

      include 'SIZE'

      open(unit=10,file=reaname_out,status='replace')

      call write_rea_top
      write(*,*) 'done write_rea_top'
      call gen_rea_xyz
      write(*,*) 'done gen_rea_xyz'
      call gen_rea_curve(2)
      write(*,*) 'done gen_rea_curve'
      call gen_rea_bc
      write(*,*) 'done gen_rea_bc'
      call write_rea_bot
      write(*,*) 'done write_rea_bot'

      close(10)

      return
      end

c-----------------------------------------------------------------------

      subroutine gen_rea_xyz

      include 'SIZE'

      parameter (lv=2**ldim)

      real       xyz(lv,ldim) 

      integer isym2pre(8)   ! Symmetric-to-prenek vertex ordering
      save    isym2pre
      data    isym2pre / 1 , 2 , 4 , 3 , 5 , 6 , 8 , 7 /

      integer     eb, eg
      character*1 letapt

      nxs = 3-1  ! always nx1=3
      nys = 3-1
      nzs = 3-1

      if (num_dim.eq.2) then
        write(10,*) 
     &  '**MESH DATA** 1st line is X of corner 1,2,3,4. 2nd line is Y.'
      else
        write(10,*)
     & '**MESH DATA** 6 lines are X,Y,Z;X,Y,Z. Columns corners 1-4;5-8'
      endif

      write(*,*) "write numbers", num_elem, num_dim, num_elf

      write(10,'(i12,i3,i12,'' NEL,NDIM,NELV'')')
     &                        num_elem, num_dim, num_elf

        do eg=1,num_elem 
            l=0
            do k=0,1
            do j=0,1
            do i=0,1
              l=l+1
              li=isym2pre(l)
              xyz(li,1) = xm1(1+i*nxs,1+j*nys,1+k*nzs,eg)
              xyz(li,2) = ym1(1+i*nxs,1+j*nys,1+k*nzs,eg)
              xyz(li,3) = zm1(1+i*nxs,1+j*nys,1+k*nzs,eg)    
            enddo
            enddo
            enddo
          letapt = 'a'
          numapt = 1
          igr    = 0 
          write(10,'(a15,i12,a2,i5,a1,a10,i6)')
     &      '      ELEMENT  ',eg,' [',numapt,letapt,']    GROUP',igr
            write(10,'(4g15.7)')(xyz(ic,1),ic=1,4)
            write(10,'(4g15.7)')(xyz(ic,2),ic=1,4)
            write(10,'(4g15.7)')(xyz(ic,3),ic=1,4)
            write(10,'(4g15.7)')(xyz(ic,1),ic=5,8)
            write(10,'(4g15.7)')(xyz(ic,2),ic=5,8)
            write(10,'(4g15.7)')(xyz(ic,3),ic=5,8)
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine gen_rea_curve

c     This routine is complex because we must first count number of
c     nontrivial curved sides.
c     A two pass strategy is used:  first count, then write


      include 'SIZE'

      integer e,eb,eg,ii
 
      do e=1,num_elem
c        write(*,*) 'element: ',e
        call gen_rea_midside_e(e)
      enddo

      nedge = 4 + 8*(num_dim-2)
      ncurvn = 0
      do e=1,num_elem
        do i=1,nedge
          if (ccurve(i,e).ne.' ') ncurvn = ncurvn+1
        enddo
      enddo

      WRITE(10,*)' ***** CURVED SIDE DATA *****'
      WRITE(10,'(I10,A20,A33)') ncurvn,' Curved sides follow',
     $   ' IEDGE,IEL,CURVE(I),I=1,5, CCURVE'

      do eg=1,num_elem
      do i=1,nedge
            if (ccurve(i,eg).ne.' ') then
              if (num_elem.lt.1000) then
                write(10,'(i3,i3,5g14.6,1x,a1)') i,eg,
     $               (curve(k,i,eg),k=1,5),ccurve(i,eg)
              elseif (num_elem.lt.1000000) then
                write(10,'(i2,i6,5g14.6,1x,a1)') i,eg,
     $               (curve(k,i,eg),k=1,5),ccurve(i,eg)
              else
                write(10,'(i2,i12,5g14.6,1x,a1)') i,eg,
     $               (curve(k,i,eg),k=1,5),ccurve(i,eg)
              endif
            endif
        enddo
        enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine gen_rea_midside_e(e)

      include 'SIZE'

      real        x3(27),y3(27),z3(27),xyz(3,3)
      integer     e,nedge,ncc1

      integer e3(3,12)
      save    e3
      data    e3 /  1, 2, 3,    3, 6, 9,    9, 8, 7,    7, 4, 1
     $           , 19,20,21,   21,24,27,   27,26,25,   25,22,19
     $           ,  1,10,19,    3,12,21,    9,18,27,    7,16,25 /

      real len


      call map2reg(x3,3,xm1(1,1,1,e),1)  ! Map to 3x3x3 array
      call map2reg(y3,3,ym1(1,1,1,e),1)
      if (num_dim.eq.3) call map2reg(z3,3,zm1(1,1,1,e),1)
 

      tol   = 1.e-4
      tol2  = tol**2
      nedge = 4 + 8*(num_dim-2)

      do i=1,nedge
            do j=1,5
               curve(j,i,e)=0.0      
            enddo 
            ccurve(i,e) = ' '
            do j=1,3
               xyz(1,j)=x3(e3(j,i))
               xyz(2,j)=y3(e3(j,i))
               xyz(3,j)=z3(e3(j,i))
            enddo
            len = 0.
            h   = 0.
            do j=1,num_dim
               xmid = .5*(xyz(j,1)+xyz(j,3))
               h    = h   + (xyz(j,2)-xmid)**2
               len  = len + (xyz(j,3)-xyz(j,1))**2
            enddo
            if (h.gt.tol2*len) ccurve(i,e) = 'm'
            if (h.gt.tol2*len) then
                curve(1,i,e)=xyz(1,2)
                curve(2,i,e)=xyz(2,2)
                curve(3,i,e)=xyz(3,2)
            endif 
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine map2reg(ur,n,u,nel)
c
c     Map scalar field u() to regular n x n x n array ur
c
      include 'SIZE'

      real    ur(1), u(3*3*3,1)
      integer e

      ldr = n**num_dim

      k=1
      do e=1,nel
         if (num_dim.eq.2) call map2reg_2di_e(ur(k),n,u(1,e),3)
         if (num_dim.eq.3) call map2reg_3di_e(ur(k),n,u(1,e),3)
         k = k + ldr
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine map2reg_2di_e(uf,n,uc,m) ! Fine, uniform pt

      real      uf(n,n),uc(m,m)
      parameter (l=3) 
c      parameter (l=50)
      integer   j(l*l),jt(l*l),w(l*l),z(l)

      call zwgll     (z,w,m)
      call zuni      (w,n)
      call gen_int_gz(j,jt,w,n,z,m)
      call mxmf2     (j,n,uc,m,w ,m)
      call mxmf2     (w,n,jt,m,uf,n)

      return
      end

c-----------------------------------------------------------------------

      subroutine map2reg_3di_e(uf,n,uc,m) ! Fine, uniform pt

      real      uf(n,n,n),uc(m,m,m)
      parameter (l=3)
c      parameter (l=50)
      integer   j(l*l),jt(l*l),v(l*l*l),w(l*l*l),z(l)

      call zwgll     (z,w,m)
      call zuni      (w,n)
      call gen_int_gz(j,jt,w,n,z,m)
      mm = m*m
      mn = m*n
      nn = n*n
      call mxmf2(j,n,uc,m,v ,mm)
      iv=1
      iw=1
      do k=1,m
         call mxmf2(v(iv),n,jt,m,w(iw),n)
         iv = iv+mn
         iw = iw+nn
      enddo
      call mxmf2(w,nn,jt,m,uf,n)

      return
      end

c-----------------------------------------------------------------------

      subroutine gen_int_gz(j,jt,g,n,z,m)

c     Generate interpolater from m z points to n g points

c        j   = interpolation matrix, mapping from z to g
c        jt  = transpose of interpolation matrix
c        m   = number of points on z grid
c        n   = number of points on g grid

      real j(n,m),jt(m,n),g(n),z(m)

      mpoly  = m-1
      do i=1,n
         call fd_weights_full(g(i),z,mpoly,0,jt(1,i))
      enddo

      call transpose(j,n,jt,m)

      return
      end

c-----------------------------------------------------------------------

      subroutine gen_rea_bc

      include 'SIZE'

      integer ii,i,eg

      parameter (lblock=500)
      real       vbc(5,6,lblock),wk(5*6*lblock)
      integer    ibc(6,lblock)

      character*1 s4(4)
      character*3 s3
      integer     i4
      equivalence(i4,s4)
      equivalence(s3,s4)

      character*1 chtemp
      save        chtemp
      data        chtemp /' '/   ! For mesh bcs


      write(10,*)' ***** BOUNDARY CONDITIONS *****'

      write(10,*) ' *****    FLUID   BOUNDARY CONDITIONS *****'

      nface = 2*num_dim

              do eg=1,num_elf
                do i=1,nface
                  call chcopy(s4,cbc(i,eg,1),3)  
                  if (num_elem.lt.1000) then
                    write(10,'(a1,a3,2i3,5g14.6)')
     &                    chtemp,s3,eg,i,(bc(ii,i,eg,1),ii=1,5)
                  elseif (num_elem.lt.100000) then
                    write(10,'(a1,a3,i5,i1,5g14.6)')
     &                    chtemp,s3,eg,i,(bc(ii,i,eg,1),ii=1,5)
                  elseif (num_elem.lt.1000000) then
                    write(10,'(a1,a3,i6,5g14.6)')
     &                    chtemp,s3,eg,(bc(ii,i,eg,1),ii=1,5)
                  else
                    write(10,'(a1,a3,i12,5g18.11)')
     &                    chtemp,s3,eg,(bc(ii,i,eg,1),ii=1,5)
                  endif
                enddo
              enddo


      write(10,*) ' ***** THERMAL BOUNDARY CONDITIONS *****'

      nface = 2*num_dim

              do eg=1,num_elem
                do i=1,nface
                  call chcopy(s4,cbc(i,eg,2),3)
                  if (num_elem.lt.1000) then
                    write(10,'(a1,a3,2i3,5g14.6)')
     &                    chtemp,s3,eg,i,(bc(ii,i,eg,2),ii=1,5)
                  elseif (num_elem.lt.100000) then
                    write(10,'(a1,a3,i5,i1,5g14.6)')
     &                    chtemp,s3,eg,i,(bc(ii,i,eg,2),ii=1,5)
                  elseif (num_elem.lt.1000000) then
                    write(10,'(a1,a3,i6,5g14.6)')
     &                    chtemp,s3,eg,(bc(ii,i,eg,2),ii=1,5)
                  else
                    write(10,'(a1,a3,i12,5g18.11)')
     &                    chtemp,s3,eg,(bc(ii,i,eg,2),ii=1,5)
                  endif
                enddo
              enddo


      return
      end

c-----------------------------------------------------------------------
