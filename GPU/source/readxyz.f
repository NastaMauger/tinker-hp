c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine readxyz  --  input of Cartesian coordinates  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "readxyz" gets a set of Cartesian coordinates from
c     an external disk file
c
c
#include "tinker_precision.h"
      subroutine readxyz (ixyz)
      use abinitio
      use atmtyp
      use atoms
      use atomsMirror,only:xm=>x,ym=>y,zm=>z
     &               ,reCast_position,atomsmirror_init
     &               ,download_position
      use bound
      use boxes
      use couple
      use domdec    ,only: ranktot
      use files
      use inform
      use iounit
      use iso_c_binding ,only: c_loc, c_f_pointer
      use nvshmem
      use titles
      use timestat  ,only:timer_io,timer_enter,timer_exit,quiet_timers
      use tinheader ,only:ti_p,re_p
      use tinMemory
      implicit none
      integer i,j,k,m
      integer ixyz,nmax
      integer next,size
      integer first,last
      integer nexttext
      integer trimtext
      integer, allocatable :: list(:)
      real*8 :: xlen,ylen,zlen
      real*8 :: aang,bang,gang
      logical exist,opened
      logical quit,reorder
      logical clash,use_wrap
      character*240 xyzfile
      character*240 record
      character*240 string
 
      integer number_qm_file
      character*1 name_first
c
      call timer_enter(timer_io)
c
c     initialize the total number of atoms in the system
c
      n = 0
c
c     open the input file if it has not already been done
c
      inquire (unit=ixyz,opened=opened)
      if (.not. opened) then
         xyzfile = filename(1:leng)//'.xyz'
         call version (xyzfile,'old')
         inquire (file=xyzfile,exist=exist)
         if (exist) then
            open (unit=ixyz,file=xyzfile,status='old')
            rewind (unit=ixyz)
         else
            write (iout,10)
   10       format (/,' READXYZ  --  Unable to Find the Cartesian',
     &                 ' Coordinates File')
            call fatal
         end if
      end if
c
c     read first line and return if already at end of file
c
      quit = .false.
      abort = .true.
      size = 0
      do while (size .eq. 0)
         read (ixyz,20,err=80,end=140)  record
   20    format (a240)
         size = trimtext (record)
      end do
      abort = .false.
      quit = .true.
c
c     parse the title line to get the number of atoms
c
      i = 0
      next = 1
      call gettext (record,string,next)
      read (string,*,err=80,end=80)  n
c
c     allocate global arrays
c
      call prmem_request(x,n)
      call prmem_request(y,n)
      call prmem_request(z,n)
      call prmem_request(xold,n)
      call prmem_request(yold,n)
      call prmem_request(zold,n)
      call prmem_request(type,n)
      call prmem_request(tag,n,config=mhostonly)
      use_wrap = app_id.eq.dynamic_a.or.app_id.eq.pimd_a
      if (use_wrap) then
         call prmem_request(pbcWrap,n)
         call c_f_pointer(c_loc(pbcWrap),pbcWrapIdx,[4*n])
      end if
      call atomsmirror_init

      if (.not.associated(name)) allocate (name(n))
      if (.not.allocated (i12 )) allocate (i12(maxvalue,n))
      if (.not.allocated (n12 )) then
         allocate (n12(n))
         s_prmem = s_prmem + int(maxvalue+1,8)*sizeof(n12)
         s_prmem = s_prmem + sizeof(name)
        sd_prmem =sd_prmem + int(maxvalue+1,8)*sizeof(n12)
      end if

!$acc parallel loop present(x,y,z,xold,yold,zold)
      do i = 1, n
         x   (i) = 0_ti_p
         y   (i) = 0_ti_p
         z   (i) = 0_ti_p
         xold(i) = 0.0_ti_p
         yold(i) = 0.0_ti_p
         zold(i) = 0.0_ti_p
      end do

!$acc update host(xold,yold,zold)
      type     = 0
      n12      = 0
      i12      = 0
      tag      = 0
      maxbnd   = 4*n
      maxang   = 6*n
      maxtors  = 18*n
      maxbitor = 54*n
c
c     extract the title and determine its length
c
      string = record(next:240)
      first = nexttext (string)
      last = trimtext (string)
      if (last .eq. 0) then
         title = ' '
         ltitle = 0
      else
         title = string(first:last)
         ltitle = trimtext (title)
      end if
c
c     check for too few or too many total atoms in the file
c
      if (n .le. 0) then
         write (iout,30)
   30    format (/,' READXYZ  --  The Coordinate File Does Not',
     &              ' Contain Any Atoms')
         call fatal
      end if
c     nloop help vectorization
c     nloop is n if n is a multiple of 16, or the first one greater
      nloop = merge(n,(int(n/16)+1)*16,(mod(n,16).eq.0))
      if (npes.gt.0) n_pe = ceiling(1.0*n/npes)

c
c     initialize coordinates and connectivities for each atom
c
      do i = 1, n
         tag(i) = 0
         name(i) = '   '
         x(i) = 0.0_ti_p
         y(i) = 0.0_ti_p
         z(i) = 0.0_ti_p
         type(i) = 0
         do j = 1, maxvalue
            i12(j,i) = 0
         end do
      end do
!$acc update device(n)
c
c     read the coordinates and connectivities for each atom
c
      do i = 1, n
         next = 1
         size = 0
         do while (size .eq. 0)
            read (ixyz,50,err=80,end=80)  record
   50       format (a240)
            size = trimtext (record)
            if (i .eq. 1) then
               next = 1
               call getword (record,name(i),next)
               name_first=name(i)
               if (name(i) .ne. '   ')  goto 60
               read (record,*,err=60,end=60)  xlen,ylen,zlen,
     &                                        aang,bang,gang
               size  = 0
               xbox  = xlen
               ybox  = ylen
               zbox  = zlen
               alpha = aang
               beta  = bang
               gamma = gang
               use_bounds = .true.
!$acc update device(xbox,ybox,zbox,alpha,beta,gamma)
c               call lattice
   60       continue
            end if
         end do
         read (record,*,err=80,end=80)  tag(i)
         call getword (record,name(i),next)
         string = record(next:240)
         read (string,*,err=70,end=70)  xm(i),ym(i),zm(i),type(i),
     &                                  (i12(j,i),j=1,maxvalue)
   70    continue
      end do

!$acc update device(xm,ym,zm)

      if (aiMD) then
        if (orca_qm) then
          number_qm_file=415
        elseif (g16_qm) then
          number_qm_file=416
        elseif (psi4_qm) then
          number_qm_file=417
        elseif (pyscf_qm) then
          number_qm_file=418
        elseif (qchem_qm) then
          number_qm_file=419
        else
          number_qm_file=-1
        endif
!$acc update host(xm, ym, zm)
        do i = 1, n
          if (i == 1) then
            write(number_qm_file,'(A,3(F10.4, 1X))') trim(name_first)
     &                                                ,xm(i),ym(i),zm(i)
          else
            write(number_qm_file,'(A,3(F10.4,1X))') trim(name(i))
     &                                                ,xm(i),ym(i),zm(i)
          endif
        end do
      endif



      ! Init Wrapping state buffer
      if (app_id.eq.dynamic_a) then
         do i = 1, n
            pbcWrap(i) = 0
         end do
!$acc update device(pbcWrap)
      end if

      quit = .false.
   80 continue
      if (.not. opened)  close (unit=ixyz)

      if (t_p.ne.r_p) then
         call reCast_position
!$acc wait
         call download_position
      end if
c
c     an error occurred in reading the coordinate file
c
      if (quit) then
         write (iout,90)  i
   90    format (/,' READXYZ  --  Error in Coordinate File at Atom',i6)
         call fatal
      end if
c
c     for each atom, count and sort its attached atoms
c
      do i = 1, n
         n12(i) = 0
         do j = maxvalue, 1, -1
            if (i12(j,i) .ne. 0) then
               n12(i) = j
               goto 100
            end if
         end do
  100    continue
         call sort (n12(i),i12(1,i))
      end do
!$acc update device(n12(:))
c
c     perform dynamic allocation of some local arrays
c
      nmax = 0
      do i = 1, n
         nmax = max(tag(i),nmax)
         do j = 1, n12(i)
            nmax = max(i12(j,i),nmax)
         end do
      end do
      allocate (list(nmax))
c
c     check for scrambled atom order and attempt to renumber
c
      reorder = .false.
      do i = 1, n
         list(tag(i)) = i
         if (tag(i) .ne. i)  reorder = .true.
      end do
      if (reorder) then
         write (iout,110)
  110    format (/,' READXYZ  --  Atom Labels not Sequential,',
     &              ' Attempting to Renumber')
         do i = 1, n
            tag(i) = i
            do j = 1, n12(i)
               i12(j,i) = list(i12(j,i))
            end do
            call sort (n12(i),i12(1,i))
         end do
      end if
!$acc update device(i12(:,:))
c
c     Determine maxn12 and maxd12
c
      block
      integer b0,idx
      maxn12_ = 0
      maxd12  = 0
      do i = 1,n
         maxn12_ = max(maxn12_,n12(i))
      end do
      do i = 1,n; do j = 1,n12(i)
         b0 = (i12(j,i))
         if (b0.eq.0) cycle
         b0 = abs(b0 - i)
         if (b0.gt.maxd12) then
            maxd12 = b0
            idx    = i
         end if
      end do; end do
#ifdef _OPENACC
      if (ranktot.eq.0.And.maxn12_.gt.4) then
 111     format(" Tinker-HP: Detected atom ",I0," with a",I2
     &   ," valence !!!",/,9x," --- Operating to Legacy subroutines"
     &   ,/)
         write(*,111) maxd12,maxn12_
      end if
#endif
      end block
c
c     perform deallocation of some local arrays
c
      deallocate (list)
c
c     check for atom pairs with identical coordinates
c
      clash = .false.
      if (n .le. 10000)  call chkxyz (clash)
c
c     make sure that all connectivities are bidirectional
c
      do i = 1, n
         do j = 1, n12(i)
            k = i12(j,i)
            do m = 1, n12(k)
               if (i12(m,k) .eq. i)  goto 130
            end do
            write (iout,120)  k,i
  120       format (/,' READXYZ  --  Check Connection of Atom',
     &                 i6,' to Atom',i6)
            call fatal
  130       continue
         end do
      end do

  140 continue
      if (.not. opened)  close (unit=ixyz)

      call timer_exit(timer_io)
      end
