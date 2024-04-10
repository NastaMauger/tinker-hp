c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine shakeup  --  setup of holonomic constraints  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "shakeup" initializes any holonomic constraints for use with
c     the RATTLE algorithm during molecular dynamics
c
c
#include "tinker_precision.h"
      subroutine shakeup(init)
      use angle
      use atmlst
      use atmtyp
      use atoms
      use bond
      use bound
      use couple
      use domdec
      use freeze
      use keys
      use inform    ,only: deb_Path
      use math
      use molcul
      use potent
      use strbnd
      use tinheader ,only: ti_p,re_p
      use tinMemory
      use usage
      implicit none
      integer i,j,k,m,nh
      integer ia,ib,ic
      integer ja,jb,jc
      integer ilist,next

      integer nprocloc,start,stop,im
      integer :: maxratmol
      real(t_p) weigh,xmid,ymid,zmid
      real(t_p) rab,rbc,rac
      real(t_p) cosine
      real(t_p) xa,xb,ya,yb,za,zb

      logical done,init,docompute
      character*9 rattyp
      character*20 keyword
      character*240 record
      character*240 string
c
      if (init) then
c
c     set defaults for constraints and convergence tolerance
c
        nrat = 0
        nratx = 0
        rateps = 0.000001_ti_p
        use_rattle = .true.
c
c       perform dynamic allocation of some global arrays
c
c
c       to move to shared memory segments
c
        if (deb_Path) print*,'shakeup'
        call alloc_shared_shake
c
c       process keywords containing holonomic constraint options
c
        do k = 1, nkey
           next = 1
           record = keyline(k)
           call gettext (record,keyword,next)
           call upcase (keyword)
           string = record(next:240)
           if (keyword(1:11) .eq. 'RATTLE-EPS ') then
              read (string,*,err=10,end=10)  rateps
           end if
   10      continue
        end do
c
c       process keywords containing various constraint types
c
        do k = 1, nkey
           next = 1
           record = keyline(k)
           call upcase (record)
           call gettext (record,keyword,next)
           if (keyword(1:7) .eq. 'RATTLE ') then
              call getword (record,rattyp,next)
c
c       constrain all bond lengths at their ideal values
c
              if (rattyp(1:5) .eq. 'BONDS') then
                 do i = 1, nbond
                    ia = ibnd(1,i)
                    ib = ibnd(2,i)
                    if (use(ia) .or. use(ib)) then
                       nrat = nrat + 1
                       irat(1,nrat) = ia
                       irat(2,nrat) = ib
                       krat(nrat) = bl(i)
                    end if
                 end do
c
c       constrain bonds and independent angles at ideal values
c
              else if (rattyp(1:6) .eq. 'ANGLES') then
                 do i = 1, nbond
                    ia = ibnd(1,i)
                    ib = ibnd(2,i)
                    if (use(ia) .or. use(ib)) then
                       nrat = nrat + 1
                       irat(1,nrat) = ia
                       irat(2,nrat) = ib
                       krat(nrat) = bl(i)
                    end if
                 end do
                 do i = 1, n
                    if (n12(i) .gt. 1) then
                       do j = 1, 2*n12(i)-3
                          ilist = anglist(j,i)
                          ia = iang(1,ilist)
                          ib = iang(2,ilist)
                          ic = iang(3,ilist)
                          if (use(ia) .or. use(ib) .or. use(ic)) then
                             do m = 1, n12(ib)
                                if (i12(m,ib) .eq. ia) then
                                   rab = bl(bndlist(m,ib))
                                else if (i12(m,ib) .eq. ic) then
                                   rbc = bl(bndlist(m,ib))
                                end if
                             end do
                             cosine = cos(anat(ilist)/radian)
                             rac = sqrt(rab*rab+rbc*rbc
     &                                     -2.0_ti_p*rab*rbc*cosine)
                             nrat = nrat + 1
                             irat(1,nrat) = ia
                             irat(2,nrat) = ic
                             krat(nrat) = rac
                             call chkangle (ia,ib,ic)
                          end if
                       end do
                    end if
                 end do
c
c       fix bond length in diatomics to give a rigid molecule
c
              else if (rattyp(1:8) .eq. 'DIATOMIC') then
                 do i = 1, nbond
                    ia = ibnd(1,i)
                    ib = ibnd(2,i)
                    if (n12(ia).eq.1 .and. n12(ib).eq.1) then
                       if (use(ia) .or. use(ib)) then
                          nrat = nrat + 1
                          irat(1,nrat) = ia
                          irat(2,nrat) = ib
                          krat(nrat) = bl(i)
                       end if
                    end if
                 end do
c
c       fix bonds and angle in triatomics to give a rigid molecule
c
              else if (rattyp(1:9) .eq. 'TRIATOMIC') then
                 do i = 1, nangle
                    ia = iang(1,i)
                    ib = iang(2,i)
                    ic = iang(3,i)
                    if (n12(ia)+n12(ib)+n12(ic) .eq. 4) then
                       rab = bl(bndlist(1,ia))
                       rbc = bl(bndlist(1,ic))
                       cosine = cos(anat(i)/radian)
                       rac = sqrt(rab**2+rbc**2-2.0_ti_p*rab*rbc*cosine)
                       if (use(ia) .or. use(ib)) then
                          nrat = nrat + 1
                          irat(1,nrat) = ia
                          irat(2,nrat) = ib
                          krat(nrat) = rab
                       end if
                       if (use(ib) .or. use(ic)) then
                          nrat = nrat + 1
                          irat(1,nrat) = ib
                          irat(2,nrat) = ic
                          krat(nrat) = rbc
                       end if
                       if (use(ia) .or. use(ib) .or. use(ic)) then
                          nrat = nrat + 1
                          irat(1,nrat) = ia
                          irat(2,nrat) = ic
                          krat(nrat) = rac
                       end if
                    end if
                 end do
c
c       fix bonds and angles of each water to give a rigid molecule
c
              else if (rattyp(1:5) .eq. 'WATER') then
                 do i = 1, n
                    nh = 0
                    if (atomic(i) .eq. 8) then
                       do j = 1, n12(i)
                          if (atomic(i12(j,i)) .eq. 1)  nh = nh + 1
                       end do
                    end if
                    if (nh .ge. 2) then
                       do j = 1, n12(i)
                          ilist = bndlist(j,i)
                          ia = ibnd(1,ilist)
                          ib = ibnd(2,ilist)
                          ja = atomic(ia)
                          jb = atomic(ib)
                          if (use(ia) .or. use(ib)) then
                             if (ja.eq.1 .or. jb.eq.1) then
                                nrat = nrat + 1
                                irat(1,nrat) = ia
                                irat(2,nrat) = ib
                                krat(nrat) = bl(ilist)
                             end if
                          end if
                       end do
                       do j = 1, 2*n12(i)-3
                          ilist = anglist(j,i)
                          ia = iang(1,ilist)
                          ib = iang(2,ilist)
                          ic = iang(3,ilist)
                          ja = atomic(ia)
                          jc = atomic(ic)
                          if (use(ia) .or. use(ib) .or. use(ic)) then
                             if (ja.eq.1 .and. jc.eq.1) then
                             do m = 1, n12(ib)
                                if (i12(m,ib) .eq. ia) then
                                   rab = bl(bndlist(m,ib))
                                else if (i12(m,ib) .eq. ic) then
                                   rbc = bl(bndlist(m,ib))
                                end if
                             end do
                             cosine = cos(anat(ilist)/radian)
                             rac = sqrt(rab*rab+rbc*rbc
     &                                     -2.0_ti_p*rab*rbc*cosine)
                             nrat = nrat + 1
                             irat(1,nrat) = ia
                             irat(2,nrat) = ic
                             krat(nrat) = rac
                             call chkangle (ia,ib,ic)
                             end if
                          end if
                       end do
                    end if
                 end do
c
c       fix all bonds to hydrogen atoms at their ideal length
c
              else
                 do i = 1, nbond
                    ia = ibnd(1,i)
                    ib = ibnd(2,i)
                    if (atomic(ia).eq.1 .or. atomic(ib).eq.1) then
                       if (use(ia) .or. use(ib)) then
                          nrat = nrat + 1
                          irat(1,nrat) = ia
                          irat(2,nrat) = ib
                          krat(nrat) = bl(i)
                       end if
                    end if
                 end do
              end if
           end if
        end do
c
c       process keywords containing specific distance constraints
c
        do k = 1, nkey
           next = 1
           record = keyline(k)
           call gettext (record,keyword,next)
           call upcase (keyword)
           if (keyword(1:16) .eq. 'RATTLE-DISTANCE ') then
              call getnumb (record,ia,next)
              call getnumb (record,ib,next)
              rab = 0.0_ti_p
              string = record(next:240)
              read (string,*,err=20,end=20)  rab
   20         continue
              if (rab .eq. 0.0_ti_p) then
                 do i = 1, n12(ia)
                    if (i12(i,ia) .eq. ib) then
                       rab = bl(bndlist(i,ia))
                    end if
                 end do
              end if
              if (rab .eq. 0.0_ti_p) then
                 rab = sqrt((x(ia)-x(ib))**2 + (y(ia)-y(ib))**2
     &                             + (z(ia)-z(ib))**2)
              end if
              done = .false.
              do j = 1, nrat
                 ja = irat(1,j)
                 jb = irat(2,j)
                 if ((ia.eq.ja .and. ib.eq.jb) .or.
     &               (ia.eq.jb .and. ib.eq.ja)) then
                    done = .true.
                    krat(j) = rab
                 end if
              end do
              if (.not. done) then
                 nrat = nrat + 1
                 irat(1,nrat) = ia
                 irat(2,nrat) = ib
                 krat(nrat) = rab
              end if
c
c       process keywords containing atom group spatial constraints
c
           else if (keyword(1:13) .eq. 'RATTLE-PLANE ') then
              call getnumb (record,ia,next)
              nratx = nratx + 1
              iratx(nratx) = ia
              kratx(nratx) = 1
           else if (keyword(1:12) .eq. 'RATTLE-LINE ') then
              call getnumb (record,ia,next)
              nratx = nratx + 1
              iratx(nratx) = ia
              kratx(nratx) = 2
           else if (keyword(1:14) .eq. 'RATTLE-ORIGIN ') then
              call getnumb (record,ia,next)
              nratx = nratx + 1
              iratx(nratx) = ia
              kratx(nratx) = 3
           end if
        end do
c
c       find and remove any duplicate distance constraints
c
        do i = 1, nrat-1
           ia = irat(1,i)
           ib = irat(2,i)
           do j = i+1, nrat
              ja = irat(1,j)
              jb = irat(2,j)
              if ((ia.eq.ja .and. ib.eq.jb) .or.
     &            (ia.eq.jb .and. ib.eq.ja))  krat(j) = -1.0_ti_p
           end do
        end do
        k = nrat
        do i = k, 1, -1
           if (krat(i) .lt. 0.0_ti_p) then
              do j = i, k-1
                 irat(1,j) = irat(1,j+1)
                 irat(2,j) = irat(2,j+1)
                 krat(j) = krat(j+1)
              end do
              nrat = nrat - 1
           end if
        end do
c
c       set flag to apply minimum image to intermolecular constraints
c
        do i = 1, nrat
           ia = irat(1,i)
           ib = irat(2,i)
           if (use_bounds .and. (molcule(ia).ne.molcule(ib))) then
              ratimage(i) = .true.
           else if (use_polymer) then
              ratimage(i) = .true.
           else
              ratimage(i) = .false.
           end if
        end do
c
c       if no holonomic constraints are present, turn off their use
c
        if (nrat.eq.0 .and. nratx.eq.0)  then
            use_rattle = .false.
            return
        end if

        allocate(nratmol(nmol))
        nratmol(:) = 0
        do i = 1, nrat
          ia = irat(1,i)
          ib = irat(2,i)
          if (molcule(ia).ne.molcule(ib)) then
            write(*,*) 'RATTLE only compatible with constrains 
     $    involving atoms inside the same molecule'
            call fatal
          end if
          nratmol(molcule(ia)) = nratmol(molcule(ia)) + 1
        enddo
        maxratmol = maxval(nratmol)
        allocate(iratmol(maxratmol,nmol))
        iratmol(:,:) = -1
        nratmol(:) = 0
        do i = 1, nrat
            ia = irat(1,i)
            ib = irat(2,i)
            nratmol(molcule(ia)) = nratmol(molcule(ia)) + 1
            iratmol(nratmol(molcule(ia)),molcule(ia)) = i
        enddo
!$acc enter data copyin(iratx,kratx,irat,krat,ratimage,nratmol,iratmol)

      end if

      if (nrat.eq.0) return

      if (nproc > 1) then
         write (0,*) 'Error: RATTLE is not implemented for nproc>1'
         call fatal
      end if

      return
      !!! THE REMAINING CODE IS ONLY FOR MPI => we do not need it for now !!!
      
c
c     get local constraints
c
      call prmem_request(ratglob,nrat) !,config=mhostonly)
      nratloc = 0

      if (use_pmecore) then
        nprocloc = ndir
      else
        nprocloc = nproc
      end if
c
      do i = 1, nrat
        ia = irat(1,i)
        ib = irat(2,i)
        if (molcule(ia).ne.molcule(ib)) then
          write(*,*) 'RATTLE only compatible with constrains involving
     $    atoms inside the same molecule'
          call fatal
        end if
c
c       constrain treated by the proc that has the center of mass of the
c       molecule
c
        im = molcule(ia)
        start = imol(1,im)
        stop = imol(2,im)
        xmid = 0.0_ti_p
        ymid = 0.0_ti_p
        zmid = 0.0_ti_p
        do j = start, stop
           k = kmol(j)
           weigh = mass(k)
           xmid = xmid + x(k)*weigh
           ymid = ymid + y(k)*weigh
           zmid = zmid + z(k)*weigh
        end do
        weigh = molmass(im)
        xmid = xmid / weigh
        ymid = ymid / weigh
        zmid = zmid / weigh
        call image(xmid,ymid,zmid)
        docompute = .false.
        if ((zmid.ge.zbegproc(rank+1)).and.
     $    (zmid.lt.zendproc(rank+1)).and.(ymid.ge.ybegproc(rank+1))
     $    .and.(ymid.lt.yendproc(rank+1)).and.
     $    (xmid.ge.xbegproc(rank+1)).and.(xmid.lt.xendproc(rank+1)))
     $       then
           docompute = .true.
        end if
        if (docompute) then
          nratloc = nratloc + 1
          ratglob(nratloc) = i
        end if
      end do
c
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine chkangle  --  eliminate redundant constraints  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "chkangle" tests angles to be constrained for their presence
c     in small rings and removes constraints that are redundant
c
c     note this version correctly handles isolated small rings,
c     but should remove one additional redundant constraint for
c     each ring fusion
c
c
      subroutine chkangle (ia,ib,ic)
      use couple
      use freeze
      use ring
      implicit none
      integer i,j,k
      integer ia,ib,ic
      integer id,ie,imin
      logical remove
c
c
c     all internal angles of 3-membered rings are redundant
c
      remove = .false.
      if (nring3 .ne. 0) then
         do i = 1, n12(ia)
            j = i12(i,ia)
            if (j .eq. ic)  remove = .true.
         end do
      end if
c
c     for 4-membered rings, two internal angles are redundant
c
      if (nring4 .ne. 0) then
         do i = 1, n12(ia)
            id = i12(i,ia)
            if (id .ne. ib) then
               do j = 1, n12(id)
                  k = i12(j,id)
                  if (k .eq. ic) then
                     imin = min(ia,ib,ic,id)
                     if (ib .eq. imin)  remove = .true.
                     if (id .eq. imin)  remove = .true.
                  end if
               end do
            end if
         end do
      end if
c
c     for 5-membered rings, one internal angle is redundant
c
      if (nring5 .ne. 0) then
         do i = 1, n12(ia)
            id = i12(i,ia)
            if (id.ne.ib .and. id.ne.ic) then
               do j = 1, n12(ic)
                  ie = i12(j,ic)
                  if (ie.ne.ib .and. ie.ne.ia) then
                     do k = 1, n12(id)
                        if (i12(k,id) .eq. ie) then
                           imin = min(ia,ib,ic,id,ie)
                           if (ib .eq. imin)  remove = .true.
                        end if
                     end do
                  end if
               end do
            end if
         end do
      end if
c
c     remove the constraint from the list if it is redundant
c
      if (remove)  nrat = nrat - 1
      return
      end
c
c     subroutine alloc_shared_shake : allocate shared memory pointers for rattle
c     parameter arrays
c
      subroutine alloc_shared_shake
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
      use atoms
      use domdec
      use freeze
      use mpi
      use tinMemory
      implicit none
      
    !   if (associated(iratx).and.size(iratx).eq.n) return
        if (allocated(iratx) .and. size(iratx).eq.n) return
c
    !   call shmem_request(iratx,winiratx,[n])
    !   call shmem_request(kratx,winkratx,[n])
    !   call shmem_request(irat,winirat,[2,n])
    !   call shmem_request(krat,winkrat,[n])
    !   call shmem_request(ratimage,winratimage,[n])
        ! call prmem_request(iratx,n)
        ! call prmem_request(kratx,n)
        ! call prmem_request(irat,2,n)
        ! call prmem_request(krat,n)
        ! call prmem_request(ratimage,n)
        allocate(iratx(n))
        allocate(kratx(n))
        allocate(irat(2,n))
        allocate(krat(n))
        allocate(ratimage(n))
      end
