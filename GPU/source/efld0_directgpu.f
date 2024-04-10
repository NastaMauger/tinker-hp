c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
#include "tinker_macro.h"
      module efld0_directgpu_inl
        use tinheader, only: ti_p, zeror, zerom
        use tintypes , only: real3, rpole_elt, real7
        implicit none
        include "erfcore_data.f.inc"
        contains
#include "image.f.inc"
#if defined(SINGLE) | defined(MIXED)
        include "erfcscore.f.inc"
#else
        include "erfcdcore.f.inc"
#endif
#include "pair_efld.inc.f"
      end module

c
c     Compute the direct space contribution to the permanent electric field.
c     Also compute the "p" field, which is used to
c     compute the energy according to the AMOEBA force field.
c
      subroutine efld0_directgpu(nrhs,ef)
      use atmlst   ,only: poleglobnl
      use atoms    ,only: x,y,z
      use couple
      use domdec   ,only: rank,loc,nbloc
      use ewald    ,only: aewald
      use efld0_directgpu_inl
      use inform   ,only: deb_Path
      use math     ,only: sqrtpi
      use mpole    ,only: npolebloc,ipole,rpole,poleloc,npolelocnl
      use neigh    ,only: nelst,elst,shortelst,nshortelst
      use polar    ,only: pdamp,thole
      use polgrp
      use potent   ,only: use_polarshortreal
      use polpot
      use shunt    ,only: cut2
      use utilgpu  ,only: dir_queue,rec_queue,def_queue,
     &                    warning,maxscaling,maxscaling1
#ifdef _OPENACC
     &                   ,rec_stream,dir_stream,stream_wait_async
     &                   ,rec_event
#endif
      use timestat ,only: timer_enter,timer_exit,timer_efld0_direct
      implicit none

      integer  ,intent(in)   :: nrhs
      real(t_p),intent(inout):: ef(:,:,:)

      integer i,iglob,iploc,kk
      !integer countsel
      integer ii,j,k,ksp,kd,iipole,kbis,kpole,kglob
      integer nn12,nn13,nn14,nn4,ntot,nn15
      integer nnp11,nnp12,nnp13,nnp14
      integer nnelst
      integer,pointer :: lst(:,:),nlst(:)    ! pointers to neighbor list
      real(t_p) pdi,pti
      real(t_p) xi,yi,zi,d2
      real(t_p) thole1,pgamma,damp
      real(t_p) alsq2, alsq2n
      real(t_p) half,one
      real(t_p) pscale,dscale
      real(t_p) d,bn1,bn2,sc3,sc5
      type(rpole_elt) ip,kp
      type(real3) fid,fip,fkd,fkp,pos
      integer   ipscal(maxscaling),idscal(maxscaling1)
      real(t_p) fpscal(maxscaling),fdscal(maxscaling1)
      real(t_p) dscalevec(5),pscalevec(5)
      character*10 mode

      parameter(half=0.5_ti_p)
      parameter(one =1.0_ti_p)

      if (use_polarshortreal) then
          lst =>  shortelst
         nlst => nshortelst
         mode = 'SHORTEWALD'
      else
          lst =>  elst
         nlst => nelst
         mode = 'EWALD'
      end if
      dscalevec(:) = [d1scale,d2scale,d3scale,d4scale,0.0_ti_p]
      pscalevec(:) = [0.0_ti_p,p2scale,p3scale,p4scale,p5scale]

!$acc enter data attach(lst,nlst) async(def_queue)
      call switch (mode)
c
      if (deb_Path)
     &   write(*,'(3x,a)') 'efld0_directgpu'
      call timer_enter( timer_efld0_direct )

      alsq2  = 2.0_ti_p * aewald**2
      alsq2n = 0.0_ti_p
      if (aewald .gt. 0.0_ti_p) alsq2n = 1.0_ti_p / (sqrtpi*aewald)
      def_queue = dir_queue

!$acc parallel loop gang vector_length(32)
!$acc&         copyin(dscalevec,pscalevec) present(ef)
!$acc&         default(present)
!$acc&         private(ksp,kd,idscal,fdscal,ipscal,fpscal,ip,nn4,nn15,
!$acc&  nnp11)
!$acc&         async(def_queue)
      MAINLOOP:
     &do ii = 1, npolelocnl
         iipole   = poleglobnl(ii)
         iglob    = ipole     (iipole)
         i        = loc       (iglob)
         iploc    = poleloc   (iipole)
         if ((i.eq.0).or.(i.gt.nbloc)) then
            print*, warning,'efld0'
            cycle MAINLOOP
         end if

         nnelst   = nelst(ii)
         if(nnelst.eq.0) cycle MAINLOOP
         !countsel = 0
         pdi      = pdamp(iipole)
         pti      = thole(iipole)
         xi       = x    (iglob) 
         yi       = y    (iglob) 
         zi       = z    (iglob) 

         ip%c     = rpole(01,iipole)
         ip%dx    = rpole(02,iipole)
         ip%dy    = rpole(03,iipole)
         ip%dz    = rpole(04,iipole)
         ip%qxx   = rpole(05,iipole)
         ip%qxy   = rpole(06,iipole)
         ip%qxz   = rpole(07,iipole)
         ip%qyy   = rpole(09,iipole)
         ip%qyz   = rpole(10,iipole)
         ip%qzz   = rpole(13,iipole)
c
c        default values
c
         ksp      = 0
         kd       = 0
c
c       get number of atoms  directly (1-2) or 1-3 or 1-4 bonded
c
         nn4   = 0
         nn15  = 0
         nnp11 = 0
         ntot  = numscal_n(iglob)
         nn12  = scalbeg_n(iglob)
         nnp14 = numscal_p(iglob)
         nnp12 = scalbeg_p(iglob)
         if (ntot.gt.maxscaling) 
     &      print*,'scaling array too short efldo_dir'
         if (nnp14.gt.maxscaling1) 
     &      print*,'pscaling array to short in efld0_dir',ii
c
c       fill scaling factor along kglob and interaction type
c
!$acc loop vector
         do j = 1,ntot
            ipscal(j) = allscal_n(nn12+j)
            fpscal(j) = pscalevec(typscal_n(nn12+j))
            if (int(typscal_n(nn12+j)).eq.4) then
!$acc atomic
               nn4 = nn4 +1
            end if
            if (int(typscal_n(nn12+j)).eq.5) then
!$acc atomic
               nn15 = nn15 +1
            end if
         end do

!$acc loop vector
         do j = 1,nnp14
            idscal(j) = allscal_p(nnp12+j)
            fdscal(j) = dscalevec(int(typscal_p(nnp12+j)))
            if (int(typscal_p(nnp12+j)).eq.1) then
!$acc atomic
               nnp11 = nnp11 +1
            end if
         end do
         nn13 = ntot - nn4 - nn15
         nn14 = ntot - nn15

!$acc loop vector private(fip,fid,fkp,fkd,kp)
         do k =  1, nnelst
            kpole = elst(k,ii)
            kbis  = poleloc(kpole)
            kglob = ipole  (kpole)

            if (kbis.gt.npolebloc) then 
               print*,warning,'efld0 neighbour',kpole
               cycle
            end if
            !countsel = countsel + 1
            pos%x  = x (kglob) - xi
            pos%y  = y (kglob) - yi
            pos%z  = z (kglob) - zi
            call image_inl(pos%x,pos%y,pos%z)

            d2 = pos%x**2 + pos%y**2 + pos%z**2
            if (d2.gt.cut2) cycle

            kp%c    = rpole( 1, kpole)
            kp%dx   = rpole( 2, kpole)
            kp%dy   = rpole( 3, kpole)
            kp%dz   = rpole( 4, kpole)
            kp%qxx  = rpole( 5, kpole)
            kp%qxy  = rpole( 6, kpole)
            kp%qxz  = rpole( 7, kpole)
            kp%qyy  = rpole( 9, kpole)
            kp%qyz  = rpole(10, kpole)
            kp%qzz  = rpole(13, kpole)

            thole1  = thole(kpole)
            damp    = pdi * pdamp(kpole)
            pgamma  = min( pti,thole1 )
c
c      set exclusion coefficients for connected atoms
c
            pscale = 1.0_ti_p
            dscale = 1.0_ti_p
            if (ksp<ntot) then
!$acc loop seq
               do j=1,ntot
                  if (ipscal(j)==kglob) then
                     pscale  = fpscal(j)
                     ! deal with 4-1 interaction
                     if (nn13.lt.j.and.j.le.nn14) then
                        do kk=1,nnp11
                           if (idscal(kk).eq.kglob) then
                              pscale = pscale*p41scale
                              exit
                           end if
                        end do
                     end if
!$acc atomic update
                     ksp = ksp+1
                     exit
                  end if
               end do
            end if
            if (kd<nnp14) then
!$acc loop seq
               do j=1,nnp14
                  if (idscal(j)==kglob) then
                     dscale  = fdscal(j)
!$acc atomic update
                     kd = kd+1
                     exit
                  end if
               end do
            end if

            call efld0_couple(d2,pos,ip,kp,alsq2,alsq2n,
     &                 aewald,damp,pgamma,dscale,pscale,
     &                 fid,fip,fkd,fkp,d,bn1,bn2,sc3,sc5,.false.)

!$acc atomic update
            ef(1,1,iploc) = ef(1,1,iploc) + fid%x
!$acc atomic update
            ef(2,1,iploc) = ef(2,1,iploc) + fid%y
!$acc atomic update
            ef(3,1,iploc) = ef(3,1,iploc) + fid%z
!$acc atomic update
            ef(1,2,iploc) = ef(1,2,iploc) + fip%x
!$acc atomic update
            ef(2,2,iploc) = ef(2,2,iploc) + fip%y
!$acc atomic update
            ef(3,2,iploc) = ef(3,2,iploc) + fip%z
!$acc atomic update
            ef(1,1,kbis)  = ef(1,1,kbis ) + fkd%x
!$acc atomic update       
            ef(2,1,kbis)  = ef(2,1,kbis ) + fkd%y
!$acc atomic update       
            ef(3,1,kbis)  = ef(3,1,kbis ) + fkd%z
!$acc atomic update       
            ef(1,2,kbis)  = ef(1,2,kbis ) + fkp%x
!$acc atomic update       
            ef(2,2,kbis)  = ef(2,2,kbis ) + fkp%y
!$acc atomic update       
            ef(3,2,kbis)  = ef(3,2,kbis ) + fkp%z
         enddo

      end do MAINLOOP
!$acc exit data detach(nlst,lst) async(def_queue)

      call timer_exit( timer_efld0_direct )
      end subroutine

      subroutine efld0_directgpu2(nrhs,ef)
      use atmlst  ,only: poleglobnl
      use atoms   ,only: x,y,z
      !use couple  ,only: i12,i13,i14,i15,n12,n13,n14,n15
      use domdec  ,only: rank,loc,nbloc
      use ewald   ,only: aewald
      use efld0_directgpu_inl
      use inform  ,only: deb_Path
      use interfaces ,only: efld0_direct_correct_scaling
      use math    ,only: sqrtpi
      use mpole   ,only: npolebloc,ipole,rpole,poleloc,npolelocnl
      use mutant  ,only: deflambda,elambda,mut
      use neigh   ,only: nelst,elst,shortelst,nshortelst
      use polar   ,only: pdamp,thole
      use potent  , only : use_polarshortreal,use_chgpen,use_lambdadyn
      !use polgrp  ,only: ip11,ip12,ip13,ip14,np11,np12,np13,np14
      use polpot  ,only: dpcorrect_ik,dpcorrect_scale,n_dpscale
      use shunt   ,only: cut2
      use utilgpu ,only: dir_queue,rec_queue,def_queue
#ifdef _OPENACC
     &                  ,dir_stream,stream_wait_async,
     &                   rec_stream,rec_event
#endif
      use timestat   ,only: timer_enter,timer_exit,timer_efld0_direct
      implicit none

      integer  ,intent(in)   :: nrhs
      ! shape (ef) = (/3,nrhs,npolebloc/)
      real(t_p),intent(inout):: ef(:,:,:)

      integer i,iglob,iploc,kk
      integer ii,j,k,ksp,kd,iipole,kbis,kpole,kglob
      integer nn12,nn13,nn14,ntot
      integer nnp11,nnp12,nnp13,nnp14
      integer nnelst
      integer,pointer :: lst(:,:),nlst(:)    ! pointers to neighbor list
      real(t_p) pdi,pti
      real(t_p) xi,yi,zi,d2
      real(t_p) pgamma,damp
      real(t_p) alsq2, alsq2n
      real(t_p) one
      real(t_p) pscale,dscale
      real(t_p) d,bn1,bn2,sc3,sc5
      type(rpole_elt) ip,kp
      type(real3) fid,fip,fkd,fkp,pos
      character*10 mode

      parameter(one =1.0_ti_p)

      if (use_chgpen) then
#ifdef _OPENACC
         __TINKER_FATAL__
#else
         call efld0_direct(nrhs,ef)
         return
#endif
      end if
         
      if (deb_Path)
     &   write(*,'(3x,a,L3)') 'efld0_directgpu2',use_polarshortreal
      call timer_enter( timer_efld0_direct )

      if (use_polarshortreal) then
          lst =>  shortelst
         nlst => nshortelst
         mode = 'SHORTEWALD'
      else
          lst =>  elst
         nlst => nelst
         mode = 'EWALD'
      end if
!$acc enter data attach(lst,nlst) async(def_queue)
      call switch (mode)
c
      def_queue = dir_queue
      pscale = 1.0
      dscale = 1.0
      alsq2  = 2.0_ti_p * aewald**2
      alsq2n = 0.0_ti_p
      if (aewald .gt. 0.0_ti_p) alsq2n = 1.0_ti_p / (sqrtpi*aewald)

!$acc parallel loop gang vector_length(32)
!$acc&         present(ef)
!$acc&         present(poleglobnl,ipole,loc,pdamp,
!$acc&  thole,x,y,z,rpole,nelst,elst,poleloc,mut,deflambda)
!$acc&         private(ip)
!$acc&         async(def_queue)
      MAINLOOP:
     &do ii = 1, npolelocnl
         iipole  = poleglobnl(ii)
         iglob   = ipole     (iipole)
         iploc   = poleloc   (iipole)

         nnelst  = nelst(ii)
         if(nnelst.eq.0) cycle MAINLOOP
         pdi     = pdamp(iipole)
         pti     = thole(iipole)
         xi      = x    (iglob) 
         yi      = y    (iglob) 
         zi      = z    (iglob) 

         ip%c    = rpole(01,iipole)
         ip%dx   = rpole(02,iipole)
         ip%dy   = rpole(03,iipole)
         ip%dz   = rpole(04,iipole)
         ip%qxx  = rpole(05,iipole)
         ip%qxy  = rpole(06,iipole)
         ip%qxz  = rpole(07,iipole)
         ip%qyy  = rpole(09,iipole)
         ip%qyz  = rpole(10,iipole)
         ip%qzz  = rpole(13,iipole)

!$acc loop vector private(kp,fip,fid,fkp,fkd)
         NEIGHBORS LOOP:
     &   do k =  1, nnelst
            kpole  = elst(k,ii)
            kbis   = poleloc(kpole)
            kglob  = ipole  (kpole)

            pos%x  = x (kglob) - xi
            pos%y  = y (kglob) - yi
            pos%z  = z (kglob) - zi
            call image_inl(pos%x,pos%y,pos%z)

            d2     = pos%x**2 + pos%y**2 + pos%z**2
            if (d2.gt.cut2) cycle

            kp%c   = rpole( 1, kpole)
            kp%dx  = rpole( 2, kpole)
            kp%dy  = rpole( 3, kpole)
            kp%dz  = rpole( 4, kpole)
            kp%qxx = rpole( 5, kpole)
            kp%qxy = rpole( 6, kpole)
            kp%qxz = rpole( 7, kpole)
            kp%qyy = rpole( 9, kpole)
            kp%qyz = rpole(10, kpole)
            kp%qzz = rpole(13, kpole)

            damp   = pdi * pdamp(kpole)
            pgamma = min( pti,thole(kpole) )

            call efld0_couple(d2,pos,ip,kp,alsq2,alsq2n,
     &                 aewald,damp,pgamma,1.0_ti_p,1.0_ti_p,
     &                 fid,fip,fkd,fkp,d,bn1,bn2,sc3,sc5,.false.)
c
c     get vector of derivative of direct field wrt lambda for lambda-dynamics
c
            if ((use_lambdadyn).and.(elambda.gt.0)) then
              if (mut(kglob)) then
!$acc atomic
                deflambda(1,1,iploc) = deflambda(1,1,iploc)
     $               + fid%x/elambda
!$acc atomic
                deflambda(2,1,iploc) = deflambda(2,1,iploc)
     $               + fid%y/elambda
!$acc atomic
                deflambda(3,1,iploc) = deflambda(3,1,iploc)
     $               + fid%z/elambda
!$acc atomic
                deflambda(1,2,iploc) = deflambda(1,2,iploc)
     $               + fip%x/elambda
!$acc atomic
                deflambda(2,2,iploc) = deflambda(2,2,iploc)
     $               + fip%y/elambda
!$acc atomic
                deflambda(3,2,iploc) = deflambda(3,2,iploc)
     $               + fip%z/elambda
              end if
              if (mut(iglob)) then
!$acc atomic
                deflambda(1,1,kbis) = deflambda(1,1,kbis)
     $               + fkd%x/elambda
!$acc atomic
                deflambda(2,1,kbis) = deflambda(2,1,kbis)
     $               + fkd%y/elambda
!$acc atomic
                deflambda(3,1,kbis) = deflambda(3,1,kbis)
     $               + fkd%z/elambda
!$acc atomic
                deflambda(1,2,kbis) = deflambda(1,2,kbis)
     $               + fkp%x/elambda
!$acc atomic
                deflambda(2,2,kbis) = deflambda(2,2,kbis)
     $               + fkp%y/elambda
!$acc atomic
                deflambda(3,2,kbis) = deflambda(3,2,kbis)
     $               + fkp%z/elambda
              end if
            end if

!$acc atomic
            ef(1,1,iploc) = ef(1,1,iploc) + fid%x
!$acc atomic
            ef(2,1,iploc) = ef(2,1,iploc) + fid%y
!$acc atomic
            ef(3,1,iploc) = ef(3,1,iploc) + fid%z
!$acc atomic
            ef(1,2,iploc) = ef(1,2,iploc) + fip%x
!$acc atomic
            ef(2,2,iploc) = ef(2,2,iploc) + fip%y
!$acc atomic
            ef(3,2,iploc) = ef(3,2,iploc) + fip%z
!$acc atomic
            ef(1,1,kbis)  = ef(1,1,kbis ) + fkd%x
!$acc atomic
            ef(2,1,kbis)  = ef(2,1,kbis ) + fkd%y
!$acc atomic
            ef(3,1,kbis)  = ef(3,1,kbis ) + fkd%z
!$acc atomic
            ef(1,2,kbis)  = ef(1,2,kbis ) + fkp%x
!$acc atomic
            ef(2,2,kbis)  = ef(2,2,kbis ) + fkp%y
!$acc atomic
            ef(3,2,kbis)  = ef(3,2,kbis ) + fkp%z
         end do NEIGHBORS LOOP

      end do MAINLOOP
!$acc exit data detach(nlst,lst) async(def_queue)

      call efld0_direct_correct_scaling(ef)

      call timer_exit( timer_efld0_direct )
      end subroutine

      subroutine efld0_directgpu3(nrhs,ef)
#ifdef _OPENACC
      use atmlst  ,only: poleglobnl
      use atoms   ,only: x,y,z,n
      use cell
      use chgpen  ,only: pcore,pval,palpha
      use domdec  ,only: xbegproc,ybegproc,zbegproc
     &            ,nproc,rank,xendproc,yendproc,zendproc
     &            ,nbloc,loc
      use ewald   ,only: aewald
      use efld0_directgpu_inl
      use efld0_cpencu
      use inform  ,only: deb_Path,minmaxone
      use interfaces ,only: efld0_direct_correct_scaling
     &               , cu_efld0_direct
      use math    ,only: sqrtpi
      use mplpot  ,only: pentyp_i
      use mpole   ,only: npolebloc,ipole,rpole,poleloc,npolelocnl
     &            , npolelocnlb_pair,npolelocnlb
     &            ,  npnlb2=>npolelocnlb2_pair
     &            , nspnlb2=>nshortpolelocnlb2_pair
      use mutant  ,only: deflambda,elambda,mut,mutInt
      use neigh   , only : ipole_s=>celle_glob,pglob_s=>celle_pole
     &            , ploc_s=>celle_ploc, ieblst_s=>ieblst
     &            , iseblst_s=>ishorteblst, seblst_s=>shorteblst
     &            , eblst_s=>eblst, x_s=>celle_x, y_s=>celle_y
     &            , z_s=>celle_z
      use polar   ,only: pdamp,thole,polarity,dirdamp
      use polpot  ,only: dpcorrect_ik,dpcorrect_scale,n_dpscale
      use potent  , only : use_polarshortreal,use_chgpen,use_lambdadyn
      !use polgrp  ,only: ip11,ip12,ip13,ip14,np11,np12,np13,np14
      use polpot  ,only: dpcorrect_ik,dpcorrect_scale,n_dpscale
     &            ,use_thole,use_dirdamp
      use shunt   ,only: cut2
      use utilcu     ,only: check_launch_kernel,BLOCK_DIM
      use utilcomm,only: no_commdir
      use utilgpu ,only: dir_queue,rec_queue,def_queue,get_GridDim
     &                  ,dir_stream,def_stream,stream_wait_async,
     &                   rec_stream,rec_event
      use utils   ,only: associate_ptr
      use timestat,only: timer_enter,timer_exit,timer_efld0_direct
      use tinMemory    ,only: mipk
      use tmatxb_pmecu ,only: efld0_direct_scaling_cu
      implicit none

      integer  ,intent(in)   :: nrhs
      ! shape (ef) = (/3,nrhs,npolebloc/)
      real(t_p),intent(inout):: ef(:,:,:)

      integer i,j,k,iglob,iploc,kk,start_lst,gS,idirdamp
      integer(mipk) n_
      real(t_p) alsq2, alsq2n
      real(t_p) p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
      real(t_p),pointer:: gamma(:)
      character*11 mode
      character*64 rtami
c
      if (deb_Path) then
         rtami = 'efld0_directgpu3'//merge(' THOLE','',use_thole)
     &         //merge(' DIRDAMP','',use_dirdamp)
     &         //merge(' CHGPEN','',use_chgpen.and..not.use_thole)
     &         //merge(' SHORT','',use_polarshortreal)
         write(*,'(3x,a)') rtami
      end if
      call timer_enter( timer_efld0_direct )

      mode      = merge('SHORTEWALD','EWALD',use_polarshortreal)
      call switch (mode)

      p_xbeg    = xbegproc(rank+1)
      p_xend    = xendproc(rank+1)
      p_ybeg    = ybegproc(rank+1)
      p_yend    = yendproc(rank+1)
      p_zbeg    = zbegproc(rank+1)
      p_zend    = zendproc(rank+1)
      start_lst = 2*npolelocnlb_pair + 1

      alsq2     = 2.0_ti_p * aewald**2
      alsq2n    = merge(1.0_re_p/real(sqrtpi*aewald,r_p),zerom
     &                 ,aewald.gt.zeror)
      def_queue = dir_queue
      def_stream= dir_stream

      thole_cpen: if (use_thole) then
         n_ = size(dirdamp)
         if (use_dirdamp) then
            idirdamp = 1
            call associate_ptr(gamma,dirdamp,n_)
         else
            idirdamp = 0
            call associate_ptr(gamma,thole,n_)
         end if

      range0: if (use_polarshortreal) then
!$acc host_data use_device(ipole_s,pglob_s,ploc_s,iseblst_s,seblst_s,
!$acc&    mutInt,x_s,y_s,z_s,pdamp,rpole,gamma,polarity,ef,deflambda)

      call cu_efld0_direct    !Find his interface inside MOD_inteface
     &     (ipole_s,pglob_s,ploc_s,iseblst_s,seblst_s(start_lst),mutInt
     &     ,x_s,y_s,z_s,pdamp,gamma,polarity,rpole,ef,deflambda
     &     ,npolelocnlb,nspnlb2,npolebloc,n,nproc,idirdamp
     &     ,use_lambdadyn,.not.no_commdir
     &     ,cut2,alsq2,alsq2n,aewald,elambda
     &     ,xcell,ycell,zcell,xcell2,ycell2,zcell2
     &     ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend,def_stream)

!$acc end host_data
      else
!$acc host_data use_device(ipole_s,pglob_s,ploc_s,ieblst_s,eblst_s,
!$acc&    mutInt,x_s,y_s,z_s,pdamp,rpole,gamma,polarity,ef,deflambda)

      call cu_efld0_direct    !Find his interface inside MOD_inteface
     &     (ipole_s,pglob_s,ploc_s,ieblst_s,eblst_s(start_lst),mutInt
     &     ,x_s,y_s,z_s,pdamp,gamma,polarity,rpole,ef,deflambda
     &     ,npolelocnlb,npnlb2,npolebloc,n,nproc,idirdamp
     &     ,use_lambdadyn,.not.no_commdir
     &     ,cut2,alsq2,alsq2n,aewald,elambda
     &     ,xcell,ycell,zcell,xcell2,ycell2,zcell2
     &     ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend,def_stream)

!$acc end host_data
      end if range0
      
      if (n_dpscale.gt.0) then
         gS = get_GridDim(n_dpscale,BLOCK_DIM,1)
!$acc host_data use_device(dpcorrect_ik,dpcorrect_scale,poleloc,ipole
!$acc&         ,mut,pdamp,gamma,rpole,x,y,z,ef,deflambda)
         call efld0_direct_scaling_cu<<<gS,BLOCK_DIM,0,def_stream>>>
     &       (dpcorrect_ik,dpcorrect_scale,poleloc,ipole,mut,pdamp,gamma
     &       ,x,y,z,rpole,ef,deflambda,n_dpscale,n,npolebloc,use_dirdamp
     &       ,use_lambdadyn,elambda,cut2,aewald,alsq2,alsq2n)
!$acc end host_data
         call check_launch_kernel( 'efld0_direct_scaling_cu' )
      end if

      else if (use_chgpen) then

      range1: if (use_polarshortreal) then
         gS = get_GridDim(nspnlb2,BLOCK_DIM)
!$acc host_data use_device(ipole_s,pglob_s,ploc_s,iseblst_s,seblst_s
!$acc&    ,x_s,y_s,z_s,pcore,pval,palpha,rpole,ef
!$acc&    ,dpcorrect_ik,dpcorrect_scale,poleloc,ipole,x,y,z)
         call efld0_cpen_kcu<<<gS,BLOCK_DIM,0,def_stream>>>
     &       (ipole_s,pglob_s,ploc_s,iseblst_s,seblst_s(start_lst)
     &       ,x_s,y_s,z_s,pcore,pval,palpha,rpole,npolelocnl
     &       ,npolelocnlb,nspnlb2,npolebloc,n,pentyp_i,cut2,aewald
     &       ,ef
     &       ,xcell,ycell,zcell,xcell2,ycell2,zcell2
     &       ,dpcorrect_ik,dpcorrect_scale,poleloc,ipole,x,y,z,n_dpscale
     &       )
!$acc end host_data
      else
         gS = get_GridDim(npnlb2,BLOCK_DIM)
!$acc host_data use_device(ipole_s,pglob_s,ploc_s,ieblst_s,eblst_s
!$acc&    ,x_s,y_s,z_s,pcore,pval,palpha,rpole,ef
!$acc&    ,dpcorrect_ik,dpcorrect_scale,poleloc,ipole,x,y,z)
         call efld0_cpen_kcu<<<gS,BLOCK_DIM,0,def_stream>>>
     &       (ipole_s,pglob_s,ploc_s,ieblst_s,eblst_s(start_lst)
     &       ,x_s,y_s,z_s,pcore,pval,palpha,rpole,npolelocnl
     &       ,npolelocnlb,npnlb2,npolebloc,n,pentyp_i,cut2,aewald
     &       ,ef
     &       ,xcell,ycell,zcell,xcell2,ycell2,zcell2
     &       ,dpcorrect_ik,dpcorrect_scale,poleloc,ipole,x,y,z,n_dpscale
     &       )
!$acc end host_data
      end if range1

      end if thole_cpen

      if (use_lambdadyn.and.elambda.gt.0.0) then
!$acc parallel loop async(def_queue) collapse(3) present(deflambda)
         do i =1,npolebloc; do j=1,2; do k=1,3;
         if (deflambda(k,j,i).ne.0.0)
     &      deflambda(k,j,i) = deflambda(k,j,i)/elambda
         end do; end do; end do
      end if

      call timer_exit( timer_efld0_direct )
#else
      print 100
 100  format('eld0_directgpu3 is a specific device routine',/,
     &       'you are not supposed to get inside with your compile',
     &       'type.')
      __TINKER_FATAL__
#endif
      end subroutine

c""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
c---------------------------
c     Otf routines for efld0
c---------------------------
c""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

      subroutine otf_dc_efld0_directgpu2(nrhs,ef)
      use atmlst  ,only: poleglobnl
      use atoms   ,only: x,y,z
      !use couple  ,only: i12,i13,i14,i15,n12,n13,n14,n15
      use domdec  ,only: rank,loc,nbloc
      use divcon
      use ewald   ,only: aewald
      use efld0_directgpu_inl
      use inform  ,only: deb_Path
      use interfaces ,only: efld0_direct_correct_scaling
      use math    ,only: sqrtpi
      use mpole   ,only: npolebloc,ipole,rpole,poleloc,npolelocnl
      use neigh   ,only: nelst,elst,shortelst,nshortelst
      use polar   ,only: pdamp,thole,polarity
      use potent  , only : use_polarshortreal
      !use polgrp  ,only: ip11,ip12,ip13,ip14,np11,np12,np13,np14
      use polpot  ,only: dpcorrect_ik,dpcorrect_scale,n_dpscale
      use shunt   ,only: cut2
      use utilgpu ,only: dir_queue,rec_queue,def_queue
#ifdef _OPENACC
     &                  ,dir_stream,stream_wait_async,
     &                   rec_stream,rec_event
#endif
      use timestat   ,only: timer_enter,timer_exit,timer_efld0_direct
      implicit none

      integer  ,intent(in)   :: nrhs
      ! shape (ef) = (/3,nrhs,npolebloc/)
      real(t_p),intent(inout):: ef(:,:,:)

      integer i,iglob,iploc,kk
      integer atii,atkk,l,cofst1,cofst2,cofst3,rofst
      integer maxrow,ikof
      integer ii,j,k,ksp,kd,iipole,kbis,kpole,kglob
      integer nn12,nn13,nn14,ntot
      integer nnp11,nnp12,nnp13,nnp14
      integer nnelst
      integer,pointer :: lst(:,:),nlst(:)    ! pointers to neighbor list
      real(t_p) pdi,pti
      real(t_p) xi,yi,zi,d2
      real(t_p) pgamma,damp
      real(t_p) alsq2, alsq2n
      real(t_p) one,invpol
      real(t_p) pscale,dscale
      real(t_p) d,bn1,bn2,sc3,sc5
      type(rpole_elt) ip,kp
      type(real3) fid,fip,fkd,fkp,pos
      character*10 mode

      parameter(one =1.0_ti_p)

      if (use_polarshortreal) then
          lst =>  shortelst
         nlst => nshortelst
         mode = 'SHORTEWALD'
      else
          lst =>  elst
         nlst => nelst
         mode = 'EWALD'
      end if
      call switch (mode)
c
      if (deb_Path)
     &   write(*,'(3x,a)') 'oft_dc_efld0_directgpu2'
      call timer_enter( timer_efld0_direct )

      def_queue = dir_queue
      pscale = 1.0
      dscale = 1.0
      alsq2  = 2.0_ti_p * aewald**2
      alsq2n = 0.0_ti_p
      if (aewald .gt. 0.0_ti_p) alsq2n = 1.0_ti_p / (sqrtpi*aewald)

!$acc enter data attach(lst,nlst) async(def_queue)
#ifdef _OPENACC
      if (dir_queue.ne.rec_queue) then
!!$acc wait(rec_queue) async(rec_queue)
         call stream_wait_async(rec_stream,dir_stream,rec_event)
      end if
#endif

c
c     diagonal of Z mat.
c
!$acc parallel loop default(present) async(def_queue)
      do ii = 1, npolelocnl
        iipole = poleglobnl(ii)
        iglob  = ipole(iipole)
        l      = grplst(iglob)
        if(l .gt. 0) then
          rofst = (atmofst(iglob) - 1)*3
          if (polarity(iipole).ne.0_ti_p) then
            invpol = 1.0_ti_p/polarity(iipole)
          else
            invpol = 1000.0_ti_p
          end if
          maxrow=npergrp(l)*3
          cofst1 = rofst + 1
          cofst2 = rofst + 2
          cofst3 = rofst + 3
          cofst1 = (cofst1-1)*(2*maxrow-cofst1)/2
          cofst2 = (cofst2-1)*(2*maxrow-cofst2)/2
          cofst3 = (cofst3-1)*(2*maxrow-cofst3)/2 
          zmat(rofst+1+cofst1+kofst(l)) = invpol
          zmat(rofst+2+cofst2+kofst(l)) = invpol
          zmat(rofst+3+cofst3+kofst(l)) = invpol
        end if
      end do

!$acc parallel loop gang vector_length(32)
!$acc&         present(ef)
!$acc&         present(poleglobnl,ipole,loc,pdamp,
!$acc&  thole,x,y,z,rpole,nlst,lst,poleloc,zmat,
!$acc&  grplst,atmofst,kofst,npergrp)
!$acc&         async(def_queue)
      MAINLOOP:
     &do ii = 1, npolelocnl
         iipole  = poleglobnl(ii)
         iglob   = ipole     (iipole)
         iploc   = poleloc   (iipole)
         l       = grplst(iglob)

         nnelst  = nlst(ii)
         if(nnelst.eq.0) cycle MAINLOOP
         pdi     = pdamp(iipole)
         pti     = thole(iipole)
         xi      = x    (iglob) 
         yi      = y    (iglob) 
         zi      = z    (iglob) 

         ip%c    = rpole(01,iipole)
         ip%dx   = rpole(02,iipole)
         ip%dy   = rpole(03,iipole)
         ip%dz   = rpole(04,iipole)
         ip%qxx  = rpole(05,iipole)
         ip%qxy  = rpole(06,iipole)
         ip%qxz  = rpole(07,iipole)
         ip%qyy  = rpole(09,iipole)
         ip%qyz  = rpole(10,iipole)
         ip%qzz  = rpole(13,iipole)

!$acc loop vector private(pos,kp,fid,fip,fkd,fkp)
         NEIGHBORS LOOP:
     &   do k =  1, nnelst
            kpole  = lst(k,ii)
            kbis   = poleloc(kpole)
            kglob  = ipole  (kpole)

            pos%x  = x (kglob) - xi
            pos%y  = y (kglob) - yi
            pos%z  = z (kglob) - zi
            call image_inl(pos%x,pos%y,pos%z)

            d2     = pos%x**2 + pos%y**2 + pos%z**2
            if (d2.gt.cut2) cycle

            kp%c   = rpole( 1, kpole)
            kp%dx  = rpole( 2, kpole)
            kp%dy  = rpole( 3, kpole)
            kp%dz  = rpole( 4, kpole)
            kp%qxx = rpole( 5, kpole)
            kp%qxy = rpole( 6, kpole)
            kp%qxz = rpole( 7, kpole)
            kp%qyy = rpole( 9, kpole)
            kp%qyz = rpole(10, kpole)
            kp%qzz = rpole(13, kpole)

            damp   = pdi * pdamp(kpole)
            pgamma = min( pti,thole(kpole) )

            call efld0_couple(d2,pos,ip,kp,alsq2,alsq2n,
     &                 aewald,damp,pgamma,1.0_ti_p,1.0_ti_p,
     &                 fid,fip,fkd,fkp,d,bn1,bn2,sc3,sc5,.false.)

!$acc atomic update
            ef(1,1,iploc) = ef(1,1,iploc) + fid%x
!$acc atomic update
            ef(2,1,iploc) = ef(2,1,iploc) + fid%y
!$acc atomic update
            ef(3,1,iploc) = ef(3,1,iploc) + fid%z
!$acc atomic update
            ef(1,2,iploc) = ef(1,2,iploc) + fip%x
!$acc atomic update
            ef(2,2,iploc) = ef(2,2,iploc) + fip%y
!$acc atomic update
            ef(3,2,iploc) = ef(3,2,iploc) + fip%z
!$acc atomic update
            ef(1,1,kbis)  = ef(1,1,kbis ) + fkd%x
!$acc atomic update       
            ef(2,1,kbis)  = ef(2,1,kbis ) + fkd%y
!$acc atomic update       
            ef(3,1,kbis)  = ef(3,1,kbis ) + fkd%z
!$acc atomic update       
            ef(1,2,kbis)  = ef(1,2,kbis ) + fkp%x
!$acc atomic update       
            ef(2,2,kbis)  = ef(2,2,kbis ) + fkp%y
!$acc atomic update       
            ef(3,2,kbis)  = ef(3,2,kbis ) + fkp%z

            if (l.eq.grplst(kglob) .and. l.ne.-1) then

              bn1    = bn1 -          (1.0_ti_p - sc3)/ (d*d2)
              bn2    = bn2 - 3.0_ti_p*(1.0_ti_p - sc5)/ (d*d2*d2)
              atii   = (atmofst(iglob) - 1)*3
              atkk   = (atmofst(kglob) - 1)*3
              if (l.eq.0) print*,kglob
              maxrow = npergrp(l)*3
              ikof   = kofst(l)
              if(atii .lt. atkk) then
                cofst1 = atii + 1
                cofst2 = atii + 2
                cofst3 = atii + 3
                rofst  = atkk
              else
                cofst1 = atkk + 1
                cofst2 = atkk + 2
                cofst3 = atkk + 3
                rofst  = atii
              end if

              cofst1 = (cofst1-1)*(2*maxrow-cofst1)/2
              cofst2 = (cofst2-1)*(2*maxrow-cofst2)/2
              cofst3 = (cofst3-1)*(2*maxrow-cofst3)/2 
              zmat(rofst+1+cofst1+ikof) =  bn1 - bn2*pos%x*pos%x
              zmat(rofst+2+cofst1+ikof) =      - bn2*pos%x*pos%y
              zmat(rofst+3+cofst1+ikof) =      - bn2*pos%x*pos%z
              zmat(rofst+1+cofst2+ikof) =      - bn2*pos%x*pos%y
              zmat(rofst+2+cofst2+ikof) =  bn1 - bn2*pos%y*pos%y
              zmat(rofst+3+cofst2+ikof) =      - bn2*pos%y*pos%z
              zmat(rofst+1+cofst3+ikof) =      - bn2*pos%x*pos%z
              zmat(rofst+2+cofst3+ikof) =      - bn2*pos%y*pos%z
              zmat(rofst+3+cofst3+ikof) =  bn1 - bn2*pos%z*pos%z

            end if
         end do NEIGHBORS LOOP

      end do MAINLOOP

!$acc exit data detach(nlst,lst) async(def_queue)

      call efld0_direct_correct_scaling(ef)

      call timer_exit( timer_efld0_direct )
      end subroutine

      subroutine otf_dc_efld0_directgpu3(nrhs,ef)
      use atmlst  ,only: poleglobnl
      use atoms   ,only: x,y,z,n
      use cell
      use divcon
      use domdec  ,only: xbegproc,ybegproc,zbegproc
     &            ,nproc,rank,xendproc,yendproc,zendproc
     &            ,nbloc,loc
      use ewald   ,only: aewald
      use efld0_directgpu_inl
      use inform  ,only: deb_Path
      use interfaces ,only: efld0_direct_correct_scaling
#ifdef _CUDA
     &               , cu_otfdc_efld0_direct
#endif
      use math    ,only: sqrtpi
      use mpole   ,only: npolebloc,ipole,rpole,poleloc,npolelocnl
     &            , npolelocnlb_pair,npolelocnlb2_pair,npolelocnlb
     &            , nspnlb2=>nshortpolelocnlb2_pair
      use neigh   , only : ipole_s=>celle_glob,pglob_s=>celle_pole
     &            , ploc_s=>celle_ploc, ieblst_s=>ieblst
     &            , iseblst_s=>ishorteblst, seblst_s=>shorteblst
     &            , eblst_s=>eblst, x_s=>celle_x, y_s=>celle_y
     &            , z_s=>celle_z
      use polar   ,only: pdamp,thole,polarity
      use potent  , only : use_polarshortreal
      !use polgrp  ,only: ip11,ip12,ip13,ip14,np11,np12,np13,np14
      use polpot  ,only: dpcorrect_ik,dpcorrect_scale,n_dpscale
      use utilcomm,only: no_commdir
      use shunt   ,only: cut2
      use utilgpu ,only: dir_queue,rec_queue,def_queue
#ifdef _OPENACC
     &                  ,dir_stream,def_stream,stream_wait_async,
     &                   rec_stream,rec_event
#endif
      use timestat   ,only: timer_enter,timer_exit,timer_efld0_direct
      implicit none

      integer  ,intent(in)   :: nrhs
      ! shape (ef) = (/3,nrhs,npolebloc/)
      real(t_p),intent(inout):: ef(:,:,:)

      integer i,ii,iipole,l,iglob,iploc,kk,start_lst
      integer rofst,cofst1,cofst2,cofst3, maxrow
      real(t_p) alsq2, alsq2n, invpol
      real(t_p) p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
      character*10 mode

      if (use_polarshortreal) then
         mode = 'SHORTEWALD'
      else
         mode = 'EWALD'
      end if
      call switch (mode)
c
      if (deb_Path) write(*,'(3x,a)') 'otf_dc_efld0_directgpu3'
      call timer_enter( timer_efld0_direct )

      p_xbeg = xbegproc(rank+1)
      p_xend = xendproc(rank+1)
      p_ybeg = ybegproc(rank+1)
      p_yend = yendproc(rank+1)
      p_zbeg = zbegproc(rank+1)
      p_zend = zendproc(rank+1)
      start_lst = 2*npolelocnlb_pair + 1

      alsq2  = 2.0_ti_p * aewald**2
      alsq2n = 0.0_ti_p
      if (aewald .gt. 0.0_ti_p) alsq2n = 1.0_ti_p / (sqrtpi*aewald)
      def_queue = dir_queue

#ifdef _OPENACC
      def_stream = dir_stream
      def_queue  = dir_queue
      if (dir_queue.ne.rec_queue) then
!!$acc wait(rec_queue) async(rec_queue)
         call stream_wait_async(rec_stream,dir_stream,rec_event)
      end if
#endif
c
c     diagonal of Z mat.
c
!$acc parallel loop default(present) async(def_queue)
      do ii = 1, npolelocnl
         iipole = poleglobnl(ii)
         iglob  = ipole(iipole)
         l      = grplst(iglob)
         if (l .gt. 0) then
            rofst = (atmofst(iglob) - 1)*3
            maxrow=npergrp(l)*3
            if (polarity(iipole).ne.0_ti_p) then
              invpol = 1.0_ti_p/polarity(iipole)
            else
              invpol = 1000.0_ti_p
            end if
            cofst1 = (rofst+0)*(2*maxrow-rofst-1)/2
            cofst2 = (rofst+1)*(2*maxrow-rofst-2)/2
            cofst3 = (rofst+2)*(2*maxrow-rofst-3)/2
            zmat(rofst+1+cofst1+kofst(l)) = invpol
            zmat(rofst+2+cofst2+kofst(l)) = invpol
            zmat(rofst+3+cofst3+kofst(l)) = invpol
         end if
      end do

#ifdef _CUDA
      if (use_polarshortreal) then
!$acc host_data use_device(ipole_s,pglob_s,ploc_s,grplst,atmofst
!$acc&    ,npergrp,kofst,iseblst_s,seblst_s,x_s,y_s,z_s,pdamp
!$acc&    ,rpole,thole,polarity,ef,zmat)

      call cu_otfdc_efld0_direct    !Find his interface inside MOD_inteface
     &     (ipole_s,pglob_s,ploc_s,grplst,atmofst,npergrp,kofst
     &     ,iseblst_s,seblst_s(start_lst)
     &     ,x_s,y_s,z_s,pdamp,thole,polarity,rpole
     &     ,ef,zmat
     &     ,npolelocnlb,nspnlb2,npolebloc,n,nproc
     &     ,cut2,alsq2,alsq2n,aewald
     &     , xcell, ycell, zcell,xcell2,ycell2,zcell2
     &     ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend,def_stream)

!$acc end host_data
      else
!$acc host_data use_device(ipole_s,pglob_s,ploc_s,grplst,atmofst
!$acc&    ,npergrp,kofst,ieblst_s,eblst_s,x_s,y_s,z_s,pdamp,rpole
!$acc&    ,thole,polarity,ef,zmat)

      call cu_otfdc_efld0_direct    !Find his interface inside MOD_inteface
     &     (ipole_s,pglob_s,ploc_s,grplst,atmofst,npergrp,kofst
     &     ,ieblst_s,eblst_s(start_lst)
     &     ,x_s,y_s,z_s,pdamp,thole,polarity,rpole
     &     ,ef,zmat
     &     ,npolelocnlb,npolelocnlb2_pair,npolebloc,n,nproc
     &     ,cut2,alsq2,alsq2n,aewald
     &     , xcell, ycell, zcell,xcell2,ycell2,zcell2
     &     ,p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend,def_stream)


!$acc end host_data
      end if
#else
      print 100
 100  format('eld0_directgpu3 is a specific device routine',/,
     &       'you are not supposed to get inside with your compile',
     &       'type.')
      call fatal
#endif

      call efld0_direct_correct_scaling(ef)

      call timer_exit( timer_efld0_direct )
      end subroutine


c""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
c     --------------------------
c     Correcting scaling interactions routines
c     --------------------------
c""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

      subroutine efld0_direct_correct_scaling(ef)
      use atmlst   ,only: poleglobnl
      use atoms    ,only: x,y,z
      use domdec   ,only: rank,loc,nbloc
      use ewald    ,only: aewald
      use efld0_directgpu_inl
      use inform   ,only: deb_Path
      use math     ,only: sqrtpi
      use mpole    ,only: npolebloc,ipole,rpole,poleloc,npolelocnl
      use mutant   ,only: deflambda,elambda,mut
      use neigh    ,only: nelst,elst
      use potent  , only : use_polarshortreal,use_lambdadyn
      use polar    ,only: pdamp,thole
      use polpot   ,only: dpcorrect_ik,dpcorrect_scale,n_dpscale
      use shunt    ,only: cut2
      use utilgpu  ,only: dir_queue,def_queue
      use timestat ,only: timer_enter,timer_exit,timer_efld0_direct
      implicit none

      ! shape(ef) = (/3,nrhs,npolebloc/)
      real(t_p),intent(inout):: ef(:,:,:)

      integer i,iglob,iploc,kk
      integer ii,j,k,ksp,kd,iipole,kbis,kpole,kglob
      real(t_p) pdi,pti
      real(t_p) xi,yi,zi,d2
      real(t_p) thole1,pgamma,damp
      real(t_p) alsq2, alsq2n
      real(t_p) one
      real(t_p) pscale,dscale
      real(t_p) d,bn1,bn2,sc3,sc5
      type(real3) fid,fip,fkd,fkp,pos
      type(rpole_elt) ip,kp

      parameter(one =1.0_ti_p)
c
      if (deb_Path) write(*,'(4x,a)') 'efld0_direct_correct_scaling'

      alsq2  = 2.0_ti_p * aewald**2
      alsq2n = 0.0_ti_p
      if (aewald .gt. 0.0_ti_p) alsq2n = one / (sqrtpi*aewald)

!$acc parallel loop gang vector
!$acc&         present(ef)
!$acc&         present(poleglobnl,ipole,loc,pdamp,
!$acc&  thole,x,y,z,rpole,nelst,elst,poleloc,
!$acc&  dpcorrect_ik,dpcorrect_scale,mut,deflambda)
!$acc&         private(fid,fip,fkd,fkp,pos,ip,kp)
!$acc&         async(def_queue)
      do ii = 1, n_dpscale
         iipole = dpcorrect_ik(2*(ii-1)+1)
         kpole  = dpcorrect_ik(2*(ii-1)+2)

         dscale = dpcorrect_scale(2*(ii-1)+1)
         pscale = dpcorrect_scale(2*(ii-1)+2)

         iploc  = poleloc   (iipole)
         iglob  = ipole     (iipole)

         kbis   = poleloc(kpole)
         kglob  = ipole  (kpole)

         if (iploc.lt.1.or.iploc.gt.npolebloc.or.
     &        kbis.lt.1.or. kbis.gt.npolebloc) cycle

         pdi    = pdamp(iipole)
         pti    = thole(iipole)

         pos%x  = x(kglob) - x(iglob)
         pos%y  = y(kglob) - y(iglob)
         pos%z  = z(kglob) - z(iglob)
         call image_inl(pos%x,pos%y,pos%z)
         d2     = pos%x**2 + pos%y**2 + pos%z**2
         if (d2.gt.cut2) cycle

         ip%c   = rpole(01,iipole)
         ip%dx  = rpole(02,iipole)
         ip%dy  = rpole(03,iipole)
         ip%dz  = rpole(04,iipole)
         ip%qxx = rpole(05,iipole)
         ip%qxy = rpole(06,iipole)
         ip%qxz = rpole(07,iipole)
         ip%qyy = rpole(09,iipole)
         ip%qyz = rpole(10,iipole)
         ip%qzz = rpole(13,iipole)

         kp%c   = rpole(01, kpole)
         kp%dx  = rpole(02, kpole)
         kp%dy  = rpole(03, kpole)
         kp%dz  = rpole(04, kpole)
         kp%qxx = rpole(05, kpole)
         kp%qxy = rpole(06, kpole)
         kp%qxz = rpole(07, kpole)
         kp%qyy = rpole(09, kpole)
         kp%qyz = rpole(10, kpole)
         kp%qzz = rpole(13, kpole)

         thole1 = thole(kpole)
         damp   = pdi * pdamp(kpole)
         pgamma = min( pti,thole1 )

         call efld0_couple(d2,pos,ip,kp,alsq2,alsq2n,
     &              aewald,damp,pgamma,dscale,pscale,
     &              fid,fip,fkd,fkp,d,bn1,bn2,sc3,sc5,.true.)
c
c     get vector of derivative of direct field wrt lambda for lambda-dynamics         
c
         if ((use_lambdadyn).and.(elambda.gt.0)) then
              if (mut(kglob)) then
                if (dscale.ne.0.0_ti_p) then
!$acc atomic
                  deflambda(1,1,iploc) = deflambda(1,1,iploc)
     $               + fid%x/elambda
!$acc atomic
                  deflambda(2,1,iploc) = deflambda(2,1,iploc)
     $               + fid%y/elambda
!$acc atomic
                  deflambda(3,1,iploc) = deflambda(3,1,iploc)
     $               + fid%z/elambda
                end if
                if (pscale.ne.0.0_ti_p) then
!$acc atomic
                  deflambda(1,2,iploc) = deflambda(1,2,iploc)
     $               + fip%x/elambda
!$acc atomic
                  deflambda(2,2,iploc) = deflambda(2,2,iploc)
     $               + fip%y/elambda
!$acc atomic
                  deflambda(3,2,iploc) = deflambda(3,2,iploc)
     $               + fip%z/elambda
                end if
              end if
              if (mut(iglob)) then
                if (dscale.ne.0.0_ti_p) then
!$acc atomic
                  deflambda(1,1,kbis) = deflambda(1,1,kbis)
     $               + fkd%x/elambda
!$acc atomic
                  deflambda(2,1,kbis) = deflambda(2,1,kbis)
     $               + fkd%y/elambda
!$acc atomic
                  deflambda(3,1,kbis) = deflambda(3,1,kbis)
     $               + fkd%z/elambda
                end if
                if (pscale.ne.0.0_ti_p) then
!$acc atomic
                  deflambda(1,2,kbis) = deflambda(1,2,kbis)
     $               + fkp%x/elambda
!$acc atomic
                  deflambda(2,2,kbis) = deflambda(2,2,kbis)
     $               + fkp%y/elambda
!$acc atomic
                  deflambda(3,2,kbis) = deflambda(3,2,kbis)
     $               + fkp%z/elambda
                end if
              end if
         end if

         if (dscale.ne.0.0_ti_p) then
!$acc atomic
            ef(1,1,iploc) = ef(1,1,iploc) + fid%x
!$acc atomic
            ef(2,1,iploc) = ef(2,1,iploc) + fid%y
!$acc atomic
            ef(3,1,iploc) = ef(3,1,iploc) + fid%z
!$acc atomic
            ef(1,1,kbis)  = ef(1,1,kbis ) + fkd%x
!$acc atomic
            ef(2,1,kbis)  = ef(2,1,kbis ) + fkd%y
!$acc atomic
            ef(3,1,kbis)  = ef(3,1,kbis ) + fkd%z
         end if
         if (pscale.ne.0.0_ti_p) then
!$acc atomic
            ef(1,2,iploc) = ef(1,2,iploc) + fip%x
!$acc atomic
            ef(2,2,iploc) = ef(2,2,iploc) + fip%y
!$acc atomic
            ef(3,2,iploc) = ef(3,2,iploc) + fip%z
!$acc atomic
            ef(1,2,kbis)  = ef(1,2,kbis ) + fkp%x
!$acc atomic
            ef(2,2,kbis)  = ef(2,2,kbis ) + fkp%y
!$acc atomic
            ef(3,2,kbis)  = ef(3,2,kbis ) + fkp%z
         end if
         end do
      end subroutine
