!
!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!
!     "bspline_fill_site" finds B-spline coefficients and derivatives
!     for PME i-th atomic sites along the fractional coordinate axes
!
!
#include "tinker_macro.h"
!
!
!     #################################################################
!     ##                                                             ##
!     ##  subroutine bsplgen  --  B-spline coefficients for an atom  ##
!     ##                                                             ##
!     #################################################################
!
!
!     "bsplgen" gets B-spline coefficients and derivatives for
!     a single PME atomic site along a particular direction
!
!
subroutine bsplgen (w,isite,thetai)
!$acc routine
   use tinheader
   use pme,only:bsorder,maxorder
   use potent,only:use_mpole,use_polar
   implicit none
   integer i,j,k
   integer isite
   integer level
   real(t_p) w,denom
   real(t_p) thetai(4,bsorder)
   real(t_p) temp(maxorder,maxorder)
!
!     set B-spline depth for partial charges or multipoles
!
   level = 2
   if (use_mpole .or. use_polar)  level = 4
!
!     initialization to get to 2nd order recursion
!
   temp(2,2) = w
   temp(2,1) = 1.0_ti_p - w
!
!     perform one pass to get to 3rd order recursion
!
   temp(3,3) = 0.5_ti_p * w * temp(2,2)
   temp(3,2) = 0.5_ti_p * ((1.0_ti_p+w)*temp(2,1)+&
      &(2.0_ti_p-w)*temp(2,2))
   temp(3,1) = 0.5_ti_p * (1.0_ti_p-w) * temp(2,1)
!
!     compute standard B-spline recursion to desired order
!
!$acc loop seq
   do i = 4, bsorder
      k = i - 1
      denom = 1.0_ti_p / real(k,t_p)
      temp(i,i) = denom * w * temp(k,k)
      do j = 1, i-2
         temp(i,i-j) = denom * ((w+real(j,t_p))*temp(k,i-j-1)&
            &+(real(i-j,t_p)-w)*temp(k,i-j))
      end do
      temp(i,1) = denom * (1.0_ti_p-w) * temp(k,1)
   end do
!
!     get coefficients for the B-spline first derivative
!
   k = bsorder - 1
   temp(k,bsorder) = temp(k,bsorder-1)
!$acc loop seq
   do i = bsorder-1, 2, -1
      temp(k,i) = temp(k,i-1) - temp(k,i)
   end do
   temp(k,1) = -temp(k,1)
!
!     get coefficients for the B-spline second derivative
!
   if (level .eq. 4) then
      k = bsorder - 2
      temp(k,bsorder-1) = temp(k,bsorder-2)
      do i = bsorder-2, 2, -1
         temp(k,i) = temp(k,i-1) - temp(k,i)
      end do
      temp(k,1) = -temp(k,1)
      temp(k,bsorder) = temp(k,bsorder-1)
      do i = bsorder-1, 2, -1
         temp(k,i) = temp(k,i-1) - temp(k,i)
      end do
      temp(k,1) = -temp(k,1)
!
!     get coefficients for the B-spline third derivative
!
      k = bsorder - 3
      temp(k,bsorder-2) = temp(k,bsorder-3)
      do i = bsorder-3, 2, -1
         temp(k,i) = temp(k,i-1) - temp(k,i)
      end do
      temp(k,1) = -temp(k,1)
      temp(k,bsorder-1) = temp(k,bsorder-2)
      do i = bsorder-2, 2, -1
         temp(k,i) = temp(k,i-1) - temp(k,i)
      end do
      temp(k,1) = -temp(k,1)
      temp(k,bsorder) = temp(k,bsorder-1)
      do i = bsorder-1, 2, -1
         temp(k,i) = temp(k,i-1) - temp(k,i)
      end do
      temp(k,1) = -temp(k,1)
   end if
!
!     copy coefficients from temporary to permanent storage
!
!$acc loop seq
   do i = 1, bsorder
!$acc loop seq
      do j = 1, level
         thetai(j,i) = temp(bsorder-j+1,i)
      end do
   end do
   return
end

 ! bslgen special version
 ! To be used only with point charge forcefield (use_charge)
subroutine bsplgen_chg (w,isite,thetai)
!$acc routine
   use tinheader,only: ti_p
   use pme      ,only: bsorder,maxorder
   use potent   ,only: use_charge
   implicit none
   integer i,j,k
   integer isite
   integer,parameter:: level=2
   real(t_p) w,denom
   real(t_p) thetai(level,bsorder,*)
   real(t_p) temp(maxorder,maxorder)

#ifdef TINKER_DEBUG
!$acc routine(fatal_acc)
   if (.not.use_charge) then
      print*,"FATAL ERROR !! bsplgen routine is specific to point",&
         &"charge"
      call fatal_acc
   end if
#endif
!
!     initialization to get to 2nd order recursion
!
   temp(2,2) = w
   temp(2,1) = 1.0_ti_p - w
!
!     perform one pass to get to 3rd order recursion
!
   temp(3,3) = 0.5_ti_p * w * temp(2,2)
   temp(3,2) = 0.5_ti_p * ((1.0_ti_p+w)*temp(2,1)+&
      &(2.0_ti_p-w)*temp(2,2))
   temp(3,1) = 0.5_ti_p * (1.0_ti_p-w) * temp(2,1)
!
!     compute standard B-spline recursion to desired order
!
!$acc loop seq
   do i = 4, bsorder
      k = i - 1
      denom = 1.0_ti_p / real(k,t_p)
      temp(i,i) = denom * w * temp(k,k)
      do j = 1, i-2
         temp(i,i-j) = denom * ((w+real(j,t_p))*temp(k,i-j-1)&
            &+(real(i-j,t_p)-w)*temp(k,i-j))
      end do
      temp(i,1) = denom * (1.0_ti_p-w) * temp(k,1)
   end do
!
!     get coefficients for the B-spline first derivative
!
   k = bsorder - 1
   temp(k,bsorder) = temp(k,bsorder-1)
!$acc loop seq
   do i = bsorder-1, 2, -1
      temp(k,i) = temp(k,i-1) - temp(k,i)
   end do
   temp(k,1) = -temp(k,1)
!
!     copy coefficients from temporary to permanent storage
!
!$acc loop seq
   do i = 1, bsorder
!$acc loop seq
      do j = 1, level
         thetai(j,i,isite) = temp(bsorder-j+1,i)
      end do
   end do
end


subroutine bspline_fill_sitegpu(config)
   use atmlst
   use atoms
   use boxes
   use charge   ,only:iion,nion,nionrecloc
   use domdec   ,only:nproc,nrec
   use disp     ,only:idisp,ndisprecloc
   use inform   ,only:deb_Path,minmaxone
   use interfaces,only:grid_pchg_site_p,grid_pchg_sitecu&
      &,grid_mpole_site_p,grid_mpole_sitecu&
      &,grid_disp_site_p,grid_disp_sitecu
   use mpole
   use neigh    ,only:celle_pole,celle_chg,disp_nbl
   use pme
   use tinheader,only:ti_p,prec1_eps
   use utilgpu  ,only:rec_queue,nSMP
   use potent
   implicit none
   integer,intent(in),optional::config
   integer c_mpole,c_charge,c_disp,c_n
   integer i,ifr,k,cfg,iipole,iichg
   integer,pointer,save:: glob_p(:),type_p(:)
   real(t_p) xi,yi,zi
   real(t_p) w,fr
   parameter(c_mpole=0,c_charge=1,c_disp=2)
!$acc routine(bsplgen) seq
!$acc routine(bsplgen_chg) seq

   if (present(config)) then
      cfg = config
   else
      cfg = c_mpole
   end if

   if (cfg.eq.c_mpole) then
      glob_p => ipole
      type_p => polerecglob
      c_n    =  npolerecloc
#ifdef _OPENACC
      if (associated(grid_mpole_site_p,grid_mpole_sitecu)) then
         type_p => celle_pole
      end if
#endif
   else if (cfg.eq.c_charge) then
      glob_p => iion
      type_p => chgrecglob
      c_n    =  nionrecloc
#ifdef _OPENACC
      if (associated(grid_pchg_site_p,grid_pchg_sitecu)) then
         type_p => celle_chg
      end if
#endif
   else if (cfg.eq.c_disp) then
      glob_p => idisp
      type_p => disprecglob
      c_n    =  ndisprecloc
#ifdef _OPENACC
      if (associated(grid_disp_site_p,grid_disp_sitecu)) then
         type_p => disp_nbl%s_kind
      end if
#endif
   end if

   if (cfg.eq.c_mpole) then

!$acc parallel loop num_gangs(4*nSMP) &
!$acc         present(x,y,z,recip,igrid,thetai1,thetai2, &
!$acc  thetai3,type_p,glob_p) &
!$acc         async(rec_queue)
      do k=1,c_n
         iipole = type_p(k)
         i      = glob_p(iipole)
!
!     get the b-spline coefficients for the i-th atomic site
!
         xi   = x(i)
         yi   = y(i)
         zi   = z(i)
         w    = xi*recip(1,1) + yi*recip(2,1) + zi*recip(3,1)
         fr   = real(nfft1,t_p) * (w-anint(w)+0.5_ti_p)
         ifr  = int(fr-pme_eps)
         w    = fr - real(ifr,t_p)
         igrid(1,i) = ifr - bsorder
         call bsplgen (w,k,thetai1(1,1,k))
         w    = xi*recip(1,2) + yi*recip(2,2) + zi*recip(3,2)
         fr   = real(nfft2,t_p) * (w-anint(w)+0.5_ti_p)
         ifr  = int(fr-pme_eps)
         w    = fr - real(ifr,t_p)
         igrid(2,i) = ifr - bsorder
         call bsplgen (w,k,thetai2(1,1,k))
         w    = xi*recip(1,3) + yi*recip(2,3) + zi*recip(3,3)
         fr   = real(nfft3,t_p) * (w-anint(w)+0.5_ti_p)
         ifr  = int(fr-pme_eps)
         w    = fr - real(ifr,t_p)
         igrid(3,i) = ifr - bsorder
         call bsplgen (w,k,thetai3(1,1,k))
      end do

   else if (cfg.eq.c_charge.or.cfg.eq.c_disp) then
      ! FIXME GnU thetai1_p(1,1,1) instead of thetai1_p

!$acc parallel loop num_gangs(4*nSMP) &
!$acc         present(x,y,z,recip,igrid,thetai1_p &
!$acc                ,thetai2_p,thetai3_p,type_p,glob_p) &
!$acc         async(rec_queue)
      do k=1,c_n
         iichg  = type_p(k)
         i      = glob_p(iichg)
!
!     get the b-spline coefficients for the i-th atomic site
!
         xi   = x(i)
         yi   = y(i)
         zi   = z(i)
         w    = xi*recip(1,1) + yi*recip(2,1) + zi*recip(3,1)
         fr   = real(nfft1,t_p) * (w-anint(w)+0.5_ti_p)
         ifr  = int(fr-pme_eps)
         w    = fr - real(ifr,t_p)
         igrid(1,i) = ifr - bsorder
         call bsplgen_chg (w,k,thetai1_p)
         w    = xi*recip(1,2) + yi*recip(2,2) + zi*recip(3,2)
         fr   = real(nfft2,t_p) * (w-anint(w)+0.5_ti_p)
         ifr  = int(fr-pme_eps)
         w    = fr - real(ifr,t_p)
         igrid(2,i) = ifr - bsorder
         call bsplgen_chg (w,k,thetai2_p)
         w    = xi*recip(1,3) + yi*recip(2,3) + zi*recip(3,3)
         fr   = real(nfft3,t_p) * (w-anint(w)+0.5_ti_p)
         ifr  = int(fr-pme_eps)
         w    = fr - real(ifr,t_p)
         igrid(3,i) = ifr - bsorder
         call bsplgen_chg (w,k,thetai3_p)
      end do

   end if
end
!
!       "grid_mpole_site" places the i-th fractional atomic multipole onto
!       the particle mesh Ewald grid
!
!
subroutine grid_mpole_sitegpu(fmpvec)
   use atmlst
   use chunks
   use domdec
   use fft
   use inform ,only:deb_Path
   use mpole
   use pme
   use potent
   use sizes
   use utils
   use utilgpu,only:rec_queue,ngangs_rec
   implicit none
   integer istart,iend,jstart,jend,kstart,kend
   integer istat,ied,jstat,jed,kstat,ked
   integer iipole
   integer i,j,k,m,impi,rankloc
   integer mk,mj
   integer iproc,proc
   integer ii,jj,kk
   integer isite,iatm
   integer offsetx,offsety
   integer offsetz
   integer location,twonlpts_1,twonlpts_12
   integer nlptsit
   real(t_p) v0,u0,t0
   real(t_p) v1,u1,t1
   real(t_p) v2,u2,t2
   real(t_p) vut(6)
   real(t_p) term,term0,term1,term2
   real(t_p) fmp(10)
   real(t_p) fmpvec(10,max(npolerecloc,1))
!
   if (deb_Path) write(*,'(4X,A)') "grid_mpole_sitegpu"
   if (use_pmecore) then
      rankloc  = rank_bis
   else
      rankloc  = rank
   end if
   kstat = kstart1(rankloc+1)
   ked   = kend1  (rankloc+1)
   jstat = jstart1(rankloc+1)
   jed   = jend1  (rankloc+1)
   istat = istart1(rankloc+1)
   ied   = iend1  (rankloc+1)
   twonlpts_1  = 2*nlpts+1
   twonlpts_12 = twonlpts_1**2
   nlptsit     = (2*nlpts+1)**3 - 1

!
!     put the permanent multipole moments onto the grid
!
!$acc parallel loop gang worker vector_length(32) &
!$acc         present(polerecglob,igrid,ipole,thetai1,thetai2) &
!$acc         present(thetai3,fmpvec,qgridin_2d) &
!$acc         async(rec_queue)
   do impi = 1,npolerecloc
      isite   = ipole(polerecglob(impi))
      offsetx = 1 - (igrid(1,isite) + grdoff - nlpts)
      offsety = 1 - (igrid(2,isite) + grdoff - nlpts)
      offsetz = 1 - (igrid(3,isite) + grdoff - nlpts)
!
!       put the induced dipole moments onto the grid
!
!$acc loop vector private(vut)
      do kk = 0, nlptsit
         k  = igrid(3,isite) + grdoff + kk/twonlpts_12 - nlpts
         j  = igrid(2,isite) + grdoff&
            &+ mod(kk/twonlpts_1,twonlpts_1) - nlpts
         i  = igrid(1,isite) + grdoff&
            &+ mod(kk,twonlpts_1)-nlpts
         mk = k + offsetz
         if (k .lt. 1) k = k + nfft3
         mj = j + offsety
         if (j .lt. 1) j = j + nfft2
         m  = i + offsetx
         if (i .lt. 1) i = i + nfft1
!
         vut(1)= thetai3(1,mk,impi)
         vut(2)= thetai3(2,mk,impi)
         vut(3)= thetai3(3,mk,impi)
         vut(4)= thetai2(1,mj,impi)
         vut(5)= thetai2(2,mj,impi)
         vut(6)= thetai2(3,mj,impi)
         term0 = fmpvec(1,impi)*vut(4)*vut(1)&
            &+ fmpvec(3,impi)*vut(5)*vut(1)&
            &+ fmpvec(4,impi)*vut(4)*vut(2)&
            &+ fmpvec(6,impi)*vut(6)*vut(1)&
            &+ fmpvec(7,impi)*vut(4)*vut(3)&
            &+ fmpvec(10,impi)*vut(5)*vut(2)
         term1 = fmpvec(2,impi)*vut(4)*vut(1)&
            &+ fmpvec(8,impi)*vut(5)*vut(1)&
            &+ fmpvec(9,impi)*vut(4)*vut(2)
         term2 = fmpvec(5,impi)*vut(4)*vut(1)
         vut(1)= thetai1(1,m,impi)
         vut(2)= thetai1(2,m,impi)
         vut(3)= thetai1(3,m,impi)
         term  = term0*vut(1) + term1*vut(2) + term2*vut(3)

         if (((k.ge.kstat).and.(k.le.ked)).and.&
            &((j.ge.jstat).and.(j.le.jed)).and.&
            &((i.ge.istat).and.(i.le.ied))) then
!            location = 1+(i-istat)*2+(j-jstat)*2*n1mpimax
!    &                   +(k-kstat)*2*n1mpimax*n2mpimax
!            call atomic_adds(location,term,qgridin_p(1))
!$acc atomic update
            qgridin_2d(1,i-istat+1,j-jstat+1,k-kstat+1,1) =&
               &qgridin_2d(1,i-istat+1,j-jstat+1,k-kstat+1,1)&
               &+ term
            cycle
         end if

!$acc loop seq
         do iproc = 1, nrec_send
            proc   = prec_send(iproc)
            kstart = kstart1(proc+1)
            kend   = kend1  (proc+1)
            jstart = jstart1(proc+1)
            jend   = jend1  (proc+1)
            istart = istart1(proc+1)
            iend   = iend1  (proc+1)
            if (((k.ge.kstart).and.(k.le.kend)).and.&
               &((j.ge.jstart).and.(j.le.jend)).and.&
               &((i.ge.istart).and.(i.le.iend))) then
!             location = 1+(i-istart)*2+(j-jstart)*2*n1mpimax
!    &                    +(k-kstart)*2*n1mpimax*n2mpimax
!    &                    +iproc*2*n1mpimax*n2mpimax*n3mpimax
!             call atomic_adds(location,term,qgridin_p(1))
!$acc atomic update
               qgridin_2d(1,i-istart+1,j-jstart+1,k-kstart+1,iproc+1)&
                  &= qgridin_2d(1,i-istart+1,j-jstart+1,k-kstart+1,iproc+1)&
                  &+ term
               exit
            end if
         end do
      end do
   end do
end
!
!     "grid_pchg_sitegpu" places the i-th fractional atomic charge onto
!     the particle mesh Ewald grid
!
subroutine grid_pchg_sitegpu
   use atmlst ,only: chgrecglob
   use charge ,only: pchg,iion,nionrecloc
   implicit none
   call grid_put_site( chgrecglob,iion,pchg,nionrecloc,0)
end subroutine
!
!
!     "grid_disp_sitegpu" places the i-th damped dispersion coefficients onto
!     the particle mesh Ewald grid
!
!
subroutine grid_disp_sitegpu
   use atmlst
   use disp
   implicit none
   call grid_put_site( disprecglob,idisp,csix,ndisprecloc,1 )
end subroutine

subroutine grid_put_site( kind_id,atom_id,attr,nk,prm)
   use atoms  ,only: n
   use chunks
   use domdec
   use fft
   use inform ,only: deb_Path
   use mpole
   use pme
   use potent ,only: use_pmecore
   use sizes
   use utils
   use utilgpu,only: rec_queue
   implicit none
   integer  ,intent(in):: nk, prm, kind_id(nk), atom_id(n)
   real(t_p),intent(in):: attr(n)
   integer istart,iend,jstart,jend,kstart,kend
   integer istat,ied,jstat,jed,kstat,ked
   integer iichg
   integer i,j,k,m,impi,rankloc
   integer mk,mj
   integer iproc,proc,tsize
   integer ii,jj,kk
   integer isite,iatm
   integer offsetx,offsety
   integer offsetz
   integer location,twonlpts_1,twonlpts_12
   integer nlptsit
   real(t_p) q
   real(t_p) v0,u0,t0
   real(t_p) term,term0
   logical,save:: f_in=.true.
!
   if (deb_Path) write(*,'(3x,a,I3)') "grid_put_site",prm
   if (use_pmecore) then
      rankloc  = rank_bis
   else
      rankloc  = rank
   end if

   kstat = kstart1(rankloc+1)
   ked   = kend1  (rankloc+1)
   jstat = jstart1(rankloc+1)
   jed   = jend1  (rankloc+1)
   istat = istart1(rankloc+1)
   ied   = iend1  (rankloc+1)
   twonlpts_1  = 2*nlpts+1
   twonlpts_12 = twonlpts_1**2
   nlptsit     = (2*nlpts+1)**3 - 1
!
!     put the permanent multipole moments onto the grid
!
!$acc parallel loop gang worker vector_length(32) &
!$acc         present(kind_id,atom_id,attr,igrid,thetai1_p,thetai2_p &
!$acc   ,thetai3_p,qgridin_2d) async(rec_queue)
   do impi = 1,nk
      iichg   = kind_id(impi)
      isite   = atom_id(iichg)
      q       = attr   (iichg)
      offsetx = 1 - (igrid(1,isite) + grdoff - nlpts)
      offsety = 1 - (igrid(2,isite) + grdoff - nlpts)
      offsetz = 1 - (igrid(3,isite) + grdoff - nlpts)
!
!       put the induced dipole moments onto the grid
!
!$acc loop vector
      do kk = 0, nlptsit
         k  = igrid(3,isite) + grdoff + kk/twonlpts_12 - nlpts
         j  = igrid(2,isite) + grdoff&
            &+ mod(kk/twonlpts_1,twonlpts_1) - nlpts
         i  = igrid(1,isite) + grdoff&
            &+ mod(kk,twonlpts_1)-nlpts
         mk = k + offsetz
         if (k .lt. 1) k = k + nfft3
         mj = j + offsety
         if (j .lt. 1) j = j + nfft2
         m  = i + offsetx
         if (i .lt. 1) i = i + nfft1
!
         v0    = thetai3_p(1,mk,impi)
         u0    = thetai2_p(1,mj,impi)
         t0    = thetai1_p(1,m ,impi)
         term  = q*u0*v0*t0

         if (((k.ge.kstat).and.(k.le.ked)).and.&
            &((j.ge.jstat).and.(j.le.jed)).and.&
            &((i.ge.istat).and.(i.le.ied))) then
!            location = 1+(i-istat)*2+(j-jstat)*2*n1mpimax
!    &                   +(k-kstat)*2*n1mpimax*n2mpimax
!            call atomic_adds(location,term,qgridin_p(1))
!$acc atomic update
            qgridin_2d(1,i-istat+1,j-jstat+1,k-kstat+1,1) =&
               &qgridin_2d(1,i-istat+1,j-jstat+1,k-kstat+1,1)&
               &+ term
            cycle
         end if

!$acc loop seq
         do iproc = 1, nrec_send
            proc   = prec_send(iproc)
            kstart = kstart1(proc+1)
            kend   = kend1  (proc+1)
            jstart = jstart1(proc+1)
            jend   = jend1  (proc+1)
            istart = istart1(proc+1)
            iend   = iend1  (proc+1)
            if (((k.ge.kstart).and.(k.le.kend)).and.&
               &((j.ge.jstart).and.(j.le.jend)).and.&
               &((i.ge.istart).and.(i.le.iend))) then
!             location = 1+(i-istart)*2+(j-jstart)*2*n1mpimax
!    &                    +(k-kstart)*2*n1mpimax*n2mpimax
!    &                    +iproc*2*n1mpimax*n2mpimax*n3mpimax
!             call atomic_adds(location,term,qgridin_p(1))
!$acc atomic update
               qgridin_2d(1,i-istart+1,j-jstart+1,k-kstart+1,iproc+1)&
                  &= qgridin_2d(1,i-istart+1,j-jstart+1,k-kstart+1,iproc+1)&
                  &+ term
               exit
            end if
         end do
      end do
   end do
end
!
!
!       "grid_uind_site" places the fractional induced dipoles onto the
!       particle mesh Ewald grid
!
!
subroutine grid_uind_sitegpu(fuindvec,fuinpvec,qgrid2loc)
   use atmlst,only:polerecglob
   use chunks
   use domdec
   use fft
   use inform,only:deb_Path
   use mpole
   use pme
   use potent
   use utils
   use utilgpu,only:rec_queue,ngangs_rec
   implicit none
   integer istart,iend,jstart,jend,kstart,kend
   integer istat,ied,jstat,jed,kstat,ked
   integer i,j,k,impi
   integer mk,mj,mi
   integer iproc,proc,rankloc
   integer ii,jj,kk
   integer iipole,iglob,isite,iatm
   integer offsetx,offsety
   integer offsetz
   integer twonlpts_1,twonlpts_12,nlptsit
   integer location
   real(t_p) v0,u0,t0
   real(t_p) v1,u1,t1
   real(t_p) vut(6)
   real(t_p) term
   real(t_p) term01,term11
   real(t_p) term02,term12
   integer igrid1,igrid2,igrid3
   integer nearpt(3),abound(6)
   real(t_p)  fuinp(3),fuind(3)
   real(t_p),intent(in),dimension(3,npolerecloc)::fuindvec,fuinpvec
   real(t_p) qgrid2loc(:,:,:,:,:)
   logical,save::first_in=.true.

   if (deb_Path)  print '(5x,a)','grid_uind_sitegpu'
   rankloc = merge(rank_bis,rank,use_pmecore)
   kstat = kstart1(rankloc+1)
   ked   = kend1  (rankloc+1)
   jstat = jstart1(rankloc+1)
   jed   = jend1  (rankloc+1)
   istat = istart1(rankloc+1)
   ied   = iend1  (rankloc+1)
   twonlpts_1  = 2*nlpts+1
   twonlpts_12 = twonlpts_1**2
   nlptsit     = (2*nlpts+1)**3-1

!
!       put the induced dipole moments onto the grid
!
!$acc parallel loop gang worker vector_length(32) &
!$acc         present(polerecglob,igrid,ipole,thetai1,thetai2) &
!$acc         present(fuindvec,fuinpvec,qgrid2loc) &
!$acc         async(rec_queue)
   do impi = 1,npolerecloc
      isite     = ipole(polerecglob(impi))
      offsetx   = 1 - (igrid(1,isite) + grdoff - nlpts)
      offsety   = 1 - (igrid(2,isite) + grdoff - nlpts)
      offsetz   = 1 - (igrid(3,isite) + grdoff - nlpts)
!
!     Three dimensional loop on the grid collapse by hand
!
!$acc loop vector private(vut)
      do kk = 0, nlptsit
         k  = igrid(3,isite) + grdoff + kk/twonlpts_12 - nlpts
         j  = igrid(2,isite) + grdoff&
            &+ mod(kk/twonlpts_1,twonlpts_1) - nlpts
         i  = igrid(1,isite) + grdoff&
            &+ mod(kk,twonlpts_1)-nlpts
         mk     = k + offsetz
         if (k .lt. 1) k = k + nfft3
         mj     = j + offsety
         if (j .lt. 1) j = j + nfft2
         mi     = i + offsetx
         if (i .lt. 1) i = i + nfft1

         vut(1) = thetai3(1,mk,impi)
         vut(2) = thetai3(2,mk,impi)
         vut(3) = thetai2(1,mj,impi)
         vut(4) = thetai2(2,mj,impi)
         vut(5) = thetai1(1,mi,impi)
         vut(6) = thetai1(2,mi,impi)
         term01 = fuindvec(2,impi)*vut(4)*vut(1)&
            &+ fuindvec(3,impi)*vut(3)*vut(2)
         term11 = fuindvec(1,impi)*vut(3)*vut(1)
         term02 = fuinpvec(2,impi)*vut(4)*vut(1)&
            &+ fuinpvec(3,impi)*vut(3)*vut(2)
         term12 = fuinpvec(1,impi)*vut(3)*vut(1)
         term01 = term01*vut(5) + term11*vut(6)
         term02 = term02*vut(5) + term12*vut(6)
!
         if (((k.ge.kstat).and.(k.le.ked)).and.&
            &((j.ge.jstat).and.(j.le.jed)).and.&
            &((i.ge.istat).and.(i.le.ied)))    then
!              location = 1+(i-istat)*2+(j-jstat)*2*n1mpimax
!    &                     +(k-kstat)*2*n1mpimax*n2mpimax
!              call atomic_Adds2(location,term01,term02,qgrid2in_p(1))
!$acc atomic update
            qgrid2loc(1,i-istat+1,j-jstat+1,k-kstat+1,1)&
               &= qgrid2loc(1,i-istat+1,j-jstat+1,k-kstat+1,1)&
               &+ term01
!$acc atomic update
            qgrid2loc(2,i-istat+1,j-jstat+1,k-kstat+1,1)&
               &= qgrid2loc(2,i-istat+1,j-jstat+1,k-kstat+1,1)&
               &+ term02
            cycle
         end if

!$acc loop seq
         do iproc = 1, nrec_send
            proc   = prec_send(iproc)
            kstart = kstart1(proc+1)
            kend   = kend1  (proc+1)
            jstart = jstart1(proc+1)
            jend   = jend1  (proc+1)
            istart = istart1(proc+1)
            iend   = iend1  (proc+1)
            if (((k.ge.kstart).and.(k.le.kend)).and.&
               &((j.ge.jstart).and.(j.le.jend)).and.&
               &((i.ge.istart).and.(i.le.iend))) then
!              location = 1+(i-istart)*2+(j-jstart)*2*n1mpimax
!    &                     +(k-kstart)*2*n1mpimax*n2mpimax
!              call atomic_Adds2(location,term01,term02,qgrid2in_p(1))
!$acc atomic update
               qgrid2loc(1,i-istart+1,j-jstart+1,k-kstart+1,iproc+1)&
                  &= qgrid2loc(1,i-istart+1,j-jstart+1,k-kstart+1,iproc+1)&
                  &+ term01
!$acc atomic update
               qgrid2loc(2,i-istart+1,j-jstart+1,k-kstart+1,iproc+1)&
                  &= qgrid2loc(2,i-istart+1,j-jstart+1,k-kstart+1,iproc+1)&
                  &+ term02
               exit
            end if
         end do
      end do
   end do
   first_in=.false.
end
!
!
!       "fphi_mpole_site" extracts the permanent multipole potential on the i-th site from
!       the particle mesh Ewald grid
!
!
subroutine fphi_mpole_sitegpu
   use atmlst
   use domdec
   use fft
   use inform ,only:deb_Path
   use mpole
   use pme
   use potent
   use utilgpu
   implicit none
   integer istart,iend,jstart,jend,kstart,kend
   integer istarti,iendi,jstarti,jendi,kstarti,kendi
   integer i,j,k,impi,rankloc
   integer iproc,proc
   integer isite,iatm,iipole
   integer i0,j0,k0
   integer it1,it2,it3
   integer igrd0,jgrd0,kgrd0
   real(t_p) v0,v1,v2,v3
   real(t_p) u0,u1,u2,u3
   real(t_p) t0,t1,t2,t3,tq
   real(t_p) tu00,tu10,tu01,tu20,tu11
   real(t_p) tu02,tu21,tu12,tu30,tu03
   real(t_p) tuv000,tuv100,tuv010,tuv001
   real(t_p) tuv200,tuv020,tuv002,tuv110
   real(t_p) tuv101,tuv011,tuv300,tuv030
   real(t_p) tuv003,tuv210,tuv201,tuv120
   real(t_p) tuv021,tuv102,tuv012,tuv111
!
   if (deb_Path) write(*,'(5x,a)') "fphi_mpole_sitegpu"
   if (use_pmecore) then
      rankloc  = rank_bis
   else
      rankloc  = rank
   end if
   kstarti = kstart1(rankloc+1)
   kendi   = kend1  (rankloc+1)
   jstarti = jstart1(rankloc+1)
   jendi   = jend1  (rankloc+1)
   istarti = istart1(rankloc+1)
   iendi   = iend1  (rankloc+1)
!
!$acc parallel loop num_gangs(ngangs_rec) &
!$acc         present(polerecglob,ipole,qgridin_2d,igrid,thetai1, &
!$acc  thetai2,thetai3,kstart1,kend1,jstart1,jend1,kstart1,kend1, &
!$acc  fphirec) &
!$acc         async(rec_queue)
   do impi = 1,npolerecloc
      iipole = polerecglob(impi)
      isite  = ipole(iipole)
      igrd0  = igrid(1,isite)
      jgrd0  = igrid(2,isite)
      kgrd0  = igrid(3,isite)
      tuv000 = 0.0_ti_p
      tuv001 = 0.0_ti_p
      tuv010 = 0.0_ti_p
      tuv100 = 0.0_ti_p
      tuv200 = 0.0_ti_p
      tuv020 = 0.0_ti_p
      tuv002 = 0.0_ti_p
      tuv110 = 0.0_ti_p
      tuv101 = 0.0_ti_p
      tuv011 = 0.0_ti_p
      tuv300 = 0.0_ti_p
      tuv030 = 0.0_ti_p
      tuv003 = 0.0_ti_p
      tuv210 = 0.0_ti_p
      tuv201 = 0.0_ti_p
      tuv120 = 0.0_ti_p
      tuv021 = 0.0_ti_p
      tuv102 = 0.0_ti_p
      tuv012 = 0.0_ti_p
      tuv111 = 0.0_ti_p
      k0     = kgrd0
!$acc loop seq
      do it3 = 1, bsorder
         k0   = k0 + 1
!           k   = k0 + 1 + (nfft3-isign(nfft3,k0))/2
         k   = k0 + 1 + ishft(nfft3-isign(nfft3,k0),-1)
         v0   = thetai3(1,it3,impi)
         v1   = thetai3(2,it3,impi)
         v2   = thetai3(3,it3,impi)
         v3   = thetai3(4,it3,impi)
         tu00 = 0.0_ti_p
         tu10 = 0.0_ti_p
         tu01 = 0.0_ti_p
         tu20 = 0.0_ti_p
         tu11 = 0.0_ti_p
         tu02 = 0.0_ti_p
         tu30 = 0.0_ti_p
         tu21 = 0.0_ti_p
         tu12 = 0.0_ti_p
         tu03 = 0.0_ti_p
         j0   = jgrd0
!$acc loop seq
         do it2 = 1, bsorder
            j0 = j0 + 1
!              j  = j0 + 1 + (nfft2-isign(nfft2,j0))/2
            j  = j0 + 1 + ishft(nfft2-isign(nfft2,j0),-1)
            u0 = thetai2(1,it2,impi)
            u1 = thetai2(2,it2,impi)
            u2 = thetai2(3,it2,impi)
            u3 = thetai2(4,it2,impi)
            t0 = 0.0_ti_p
            t1 = 0.0_ti_p
            t2 = 0.0_ti_p
            t3 = 0.0_ti_p
            i0 = igrd0
!$acc loop seq
            do it1 = 1, bsorder
               i0 = i0 + 1
!                 i = i0 + 1 + (nfft1-isign(nfft1,i0))/2
               i = i0 + 1 + ishft(nfft1-isign(nfft1,i0),-1)
!
               tq     = 0.0_ti_p
               if (((k.ge.kstarti).and.(k.le.kendi)).and.&
                  &((j.ge.jstarti).and.(j.le.jendi)).and.&
                  &((i.ge.istarti).and.(i.le.iendi))) then
                  tq = qgridin_2d(1,i-istarti+1,j-jstarti+1,&
                     &k-kstarti+1,1)
                  goto 10
               end if
!$acc loop seq
               do iproc = 1, nrec_send
                  proc   = prec_send(iproc)
                  kstart = kstart1(proc+1)
                  kend   = kend1  (proc+1)
                  jstart = jstart1(proc+1)
                  jend   = jend1  (proc+1)
                  istart = istart1(proc+1)
                  iend   = iend1  (proc+1)
                  if (((k.ge.kstart).and.(k.le.kend)).and.&
                     &((j.ge.jstart).and.(j.le.jend)).and.&
                     &((i.ge.istart).and.(i.le.iend))) then
                     tq  = qgridin_2d(1,i-istart+1,j-jstart+1,k-kstart+1,&
                        &iproc+1)
                     goto 10
                  end if
               end do
               cycle
!
10             continue
               t0 = t0 + tq*thetai1(1,it1,impi)
               t1 = t1 + tq*thetai1(2,it1,impi)
               t2 = t2 + tq*thetai1(3,it1,impi)
               t3 = t3 + tq*thetai1(4,it1,impi)
            end do
            tu00 = tu00 + t0*u0
            tu10 = tu10 + t1*u0
            tu01 = tu01 + t0*u1
            tu20 = tu20 + t2*u0
            tu11 = tu11 + t1*u1
            tu02 = tu02 + t0*u2
            tu30 = tu30 + t3*u0
            tu21 = tu21 + t2*u1
            tu12 = tu12 + t1*u2
            tu03 = tu03 + t0*u3
         end do
         tuv000 = tuv000 + tu00*v0
         tuv100 = tuv100 + tu10*v0
         tuv010 = tuv010 + tu01*v0
         tuv001 = tuv001 + tu00*v1
         tuv200 = tuv200 + tu20*v0
         tuv020 = tuv020 + tu02*v0
         tuv002 = tuv002 + tu00*v2
         tuv110 = tuv110 + tu11*v0
         tuv101 = tuv101 + tu10*v1
         tuv011 = tuv011 + tu01*v1
         tuv300 = tuv300 + tu30*v0
         tuv030 = tuv030 + tu03*v0
         tuv003 = tuv003 + tu00*v3
         tuv210 = tuv210 + tu21*v0
         tuv201 = tuv201 + tu20*v1
         tuv120 = tuv120 + tu12*v0
         tuv021 = tuv021 + tu02*v1
         tuv102 = tuv102 + tu10*v2
         tuv012 = tuv012 + tu01*v2
         tuv111 = tuv111 + tu11*v1
      end do
      fphirec(1,impi) = tuv000
      fphirec(2,impi) = tuv100
      fphirec(3,impi) = tuv010
      fphirec(4,impi) = tuv001
      fphirec(5,impi) = tuv200
      fphirec(6,impi) = tuv020
      fphirec(7,impi) = tuv002
      fphirec(8,impi) = tuv110
      fphirec(9,impi) = tuv101
      fphirec(10,impi) = tuv011
      fphirec(11,impi) = tuv300
      fphirec(12,impi) = tuv030
      fphirec(13,impi) = tuv003
      fphirec(14,impi) = tuv210
      fphirec(15,impi) = tuv201
      fphirec(16,impi) = tuv120
      fphirec(17,impi) = tuv021
      fphirec(18,impi) = tuv102
      fphirec(19,impi) = tuv012
      fphirec(20,impi) = tuv111
   end do
end
!
!
!       "fphi_chg_site" extracts the permanent charge potential on the i-th site from
!       the particle mesh Ewald grid
!
!
subroutine fphi_chg_sitegpu
   use atmlst
   use domdec
   use fft
   use inform ,only:deb_Path
   use charge
   use pme
   use potent
   use utilgpu
   implicit none
   integer istart,iend,jstart,jend,kstart,kend
   integer istarti,iendi,jstarti,jendi,kstarti,kendi
   integer i,j,k,impi,rankloc
   integer iproc,proc
   integer isite,iatm,iichg
   integer i0,j0,k0
   integer it1,it2,it3
   integer igrd0,jgrd0,kgrd0
   real(t_p) v0,v1
   real(t_p) u0,u1
   real(t_p) t0,t1,tq
   real(t_p) tu00,tu10,tu01
   real(t_p) tuv000,tuv100,tuv010,tuv001
!
   if (deb_Path) write(*,'(5x,a)') "fphi_chg_sitegpu"
   if (use_pmecore) then
      rankloc  = rank_bis
   else
      rankloc  = rank
   end if
   kstarti = kstart1(rankloc+1)
   kendi   = kend1  (rankloc+1)
   jstarti = jstart1(rankloc+1)
   jendi   = jend1  (rankloc+1)
   istarti = istart1(rankloc+1)
   iendi   = iend1  (rankloc+1)
!
!$acc parallel loop gang vector async(rec_queue) &
!$acc         present(chgrecglob,iion,qgridin_2d,igrid,thetai1_p, &
!$acc  thetai2_p,thetai3_p,kstart1,kend1,jstart1,jend1,kstart1,kend1, &
!$acc  fphirec)
   do impi = 1,nionrecloc
      iichg = chgrecglob(impi)
      isite  = iion(iichg)
      igrd0  = igrid(1,isite)
      jgrd0  = igrid(2,isite)
      kgrd0  = igrid(3,isite)
      tuv000 = 0.0_ti_p
      tuv001 = 0.0_ti_p
      tuv010 = 0.0_ti_p
      tuv100 = 0.0_ti_p
      k0     = kgrd0
!$acc loop seq
      do it3 = 1, bsorder
         k0   = k0 + 1
!           k   = k0 + 1 + (nfft3-isign(nfft3,k0))/2
         k   = k0 + 1 + ishft(nfft3-isign(nfft3,k0),-1)
         v0   = thetai3_p(1,it3,impi)
         v1   = thetai3_p(2,it3,impi)
         tu00 = 0.0_ti_p
         tu10 = 0.0_ti_p
         tu01 = 0.0_ti_p
         j0   = jgrd0
!$acc loop seq
         do it2 = 1, bsorder
            j0 = j0 + 1
!              j  = j0 + 1 + (nfft2-isign(nfft2,j0))/2
            j  = j0 + 1 + ishft(nfft2-isign(nfft2,j0),-1)
            u0 = thetai2_p(1,it2,impi)
            u1 = thetai2_p(2,it2,impi)
            t0 = 0.0_ti_p
            t1 = 0.0_ti_p
            i0 = igrd0
!$acc loop seq
            do it1 = 1, bsorder
               i0 = i0 + 1
!                 i = i0 + 1 + (nfft1-isign(nfft1,i0))/2
               i = i0 + 1 + ishft(nfft1-isign(nfft1,i0),-1)
!
               tq     = 0.0_ti_p
               if (((k.ge.kstarti).and.(k.le.kendi)).and.&
                  &((j.ge.jstarti).and.(j.le.jendi)).and.&
                  &((i.ge.istarti).and.(i.le.iendi))) then
                  tq = qgridin_2d(1,i-istarti+1,j-jstarti+1,&
                     &k-kstarti+1,1)
                  goto 10
               end if
!$acc loop seq
               do iproc = 1, nrec_send
                  proc   = prec_send(iproc)
                  kstart = kstart1(proc+1)
                  kend   = kend1  (proc+1)
                  jstart = jstart1(proc+1)
                  jend   = jend1  (proc+1)
                  istart = istart1(proc+1)
                  iend   = iend1  (proc+1)
                  if (((k.ge.kstart).and.(k.le.kend)).and.&
                     &((j.ge.jstart).and.(j.le.jend)).and.&
                     &((i.ge.istart).and.(i.le.iend))) then
                     tq  = qgridin_2d(1,i-istart+1,j-jstart+1,k-kstart+1,&
                        &iproc+1)
                     goto 10
                  end if
               end do
               cycle
!
10             continue
               t0 = t0 + tq*thetai1_p(1,it1,impi)
               t1 = t1 + tq*thetai1_p(2,it1,impi)
            end do
            tu00 = tu00 + t0*u0
            tu10 = tu10 + t1*u0
            tu01 = tu01 + t0*u1
         end do
         tuv000 = tuv000 + tu00*v0
         tuv100 = tuv100 + tu10*v0
         tuv010 = tuv010 + tu01*v0
         tuv001 = tuv001 + tu00*v1
      end do
      fphirec(1,impi) = tuv000
      fphirec(2,impi) = tuv100
      fphirec(3,impi) = tuv010
      fphirec(4,impi) = tuv001
   end do
end

subroutine grid_disp_force
   write(0,*) " TINKER_ERROR !! Code Source Alteration !! "
   write(0,*) " grid_disp_force is not supposed to be called "
   __TINKER_FATAL__
end subroutine
!
!     "grid_pchg_force"
!     get first derivatives of the reciprocal space energy
!
subroutine grid_pchg_force
   use atmlst ,only: chgrecglob
   use boxes  ,only: recip
   use charge
   use chgpot
   use domdec ,only: rank,rank_bis,nproc,nrec,ndir,comm_rec&
      &,COMM_TINKER,nrec_send,locrec,prec_send,nbloc
   use deriv  ,only: decrec
   use fft
   use inform ,only: deb_Path
   use pme    ,only: igrid,nfft1,nfft2,nfft3,bsorder&
      &,thetai1_p,thetai2_p,thetai3_p,qgridin_2d
   use potent ,only: use_pmecore
   use utilgpu,only: ngangs_rec,rec_queue
   use tinheader ,only: ti_p
   implicit none
   integer iichg,iglob,iloc,iatm,isite,rankloc&
      &,igrd0,jgrd0,kgrd0,i,j,k,i0,j0,k0,it1,it2,it3&
      &,istart,iend,jstart,jend,kstart,kend&
      &,istartl,iendl,jstartl,jendl,kstartl,kendl&
      &,iproc,proc,commloc,nprocloc

   real(t_p) f,fi,dn1,dn2,dn3,de1,de2,de3,t1,t2,t3&
      &,dt1,dt2,dt3,term

   if (use_pmecore) then
      nprocloc = nrec
      commloc  = comm_rec
      rankloc = rank_bis
   else
      nprocloc = nproc
      commloc  = COMM_TINKER
      rankloc = rank
   end if
   if (deb_Path) write(*,'(3x,a)') "grid_pchg_force"

   istartl = istart1(rankloc+1)
   jstartl = jstart1(rankloc+1)
   kstartl = kstart1(rankloc+1)
   iendl   =   iend1(rankloc+1)
   jendl   =   jend1(rankloc+1)
   kendl   =   kend1(rankloc+1)
   f       = electric / dielec
   dn1     = real(nfft1,t_p)
   dn2     = real(nfft2,t_p)
   dn3     = real(nfft3,t_p)

!$acc parallel loop gang vector async(rec_queue) &
!$acc         present(chgrecglob,iion,locrec,igrid, &
!$acc     pchg,thetai1_p,thetai2_p,thetai3_p,qgridin_2d, &
!$acc     decrec)
   do isite = 1, nionrecloc
      iichg = chgrecglob(isite)
      iglob = iion(iichg)
      iloc  = locrec(iglob)
      iatm  = iglob
      igrd0 = igrid(1,iatm)
      jgrd0 = igrid(2,iatm)
      kgrd0 = igrid(3,iatm)
      fi    = f * pchg(iichg)
      de1   = 0.0_ti_p
      de2   = 0.0_ti_p
      de3   = 0.0_ti_p
      k0    = kgrd0
!$acc loop seq
      do it3 = 1, bsorder
         k0  = k0 + 1
         k   = k0 + 1 + (nfft3-sign(nfft3,k0))/2
         t3  = thetai3_p(1,it3,isite)
         dt3 = dn3 * thetai3_p(2,it3,isite)
         j0  = jgrd0
!$acc loop seq
         do it2 = 1, bsorder
            j0  = j0 + 1
            j   = j0 + 1 + (nfft2-sign(nfft2,j0))/2
            t2  = thetai2_p(1,it2,isite)
            dt2 = dn2 * thetai2_p(2,it2,isite)
            i0  = igrd0
!$acc loop seq
            do it1 = 1, bsorder
               i0  = i0 + 1
               i   = i0 + 1 + (nfft1-sign(nfft1,i0))/2
               t1  = thetai1_p(1,it1,isite)
               dt1 = dn1 * thetai1_p(2,it1,isite)
!
               if (((k.ge.kstartl).and.(k.le.kendl)).and.&
                  &((j.ge.jstartl).and.(j.le.jendl)).and.&
                  &((i.ge.istartl).and.(i.le.iendl)))  then
                  term = qgridin_2d(1,i-istartl+1,j-jstartl+1,&
                     &k-kstartl+1,1)
                  goto 100
               end if
!
!$acc loop seq
               do iproc = 1, nrec_send
                  proc   = prec_send(iproc)
                  kstart = kstart1(proc+1)
                  kend   = kend1(proc+1)
                  jstart = jstart1(proc+1)
                  jend   = jend1(proc+1)
                  istart = istart1(proc+1)
                  iend   = iend1(proc+1)
                  if (((k.ge.kstart).and.(k.le.kend)).and.&
                     &((j.ge.jstart).and.(j.le.jend)).and.&
                     &((i.ge.istart).and.(i.le.iend))) then
                     term=qgridin_2d(1,i-istart+1,j-jstart+1&
                        &,k-kstart+1,iproc+1)
                     goto 100
                  end if
               end do
100            continue
!
               de1 = de1 + term*dt1*t2*t3
               de2 = de2 + term*dt2*t1*t3
               de3 = de3 + term*dt3*t1*t2
            end do
         end do
      end do
      decrec(1,iloc) =decrec(1,iloc)+fi*(recip(1,1)*de1+recip(1,2)*de2&
         &+recip(1,3)*de3)
      decrec(2,iloc) =decrec(2,iloc)+fi*(recip(2,1)*de1+recip(2,2)*de2&
         &+recip(2,3)*de3)
      decrec(3,iloc) =decrec(3,iloc)+fi*(recip(3,1)*de1+recip(3,2)*de2&
         &+recip(3,3)*de3)
   end do

end subroutine
!
!     "fphi_uind_site" extracts the induced dipole potential at the i-th site from
!     the particle mesh Ewald grid
!
!
subroutine fphi_uind_sitegpu1(fdip_phi1,fdip_phi2,fdip_sum_phi)
   use atmlst
   use domdec
   use fft
   use inform ,only: deb_Path
   use mpole  ,only: npolerecloc, ipole
   use pme
   use potent
   use utilgpu,only: rec_queue,ngangs_rec
   use tinheader ,only: ti_p
   implicit none
   integer istart,iend,jstart,jend,kstart,kend
   integer istat,ied,jstat,jed,kstat,ked
   integer i,j,k,impi
   integer iproc,proc
   integer isite,iatm,iipole
   integer i0,j0,k0
   integer it1,it2,it3
   integer igrd0,jgrd0,kgrd0
   real(t_p) v0,v1,v2,v3
   real(t_p) u0,u1,u2,u3
   real(t_p) t0,t1,t2,t3
   real(t_p) t0_1,t0_2,t1_1,t1_2
   real(t_p) t2_1,t2_2,tq_1,tq_2
   real(t_p) tu00,tu10,tu01,tu20,tu11
   real(t_p) tu02,tu30,tu21,tu12,tu03
   real(t_p) tu00_1,tu01_1,tu10_1
   real(t_p) tu00_2,tu01_2,tu10_2
   real(t_p) tu20_1,tu11_1,tu02_1
   real(t_p) tu20_2,tu11_2,tu02_2
   real(t_p) tuv100_1,tuv010_1,tuv001_1
   real(t_p) tuv100_2,tuv010_2,tuv001_2
   real(t_p) tuv200_1,tuv020_1,tuv002_1
   real(t_p) tuv110_1,tuv101_1,tuv011_1
   real(t_p) tuv200_2,tuv020_2,tuv002_2
   real(t_p) tuv110_2,tuv101_2,tuv011_2
   real(t_p) tuv000,tuv100,tuv010,tuv001
   real(t_p) tuv200,tuv020,tuv002,tuv110
   real(t_p) tuv101,tuv011,tuv300,tuv030
   real(t_p) tuv003,tuv210,tuv201,tuv120
   real(t_p) tuv021,tuv102,tuv012,tuv111
   real(t_p),dimension(:,:)::fdip_phi1,fdip_phi2
   real(t_p) fdip_sum_phi(:,:)

   if (deb_Path) write(*,'(4X,A)') "fphi_uind_sitegpu1"
   if (use_pmecore) then
      kstat = kstart1(rank_bis+1)
      ked   = kend1  (rank_bis+1)
      jstat = jstart1(rank_bis+1)
      jed   = jend1  (rank_bis+1)
      istat = istart1(rank_bis+1)
      ied   = iend1  (rank_bis+1)
   else
      kstat = kstart1(rank+1)
      ked   = kend1  (rank+1)
      jstat = jstart1(rank+1)
      jed   = jend1  (rank+1)
      istat = istart1(rank+1)
      ied   = iend1  (rank+1)
   end if
!
!$acc parallel loop num_gangs(ngangs_rec) &
!$acc         present(polerecglob,ipole,igrid,thetai1,thetai2, &
!$acc   kstart1,kend1,jstart1,jend1,istart1,iend1,fdip_phi1, &
!$acc   fdip_phi2,fdip_sum_phi) &
!$acc         async(rec_queue)
   do impi = 1,npolerecloc
      iipole   = polerecglob(impi)
      iatm     = ipole(iipole)
      igrd0    = igrid(1,iatm)
      jgrd0    = igrid(2,iatm)
      kgrd0    = igrid(3,iatm)
      tuv100_1 = 0.0_ti_p
      tuv010_1 = 0.0_ti_p
      tuv001_1 = 0.0_ti_p
      tuv200_1 = 0.0_ti_p
      tuv020_1 = 0.0_ti_p
      tuv002_1 = 0.0_ti_p
      tuv110_1 = 0.0_ti_p
      tuv101_1 = 0.0_ti_p
      tuv011_1 = 0.0_ti_p
      tuv100_2 = 0.0_ti_p
      tuv010_2 = 0.0_ti_p
      tuv001_2 = 0.0_ti_p
      tuv200_2 = 0.0_ti_p
      tuv020_2 = 0.0_ti_p
      tuv002_2 = 0.0_ti_p
      tuv110_2 = 0.0_ti_p
      tuv101_2 = 0.0_ti_p
      tuv011_2 = 0.0_ti_p
      tuv000   = 0.0_ti_p
      tuv001   = 0.0_ti_p
      tuv010   = 0.0_ti_p
      tuv100   = 0.0_ti_p
      tuv200   = 0.0_ti_p
      tuv020   = 0.0_ti_p
      tuv002   = 0.0_ti_p
      tuv110   = 0.0_ti_p
      tuv101   = 0.0_ti_p
      tuv011   = 0.0_ti_p
      tuv300   = 0.0_ti_p
      tuv030   = 0.0_ti_p
      tuv003   = 0.0_ti_p
      tuv210   = 0.0_ti_p
      tuv201   = 0.0_ti_p
      tuv120   = 0.0_ti_p
      tuv021   = 0.0_ti_p
      tuv102   = 0.0_ti_p
      tuv012   = 0.0_ti_p
      tuv111   = 0.0_ti_p
      k0       = kgrd0
!$acc    loop seq
      do it3 = 1, bsorder
         k0 = k0 + 1
!           k      = k0 + 1 + (nfft3-isign(nfft3,k0))/2
         k      = k0 + 1 + ishft(nfft3-isign(nfft3,k0),-1)
         v0     = thetai3(1,it3,impi)
         v1     = thetai3(2,it3,impi)
         v2     = thetai3(3,it3,impi)
         v3     = thetai3(4,it3,impi)
         tu00_1 = 0.0_ti_p
         tu01_1 = 0.0_ti_p
         tu10_1 = 0.0_ti_p
         tu20_1 = 0.0_ti_p
         tu11_1 = 0.0_ti_p
         tu02_1 = 0.0_ti_p
         tu00_2 = 0.0_ti_p
         tu01_2 = 0.0_ti_p
         tu10_2 = 0.0_ti_p
         tu20_2 = 0.0_ti_p
         tu11_2 = 0.0_ti_p
         tu02_2 = 0.0_ti_p
         tu00   = 0.0_ti_p
         tu10   = 0.0_ti_p
         tu01   = 0.0_ti_p
         tu20   = 0.0_ti_p
         tu11   = 0.0_ti_p
         tu02   = 0.0_ti_p
         tu30   = 0.0_ti_p
         tu21   = 0.0_ti_p
         tu12   = 0.0_ti_p
         tu03   = 0.0_ti_p
         j0 = jgrd0
!$acc    loop seq
         do it2 = 1, bsorder
            j0 = j0 + 1
!              j    = j0 + 1 + (nfft2-isign(nfft2,j0))/2
            j    = j0 + 1 + ishft(nfft2-isign(nfft2,j0),-1)
            u0   = thetai2(1,it2,impi)
            u1   = thetai2(2,it2,impi)
            u2   = thetai2(3,it2,impi)
            u3   = thetai2(4,it2,impi)
            t0_1 = 0.0_ti_p
            t1_1 = 0.0_ti_p
            t2_1 = 0.0_ti_p
            t0_2 = 0.0_ti_p
            t1_2 = 0.0_ti_p
            t2_2 = 0.0_ti_p
            t3   = 0.0_ti_p
            i0   = igrd0
!$acc    loop seq
            do it1 = 1, bsorder
               i0 = i0 + 1
!                 i    = i0 + 1 + (nfft1-isign(nfft1,i0))/2
               i    = i0 + 1 + ishft(nfft1-isign(nfft1,i0),-1)
!
               tq_1 = 0_ti_p
               tq_2 = 0_ti_p
               if (((k.ge.kstat).and.(k.le.ked)).and.&
                  &((j.ge.jstat).and.(j.le.jed)).and.&
                  &((i.ge.istat).and.(i.le.ied))) then
                  tq_1 = qgrid2in_2d(1,i-istat+1,j-jstat+1,k-kstat+1,1)
                  tq_2 = qgrid2in_2d(2,i-istat+1,j-jstat+1,k-kstat+1,1)
                  goto 10
               end if
!$acc    loop seq
               do iproc = 1, nrec_send
                  proc   = prec_send(iproc)
                  kstart = kstart1(proc+1)
                  kend   = kend1  (proc+1)
                  jstart = jstart1(proc+1)
                  jend   = jend1  (proc+1)
                  istart = istart1(proc+1)
                  iend   = iend1  (proc+1)
                  if (((k.ge.kstart).and.(k.le.kend)).and.&
                     &((j.ge.jstart).and.(j.le.jend)).and.&
                     &((i.ge.istart).and.(i.le.iend))) then
                     tq_1 = qgrid2in_2d(1,i-istart+1,j-jstart+1,k-kstart+1,&
                        &iproc+1)
                     tq_2 = qgrid2in_2d(2,i-istart+1,j-jstart+1,k-kstart+1,&
                        &iproc+1)
                     goto 10
                  end if
               end do
               cycle
10             continue
!
               t0_1 = t0_1 +  tq_1      *thetai1(1,it1,impi)
               t1_1 = t1_1 +  tq_1      *thetai1(2,it1,impi)
               t2_1 = t2_1 +  tq_1      *thetai1(3,it1,impi)
               t0_2 = t0_2 +        tq_2*thetai1(1,it1,impi)
               t1_2 = t1_2 +        tq_2*thetai1(2,it1,impi)
               t2_2 = t2_2 +        tq_2*thetai1(3,it1,impi)
               t3   = t3   + (tq_1+tq_2)*thetai1(4,it1,impi)
            end do
            tu00_1 = tu00_1 + t0_1*u0
            tu10_1 = tu10_1 + t1_1*u0
            tu01_1 = tu01_1 + t0_1*u1
            tu20_1 = tu20_1 + t2_1*u0
            tu11_1 = tu11_1 + t1_1*u1
            tu02_1 = tu02_1 + t0_1*u2
            tu00_2 = tu00_2 + t0_2*u0
            tu10_2 = tu10_2 + t1_2*u0
            tu01_2 = tu01_2 + t0_2*u1
            tu20_2 = tu20_2 + t2_2*u0
            tu11_2 = tu11_2 + t1_2*u1
            tu02_2 = tu02_2 + t0_2*u2
            t0     = t0_1 + t0_2
            t1     = t1_1 + t1_2
            t2     = t2_1 + t2_2
            tu00   = tu00 + t0*u0
            tu10   = tu10 + t1*u0
            tu01   = tu01 + t0*u1
            tu20   = tu20 + t2*u0
            tu11   = tu11 + t1*u1
            tu02   = tu02 + t0*u2
            tu30   = tu30 + t3*u0
            tu21   = tu21 + t2*u1
            tu12   = tu12 + t1*u2
            tu03   = tu03 + t0*u3
         end do
         tuv100_1 = tuv100_1 + tu10_1*v0
         tuv010_1 = tuv010_1 + tu01_1*v0
         tuv001_1 = tuv001_1 + tu00_1*v1
         tuv200_1 = tuv200_1 + tu20_1*v0
         tuv020_1 = tuv020_1 + tu02_1*v0
         tuv002_1 = tuv002_1 + tu00_1*v2
         tuv110_1 = tuv110_1 + tu11_1*v0
         tuv101_1 = tuv101_1 + tu10_1*v1
         tuv011_1 = tuv011_1 + tu01_1*v1
         tuv100_2 = tuv100_2 + tu10_2*v0
         tuv010_2 = tuv010_2 + tu01_2*v0
         tuv001_2 = tuv001_2 + tu00_2*v1
         tuv200_2 = tuv200_2 + tu20_2*v0
         tuv020_2 = tuv020_2 + tu02_2*v0
         tuv002_2 = tuv002_2 + tu00_2*v2
         tuv110_2 = tuv110_2 + tu11_2*v0
         tuv101_2 = tuv101_2 + tu10_2*v1
         tuv011_2 = tuv011_2 + tu01_2*v1
         tuv000   = tuv000 + tu00*v0
         tuv100   = tuv100 + tu10*v0
         tuv010   = tuv010 + tu01*v0
         tuv001   = tuv001 + tu00*v1
         tuv200   = tuv200 + tu20*v0
         tuv020   = tuv020 + tu02*v0
         tuv002   = tuv002 + tu00*v2
         tuv110   = tuv110 + tu11*v0
         tuv101   = tuv101 + tu10*v1
         tuv011   = tuv011 + tu01*v1
         tuv300   = tuv300 + tu30*v0
         tuv030   = tuv030 + tu03*v0
         tuv003   = tuv003 + tu00*v3
         tuv210   = tuv210 + tu21*v0
         tuv201   = tuv201 + tu20*v1
         tuv120   = tuv120 + tu12*v0
         tuv021   = tuv021 + tu02*v1
         tuv102   = tuv102 + tu10*v2
         tuv012   = tuv012 + tu01*v2
         tuv111   = tuv111 + tu11*v1
      end do
      fdip_phi1   ( 2,impi) = tuv100_1
      fdip_phi1   ( 3,impi) = tuv010_1
      fdip_phi1   ( 4,impi) = tuv001_1
      fdip_phi1   ( 5,impi) = tuv200_1
      fdip_phi1   ( 6,impi) = tuv020_1
      fdip_phi1   ( 7,impi) = tuv002_1
      fdip_phi1   ( 8,impi) = tuv110_1
      fdip_phi1   ( 9,impi) = tuv101_1
      fdip_phi1   (10,impi) = tuv011_1
      fdip_phi2   ( 2,impi) = tuv100_2
      fdip_phi2   ( 3,impi) = tuv010_2
      fdip_phi2   ( 4,impi) = tuv001_2
      fdip_phi2   ( 5,impi) = tuv200_2
      fdip_phi2   ( 6,impi) = tuv020_2
      fdip_phi2   ( 7,impi) = tuv002_2
      fdip_phi2   ( 8,impi) = tuv110_2
      fdip_phi2   ( 9,impi) = tuv101_2
      fdip_phi2   (10,impi) = tuv011_2
      fdip_sum_phi( 1,impi) = tuv000
      fdip_sum_phi( 2,impi) = tuv100
      fdip_sum_phi( 3,impi) = tuv010
      fdip_sum_phi( 4,impi) = tuv001
      fdip_sum_phi( 5,impi) = tuv200
      fdip_sum_phi( 6,impi) = tuv020
      fdip_sum_phi( 7,impi) = tuv002
      fdip_sum_phi( 8,impi) = tuv110
      fdip_sum_phi( 9,impi) = tuv101
      fdip_sum_phi(10,impi) = tuv011
      fdip_sum_phi(11,impi) = tuv300
      fdip_sum_phi(12,impi) = tuv030
      fdip_sum_phi(13,impi) = tuv003
      fdip_sum_phi(14,impi) = tuv210
      fdip_sum_phi(15,impi) = tuv201
      fdip_sum_phi(16,impi) = tuv120
      fdip_sum_phi(17,impi) = tuv021
      fdip_sum_phi(18,impi) = tuv102
      fdip_sum_phi(19,impi) = tuv012
      fdip_sum_phi(20,impi) = tuv111
   end do
!
   return
end
!
!     "fphi_uind_sitegpu2" extracts the induced dipole potential at the i-th site from
!     the particle mesh Ewald grid second version
!
!
subroutine fphi_uind_sitegpu2(fdip_phi1,fdip_phi2)
   use atoms
   use atmlst
   use boxes
   use domdec
   use fft
   use inform ,only: deb_Path
   use mpole  ,only: npolerecloc, ipole
   use pme
   use potent
   use utils
   use utilgpu,only: rec_queue,ngangs_rec
   use tinheader ,only: ti_p
   implicit none
   real(t_p),dimension(:,:)::fdip_phi1,fdip_phi2

   integer istart,iend,jstart,jend,kstart,kend
   integer istat,ied,jstat,jed,kstat,ked
   integer i,j,k,impi,ifr
   integer iproc,proc
   integer isite,iatm,iipole
   integer i0,j0,k0
   integer it1,it2,it3
   integer igrd0,jgrd0,kgrd0
   integer ind
   real(t_p) xi,yi,zi
   real(t_p) w,fr,eps
   real(t_p) v0,v1,v2,v3
   real(t_p) u0,u1,u2,u3
   real(t_p) t0,t1,t2,t3
   real(t_p) t0_1,t0_2,t1_1,t1_2
   real(t_p) t2_1,t2_2,tq_1,tq_2
   real(t_p) tu00,tu10,tu01,tu20,tu11
   real(t_p) tu02,tu30,tu21,tu12,tu03
   real(t_p) tu00_1,tu01_1,tu10_1
   real(t_p) tu00_2,tu01_2,tu10_2
   real(t_p) tu20_1,tu11_1,tu02_1
   real(t_p) tu20_2,tu11_2,tu02_2
   real(t_p) tuv100_1,tuv010_1,tuv001_1
   real(t_p) tuv100_2,tuv010_2,tuv001_2
   real(t_p) tuv200_1,tuv020_1,tuv002_1
   real(t_p) tuv110_1,tuv101_1,tuv011_1
   real(t_p) tuv200_2,tuv020_2,tuv002_2
   real(t_p) tuv110_2,tuv101_2,tuv011_2
   real(t_p) tuv000,tuv100,tuv010,tuv001
   real(t_p) tuv200,tuv020,tuv002,tuv110
   real(t_p) tuv101,tuv011,tuv300,tuv030
   real(t_p) tuv003,tuv210,tuv201,tuv120
   real(t_p) tuv021,tuv102,tuv012,tuv111
   real(t_p) tuv1(9),tuv2(9)

   if (deb_Path) write(*,'(4X,A)') "fphi_uind_sitegpu2"
   if (use_pmecore) then
      kstat = kstart1(rank_bis+1)
      ked   = kend1  (rank_bis+1)
      jstat = jstart1(rank_bis+1)
      jed   = jend1  (rank_bis+1)
      istat = istart1(rank_bis+1)
      ied   = iend1  (rank_bis+1)
   else
      kstat = kstart1(rank+1)
      ked   = kend1  (rank+1)
      jstat = jstart1(rank+1)
      jed   = jend1  (rank+1)
      istat = istart1(rank+1)
      ied   = iend1  (rank+1)
   end if
!
!$acc parallel loop num_gangs(ngangs_rec) &
!$acc         present(polerecglob,ipole,igrid,thetai1,thetai2, &
!$acc   kstart1,kend1,jstart1,jend1,istart1,iend1) &
!$acc         present(fdip_phi1,fdip_phi2,qgrid2in_p) &
!$acc         private(tuv1,tuv2) &
!$acc         async(rec_queue)
   do impi = 1,npolerecloc
      iipole = polerecglob(impi)
      iatm   = ipole(iipole)
      igrd0  = igrid(1,iatm)
      jgrd0  = igrid(2,iatm)
      kgrd0  = igrid(3,iatm)

!$acc loop seq
      do it1 = 1,9
         tuv1(it1)  = 0.0_ti_p
         tuv2(it1)  = 0.0_ti_p
      end do
!        k0       = kgrd0
!$acc    loop seq
      do it3 = 1, bsorder
!           k0     = k0 + 1
!           k      = k0 + 1 + (nfft3-isign(nfft3,kgrd0+it3))/2
         k      = kgrd0 +it3 +1 +&
            &ishft(nfft3-isign(nfft3,kgrd0+it3),-1)
         v0     = thetai3(1,it3,impi)
         v1     = thetai3(2,it3,impi)
         v2     = thetai3(3,it3,impi)
         tu00_1 = 0.0_ti_p
         tu01_1 = 0.0_ti_p
         tu10_1 = 0.0_ti_p
         tu20_1 = 0.0_ti_p
         tu11_1 = 0.0_ti_p
         tu02_1 = 0.0_ti_p
         tu00_2 = 0.0_ti_p
         tu01_2 = 0.0_ti_p
         tu10_2 = 0.0_ti_p
         tu20_2 = 0.0_ti_p
         tu11_2 = 0.0_ti_p
         tu02_2 = 0.0_ti_p
!           j0 = jgrd0
!$acc    loop seq
         do it2 = 1, bsorder
!              j0 = j0 + 1
!              j    = jgrd0 +it2 + 1 + (nfft2-isign(nfft2,jgrd0+it2))/2
            j    = jgrd0 +it2 + 1 +&
               &ishft(nfft2-isign(nfft2,jgrd0+it2),-1)
            u0   = thetai2(1,it2,impi)
            u1   = thetai2(2,it2,impi)
            u2   = thetai2(3,it2,impi)
            t0_1 = 0.0_ti_p
            t1_1 = 0.0_ti_p
            t2_1 = 0.0_ti_p
            t0_2 = 0.0_ti_p
            t1_2 = 0.0_ti_p
            t2_2 = 0.0_ti_p
            i0   = igrd0
!$acc loop seq
            do it1 = 1, bsorder
!                 i0  = i0 + 1
!                 i    = igrd0 +it1 +1 + (nfft1-isign(nfft1,igrd0 +it1))/2
               i    = igrd0 +it1 +1 +&
                  &ishft(nfft1-isign(nfft1,igrd0+it1),-1)
!
               tq_1 = 0_ti_p
               tq_2 = 0_ti_p
               if (((k.ge.kstat).and.(k.le.ked)).and.&
                  &((j.ge.jstat).and.(j.le.jed)).and.&
                  &((i.ge.istat).and.(i.le.ied))) then
!               ind = 1+(i-istat)*2+(j-jstat)*2*n1mpimax+(k-kstat)*
!    &                2*n1mpimax*n2mpimax
!               tq_1 = qgrid2in_p(ind)
!               tq_2 = qgrid2in_p(ind+1)
                  tq_1 = qgrid2in_2d(1,i-istat+1,j-jstat+1,k-kstat+1,1)
                  tq_2 = qgrid2in_2d(2,i-istat+1,j-jstat+1,k-kstat+1,1)
                  goto 10
               end if
!$acc loop seq
               do iproc = 1, nrec_send
                  proc   = prec_send(iproc)
                  kstart = kstart1(proc+1)
                  kend   = kend1  (proc+1)
                  jstart = jstart1(proc+1)
                  jend   = jend1  (proc+1)
                  istart = istart1(proc+1)
                  iend   = iend1  (proc+1)
                  if (((k.ge.kstart).and.(k.le.kend)).and.&
                     &((j.ge.jstart).and.(j.le.jend)).and.&
                     &((i.ge.istart).and.(i.le.iend))) then
                     tq_1 = qgrid2in_2d(1,i-istart+1,j-jstart+1,k-kstart+1,&
                        &iproc+1)
                     tq_2 = qgrid2in_2d(2,i-istart+1,j-jstart+1,k-kstart+1,&
                        &iproc+1)
                     goto 10
                  end if
               end do
               cycle
10             continue
!
               t0   = thetai1(1,it1,impi)
               t1   = thetai1(2,it1,impi)
               t2   = thetai1(3,it1,impi)
               t0_1 = t0_1 + tq_1*t0
               t1_1 = t1_1 + tq_1*t1
               t2_1 = t2_1 + tq_1*t2
               t0_2 = t0_2 + tq_2*t0
               t1_2 = t1_2 + tq_2*t1
               t2_2 = t2_2 + tq_2*t2
            end do
            tu00_1 = tu00_1 + t0_1*u0
            tu10_1 = tu10_1 + t1_1*u0
            tu01_1 = tu01_1 + t0_1*u1
            tu20_1 = tu20_1 + t2_1*u0
            tu11_1 = tu11_1 + t1_1*u1
            tu02_1 = tu02_1 + t0_1*u2
            tu00_2 = tu00_2 + t0_2*u0
            tu10_2 = tu10_2 + t1_2*u0
            tu01_2 = tu01_2 + t0_2*u1
            tu20_2 = tu20_2 + t2_2*u0
            tu11_2 = tu11_2 + t1_2*u1
            tu02_2 = tu02_2 + t0_2*u2
         end do
         tuv1(1) = tuv1(1) + tu10_1*v0
         tuv1(2) = tuv1(2) + tu01_1*v0
         tuv1(3) = tuv1(3) + tu00_1*v1
         tuv1(4) = tuv1(4) + tu20_1*v0
         tuv1(5) = tuv1(5) + tu02_1*v0
         tuv1(6) = tuv1(6) + tu00_1*v2
         tuv1(7) = tuv1(7) + tu11_1*v0
         tuv1(8) = tuv1(8) + tu10_1*v1
         tuv1(9) = tuv1(9) + tu01_1*v1
         tuv2(1) = tuv2(1) + tu10_2*v0
         tuv2(2) = tuv2(2) + tu01_2*v0
         tuv2(3) = tuv2(3) + tu00_2*v1
         tuv2(4) = tuv2(4) + tu20_2*v0
         tuv2(5) = tuv2(5) + tu02_2*v0
         tuv2(6) = tuv2(6) + tu00_2*v2
         tuv2(7) = tuv2(7) + tu11_2*v0
         tuv2(8) = tuv2(8) + tu10_2*v1
         tuv2(9) = tuv2(9) + tu01_2*v1
      end do
!        call set_dip(impi,tuv1,fdip_phi1(1,1))
!        call set_dip(impi,tuv2,fdip_phi2(1,1))
!$acc loop seq
      do i = 2,10
         fdip_phi1(i,impi) = tuv1(i-1)
         fdip_phi2(i,impi) = tuv2(i-1)
      end do
   end do
end

subroutine pme_conv_gpu(e,vxx,vxy,vxz,vyy,vyz,vzz)
   use atmlst
   use atoms
   use bound
   use boxes
   use charge
   use chgpot
   use deriv
   use domdec
   use energi
   use ewald
   use fft
   use inform
   use math
   use pme
   use pme1
   use potent
   use utils
   use utilgpu
   use timestat
   implicit none
   real(r_p) e,vxx,vxy,vxz,vyy,vyz,vzz
   integer k1,k2,k3,m1,m2,m3,rankloc
   integer nf1,nf2,nf3,nff
   integer ist2,jst2,kst2,ien2,jen2,ken2
   integer kd
   real(t_p) e0,term,expterm
   real(t_p) vterm,pterm
   real(t_p) volterm
   real(t_p) f,denom
   real(t_p) hsq,struc2
   real(t_p) h1,h2,h3
   real(t_p) r1,r2,r3
!
!     use scalar sum to get reciprocal space energy and virial
!
   call timer_enter( timer_scalar )
   if(deb_Path) write(*,'(4x,A)') 'pme_conv_gpu'

   rankloc = merge(rank_bis,rank,use_pmecore)
   ist2    = istart2(rankloc+1)
   jst2    = jstart2(rankloc+1)
   kst2    = kstart2(rankloc+1)
   ien2    =   iend2(rankloc+1)
   jen2    =   jend2(rankloc+1)
   ken2    =   kend2(rankloc+1)
   f       = 0.5d0 * electric / dielec
   pterm   = (pi/aewald)**2
   volterm = pi * volbox
   nff     = nfft1 * nfft2
   nf1     = (nfft1+1) / 2
   nf2     = (nfft2+1) / 2
   nf3     = (nfft3+1) / 2

   if ((istart2(rankloc+1).eq.1).and.(jstart2(rankloc+1).eq.1).and.&
      &(kstart2(rankloc+1).eq.1)) then
!$acc serial async(rec_queue) present(qgridout_2d)
      qgridout_2d(1,1,1,1) = 0.0
      qgridout_2d(2,1,1,1) = 0.0
!$acc end serial
   end if

!$acc host_data use_device( &
!$acc          qgridout_2d,e,vxx,vxy,vxz,vyy,vyz,vzz)
!$acc parallel loop collapse(3) async(rec_queue) &
!$acc         deviceptr( &
!$acc         qgridout_2d,e,vxx,vxy,vxz,vyy,vyz,vzz) &
!$acc         reduction(+:e,vxx,vxy,vxz,vyy,vyz,vzz)
   do k3 = kst2,ken2
      do k2 = jst2,jen2
         do k1 = ist2,ien2
            m1 = k1 - 1
            m2 = k2 - 1
            m3 = k3 - 1
            if (k1 .gt. nf1)  m1 = m1 - nfft1
            if (k2 .gt. nf2)  m2 = m2 - nfft2
            if (k3 .gt. nf3)  m3 = m3 - nfft3
            if ((m1.eq.0).and.(m2.eq.0).and.(m3.eq.0)) goto 10
            r1 = real(m1,t_p)
            r2 = real(m2,t_p)
            r3 = real(m3,t_p)
            h1 = recip(1,1)*r1 + recip(1,2)*r2 + recip(1,3)*r3
            h2 = recip(2,1)*r1 + recip(2,2)*r2 + recip(2,3)*r3
            h3 = recip(3,1)*r1 + recip(3,2)*r2 + recip(3,3)*r3
            hsq = h1*h1 + h2*h2 + h3*h3
            term = -pterm * hsq
            expterm = 0.0_ti_p
            kd = (k1-ist2) + (k2-jst2)*(ien2-ist2+1)
            if (term .gt. -50.0_ti_p) then
               denom = volterm*hsq*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
               expterm = exp(term) / denom
               if (.not. use_bounds) then
                  expterm = expterm * (1.0_ti_p-cos(pi*xbox*sqrt(hsq)))
               else if (octahedron) then
                  if (mod(m1+m2+m3,2).ne.0)  expterm = 0.0_ti_p
               end if
               struc2 = qgridout_2d(1,k1-ist2+1,k2-jst2+1,k3-kst2+1)**2&
                  &+ qgridout_2d(2,k1-ist2+1,k2-jst2+1,k3-kst2+1)**2
               e0 = f * expterm * struc2
               e = e + e0
               vterm = (2.0_ti_p/hsq) * (1.0_ti_p-term) * e0
               vxx   = vxx + h1*h1*vterm - e0
               vxy   = vxy + h1*h2*vterm
               vxz   = vxz + h1*h3*vterm
               vyy   = vyy + h2*h2*vterm - e0
               vyz   = vyz + h3*h2*vterm
               vzz   = vzz + h3*h3*vterm - e0
            end if
            qgridout_2d(1,k1-ist2+1,k2-jst2+1,k3-kst2+1) = expterm *&
               &qgridout_2d(1,k1-ist2+1,k2-jst2+1,k3-kst2+1)

            qgridout_2d(2,k1-ist2+1,k2-jst2+1,k3-kst2+1) = expterm *&
               &qgridout_2d(2,k1-ist2+1,k2-jst2+1,k3-kst2+1)
10          continue
         end do
      end do
   end do
!$acc end host_data
!
!     account for zeroth grid point for nonperiodic system
!
   if ((istart2(rankloc+1).eq.1).and.(jstart2(rankloc+1).eq.1).and.&
      &(kstart2(rankloc+1).eq.1)) then
      if (.not. use_bounds) then
         expterm = 0.5_ti_p * pi / xbox
!$acc serial async(rec_queue) present(ecrec,qgridout_2d)
         struc2 = qgridout_2d(1,1,1,1)**2 + qgridout_2d(2,1,1,1)**2
         e = f * expterm * struc2
         ecrec = ecrec + e
         qgridout_2d(1,1,1,1) = expterm * qgridout_2d(1,1,1,1)
         qgridout_2d(2,1,1,1) = expterm * qgridout_2d(2,1,1,1)
!$acc end serial
      end if
   end if
   call timer_exit ( timer_scalar,quiet_timers )
end subroutine

 ! Scalar sum for reciprocal space (dispersion)
subroutine pme_convd_gpu(e,vxx,vxy,vxz,vyy,vyz,vzz)
   use atmlst
   use atoms
   use bound
   use boxes
   use charge
   use chgpot
   use deriv
   use domdec
   use energi
   use ewald
   use fft
   use inform
   use math
   use pme
   use pme1
   use potent
   use utils
   use utilgpu
   use timestat
   implicit none
   real(r_p) e,vxx,vxy,vxz,vyy,vyz,vzz
   integer k1,k2,k3,m1,m2,m3,rankloc
   integer nf1,nf2,nf3,nff
   integer ist2,jst2,kst2,ien2,jen2,ken2
   real(t_p) e0,bfac,fac1,fac2,fac3,term,expterm,erfcterm
   real(t_p) vterm,pterm
   real(t_p) h,b,hhh,denom,denom0
   real(t_p) hsq,struc2
   real(t_p) h1,h2,h3
   real(t_p) r1,r2,r3
!
!     use scalar sum to get reciprocal space energy and virial
!
   call timer_enter( timer_scalar )
   if(deb_Path) write(*,'(4x,A)') 'pme_convd_gpu'

   rankloc = merge(rank_bis,rank,use_pmecore)
   ist2    = istart2(rankloc+1)
   jst2    = jstart2(rankloc+1)
   kst2    = kstart2(rankloc+1)
   ien2    =   iend2(rankloc+1)
   jen2    =   jend2(rankloc+1)
   ken2    =   kend2(rankloc+1)
   bfac    = pi / aewald
   fac1    = 2.0*pi**(3.5)
   fac2    = aewald**3
   fac3    = -2.0d0*aewald*pi**2
   denom0 = (6.0d0*volbox)/(pi**1.5d0)
   pterm   = (pi/aewald)**2
   nff     = nfft1 * nfft2
   nf1     = (nfft1+1) / 2
   nf2     = (nfft2+1) / 2
   nf3     = (nfft3+1) / 2

   if ((istart2(rankloc+1).eq.1).and.(jstart2(rankloc+1).eq.1).and.&
      &(kstart2(rankloc+1).eq.1)) then
!$acc serial async(rec_queue) present(qgridout_2d)
      qgridout_2d(1,1,1,1) = 0.0d0
      qgridout_2d(2,1,1,1) = 0.0d0
!$acc end serial
   end if

!$acc host_data use_device( &
!$acc          qgridout_2d,e,vxx,vxy,vxz,vyy,vyz,vzz)
!$acc parallel loop collapse(3) async(rec_queue) &
!$acc         deviceptr( &
!$acc         qgridout_2d,e,vxx,vxy,vxz,vyy,vyz,vzz) &
!$acc         reduction(+:e,vxx,vxy,vxz,vyy,vyz,vzz)
   do k3 = kst2,ken2
      do k2 = jst2,jen2
         do k1 = ist2,ien2
            m1 = k1 - 1
            m2 = k2 - 1
            m3 = k3 - 1
            if (k1 .gt. nf1)  m1 = m1 - nfft1
            if (k2 .gt. nf2)  m2 = m2 - nfft2
            if (k3 .gt. nf3)  m3 = m3 - nfft3
            if ((m1.eq.0).and.(m2.eq.0).and.(m3.eq.0)) goto 10
            r1 = real(m1,t_p)
            r2 = real(m2,t_p)
            r3 = real(m3,t_p)
            h1 = recip(1,1)*r1 + recip(1,2)*r2 + recip(1,3)*r3
            h2 = recip(2,1)*r1 + recip(2,2)*r2 + recip(2,3)*r3
            h3 = recip(3,1)*r1 + recip(3,2)*r2 + recip(3,3)*r3
            hsq = h1*h1 + h2*h2 + h3*h3
            h   = sqrt(hsq)
            b   = h*bfac
            hhh = h*hsq
            term = -b * b
            expterm = 0.0_ti_p
            denom = denom0*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
            if (term .gt. -50.0_ti_p) then
               expterm = exp(term)
               erfcterm = erfc(b)
               if (.not. use_bounds) then
                  expterm =  expterm * (1.0-cos(pi*xbox*h))
                  erfcterm = erfcterm * (1.0-cos(pi*xbox*h))
               else if (octahedron) then
                  if (mod(m1+m2+m3,2).ne.0) then
                     expterm = 0.0_ti_p
                     erfcterm = 0.0_ti_p
                  end if
               end if
               term = fac1*erfcterm*hhh + expterm*(fac2 + fac3*hsq)
               struc2 = qgridout_2d(1,k1-ist2+1,k2-jst2+1,k3-kst2+1)**2&
                  &+ qgridout_2d(2,k1-ist2+1,k2-jst2+1,k3-kst2+1)**2
               e0 = -(term / denom) * struc2
               e = e + e0
               vterm = 3.0*(fac1*erfcterm*h + fac3*expterm)*struc2/denom
               vxx   = vxx + h1*h1*vterm - e0
               vxy   = vxy + h1*h2*vterm
               vxz   = vxz + h1*h3*vterm
               vyy   = vyy + h2*h2*vterm - e0
               vyz   = vyz + h3*h2*vterm
               vzz   = vzz + h3*h3*vterm - e0
            end if
            qgridout_2d(1,k1-ist2+1,k2-jst2+1,k3-kst2+1) = expterm *&
               &qgridout_2d(1,k1-ist2+1,k2-jst2+1,k3-kst2+1)

            qgridout_2d(2,k1-ist2+1,k2-jst2+1,k3-kst2+1) = expterm *&
               &qgridout_2d(2,k1-ist2+1,k2-jst2+1,k3-kst2+1)
10          continue
         end do
      end do
   end do
!$acc end host_data
   call timer_exit ( timer_scalar,quiet_timers )
end subroutine

#ifdef _CUDA
subroutine grid_mpole_sitecu(fmpvec)
   use atmlst    ,only:polerecglob
   use chunks
   use domdec
   use fft
   use inform    ,only:deb_Path
   use mpole
   use neigh     ,only:celle_glob
   use pme
   use potent
   use utils
   use pmestuffcu,only:grid_mpole_sitecu_core_1p
   use utilcu    ,only:PME_BLOCK_DIM1,check_launch_kernel
   use utilgpu   ,only:rec_stream,rec_queue,nSMP
   implicit none
   integer istat,ied,jstat,jed,kstat,ked
   integer i,j,k,impi
   integer twonlpts_1,twonlpts_12,nlptsit,rankloc
   real(t_p) fmpvec(10,max(npolerecloc,1))
   integer,save::gS
   logical,save::first_in=.true.

   if (deb_Path) print '(4x,A)', 'grid_mpole_sitecu'
   rankloc = merge(rank_bis,rank,use_pmecore)
   kstat   = kstart1(rankloc+1)
   ked     = kend1  (rankloc+1)
   jstat   = jstart1(rankloc+1)
   jed     = jend1  (rankloc+1)
   istat   = istart1(rankloc+1)
   ied     = iend1  (rankloc+1)
   twonlpts_1  = 2*nlpts+1
   twonlpts_12 = twonlpts_1**2
   nlptsit     = (2*nlpts+1)**3-1

   if (first_in) then
      call cudaMaxgridsize("grid_mpole_core",gS)
   end if
   call set_pme_texture
!
!     put the permanent multipole moments onto the grid
!     Work only with 1 MPI process
!
!$acc host_data use_device(polerecglob,celle_glob,ipole,igrid, &
!$acc        thetai1,thetai2,thetai3,fmpvec,qgridin_2d)

   call grid_mpole_sitecu_core_1p<<<gS,PME_BLOCK_DIM1,0,rec_stream>>>&
      &(celle_glob,ipole,igrid,thetai1,thetai2,thetai3&
      &,kstat,ked,jstat,jed,istat,ied&
      &,nfft1,nfft2,nfft3,npolerecloc&
      &,nlpts,twonlpts_1,twonlpts_12,nlptsit,grdoff&
      &,bsorder,n1mpimax,n2mpimax,n3mpimax,nrec_send&
      &,fmpvec,qgridin_2d,first_in)
   call check_launch_kernel(" grid_mpole_sitecu_core_1p")

!$acc end host_data
   if (first_in) first_in=.false.
end subroutine
!
!     "grid_pchg_sitecu" places the i-th fractional atomic charge onto
!     the particle mesh Ewald grid (CUDA Routine)
!
subroutine grid_pchg_sitecu
   use charge ,only: pchg,nionrecloc
   use neigh  ,only: celle_chg, celle_glob, celle_x, celle_y, celle_z
   implicit none
   call grid_put_sitecu( celle_chg,celle_glob,pchg,celle_x,celle_y&
      &,celle_z,nionrecloc,0 )
end subroutine
!
!     "grid_disp_sitecu" places the i-th damped dispersion coefficients onto
!     the particle mesh Ewald grid
!
subroutine grid_disp_sitecu
   use atmlst
   use disp
   use neigh  ,only: s_kind,sgl_id,so_x,so_y,so_z
   implicit none
   call grid_put_sitecu( s_kind,sgl_id,csix,so_x,so_y,so_z&
      &, ndisprecloc,1 )
end subroutine

subroutine grid_put_sitecu( kind_id,atom_id,entity,x,y,z,nk,prm )
   use chunks
   use domdec
   use fft
   use inform  ,only: deb_Path
   use mpole
   use pme
   use pmestuffcu,only: grid_put_site_kcu1
   use potent
   use sizes
   use utils
   use utilcu  ,only: PME_GRID_BDIM,check_launch_kernel
   use utilgpu ,only: rec_stream
   implicit none
   integer  ,intent(in):: nk, prm, kind_id(*),atom_id(*)
   real(t_p),intent(in):: entity(*),x(*),y(*),z(*)

   integer istart,iend,jstart,jend,kstart,kend
   integer istat,ied,jstat,jed,kstat,ked
   integer iichg
   integer i,j,k,m,impi,rankloc
   integer mk,mj
   integer iproc,proc,tsize
   integer ii,jj,kk
   integer isite,iatm
   integer offsetx,offsety
   integer offsetz
   integer location,twonlpts_1,twonlpts_12
   integer nlptsit
   real(t_p) q
   real(t_p) v0,u0,t0
   real(t_p) term,term0
   logical,save::f_in=.true.
   integer,save::gS
!
   if (deb_Path) write(*,'(3x,a,I3)') "grid_put_sitecu",prm
   if (use_pmecore) then
      rankloc  = rank_bis
   else
      rankloc  = rank
   end if

   kstat = kstart1(rankloc+1)
   ked   = kend1  (rankloc+1)
   jstat = jstart1(rankloc+1)
   jed   = jend1  (rankloc+1)
   istat = istart1(rankloc+1)
   ied   = iend1  (rankloc+1)
   twonlpts_1  = 2*nlpts+1
   twonlpts_12 = twonlpts_1**2
   nlptsit     = (2*nlpts+1)**3 - 1

   if (f_in) then
      call cudaMaxgridsize("grid_put_site_kcu1",gS)
   end if
   call set_pme_texture
!
!     put the permanent multipole moments onto the grid
!
   if (nproc.eq.1.or.nrec.eq.1) then
!$acc host_data use_device(kind_id,atom_id,igrid,entity &
!$acc         ,thetai1_p,thetai2_p,thetai3_p,qgridin_2d &
!$acc         ,x,y,z)
      call grid_put_site_kcu1<<<gS,PME_GRID_BDIM,0,rec_stream>>>&
         &(kind_id,atom_id,igrid,entity,thetai1_p,thetai2_p,thetai3_p&
         &,x,y,z&
         &,kstat,ked,jstat,jed,istat,ied&
         &,nfft1,nfft2,nfft3,nk&
         &,nlpts,twonlpts_1,twonlpts_12,nlptsit,grdoff&
         &,bsorder,n1mpimax,n2mpimax,n3mpimax,nrec_send&
         &,qgridin_2d,f_in)
      call check_launch_kernel(" grid_put_site_kcu1")
!$acc end host_data
   else
      print*,'FATAL ERROR, grid_put_site_kcu1 is sequential specific'
      __TINKER_FATAL__
   end if
   if (f_in) f_in=.false.
end subroutine

subroutine grid_uind_sitecu(fuindvec,fuinpvec,qgrid2loc)
   use atmlst    ,only:polerecglob
   use chunks
   use domdec
   use fft
   use mpole
   use neigh  ,only:celle_glob
   use pme
   use potent
   use utils
   use pmestuffcu,only:grid_uind_sitecu_core_1p
   use utilcu    ,only:PME_BLOCK_DIM1,check_launch_kernel
   use utilgpu   ,only:rec_stream,rec_queue,nSMP
   implicit none
   integer istat,ied,jstat,jed,kstat,ked
   integer i,j,k,impi
   integer twonlpts_1,twonlpts_12,nlptsit,rankloc
   real(t_p),intent(in),dimension(3,npolerecloc)::fuindvec,fuinpvec
   real(t_p) qgrid2loc(:,:,:,:,:)
   integer,save::gS
   logical,save::first_in=.true.
!$acc routine(adjust) seq
   if (use_pmecore) then
      rankloc = rank_bis
   else
      rankloc = rank
   end if
   kstat = kstart1(rankloc+1)
   ked   = kend1  (rankloc+1)
   jstat = jstart1(rankloc+1)
   jed   = jend1  (rankloc+1)
   istat = istart1(rankloc+1)
   ied   = iend1  (rankloc+1)
   twonlpts_1  = 2*nlpts+1
   twonlpts_12 = twonlpts_1**2
   nlptsit     = (2*nlpts+1)**3-1

   if (first_in) then
      call cudaMaxgridsize("grid_uind_core",gS)
   end if
   call set_pme_texture
!
!     put the induced dipole moments onto the grid
!     Work only with 1 MPI process
!
!$acc host_data use_device(celle_glob,ipole,igrid,thetai1,thetai2 &
!$acc        ,thetai3,fuindvec,fuinpvec,qgrid2loc)

   call grid_uind_sitecu_core_1p<<<gS,PME_BLOCK_DIM1,0,rec_stream>>>&
      &(celle_glob,ipole,igrid,thetai1,thetai2,thetai3&
      &,kstat,ked,jstat,jed,istat,ied&
      &,nfft1,nfft2,nfft3,npolerecloc&
      &,nlpts,twonlpts_1,twonlpts_12,nlptsit,grdoff&
      &,bsorder,n1mpimax,n2mpimax,n3mpimax,nrec_send&
      &,fuindvec,fuinpvec,qgrid2loc,first_in)
   call check_launch_kernel(" grid_uind_sitecu_core_1p")

!$acc end host_data
   if (first_in) first_in=.false.
end subroutine

subroutine fphi_mpole_sitecu
   use atoms
   use atmlst
   use boxes
   use domdec
   use fft
   use mpole
   use pme
   use potent
   use pmestuffcu,only:fphi_mpole_core
   use utilcu    ,only:PME_BLOCK_DIM,check_launch_kernel
   use utilgpu   ,only:rec_stream,rec_queue,nSMP
   implicit none
   integer istat,ied,jstat,jed,kstat,ked
   integer i,rankloc
   integer,save::gS
   logical,save::first_in=.true.
!
   if (use_pmecore) then
      rankloc  = rank_bis
   else
      rankloc  = rank
   end if
   kstat = kstart1(rankloc+1)
   ked   = kend1  (rankloc+1)
   jstat = jstart1(rankloc+1)
   jed   = jend1  (rankloc+1)
   istat = istart1(rankloc+1)
   ied   = iend1  (rankloc+1)

   if (first_in) then
      call cudaMaxgridsize("fphi_mpole_core",gS)
      first_in=.false.
   end if
   call set_PME_texture

!$acc host_data use_device(polerecglob,ipole,igrid,fphirec)
   call fphi_mpole_core<<<gS,PME_BLOCK_DIM,0,rec_stream>>>&
      &(kstat,ked,jstat,jed,istat,ied,nfft1,nfft2,nfft3&
      &,npolerecloc,nlocrec,bsorder,nrec_send,nproc&
      &,polerecglob,ipole,igrid,fphirec)
   call check_launch_kernel("fphi_mpole_core")
!$acc end host_data

end subroutine

subroutine grid_disp_forcecu
   use atoms  ,only: n
   use atmlst ,only: disprecglob
   use boxes  ,only: recip
   use disp
   use domdec ,only: rank,rank_bis,nproc,nrec,ndir,comm_rec&
      &,COMM_TINKER,nrec_send,locrec,prec_send,nbloc
   use deriv  ,only: dedsprec
   use fft
   use inform ,only: deb_Path
   use pme    ,only: igrid,nfft1,nfft2,nfft3,bsorder&
      &,thetai1_p,thetai2_p,thetai3_p,qgridin_2d
   use pmestuffcu,only: grid_calc_frc_kcu,grid_calc_frc_kcu1
   use potent ,only: use_pmecore
   use utilcu ,only: PME_BLOCK_DIM,check_launch_kernel
   use utilgpu,only: rec_queue,rec_stream
   implicit none
   integer iichg,iglob,iloc,iatm,isite,rankloc&
      &,igrd0,jgrd0,kgrd0,i,j,k,i0,j0,k0,it1,it2,it3&
      &,istat,ied,jstat,jed,kstat,ked&
      &,iproc,proc,commloc,nprocloc
   integer,save:: gS,gS1
   real(t_p) f,fi,dn1,dn2,dn3,de1,de2,de3,t1,t2,t3&
      &,dt1,dt2,dt3,term
   logical,save:: f_in=.true.

   if (use_pmecore) then
      nprocloc = nrec
      commloc  = comm_rec
      rankloc  = rank_bis
   else
      nprocloc = nproc
      commloc  = COMM_TINKER
      rankloc  = rank
   end if
   if (deb_Path) write(*,'(3x,a)') ' grid_disp_forcecu'

   istat = istart1(rankloc+1)
   jstat = jstart1(rankloc+1)
   kstat = kstart1(rankloc+1)
   ied   =   iend1(rankloc+1)
   jed   =   jend1(rankloc+1)
   ked   =   kend1(rankloc+1)
   dn1   = real(nfft1,t_p)
   dn2   = real(nfft2,t_p)
   dn3   = real(nfft3,t_p)

   if (f_in) then
      call cudaMaxgridsize("grid_calc_frc_kcu",gS)
      call cudaMaxgridsize("grid_calc_frc_kcu1",gS1)
      !gS = gS - 2*nSMP
      f_in =.false.
   end if
   call set_pme_texture

!$acc host_data use_device(disprecglob,idisp,locrec,igrid,csix &
!$acc         ,thetai1_p,thetai2_p,thetai3_p,dedsprec)
   if (nproc.eq.1.or.nrec.eq.1) then
      call grid_calc_frc_kcu1<<<gS1,PME_BLOCK_DIM,0,rec_stream>>>&
         &(disprecglob,idisp,locrec,igrid,csix&
         &,thetai1_p,thetai2_p,thetai3_p&
         &,dedsprec&
         &,kstat,ked,jstat,jed,istat,ied&
         &,nrec_send,ndisprecloc,n,nfft1,nfft2,nfft3&
         &,1.0,2.0,dn1,dn2,dn3)
      call check_launch_kernel(" grid_calc_frc_kcu1")
   else
      call grid_calc_frc_kcu<<<gS,PME_BLOCK_DIM,0,rec_stream>>>&
         &(disprecglob,idisp,locrec,igrid,csix&
         &,thetai1_p,thetai2_p,thetai3_p&
         &,dedsprec&
         &,kstat,ked,jstat,jed,istat,ied&
         &,nrec_send,ndisprecloc,n,nfft1,nfft2,nfft3&
         &,1.0,2.0,dn1,dn2,dn3)
      call check_launch_kernel(" grid_calc_frc_kcu")
   end if
!$acc end host_data

end subroutine

subroutine grid_pchg_forcecu
   use atoms  ,only: n
   use atmlst ,only: chgrecglob
   use boxes  ,only: recip
   use charge
   use chgpot
   use domdec ,only: rank,rank_bis,nproc,nrec,ndir,comm_rec&
      &,COMM_TINKER,nrec_send,locrec,prec_send,nbloc
   use deriv  ,only: decrec
   use fft
   use inform ,only: deb_Path
   use pme    ,only: igrid,nfft1,nfft2,nfft3,bsorder&
      &,thetai1_p,thetai2_p,thetai3_p,qgridin_2d
   use pmestuffcu,only: grid_calc_frc_kcu,grid_calc_frc_kcu1
   use potent ,only: use_pmecore
   use utilcu ,only: PME_BLOCK_DIM,check_launch_kernel
   use utilgpu,only: rec_queue,rec_stream
   implicit none
   integer iichg,iglob,iloc,iatm,isite,rankloc&
      &,igrd0,jgrd0,kgrd0,i,j,k,i0,j0,k0,it1,it2,it3&
      &,istat,ied,jstat,jed,kstat,ked&
      &,iproc,proc,commloc,nprocloc
   integer,save:: gS,gS1
   real(t_p) f,fi,dn1,dn2,dn3,de1,de2,de3,t1,t2,t3&
      &,dt1,dt2,dt3,term
   logical,save:: f_in=.true.

   if (use_pmecore) then
      nprocloc = nrec
      commloc  = comm_rec
      rankloc = rank_bis
   else
      nprocloc = nproc
      commloc  = COMM_TINKER
      rankloc  = rank
   end if
   if (deb_Path) write(*,'(3x,a)') 'grid_pchg_forcecu'

   istat = istart1(rankloc+1)
   jstat = jstart1(rankloc+1)
   kstat = kstart1(rankloc+1)
   ied   =   iend1(rankloc+1)
   jed   =   jend1(rankloc+1)
   ked   =   kend1(rankloc+1)
   f     = electric / dielec
   dn1   = real(nfft1,t_p)
   dn2   = real(nfft2,t_p)
   dn3   = real(nfft3,t_p)

   if (f_in) then
      call cudaMaxgridsize("grid_calc_frc_kcu",gS)
      call cudaMaxgridsize("grid_calc_frc_kcu1",gS1)
      !gS = gS - 2*nSMP
      f_in=.false.
   end if
   call set_pme_texture

!$acc host_data use_device(chgrecglob,iion,locrec,igrid,pchg &
!$acc         ,thetai1_p,thetai2_p,thetai3_p,decrec)
   if (nproc.eq.1.or.nrec.eq.1) then
      call grid_calc_frc_kcu1<<<gS1,PME_BLOCK_DIM,0,rec_stream>>>&
         &(chgrecglob,iion,locrec,igrid,pchg&
         &,thetai1_p,thetai2_p,thetai3_p&
         &,decrec&
         &,kstat,ked,jstat,jed,istat,ied&
         &,nrec_send,nionrecloc,n,nfft1,nfft2,nfft3&
         &,f,1.0,dn1,dn2,dn3)
      call check_launch_kernel(" grid_calc_frc_kcu1")
   else
      call grid_calc_frc_kcu<<<gS,PME_BLOCK_DIM,0,rec_stream>>>&
         &(chgrecglob,iion,locrec,igrid,pchg&
         &,thetai1_p,thetai2_p,thetai3_p&
         &,decrec&
         &,kstat,ked,jstat,jed,istat,ied&
         &,nrec_send,nionrecloc,n,nfft1,nfft2,nfft3&
         &,f,1.0,dn1,dn2,dn3)
      call check_launch_kernel(" grid_calc_frc_kcu")
   end if
!$acc end host_data

end subroutine
!
!
!       "fphi_chg_sitecu" extracts the permanent charge potential on the i-th site from
!       the particle mesh Ewald grid ( CUDA Fortran wrapper )
!
!
subroutine fphi_chg_sitecu
   use atoms
   use atmlst
   use domdec
   use fft
   use inform ,only: deb_Path
   use charge
   use pme
   use potent
   use pmestuffcu ,only: fphi_chg_site_kcu
   use utilcu ,only: PME_BLOCK_DIM,check_launch_kernel
   use utilgpu
   implicit none
   integer istat,ied,jstat,jed,kstat,ked
   integer iproc,proc,rankloc
   integer,save:: gS
   logical,save:: f_in=.true.
!
   if (deb_Path) write(*,'(5x,a)') "fphi_chg_sitecu"
!
   rankloc = merge(rank_bis,rank,use_pmecore)
   kstat   = kstart1(rankloc+1)
   ked     = kend1  (rankloc+1)
   jstat   = jstart1(rankloc+1)
   jed     = jend1  (rankloc+1)
   istat   = istart1(rankloc+1)
   ied     = iend1  (rankloc+1)
!
   if (f_in) then
      call cudaMaxgridsize("fphi_chg_site_kcu",gS)
      f_in=.false.
   end if
   call set_pme_texture

!$acc host_data use_device(chgrecglob,iion,igrid,fphirec)
   call fphi_chg_site_kcu<<<gS,PME_BLOCK_DIM,0,rec_stream>>>&
      &(chgrecglob,iion,igrid,kstat,ked,jstat,jed,istat,ied&
      &,nrec_send,nfft1,nfft2,nfft3,nionrecloc,n&
      &,fphirec)
   call check_launch_kernel("fphi_chg_site_kcu")
!$acc end host_data

end

subroutine fphi_uind_sitecu2(fdip_phi1,fdip_phi2)
   use atoms
   use atmlst
   use boxes
   use domdec
   use fft
   use mpole     ,only:npolerecloc,ipole
   use pme
   use pmestuffcu,only:fphi_uind_sitecu2_core
   use potent
   use utilcu    ,only:PME_BLOCK_DIM,check_launch_kernel
   use utilgpu   ,only:rec_stream,rec_queue,nSMP
   implicit none
   integer istat,ied,jstat,jed,kstat,ked
   integer i
   integer iproc,proc
   integer ind,psize,ierr
   real(t_p),dimension(:,:)::fdip_phi1,fdip_phi2
   integer,pointer :: ipole_p(:)
   integer,save::gS
   logical,save::first_in=.true.

   if (use_pmecore) then
      kstat = kstart1(rank_bis+1)
      ked   = kend1  (rank_bis+1)
      jstat = jstart1(rank_bis+1)
      jed   = jend1  (rank_bis+1)
      istat = istart1(rank_bis+1)
      ied   = iend1  (rank_bis+1)
   else
      kstat = kstart1(rank+1)
      ked   = kend1  (rank+1)
      jstat = jstart1(rank+1)
      jed   = jend1  (rank+1)
      istat = istart1(rank+1)
      ied   = iend1  (rank+1)
   end if

   call set_pme_texture
   psize =  size(kstart1)
   if (first_in) then
      if (nproc.eq.1.or.nrec.eq.1) then
         call cudamaxgridSize("fphi_uind_sitecu2_core_1p",gS)
         gS = gS - 2*nSMP
      else
         call cudamaxgridSize("fphi_uind_sitecu2",gS)
         gS = gS - 2*nSMP
      end if
      first_in=.false.
   end if

!$acc host_data use_device(qgrid2in_2d,polerecglob,ipole,igrid &
!$acc    ,x,y,z,recip,thetai1,thetai2,thetai3,fdip_phi1,fdip_phi2 &
!$acc    ,kstart1,kend1,jstart1,jend1,istart1,iend1,prec_send)

   if (nproc.eq.1.or.nrec.eq.1) then
      call fphi_uind_sitecu2_core_1p<<<gS,PME_BLOCK_DIM,0,rec_stream>>>&
         &(kstat,ked,jstat,jed,istat,ied&
         &,npolerecloc,nlocrec,bsorder,n1mpimax,n2mpimax,n3mpimax&
         &,nfft1,nfft2,nfft3&
         &,qgrid2in_2d,polerecglob,ipole,igrid,x,y,z,recip&
         &,thetai1,thetai2,thetai3&
         &,fdip_phi1,fdip_phi2,first_in)
      call check_launch_kernel("fphi_uind_sitecu2_core_1p")
   else
      call fphi_uind_sitecu2_core<<<gS,PME_BLOCK_DIM,0,rec_stream>>>&
         &(kstat,ked,jstat,jed,istat,ied&
         &,npolerecloc,nlocrec,bsorder,n1mpimax,n2mpimax,n3mpimax&
         &,nrec_send,nproc,prec_send,psize,nfft1,nfft2,nfft3&
         &,kstart1,kend1,jstart1,jend1,istart1,iend1&
         &,qgrid2in_2d,polerecglob,ipole,igrid,x,y,z,recip&
         &,thetai1,thetai2,thetai3&
         &,fdip_phi1,fdip_phi2,first_in)
      call check_launch_kernel("fphi_uind_sitecu2_core")
   end if

!$acc end host_data

end subroutine

subroutine fphi_uind_sitecu1(fdip_phi1,fdip_phi2,fdip_sum_phi)
   use atoms
   use atmlst
   use boxes
   use domdec
   use fft
   use mpole     ,only:npolerecloc,ipole
   use pme
   use pmestuffcu,only:fphi_uind_sitecu1_core
   use potent
   use utilcu    ,only:PME_BLOCK_DIM,check_launch_kernel
   use utilgpu   ,only:rec_stream,rec_queue,nSMP
   implicit none
   integer istat,ied,jstat,jed,kstat,ked
   integer i
   integer iproc,proc
   integer isite,iatm,iipole
   real(t_p),dimension(:,:)::fdip_phi1,fdip_phi2
   real(t_p),dimension(:,:)::fdip_sum_phi
   integer,save::gS
   logical,save::first_in=.true.

   gS = 4*nSMP
   if (use_pmecore) then
      kstat = kstart1(rank_bis+1)
      ked   = kend1  (rank_bis+1)
      jstat = jstart1(rank_bis+1)
      jed   = jend1  (rank_bis+1)
      istat = istart1(rank_bis+1)
      ied   = iend1  (rank_bis+1)
   else
      kstat = kstart1(rank+1)
      ked   = kend1  (rank+1)
      jstat = jstart1(rank+1)
      jed   = jend1  (rank+1)
      istat = istart1(rank+1)
      ied   = iend1  (rank+1)
   end if

   if (first_in) then
      call cudamaxgridSize("fphi_uind_sitecu1",gS)
      first_in=.false.
   end if
   call set_pme_texture

!$acc host_data use_device(qgrid2in_2d,polerecglob,ipole,igrid &
!$acc   ,fdip_phi1,fdip_phi2,fdip_sum_phi &
!$acc   ,kstart1,kend1,jstart1,jend1,istart1,iend1,prec_send)

   call fphi_uind_sitecu1_core<<<gS,PME_BLOCK_DIM,0,rec_stream>>>&
      &(kstat,ked,jstat,jed,istat,ied,nfft1,nfft2,nfft3&
      &,npolerecloc,nlocrec,bsorder,nrec_send,nproc&
      &,polerecglob,ipole,igrid,fdip_phi1,fdip_phi2,fdip_sum_phi)
   call check_launch_kernel("fphi_uind_sitecu1_core")

!$acc end host_data

end subroutine

!
!     "pme_conv_cu" apply convolution on cuda device
!     use scalar sum to get reciprocal space energy and virial
!
subroutine pme_conv_cu(e,vxx,vxy,vxz,vyy,vyz,vzz)
   use bound     ,only: use_bounds
   use boxes     ,only: xbox,volbox
   use chgpot    ,only: dielec,electric
   use domdec    ,only: rank_bis,rank
   use energi    ,only: calc_e
   use ewald     ,only: aewald
   use fft       ,only: istart2,jstart2,kstart2,isize2,jsize2,ksize2
   use inform    ,only: deb_Path
   use math      ,only: pi
   use pme       ,only: bsmod1,bsmod2,bsmod3,qgridout_2d&
      &,nfft1,nfft2,nfft3
   use pmestuffcu,only: pme_conv_kcu
   use potent    ,only: use_pmecore
   use timestat
   use tinheader ,only: ti_p
   use utilcu    ,only: check_launch_kernel
   use utilgpu   ,only: reduce_energy_virial,ered_buff,vred_buf1&
      &,rec_queue,rec_stream
   use virial    ,only: use_virial
   implicit none
   real(r_p) e,vxx,vxy,vxz,vyy,vyz,vzz
   integer rankloc
   integer nf1,nf2,nf3,nff
   integer ist2,jst2,kst2,qsz1,qsz2,qsz12,qsz
   integer gDim,bDim,mgDim
   real(t_p) e0,term,expterm
   real(t_p) vterm,pterm
   real(t_p) volterm
   real(t_p) f,denom
   real(t_p) hsq,struc2
   real(t_p) h1,h2,h3
   real(t_p) r1,r2,r3
   real(t_p) xbox_
   parameter(bDim=128,mgDim=2**16)

   call timer_enter( timer_scalar )
   if(deb_Path) write(*,'(4x,A)') 'pme_conv_cu'

   rankloc = merge(rank_bis,rank,use_pmecore)
   ist2    = istart2(rankloc+1)
   jst2    = jstart2(rankloc+1)
   kst2    = kstart2(rankloc+1)
   qsz1    = isize2(rankloc+1)
   qsz2    = jsize2(rankloc+1)
   qsz12   = qsz1*qsz2
   qsz     = qsz12*ksize2(rankloc+1)
   f       = 0.5_ti_p * electric / dielec
   pterm   = (pi/aewald)**2
   volterm = pi * volbox
   nff     = nfft1 * nfft2
   nf1     = (nfft1+1)/2
   nf2     = (nfft2+1)/2
   nf3     = (nfft3+1)/2
   xbox_   = real(xbox,t_p)
   gDim    = min(qsz/(2*bDim),mgDim)

   if (ist2==1.and.jst2==1.and.kst2==1) then
!$acc serial async(rec_queue) present(qgridout_2d)
      qgridout_2d(1,1,1,1) = 0.0
      qgridout_2d(2,1,1,1) = 0.0
!$acc end serial
   end if

!$acc host_data use_device(bsmod1,bsmod2,bsmod3,qgridout_2d &
!$acc         ,ered_buff,vred_buf1)
   call pme_conv_kcu<<<gDim,bDim,0,rec_stream>>>&
      &(bsmod1,bsmod2,bsmod3&
      &,kst2,jst2,ist2,qsz1,qsz2,qsz12,qsz&
      &,nff,nf1,nf2,nf3,nfft1,nfft2,nfft3&
      &,f,pterm,volterm,xbox_,calc_e,use_bounds&
      &,qgridout_2d,ered_buff,vred_buf1)
   call check_launch_kernel(" pme_conv_cu")
!$acc end host_data

   if (calc_e.or.use_virial) then
      call reduce_energy_virial(e,vxx,vyz,vyz,vyy,vyz,vzz,&
         &ered_buff,vred_buf1,rec_queue)
   end if

   call timer_exit ( timer_scalar,quiet_timers )
end subroutine

 ! Scalar sum for reciprocal space (dispersion)
subroutine pme_convd_cu(e,vxx,vxy,vxz,vyy,vyz,vzz)
   use bound     ,only: use_bounds
   use boxes     ,only: xbox,volbox
   use chgpot    ,only: dielec,electric
   use domdec    ,only: rank_bis,rank
   use energi    ,only: calc_e
   use ewald     ,only: aewald
   use fft       ,only: istart2,jstart2,kstart2,isize2,jsize2,ksize2
   use inform    ,only: deb_Path
   use math      ,only: pi
   use pme       ,only: bsmod1,bsmod2,bsmod3,qgridout_2d&
      &,nfft1,nfft2,nfft3
   use pmestuffcu,only: pme_conv_kcu
   use potent    ,only: use_pmecore
   use timestat
   use tinheader ,only: ti_p
   use utilcu    ,only: check_launch_kernel
   use utilgpu   ,only: reduce_energy_virial,ered_buff,vred_buf1&
      &,rec_queue,rec_stream
   use virial    ,only: use_virial
   implicit none
   real(r_p) e,vxx,vxy,vxz,vyy,vyz,vzz
   integer k1,k2,k3,m1,m2,m3,rankloc,gDim,bDim,mgDim
   integer nf1,nf2,nf3,nff
   integer ist2,jst2,kst2,qsz1,qsz2,qsz12,qsz
   real(t_p) e0,bfac,fac1,fac2,fac3,term,expterm,erfcterm
   real(t_p) vterm,pterm
   real(t_p) h,b,hhh,denom,denom0
   real(t_p) hsq,struc2
   real(t_p) h1,h2,h3
   real(t_p) r1,r2,r3,xbox_
   parameter(bDim=128,mgDim=2**16)
!
!     use scalar sum to get reciprocal space energy and virial
!
   call timer_enter( timer_scalar )
   if(deb_Path) write(*,'(4x,A)') 'pme_convd_cu'

   rankloc = merge(rank_bis,rank,use_pmecore)
   ist2    = istart2(rankloc+1)
   jst2    = jstart2(rankloc+1)
   kst2    = kstart2(rankloc+1)
   qsz1    = isize2(rankloc+1)
   qsz2    = jsize2(rankloc+1)
   qsz12   = qsz1*qsz2
   qsz     = qsz12*ksize2(rankloc+1)
   bfac    = pi / aewald
   fac1    = 2.0*pi**(3.5)
   fac2    = aewald**3
   fac3    = -2.0d0*aewald*pi**2
   denom0 = (6.0d0*volbox)/(pi**1.5d0)
   pterm   = (pi/aewald)**2
   nff     = nfft1 * nfft2
   nf1     = (nfft1+1) / 2
   nf2     = (nfft2+1) / 2
   nf3     = (nfft3+1) / 2
   xbox_   = xbox
   gDim    = min(qsz/(2*bDim),mgDim)

   if (ist2==1.and.jst2==1.and.kst2==1) then
!$acc serial async(rec_queue) present(qgridout_2d)
      qgridout_2d(1,1,1,1) = 0.0
      qgridout_2d(2,1,1,1) = 0.0
!$acc end serial
   end if

!$acc host_data use_device(bsmod1,bsmod2,bsmod3,qgridout_2d &
!$acc         ,ered_buff,vred_buf1)
   call pme_convd_kcu<<<gDim,bDim,0,rec_stream>>>&
      &(bsmod1,bsmod2,bsmod3&
      &,kst2,jst2,ist2,qsz1,qsz2,qsz12,qsz&
      &,nff,nf1,nf2,nf3,nfft1,nfft2,nfft3&
      &,bfac,fac1,fac2,fac3,denom0,pterm,xbox_,calc_e,use_bounds&
      &,qgridout_2d,ered_buff,vred_buf1)
   call check_launch_kernel(" pme_convd_kcu")
!$acc end host_data

   if (calc_e.or.use_virial) then
      call reduce_energy_virial(e,vxx,vyz,vyz,vyy,vyz,vzz&
         &,ered_buff,vred_buf1,rec_queue)
   end if

   call timer_exit ( timer_scalar,quiet_timers )

end subroutine

subroutine cudaMaxGridSize(kernelname,gS)
   use cudafor
   use echargecu  ,only: ecreal1_kcu,ecreal3_kcu
   use eChgLjcu
   use empole1cu  ,only: emreal1_kcu,emreal3_kcu
   use epolar1cu  ,only: epreal1c_core_cu&
      &,epreal3_cu,mpreal1c_core_cu
   use pmestuffcu
   use tmatxb_pmecu
   use utilcu
   use utilgpu
   implicit none
   character(*),intent(in)::kernelname
   integer,intent(out)::gS
   integer sbytes,ierr

200 format("error",I5," returned when calling Cuda0ccupancy"//&
      &"MaxActiveBlocksPerMulti in cudamaxGridSize subroutine",&
      &/,A100)

201 format("cudaMaxGridSize argument ( ", A,&
      &" ) does not match any known Cuda Kernel", /,&
      &" Returning 0 as gridSize ")

   if      (kernelname.eq."fphi_uind_sitecu2") then
      ierr = CUDAOCCUPANCYMAXACTIVEBLOCKSPERMULTIPROCESSOR&
         &(gS,fphi_uind_sitecu2_core,PME_BLOCK_DIM,0)
      gS   = gS*nSMP
      if (ierr.ne.cudaSuccess)&
         &print 200, ierr,cudageterrorstring(ierr)
   else if (kernelname.eq."fphi_uind_sitecu2_core_1p") then
      ierr = CUDAOCCUPANCYMAXACTIVEBLOCKSPERMULTIPROCESSOR&
         &(gS,fphi_uind_sitecu2_core_1p,PME_BLOCK_DIM,0)
      if (ierr.ne.cudaSuccess)&
         &print 200, ierr,cudageterrorstring(ierr)
      gS   = gS*nSMP
   else if (kernelname.eq."fphi_uind_sitecu1") then
      ierr = CUDAOCCUPANCYMAXACTIVEBLOCKSPERMULTIPROCESSOR&
         &(gS,fphi_uind_sitecu1_core,PME_BLOCK_DIM,0)
      if (ierr.ne.cudaSuccess)&
         &print 200, ierr,cudageterrorstring(ierr)
      gS   = gS*nSMP
   else if (kernelname.eq."fphi_mpole_core") then
      ierr = CUDAOCCUPANCYMAXACTIVEBLOCKSPERMULTIPROCESSOR&
         &(gS,fphi_mpole_core,PME_BLOCK_DIM,0)
      if (ierr.ne.cudaSuccess)&
         &print 200, ierr,cudageterrorstring(ierr)
      gS   = gS*nSMP
   else if (kernelname.eq."fphi_chg_site_kcu") then
      ierr = CUDAOCCUPANCYMAXACTIVEBLOCKSPERMULTIPROCESSOR&
         &(gS,fphi_chg_site_kcu,PME_BLOCK_DIM,0)
      if (ierr.ne.cudaSuccess)&
         &print 200, ierr,cudageterrorstring(ierr)
      gS   = gS*nSMP
   else if (kernelname.eq."grid_mpole_core") then
      ierr = CUDAOCCUPANCYMAXACTIVEBLOCKSPERMULTIPROCESSOR&
         &(gS,grid_mpole_sitecu_core_1p,PME_BLOCK_DIM1,0)
      if (ierr.ne.cudaSuccess)&
         &print 200, ierr,cudageterrorstring(ierr)
      gS   = gS*nSMP
   else if (kernelname.eq."grid_put_site_kcu1") then
      ierr = CUDAOCCUPANCYMAXACTIVEBLOCKSPERMULTIPROCESSOR&
         &(gS,grid_put_site_kcu1,PME_BLOCK_DIM1,0)
      if (ierr.ne.cudaSuccess)&
         &print 200, ierr,cudageterrorstring(ierr)
      gS   = gS*nSMP
   else if (kernelname.eq."grid_calc_frc_kcu") then
      ierr = CUDAOCCUPANCYMAXACTIVEBLOCKSPERMULTIPROCESSOR&
         &(gS,grid_calc_frc_kcu,PME_BLOCK_DIM1,0)
      if (ierr.ne.cudaSuccess)&
         &print 200, ierr,cudageterrorstring(ierr)
      gS   = gS*nSMP
   else if (kernelname.eq."grid_calc_frc_kcu1") then
      ierr = CUDAOCCUPANCYMAXACTIVEBLOCKSPERMULTIPROCESSOR&
         &(gS,grid_calc_frc_kcu1,PME_BLOCK_DIM1,0)
      if (ierr.ne.cudaSuccess)&
         &print 200, ierr,cudageterrorstring(ierr)
      gS   = gS*nSMP
   else if (kernelname.eq."grid_uind_core") then
      ierr = CUDAOCCUPANCYMAXACTIVEBLOCKSPERMULTIPROCESSOR&
         &(gS,grid_uind_sitecu_core_1p,PME_BLOCK_DIM1,0)
      if (ierr.ne.cudaSuccess)&
         &print 200, ierr,cudageterrorstring(ierr)
      gS   = gS*nSMP
   else if (kernelname.eq."ecreal1_kcu") then
      ierr = CUDAOCCUPANCYMAXACTIVEBLOCKSPERMULTIPROCESSOR&
         &(gS,ecreal1_kcu,BLOCK_DIM,0)
      if (ierr.ne.cudaSuccess)&
         &print 200, ierr,cudageterrorstring(ierr)
      gS   = gS*nSMP
   else if (kernelname.eq."ecreal3_kcu") then
      ierr = CUDAOCCUPANCYMAXACTIVEBLOCKSPERMULTIPROCESSOR&
         &(gS,ecreal3_kcu,BLOCK_DIM,0)
      if (ierr.ne.cudaSuccess)&
         &print 200, ierr,cudageterrorstring(ierr)
      gS   = gS*nSMP
   else if (kernelname.eq."echg_lj1_kcu_v0") then
      ierr = CUDAOCCUPANCYMAXACTIVEBLOCKSPERMULTIPROCESSOR&
         &(gS,echg_lj1_kcu_v0,BLOCK_DIM,0)
      if (ierr.ne.cudaSuccess)&
         &print 200, ierr,cudageterrorstring(ierr)
      gS   = gS*nSMP
   else if (kernelname.eq."emreal1_kcu") then
      ierr = CUDAOCCUPANCYMAXACTIVEBLOCKSPERMULTIPROCESSOR&
         &(gS,emreal1_kcu,BLOCK_DIM,0)
      if (ierr.ne.cudaSuccess)&
         &print 200, ierr,cudageterrorstring(ierr)
      gS   = gS*nSMP
   else if (kernelname.eq."emreal3_kcu") then
      ierr = CUDAOCCUPANCYMAXACTIVEBLOCKSPERMULTIPROCESSOR&
         &(gS,emreal3_kcu,BLOCK_DIM,0)
      if (ierr.ne.cudaSuccess)&
         &print 200, ierr,cudageterrorstring(ierr)
      gS   = gS*nSMP
   else if (kernelname.eq."epreal1c_core_cu") then
      ierr = CUDAOCCUPANCYMAXACTIVEBLOCKSPERMULTIPROCESSOR&
         &(gS,epreal1c_core_cu,BLOCK_DIM,0)
      if (ierr.ne.cudaSuccess)&
         &print 200, ierr,cudageterrorstring(ierr)
      gS   = gS*nSMP
   else if (kernelname.eq."mpreal1c_core_cu") then
      ierr = CUDAOCCUPANCYMAXACTIVEBLOCKSPERMULTIPROCESSOR&
         &(gS,mpreal1c_core_cu,BLOCK_DIM,0)
      if (ierr.ne.cudaSuccess)&
         &print 200, ierr,cudageterrorstring(ierr)
      gS   = gS*nSMP
   else if (kernelname.eq."epreal3_cu") then
      ierr = CUDAOCCUPANCYMAXACTIVEBLOCKSPERMULTIPROCESSOR&
         &(gS,epreal3_cu,BLOCK_DIM,0)
      if (ierr.ne.cudaSuccess)&
         &print 200, ierr,cudageterrorstring(ierr)
      gS   = gS*nSMP
   else if (kernelname.eq."tmatxb_pme_core_cu") then
      ierr = CUDAOCCUPANCYMAXACTIVEBLOCKSPERMULTIPROCESSOR&
         &(gS,tmatxb_pme_core_cu,BLOCK_DIM,0)
      if (ierr.ne.cudaSuccess)&
         &print 200, ierr,cudageterrorstring(ierr)
      gS   = gS*nSMP
   else if (kernelname.eq."otf_dc_tmatxb_pme_core_cu") then
      ierr = CUDAOCCUPANCYMAXACTIVEBLOCKSPERMULTIPROCESSOR&
         &(gS,otfdc_tmatxb_pme_core_cu,BLOCK_DIM,0)
      if (ierr.ne.cudaSuccess)&
         &print 200, ierr,cudageterrorstring(ierr)
      gS   = gS*nSMP
   else
      print 201, kernelname
      gS = 0
   end if

end

subroutine attach_pmecu_pointer(config)
   use atoms    ,only: x,y,z
   use atmlst   ,only: polerecglob
   use boxes    ,only: recip
   use domdec
   use fft      ,only: kstart1,kend1,istart1,iend1,jstart1,jend1
   use pme      ,only: igrid,qgrid2in_2d,qgridin_2d&
      &,thetai1,thetai2,thetai3
   use pmestuffcu
   implicit none
   integer,intent(in)::config
   enum,bind(C)
      enumerator pmeD,thetaD
   end enum

   select case (config)
    case (pmeD)
!$acc host_data use_device(thetai1,thetai2,thetai3,igrid &
!$acc  ,x,y,z,qgrid2in_2d,qgridin_2d,kstart1,kend1,jstart1,jend1 &
!$acc  ,istart1,iend1,prec_send)
      qgridin_t   => qgridin_2d
      qgrid2in_t  => qgrid2in_2d
      kstart1_t   => kstart1
      kend1_t     => kend1
      istart1_t   => istart1
      iend1_t     => iend1
      jstart1_t   => jstart1
      jend1_t     => jend1
      prec_send_t => prec_send
!$acc end host_data
    case (thetaD)
!$acc host_data use_device(thetai1,thetai2,thetai3)
      thetai1_t => thetai1
      thetai2_t => thetai2
      thetai3_t => thetai3
!$acc end host_data
    case default
      print*, 'Unknown configuration for thetai'
   end select
end subroutine

subroutine set_PME_texture
   use atoms    ,only: x,y,z
   use atmlst   ,only: polerecglob
   use boxes    ,only: recip
   use domdec
   use fft      ,only: kstart1,kend1,istart1,iend1,jstart1,jend1
   use pme      ,only: igrid,qgrid2in_2d,qgridin_2d&
      &,thetai1,thetai2,thetai3
   use pmestuffcu
   implicit none
   logical,save :: first_in=.true.


   if (first_in) then
!$acc host_data use_device(thetai1,thetai2,thetai3,igrid &
!$acc  ,x,y,z,qgrid2in_2d,qgridin_2d,kstart1,kend1,jstart1,jend1 &
!$acc  ,istart1,iend1,prec_send)
      x_t         => x
      y_t         => y
      z_t         => z
      qgridin_t   => qgridin_2d
      qgrid2in_t  => qgrid2in_2d
      kstart1_t   => kstart1
      kend1_t     => kend1
      istart1_t   => istart1
      iend1_t     => iend1
      jstart1_t   => jstart1
      jend1_t     => jend1
      prec_send_t => prec_send
!$acc end host_data
   end if

   if (nproc>1.or.first_in) then
!$acc host_data use_device(thetai1,thetai2,thetai3)
      thetai1_t => thetai1
      thetai2_t => thetai2
      thetai3_t => thetai3
!$acc end host_data
      first_in = .false.
   end if

end subroutine
#endif

subroutine set_pme_eps
#if (defined(MIXED)||defined(SINGLE))
   use pme
#ifdef _CUDA
   use pmestuffcu,only:setcu_pme_eps
#endif
   use tinheader ,only:prec1_eps
   implicit none
   pme_eps = 2*max(max(nfft1,nfft2),nfft3)*prec1_eps
# ifdef _CUDA
   call setcu_pme_eps(pme_eps)
# endif
#endif
end
!
!     "cmp_to_fmp_site" transforms the ith atomic multipole from Cartesian
!     to fractional coordinates
!
!
subroutine cmp_to_fmp_sitegpu(cmp,fmp)
   use mpole
   use sizes
   use utilgpu   ,only:rec_queue,ctf
   use tinheader ,only:ti_p
   implicit none
   integer i,j,k
   real(t_p),intent(in) ::cmp(10,max(npolerecloc,1))
   real(t_p),intent(out)::fmp(10,max(npolerecloc,1))
   logical,save:: init=.true.
!     real(t_p) time0,time1
!
!     find the matrix to convert Cartesian to fractional
!
   if (init) then
      call cart_to_fracgpu
      init=.false.
   end if

!$acc parallel loop collapse(2) async(rec_queue) &
!$acc         present(cmp,fmp,ctf)
   do i= 1,npolerecloc
      do j = 1, 10
!
!        apply the transformation to get the fractional multipoles
!
         if (j.eq.1) then
            fmp(1,i) = ctf(1,1) * cmp(1,i)
         else if (j.lt.5) then
            fmp(j,i) = 0.0_ti_p
!$acc loop seq
            do k = 2, 4
               fmp(j,i) = fmp(j,i) + ctf(j,k)*cmp(k,i)
            end do
         else
            fmp(j,i) = 0.0_ti_p
!$acc loop seq
            do k = 5, 10
               fmp(j,i) = fmp(j,i) + ctf(j,k)*cmp(k,i)
            end do
         end if
      end do
   end do
end
!
!
subroutine fphi_to_cphi_sitegpu(fphi,cphi)
   use mpole
   use sizes
   use utilgpu   ,only:rec_queue,ftc
   use tinheader ,only:ti_p
   implicit none
   integer i,j,k
   real(t_p),intent(in )::fphi(20,max(npolerecloc,1))
   real(t_p),intent(out)::cphi(10,max(npolerecloc,1))
   logical,save::init=.true.
!
!     find the matrix to convert fractional to Cartesian
!
   if (init) then
      call frac_to_cartgpu
      init=.false.
   end if

!$acc parallel loop collapse(2) async(rec_queue) &
!$acc         present(fphi,cphi,ftc)
   do i= 1,npolerecloc
      do j = 1, 10
!
!        apply the transformation to get the Cartesian potential
!
         if (j.eq.1) then
            cphi(1,i) = ftc(1,1) * fphi(1,i)
         else if (j.lt.5) then
            cphi(j,i) = 0.0_ti_p
!$acc loop seq
            do k = 2, 4
               cphi(j,i) = cphi(j,i) + ftc(j,k)*fphi(k,i)
            end do
         else
            cphi(j,i) = 0.0_ti_p
!$acc loop seq
            do k = 5, 10
               cphi(j,i) = cphi(j,i) + ftc(j,k)*fphi(k,i)
            end do
         end if
      end do
   end do
end
