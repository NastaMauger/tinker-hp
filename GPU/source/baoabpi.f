c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ######################################################################
c     ##                                                                  ##
c     ##  subroutine baoabpi  --  BAOAB Langevin PIMD step                ##
c     ##                                                                  ##
c     ######################################################################
c
c
c     "baoabpi" performs a single PIMD time step
c     via the BAOAB recursion formula 
c
c     literature reference:
c
c     Efficient molecular dynamics using geodesic integration 
c     and solvent-solute splitting B. Leimkuhler and C. Matthews,
c     Proceedings of the Royal Society A, 472: 20160138, 2016
c
c
#include "tinker_precision.h"
      
      subroutine baoabpi(istep,dt)
      use abinitio
      use atomsMirror
      use atmtyp
      use beads
      use domdec
      use energi
      use inform
      use iounit
      use mdstuf
      use moldyn
      use timestat
      use mpi
      use spectra
      use deriv
      use utilgpu
      use random_mod
      use units
      use bath
      use langevin
      use cutoff
      use utils
      use potent
      use virial
      use polar
      use moldyn
      use utilbaoabpi
      use commstuffpi
      use spectra
      implicit none
      integer, intent(in) :: istep
      real(r_p), intent(in) :: dt
      integer :: ibead,i,j,k,iglob,ierr,kk
      real(8),save::time0=0d0,time1=0d0,time00=0d0
      real(8),save::time01=0d0, timebegin=0d0
      real(8),save::timefull=0d0
      logical :: skip_parameters_copy
      real(r_p) :: dt2,sqrtnu,factor1,factor2
      real(r_p) :: scale,eigx0,eigv0,a1,a2
      integer :: nnoise,iloc,iqtb
      real(r_p) :: gammak,nuratio     
      integer :: ilocbeg,ilocend,ibeadbeg
      logical :: only_long
      character*240 filename
      integer :: iqm, freeunit

      if(istep==1) then
!$acc enter data create(eslow,vir_ctr) async
        stepint=nalt
        stepfast=nalt2
      endif

      dt2=0.5d0*dt
      sqrtnu=sqrt(real(nbeads,r_p))
      

      if(PITIMER) call stopwatchpi(timebegin,.true.,.true.)

      !! ALLOCATE NOISE VECTOR !!
      call prmem_request(noise,3*nlocpi*nbeads,async=.true.)

      !! APPLY FIRST PART OF BAOAB INTEGRATION !! 
      call apply_B_PI(polymer,dt2)
      if(isobaric) then
        call plangevinpi(polymer,dt2,-1,'A')
        call plangevinpi(polymer,dt,-1,'O')
      endif
      call apply_AOA_PI(polymer,dt)
      if(isobaric) call plangevinpi(polymer,dt2,-1,'A')
      if(PITIMER) call stopwatchpi(timeaoa,.false.,.false.)

      call full_gradient_prepare_pi(polymer,polymer_ctr,istep,.false.)
      if(PITIMER) call stopwatchpi(timecom,.false.,.false.)

      !!! COMPUTE FORCES ON CENTROID!!!  
      if(centroid_recip .and. .not. aiMD) then
        call prmem_requestm(derivs_centroid,3,nbloc,async=.TRUE.)
        call load_bead(polymer,0)
        only_long = centroid_longrange .and.
     &      (.not. (contract.and.nbeads_ctr==1) )
        call compute_gradslow_centroid(eslow,derivs_centroid
     &            ,vir_ctr,only_long, polar_allbeads)
        if(PITIMER) call stopwatchpi(timecontract,.false.,.false.)
        call commforcesrespa1(derivs_centroid,2)
        if(PITIMER) call stopwatchpi(timecom,.true.,.false.)
        if (deb_Energy) call info_energy(rank)
        if (deb_Force)  then
          call info_forces(cNBond)
          if(ranktot==0) write(*,*) 'were Forces for centroid:'
        endif
        if (deb_Atom)   call info_minmax_pva
      endif

      if (.not. aiMD) then
        !!! COMPUTE FORCES ON BEADS!!!  
        call prmem_requestm(derivs,3,nbloc,nbeadsloc,async=.true.)
        call compute_grad_beads(epotpi,derivs,virpi,polymer
     &         ,.not. (centroid_longrange .or. contract)  !LONG 
     &         ,.not. contract ! INT
     &         , .TRUE.     !SHORT 
     &         ,polar_allbeads)
  
        if(PITIMER) call stopwatchpi(timegrad,.false.,.false.)
        if(.not.(centroid_recip .or. contract))
     &    call commforcesrecpi(derivs,derivsRec,nbeadsloc)
        call commforcesblocpi(derivs,nbeadsloc,.FALSE.)  
  
        ibeadbeg=polymer%ibead_beg(rank_polymer+1)
!$acc parallel loop collapse(3) async default(present)
        DO k=1,nbeadsloc; DO iloc=1,nloc ; DO j=1,3
          i=glob(iloc)
          polymer%forces(j,i,k+ibeadbeg-1) = -derivs(j,iloc,k)
        ENDDO ; ENDDO ; ENDDO
      endif



      if (aiMD) then 
!$acc update host(polymer%pos)
        call write_qm_inputs(polymer,istep,nbeadsloc,nloc)
        call launch_qm_software(nbeadsloc, nloc) 
        call get_gradient_from_qm(nbeadsloc,nloc)
!$acc wait
!$acc serial
        ibeadbeg=polymer%ibead_beg(rank_polymer+1)
        DO k=1,nbeadsloc  
          DO iloc=1,nloc 
            DO j=1,3
              i=glob(iloc)
              polymer%forces(j,i,k+ibeadbeg-1) = gradient_qm(i,j)
            ENDDO  
          ENDDO 
        ENDDO
!$acc end serial
!$acc update host(polymer%forces)
        call organized_qm_files
      endif

      write(*,*) 'IN BAOABPI'
        DO k=1,nbeadsloc; DO iloc=1,nloc ; DO j=1,3
          i=glob(iloc)
          write(*,*) polymer%forces(j,i,k+ibeadbeg-1) 
        ENDDO ; ENDDO ; ENDDO

cc
cc      !! ADD CENTROID CONTRIBUTION TO ENERGY AND FORCES !!
cc      if(centroid_recip .and. .not. aiMD) then 
cc        if(rank_polymer==0) then
cc!$acc serial async present(epotpi,eslow,vir_ctr,virpi)
cc          epotpi = epotpi + eslow
cc          do j=1,3; do i=1,3
cc            virpi(i,j) = virpi(i,j) + vir_ctr(i,j)
cc          enddo; enddo
cc!$acc end serial  
cc        endif  
cc!$acc parallel loop collapse(3) async default(present)
cc        DO k=1,nbeadsloc; DO iloc=1,nloc ; DO j=1,3
cc          i=glob(iloc)
cc          polymer%forces(j,i,k+ibeadbeg-1) = 
cc     &       polymer%forces(j,i,k+ibeadbeg-1)
cc     &                   -derivs_centroid(j,iloc) 
cc        ENDDO ; ENDDO ; ENDDO
cc      endif
cc
      call comm_for_normal_modes(polymer,polymer%forces)
      if(PITIMER) call stopwatchpi(timecom,.true.,.false.)
cc
cc      !!! COMPUTE FORCES ON  CONTRACTED BEADS!!!  
cc      if(contract .and. nbeads_ctr>1 .and. .not. aiMD) then
cc
cc        call prmem_requestm(derivs_ctr,3,nbloc
cc     &      ,nbeadsloc_ctr,async=.true.)
cc        call compute_grad_beads(eslow,derivs_ctr,vir_ctr,polymer_ctr
cc     &       ,.not.centroid_longrange, .TRUE., .FALSE.  ![LONG, INT  , SHORT]
cc     &       ,polar_allbeads)
cc!$acc serial async present(epotpi,eslow,vir_ctr,virpi)
cc        epotpi = epotpi + eslow
cc        do j=1,3; do i=1,3
cc          virpi(i,j) = virpi(i,j) + vir_ctr(i,j)
cc        enddo; enddo
cc!$acc end serial
cc      if(PITIMER) call stopwatchpi(timecontract,.false.,.false.)
cc      
cc        if(.not.centroid_recip)
cc     &     call commforcesrecpi(derivs_ctr,derivsRec,nbeadsloc_ctr)
cc        call commforcesblocpi(derivs_ctr,nbeadsloc_ctr,.FALSE.)
cc        call commforces_ctr(polymer_ctr,derivs_ctr)
cc        call project_forces_ctr(polymer,polymer_ctr)
cc        if(PITIMER) call stopwatchpi(timecom,.true.,.false.)
cc      endif

cc      if (aiMD) then 
cc        if(register_coord) then
cc!$acc update host(polymer%pos)
cc          do k=1,nbeadsloc  
cc            iqm = freeunit()
cc            write(filename, '(A,I0.3,A)') 'coordinates_beads'
cc     $              ,k ,'.xyz'
cc            open(iqm, file=filename, status='unknown',
cc     $            action='write', position='append')
cc            write(iqm,'(A,I0)') 'Coordinates for istep = ', istep
cc            do iloc=1,nloc 
cc              write(iqm, '(F18.12,1X, F18.12,1X,F18.12)') 
cc     $                             polymer%pos(1,iloc,k) 
cc     $                            ,polymer%pos(2,iloc,k) 
cc     $                            ,polymer%pos(3,iloc,k) 
cc            enddo
cc          close(iqm)
cc          enddo
cc        endif
cc        call write_qm_inputs(polymer,istep,nbeadsloc,nloc)
cc        call launch_qm_software(nbeadsloc, nloc) 
cc        call get_gradient_from_qm(nbeadsloc,nloc)
cc
cc
cc!$acc wait
cc!$acc update self(polymer%forces)
cc        ibeadbeg=polymer%ibead_beg(rank_polymer+1)
cc        DO k=1,nbeadsloc  
cc          DO iloc=1,nloc 
cc            DO j=1,3
cc              i=glob(iloc)
cc              polymer%forces(j,i,k+ibeadbeg-1) = gradient_qm(i,j)
cc            ENDDO  
cc          ENDDO 
cc!$acc update device(polymer%forces)
cc        ENDDO
cc      endif

      call set_eigforces_pi(polymer,polymer%forces)
      call apply_B_PI(polymer,dt2)
      if(PITIMER) call stopwatchpi(timeaoa,.false.,.false.)

      if (calc_e.or.use_virial) then
        call kinetic_pi(polymer)
        call reduce_observables_pi(polymer)
      endif
      if(ir) then
!$acc update host(dippi,dipindpi) async
!$acc wait
       call save_dipole_traj(dippi,dipindpi) 
      endif
      if(isobaric) then
        call plangevinpi(polymer,dt,-1,'B')
        call ddpme3dnpt(1.0_re_p,istep)
      endif

      call mdstatpi(istep,dt)
      call mdsavebeads (istep,dt,polymer)
      if(PITIMER) call stopwatchpi(timeobs,.false.,.false.)

      if(PITIMER) then
        timefull=timefull+mpi_wtime()-timebegin

        if (mod(istep,iprint).eq.0) then
          if(verbose .and. ranktot.eq.0) then
              write(*,*)
              write(*,*) "### PIMD TIMERS ###"
c              write(*,'(A,f10.7,A,f7.3,A)') '      Time prepare: '     
c     &    , timereass/iprint , ' ( ',100*timereass/timefull,' %)'
c              write(*,'(A,f10.7,A,f7.3,A)') ' Time normal modes: '
c     &         , timepush/iprint, ' ( ',100*timepush/timefull,' %)'
              write(*,'(A,f10.7,A,f7.3,A)') 'Time grad contract: '
     & ,timecontract/iprint , ' ( ',100*timecontract/timefull,' %)'
              write(*,'(A,f10.7,A,f7.3,A)') '   Time grad beads: '
     &         , timegrad/iprint, ' ( ',100*timegrad/timefull,' %)'
              write(*,'(A,f10.7,A,f7.3,A)') '  Time observables: '
     &         , timeobs/iprint, ' ( ',100*timeobs/timefull,' %)'
              write(*,'(A,f10.7,A,f7.3,A)') '        Time comms: '
     &         ,timecom/iprint, ' ( ',100*timecom/timefull,' %)'
              write(*,'(A,f10.7,A,f7.3,A)') '        Time baoab: '
     &         , timeaoa/iprint, ' ( ',100*timeaoa/timefull,' %)'
              write(*,'(A,f10.7)')          '         Time step: '
     &         , timefull/iprint
              write(*,*) "###################"
              write(*,*)
          endif
          timefull=0.d0
          timepush=0d0
          timereass=0d0
          timeinit=0d0
          timecom=0.d0
          timeaoa=0.d0
          timegrad=0.d0
          timeobs=0.d0
          timecontract=0.d0
          timesavenl=0.d0
        endif
      endif

      end subroutine baoabpi


