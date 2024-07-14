c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  module abintio   -- AIMD MODULE                            ##
c     ##                                                             ##
c     #################################################################
c
c
c    Flemme a faire plus tard 
c
c
#include "tinker_macro.h"
      module abinitio
      implicit none

      logical :: aiMD=.false.
      logical :: register_coord
      logical :: software_found
      logical qm_input_file_found
      logical qm_grad_file_found
      real(r_p), allocatable :: gradient_qm(:,:)
      real :: energy_qm
      logical :: orca_qm, g16_qm, psi4_qm, pyscf_qm, qchem_qm

      integer :: nprocs

      integer :: compteur_aimd = 0

      character*40 filename_orca
      character*40 filename_g16 
      character*40 filename_psi4
      character*40 filename_pyscf
      character*40 filename_qchem
      character*40 qm_input_filename
      character*40 qm_output_filename
      character*40 qm_grad_filename

      character*40 compteur_aimd_str
      character*1, allocatable :: atoms_name(:)


      contains

      subroutine read_aimd_keys()
      use domdec
      use keys
      implicit none
      character*10 value
      character*20  :: keyword
      character*240 :: record
      character*240 :: text
      character*240 :: string
      character*240 :: xyzfile
      integer :: next,i, ios



c
c     set default values for AIMD simulations
c

      orca_qm=.false.
      g16_qm=.false.
      psi4_qm=.false.
      pyscf_qm=.false.
      qchem_qm=.false.
      software_found=.false.


      nprocs          = 1         !Number of CPUs core that the QM software need
      register_coord  =.FALSE.

c
c     check for keywords containing any altered parameters
c

      do i = 1, nkey
        ios = 0
        next = 1
        record = keyline(i)
        call gettext (record,keyword,next)
        call upcase (keyword)
        string = record(next:240)

        select case (trim(keyword))
        case('REGISTER_XYZ')
          register_coord = .true.
        case('NPROCS')
          read (string,*,iostat=ios) nprocs 
        case('SOFTWARE')
          software_found = .true.
          call getword(record, value, next)
          call upcase(value)
          if (value .eq. 'ORCA') then
            orca_qm = .true.
            write(*, *) 'AIMD with ORCA software.'
          elseif (value .eq. 'GAUSSIAN') then
            g16_qm = .true.
            write(*, *) 'AIMD with GAUSSIAN software.'
          elseif (value .eq. 'PSI4') then
            psi4_qm = .true.
            write(*, *) 'AIMD with PSI4 software.'
          elseif (value .eq. 'PYSCF') then
            pyscf_qm = .true.
            write(*, *) 'AIMD with PYSCF software.'
          elseif (value .eq. 'QCHEM') then
            qchem_qm = .true.
            write(*, *) 'AIMD with QCHEM software.'
          else
            write(*, *) 'Unknown QM software: ', trim(value)
            call fatal
          endif
        end select

        if (ios /= 0) then
          write(*,*) "Warning: keyword ",trim(keyword)
     &         ," not correctly read!"
        endif

      end do

      if(.not. software_found) then
        write(*,*) 'Need a SOFTWARE keyword to run AIMD'
        call fatal
      endif 

      end subroutine read_aimd_keys


      subroutine write_qm_inputs(polymer,istep,nbeadsloc,nloc) 
      use beads
      implicit none
      type(POLYMER_COMM_TYPE), intent(inout) :: polymer
      integer :: nbeadsloc, iloc, nloc
      integer, intent(in) :: istep
      real(r_p), allocatable :: pos(:,:,:)
      integer :: i,k
      integer :: ios, iostat
      integer :: iqm, iqm_xyz, freeunit
      logical :: copy_file
      character*400 line
      character*240 filename
      character*240 filename_xyz
      character*4 :: atom_name

      compteur_aimd = istep

      write(compteur_aimd_str, '(I10)') compteur_aimd

      if (qm_input_file_found  .and.
     $        orca_qm) then
        do k=1, nbeadsloc
          copy_file=.true.
          iqm = freeunit()
          open(unit=415, file=qm_input_filename, action='read'
     $              ,status='old')
          write(filename, '(A,I0.3,A,I0.3,A)') 'orca_'
     $        // trim(adjustl(compteur_aimd_str)) // '_beads'
     $              ,k ,'.inp'
          open(iqm, file=filename, status='unknown',
     $            action='write')
          if(register_coord) then
            iqm_xyz = freeunit()
            write(filename_xyz, '(A,I0.3,A)') 'coordinates_beads'
     $              ,k ,'.xyz'
            open(iqm_xyz, file=filename_xyz, status='unknown',
     $      action='write', position='append')
            write(iqm_xyz,'(I0)') nloc
            write(iqm_xyz,'(A,I0)') 'Coordinates for istep = ', istep
          endif
          
          do
            read(415,'(A)', iostat=ios) line
            if (ios /= 0) exit

            if(index(line, '*xyz') /=0 ) then
              write(iqm,'(A)') trim(line)
              if(allocated(atoms_name)) then
                deallocate(atoms_name)
              endif
              allocate(atoms_name(nloc))
              do iloc=1,nloc
                read(415,'(A)', iostat=ios) line
                line = adjustl(line)
                read(line, '(A1)', iostat=ios) atoms_name(iloc)
              enddo
              copy_file = .false.
              do iloc=1,nloc
                write(iqm, '(A,F18.12,1X, F18.12,1X,F18.12)') 
     $                               atoms_name(iloc)
     $                              ,polymer%pos(1,iloc,k) 
     $                              ,polymer%pos(2,iloc,k) 
     $                              ,polymer%pos(3,iloc,k) 
                if(register_coord) then
                  write(iqm_xyz, '(A,F18.12,1X, F18.12,1X,F18.12)') 
     $                               atoms_name(iloc)
     $                              ,polymer%pos(1,iloc,k) 
     $                              ,polymer%pos(2,iloc,k) 
     $                              ,polymer%pos(3,iloc,k) 
                endif
              enddo
              write(iqm,'(A)') '*'
              close(415)
              close(iqm)
              if (register_coord) then
                close(iqm_xyz)
              endif
            endif
  
            if(copy_file) then
              write(iqm,'(A)') trim(line)
            endif
          enddo

        enddo
      endif
        
      if (qm_input_file_found  .and.
     $        g16_qm) then
        open(unit=415, file=qm_input_filename, action='read')
        open(unit=416 + compteur_aimd, file='g16_1_beads.inp'
     $                , action='write', status='replace')
        do
          read(415,'(A)', iostat=ios) line
          if (ios /= 0) exit

          if(index(line, ',') /=0 ) then
            write(416 + compteur_aimd ,'(A)') trim(line)
            copy_file = .false.
          endif
  
          if(copy_file) then
            write(416 + compteur_aimd ,'(A)') trim(line)
          endif

        enddo
        close(415)
        close(416 + compteur_aimd)
      endif
        
      if (qm_input_file_found  .and.
     $        psi4_qm) then
        open(unit=415, file=qm_input_filename, action='read')
        open(unit=416 + compteur_aimd, file='psi4_1_beads.inp'
     $                , action='write', status='replace')
        do
          read(415,'(A)', iostat=ios) line
          if (ios /= 0) exit

          if(index(line, 'molecule') /=0 ) then
            write(416 + compteur_aimd ,'(A)') trim(line)
            copy_file = .false.
          endif
  
          if(copy_file) then
            write(416 + compteur_aimd ,'(A)') trim(line)
          endif

          if (index(line, '}') /= 0) then
            copy_file = .true.
          endif

        enddo
        close(415)
        close(416 + compteur_aimd)
      endif
    
        
      if (qm_input_file_found  .and.
     $        pyscf_qm) then
        open(unit=415, file=qm_input_filename, action='read')
        open(unit=416 + compteur_aimd, file='pyscf_1_beads.inp'
     $                , action='write', status='replace')
        write(*,*) 'NEED QM WRITER FOR PYSCF - NOT IMPLEMENTED YET'
        call fatal
        do
          read(415,'(A)', iostat=ios) line
          if (ios /= 0) exit

          if(index(line, 'molecule') /=0 ) then
            write(416 + compteur_aimd ,'(A)') trim(line)
            copy_file = .false.
          endif
  
          if(copy_file) then
            write(416 + compteur_aimd ,'(A)') trim(line)
          endif

          if (index(line, '$end') /= 0) then
            copy_file = .true.
          endif

        enddo
        close(415)
        close(416 + compteur_aimd)
      endif

      if (qm_input_file_found  .and.
     $        qchem_qm) then
        open(unit=415, file=qm_input_filename, action='read')
        open(unit=416 + compteur_aimd, file='qchem_1_beads.inp'
     $                , action='write', status='replace')
        do
          read(415,'(A)', iostat=ios) line
          if (ios /= 0) exit

          if(index(line, 'molecule') /=0 ) then
            write(416 + compteur_aimd ,'(A)') trim(line)
            copy_file = .false.
          endif
  
          if(copy_file) then
            write(416 + compteur_aimd ,'(A)') trim(line)
          endif

          if (index(line, '$end') /= 0) then
            copy_file = .true.
          endif

        enddo
        close(415)
        close(416 + compteur_aimd)
      endif

      end subroutine write_qm_inputs


      subroutine launch_qm_software(nbeadsloc, nloc)
      implicit none
      integer :: k
      integer iqm, freeunit
      integer :: nbeadsloc, iloc, nloc
      character*3 :: k_str
      character*240  command

      write(compteur_aimd_str, '(I10)') compteur_aimd

      if (compteur_aimd == 0) then
          if (orca_qm == .true.) then
            qm_input_filename = 'orca_0.inp'
            qm_output_filename = 'orca_0.out'
          
          else if(g16_qm == .true.) then
            qm_input_filename = 'g16_0.inp'
            qm_output_filename = 'g16_0.out'

          else if(psi4_qm == .true.) then
            qm_input_filename = 'psi4_0.inp'
            qm_output_filename = 'psi4_0.out'

          else if(pyscf_qm == .true.) then
            qm_input_filename = 'pyscf_0.inp'
            qm_output_filename = 'pyscf_0.out'

          else if(qchem_qm == .true.) then
            qm_input_filename = 'qchem_0.inp'
            qm_output_filename = 'qchem_0.out'
          endif
      endif

      if(compteur_aimd  .gt. 0) then
        if (orca_qm == .true.) then
          do k=1, nbeadsloc
            write(k_str, '(I3.3)') k
            qm_input_filename = 'orca_'
     $        // trim(adjustl(compteur_aimd_str)) // '_beads'
     $        // trim(adjustl(k_str)) // '.inp'
            qm_output_filename = 'orca_'
     $        // trim(adjustl(compteur_aimd_str)) // '_beads'
     $        // trim(adjustl(k_str)) // '.out'
          enddo
        endif
      endif

      inquire(file=qm_input_filename, exist=qm_input_file_found)
      if(.NOT. qm_input_file_found) then
        write(*,*) 'NEED AN INITIAL QM INPUT FILE'
      endif

      if (orca_qm) then
        command='$(which orca) ' // qm_input_filename // '> '
     $            // qm_output_filename 
      endif

      call execute_command_line(command)
      
      end subroutine launch_qm_software


      subroutine get_gradient_from_qm(nbeadsloc, nloc)
      use units
      implicit none
      integer :: i,j,k,iostat
      integer :: line_length, pos, pos_start, pos_end
      integer :: number_atoms
      logical :: found_energy, found_gradient, found_atoms
      integer iqm, freeunit
      real(r_p), allocatable :: grad(:,:)
      integer :: nbeadsloc, iloc, nloc
      character*40 energy_str 
      character*3 :: k_str
      character*400 line
      character*40 method

      write(compteur_aimd_str, '(I10)') compteur_aimd

      found_energy = .false.
      found_atoms = .false.
      found_gradient = .false.

      if (compteur_aimd == 0 ) then
        if (orca_qm == .true.) then
          qm_grad_filename = 'orca_' 
     $           // trim(adjustl(compteur_aimd_str)) // '.engrad'
        
        else if(g16_qm == .true.) then
          qm_grad_filename = 'g16_' 
     $           // trim(adjustl(compteur_aimd_str)) // '.log'

        else if(psi4_qm == .true.) then
          qm_grad_filename = 'psi4_' 
     $           // trim(adjustl(compteur_aimd_str)) // '.inp.dat'

        else if(pyscf_qm == .true.) then
          qm_grad_filename = 'pyscf_' 
     $           // trim(adjustl(compteur_aimd_str)) // '.out'

        else if(qchem_qm == .true.) then
          qm_grad_filename = 'qchem_' 
     $           // trim(adjustl(compteur_aimd_str)) // '.out'
          
        endif
      endif

      if (compteur_aimd .gt. 0 ) then
        if (orca_qm == .true.) then
          do k=1, nbeadsloc
            write(k_str, '(I3.3)') k
c            iqm = freeunit()
            qm_grad_filename = 'orca_'
     $        // trim(adjustl(compteur_aimd_str)) // '_beads'
     $        // trim(adjustl(k_str)) // '.engrad'
          enddo
        endif
      endif

      iqm= freeunit()

      inquire(file=qm_grad_filename, exist=qm_grad_file_found)
      if(qm_grad_file_found) then 
        if(orca_qm) then
          open(iqm, file=qm_grad_filename, status='old')
          do 
            read(iqm, '(A)', iostat=iostat) line
            if (iostat /= 0) exit
            if (index(line, 'current total energy') /= 0 .and.
     $            .not. found_energy) then
              read(iqm,'(A)', iostat=iostat) line 
              read(iqm,'(A)', iostat=iostat) line 
              read(line,*) energy_qm 
              found_energy=.true.

            else if (index(line, 'Number of atoms') /= 0 .and. 
     $            .not. found_atoms) then
              read(iqm,'(A)', iostat=iostat) line
              read(iqm,'(A)', iostat=iostat) line
              read(line,*) number_atoms
              found_atoms=.true.  

            else if (index(line, 'current gradient in') /= 0 .and. 
     $            .not. found_gradient) then
              if(allocated(gradient_qm)) then
                deallocate(gradient_qm)
              endif
              allocate(gradient_qm(number_atoms,3))
              read(iqm,'(A)', iostat=iostat) line
              read(iqm,'(A)', iostat=iostat) line
              do i=1, number_atoms
                read(line,*) gradient_qm(i,1) 
                gradient_qm(i,1) = gradient_qm(i,1) * hartree / bohr
                read(iqm,'(A)', iostat=iostat) line
                read(line,*) gradient_qm(i,2)
                gradient_qm(i,2) = gradient_qm(i,2) * hartree / bohr
                read(iqm,'(A)', iostat=iostat) line
                read(line,*) gradient_qm(i,3) 
                gradient_qm(i,3) = gradient_qm(i,3) * hartree / bohr
              enddo
              found_gradient=.true.

            endif
          enddo
          close(iqm)

c         else if(g16_qm) then
c           open(426, file=qm_grad_filename, status='old')
c           do 
c             read(426, '(A)', iostat=iostat) line 
c             if (iostat /= 0) exit
c             if (index(line, '\MP2=') /= 0 .and. 
c      $                  .not. found_energy) then
c               method='MP2'
c             endif
c             if (index(line, '\MP2=') == 0 .and. 
c      $         index(line, '\HF=') /=0 .and.
c      $                  .not. found_energy) then
c               method='HF'
c             endif
c             if (index(line, '\\' // trim(method) // '=') /= 0 .and.
c      $                  .not. found_energy) then
c               pos_start = index(line, trim(method) // '=')
c               pos_start = pos_start + len(trim(method)) + 1 
c               line_length=len_trim(line)
c               do i = pos_start, line_length
c                 if (line(i:i) == '\\') then
c                   pos_end = i-1 
c                   write(*,*) pos_end
c                   exit
c                 endif
c               enddo
c               energy_str = line(pos_start:pos_end)
c               read(energy_str,*) energy_qm
c               found_energy=.true.
c 
c             else if (index(line, 'NAtoms=') /= 0 .and. 
c      $            .not. found_atoms) then
c               read(line(10:),*) number_atoms
c               found_atoms=.true.
c 
c             else if (index(line, 'Forces (Hartrees/Bohr)') /= 0 .and. 
c      $            .not. found_gradient) then
c               allocate(gradient_qm(number_atoms,3))
c               read(426,'(A)', iostat=iostat) line 
c               read(426,'(A)', iostat=iostat) line 
c               read(426,'(A)', iostat=iostat) line 
c               do i=1, number_atoms
c                 read(line,*) dummy, dummy, gradient_qm(i,1) 
c      $                  , gradient_qm(i,2), gradient_qm(i,3)
c                 read(426,'(A)', iostat=iostat) line
c                 gradient_qm(i,1) = gradient_qm(i,1) * hartree / bohr
c                 gradient_qm(i,2) = gradient_qm(i,2) * hartree / bohr
c                 gradient_qm(i,3) = gradient_qm(i,3) * hartree / bohr
c               enddo
c               found_gradient=.true. 
c 
c            endif
c           enddo
c         close(426)
c
c        else if(psi4_qm) then
c          open(427, file=qm_grad_filename, status='old')
c          do 
c            read(427, '(A)', iostat=iostat) line
c            if (iostat /= 0) exit
c
c            if ((index(line, '= energy(''MP2'')') /= 0 .or. 
c     $           index(line, '= energy(''mp2'')') /= 0) .and. 
c     $                  .not. found_energy) then
c              method='MP2'
c            endif
c            if (index(line, 'energy(''scf'')') /= 0 .and. 
c     $                  .not. found_energy) then
c              method='HF'
c            endif
c 
c            if (index(trim(line), 'MP2 Total Energy (a.u.)') /= 0 
c     $            .and. .not. found_energy) then
c                do i=1,6
c                  read(427,'(A)', iostat=iostat) line 
c                enddo
c                read(line(40:), *) energy_str 
c                read(energy_str,*) energy_qm
c                found_energy = .true.
c
c            else if (index(line, 'Number of atoms') /= 0 .and. 
c     $            .not. found_atoms) then
c               read(line(40:),*) number_atoms
c               found_atoms=.true.
c   
c            else if (index(line, '-Total gradient:') /= 0 .and. 
c     $            .not. found_gradient) then
c               allocate(gradient_qm(number_atoms,3))
c               read(427,'(A)', iostat=iostat) line 
c               read(427,'(A)', iostat=iostat) line 
c               read(427,'(A)', iostat=iostat) line 
c               do i=1, number_atoms
c                 read(line,*) dummy, gradient_qm(i,1) 
c     $                  , gradient_qm(i,2), gradient_qm(i,3)
c                 read(427,'(A)', iostat=iostat) line
c                 gradient_qm(i,1) = gradient_qm(i,1) * hartree / bohr
c                 gradient_qm(i,2) = gradient_qm(i,2) * hartree / bohr
c                 gradient_qm(i,3) = gradient_qm(i,3) * hartree / bohr
c               enddo
c               found_gradient=.true. 
c
c             endif
c           enddo
c         close(427)
c
c        else if(pyscf_qm) then
c          open(428, file=qm_grad_filename, status='old')
c          do 
c            read(428, '(A)', iostat=iostat) line
c            if (iostat /= 0) exit
c            if (index(line, 'E(' // trim(method) // ') =') /= 0 .and.
c     $                  .not. found_energy) then
c              read(line(10:),*) energy
c              found_energy=.true.
c
c            elseif (index(line,trim(method) // ' gradients') /= 0 .and. 
c     $                  .not. found_atoms) then
c              start_count = .true.
c              count_lines = 1
c              do while (start_count)
c                read(428, '(A)', iostat=iostat) line
c                  if (iostat /= 0 .or. index(line, '----') /= 0) then
c                      start_count = .false.
c                  else
c                      count_lines = count_lines + 1
c                  endif
c              enddo
c              number_atoms = count_lines - 2
c              found_atoms = .true.
c      
c            rewind(428)
c            elseif (index(line, 'gradients') /= 0 .and. 
c     $            .not. found_gradient) then
c              allocate(gradient_qm(number_atoms,3))
c              read(428,'(A)', iostat=iostat) line 
c              read(428,'(A)', iostat=iostat) line 
c              do i=1, number_atoms
c                read(line,*) dummy, dummy_char, gradient_qm(i,1) 
c     $                  , gradient_qm(i,2), gradient_qm(i,3)
c                read(428,'(A)', iostat=iostat) line
c                gradient_qm(i,1) = gradient_qm(i,1) * hartree / bohr
c                gradient_qm(i,2) = gradient_qm(i,2) * hartree / bohr
c                gradient_qm(i,3) = gradient_qm(i,3) * hartree / bohr
c              enddo
c              found_gradient=.true. 
c              
c            endif
c          enddo
c        close(428)
c
c        else if(qchem_qm) then
c          open(429, file=qm_grad_filename, status='old')
c          do 
c            read(429, '(A)', iostat=iostat) line
c            if (iostat /= 0) exit
c            if (index(line, 'total energy =') /= 0 .and.
c     $                  .not. found_energy) then
c              read(line(40:),*) energy
c              found_energy=.true.
c
c            elseif (index(line,'Atom')  /= 0 .and. 
c     $                  .not. found_atoms) then
c              read(429,'(A)', iostat=iostat) line
c              read(429,'(A)', iostat=iostat) line
c              start_count = .true.
c              count_lines = 1
c              do while (start_count)
c                read(429, '(A)', iostat=iostat) line 
c                  if (iostat /= 0 .or. index(line, '----') /= 0) then
c                      start_count = .false.
c                  else
c                      count_lines = count_lines + 1
c                  endif
c              enddo
c              number_atoms = count_lines 
c              found_atoms = .true.
c
c            write(*,*)'QCHEM GRADIENT READER NOT FIXED YET'  
c     $            //  '- CHANGE SOFTWARE'
c            call fatal
c!               number_atoms=11
c!               elseif (index(line, 'Full Analytical Gradient') /= 0 .and. 
c!        $            .not. found_gradient) then
c!                 read(429, '(A)', iostat=iostat) line
c!                 read(429, '(A)', iostat=iostat) line 
c!                 allocate(gradient_qm(number_atoms,3))
c!                 total = int(number_atoms/5)
c!                 compteur = 1
c!                 write(*,*) total 
c!   !
c!                 write(*,*) number_atoms, compteur
c!                 write(*,*) number_atoms/5
c!                 do compteur=1, total
c!                   do i = 1, 3
c!                     read(line,*) dummy
c!        $                , gradient_qm(1, i), gradient_qm(2,i)
c!        $                , gradient_qm(3, i), gradient_qm(4,i)
c!        $                , gradient_qm(5, i)                             
c!                     read(429, '(A)', iostat=iostat) line 
c!                     write(*,*) 
c!        $                  gradient_qm(1, i), gradient_qm(2,i)
c!        $                , gradient_qm(3, i), gradient_qm(4,i)
c!        $                , gradient_qm(5, i)                             
c!                   enddo
c!                   read(429, '(A)', iostat=iostat) line 
c!              enddo
c                
c
cC                read(line, *) dummy, 
cCc   $                    ((gradient_qm(compteur + j  , i), j = 1, 5))
c
cC              found_gradient=.true.
c            
c
c            endif
c          enddo
c        close(429)
c        
        endif
      else
        write(*,*) 'NO QM GRADIENT FILE FOUND'
        call fatal()
      endif

c      if (found_energy) then
c      endif
c      if (found_atoms) then
c        write(*,*) 'NUMBER OF ATOMS = ', number_atoms
c      endif
      if (found_gradient) then
           write(*,*) 'Gradient values:'
           do i = 1, number_atoms 
               write(*,'(A,I3,3F16.9)') 'Atom', i, gradient_qm(i, 1), 
     $                gradient_qm(i, 2), gradient_qm(i, 3)
           enddo
      endif

      end subroutine get_gradient_from_qm


      subroutine organized_qm_files
      implicit none
      integer :: i
      integer :: compteur_save
      character*240  command_1
      character*240  command_2
      character*240  command_2b
      character*240  command_3
      character*240  command_3b
      character*40 compteur_save_str

      if(compteur_aimd .ge. 2) then
        compteur_save = compteur_aimd - 1
      else if(compteur_aimd  == 1) then
        compteur_save = 0
      endif

      write(compteur_save_str, '(I10)') compteur_save

      command_1= 'mkdir -p  timestep_' // 
     $                trim(adjustl(compteur_save_str))
      call execute_command_line(command_1)

      if (compteur_aimd .ge. 2) then
        write(command_2, '(A, A, A)') 'mv orca_', 
     $      trim(adjustl(compteur_save_str)), '_beads*'
      elseif (compteur_aimd == 1) then
        write(command_2, '(A, A, A)') 'mv orca_', 
     $      trim(adjustl(compteur_save_str)), '*'
      endif

      command_3 = trim(adjustl(command_2)) // ' timestep_' //
     $      trim(adjustl(compteur_save_str))
      call execute_command_line(command_3)
      
      end subroutine organized_qm_files
      
      end module abinitio

