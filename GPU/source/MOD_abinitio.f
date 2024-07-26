
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
      logical qm_init_file_found
      logical qm_input_file_found
      logical qm_grad_file_found
      real(r_p), allocatable :: gradient_qm_t(:,:)
      real :: energy_qm
      logical :: orca_qm, g16_qm, psi4_qm, pyscf_qm, qchem_qm

      integer :: nprocs

      integer :: compteur_aimd = 0

      character*40 filename_orca
      character*40 filename_g16 
      character*40 filename_psi4
      character*40 filename_pyscf
      character*40 filename_qchem
      character*40 qm_init_filename
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


      subroutine qm_filename()
      implicit none
      
      write(compteur_aimd_str, '(I10)') compteur_aimd

      if (orca_qm == .true.) then
        qm_input_filename = 'orca_'
     $        // trim(adjustl(compteur_aimd_str)) // '.inp'
        if (compteur_aimd .eq. 0) then
          qm_init_filename = qm_input_filename
        endif
        qm_output_filename = 'orca_'
     $        // trim(adjustl(compteur_aimd_str)) // '.out'
        qm_grad_filename = 'orca_'
     $      // trim(adjustl(compteur_aimd_str)) // '.engrad'
      endif

      end subroutine qm_filename
      


      subroutine write_qm_inputs() 
      use atomsMirror
      use atmtyp
      use beads
      use domdec
      implicit none
      integer :: i,k,iglob
      integer :: ios, iostat
      integer :: iqm, iqm_xyz, freeunit
      logical :: copy_file
      character*400 line
      character*240 filename
      character*240 filename_xyz

      write(compteur_aimd_str, '(I10)') compteur_aimd

      inquire(file=qm_init_filename, exist=qm_init_file_found)
      if(.NOT. qm_init_file_found) then
        write(*,*) 'NEED INITAL QM INPUT FILE'
      endif
      if (qm_init_file_found  .and. orca_qm) then
        copy_file=.true.
        iqm = freeunit()
        open(unit=415, file=qm_init_filename, action='read'
     $              ,status='old')
        write(qm_input_filename, '(A,I0.3,A,I0.3,A)') 'orca_'
     $        // trim(adjustl(compteur_aimd_str)) // '.inp'
        open(iqm, file=qm_input_filename, status='unknown',
     $            action='write')
        if(register_coord) then
          iqm_xyz = freeunit()
           write(filename_xyz, '(A,I0.3,A)') 'coordinates.xyz'
           open(iqm_xyz, file=filename_xyz, status='unknown',
     $            action='write', position='append')
           write(iqm_xyz,'(I0)') n
           write(iqm_xyz,'(A,I0)') 'Coordinates for istep = ', 
     $                compteur_aimd
        endif
          
        do
          read(415,'(A)', iostat=ios) line
          if (ios /= 0) exit

          if(index(line, '*xyz') /=0 ) then
            write(iqm,'(A)') trim(line)
            copy_file = .false.
            do i=1,n
              iglob=glob(i)
              write(iqm, '(A,F18.12,1X, F18.12,1X,F18.12)') 
     $                 name(iglob)
     $                 ,x(iglob), y(iglob), z(iglob)
              if(register_coord) then
                write(iqm_xyz, '(A,F18.12,1X, F18.12,1X,F18.12)') 
     $                 name(iglob)
     $                 ,x(iglob), y(iglob), z(iglob)
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
      endif
        
c      if (qm_input_file_found  .and.
c     $        g16_qm) then
c        open(unit=415, file=qm_input_filename, action='read')
c        open(unit=416 + compteur_aimd, file='g16_1_beads.inp'
c     $                , action='write', status='replace')
c        do
c          read(415,'(A)', iostat=ios) line
c          if (ios /= 0) exit
c
c          if(index(line, ',') /=0 ) then
c            write(416 + compteur_aimd ,'(A)') trim(line)
c            copy_file = .false.
c          endif
c  
c          if(copy_file) then
c            write(416 + compteur_aimd ,'(A)') trim(line)
c          endif
c
c        enddo
c        close(415)
c        close(416 + compteur_aimd)
c      endif
c        
c      if (qm_input_file_found  .and.
c     $        psi4_qm) then
c        open(unit=415, file=qm_input_filename, action='read')
c        open(unit=416 + compteur_aimd, file='psi4_1_beads.inp'
c     $                , action='write', status='replace')
c        do
c          read(415,'(A)', iostat=ios) line
c          if (ios /= 0) exit
c
c          if(index(line, 'molecule') /=0 ) then
c            write(416 + compteur_aimd ,'(A)') trim(line)
c            copy_file = .false.
c          endif
c  
c          if(copy_file) then
c            write(416 + compteur_aimd ,'(A)') trim(line)
c          endif
c
c          if (index(line, '}') /= 0) then
c            copy_file = .true.
c          endif
c
c        enddo
c        close(415)
c        close(416 + compteur_aimd)
c      endif
c    
c        
c      if (qm_input_file_found  .and.
c     $        pyscf_qm) then
c        open(unit=415, file=qm_input_filename, action='read')
c        open(unit=416 + compteur_aimd, file='pyscf_1_beads.inp'
c     $                , action='write', status='replace')
c        write(*,*) 'NEED QM WRITER FOR PYSCF - NOT IMPLEMENTED YET'
c        call fatal
c        do
c          read(415,'(A)', iostat=ios) line
c          if (ios /= 0) exit
c
c          if(index(line, 'molecule') /=0 ) then
c            write(416 + compteur_aimd ,'(A)') trim(line)
c            copy_file = .false.
c          endif
c  
c          if(copy_file) then
c            write(416 + compteur_aimd ,'(A)') trim(line)
c          endif
c
c          if (index(line, '$end') /= 0) then
c            copy_file = .true.
c          endif
c
c        enddo
c        close(415)
c        close(416 + compteur_aimd)
c      endif
c
c      if (qm_input_file_found  .and.
c     $        qchem_qm) then
c        open(unit=415, file=qm_input_filename, action='read')
c        open(unit=416 + compteur_aimd, file='qchem_1_beads.inp'
c     $                , action='write', status='replace')
c        do
c          read(415,'(A)', iostat=ios) line
c          if (ios /= 0) exit
c
c          if(index(line, 'molecule') /=0 ) then
c            write(416 + compteur_aimd ,'(A)') trim(line)
c            copy_file = .false.
c          endif
c  
c          if(copy_file) then
c            write(416 + compteur_aimd ,'(A)') trim(line)
c          endif
c
c          if (index(line, '$end') /= 0) then
c            copy_file = .true.
c          endif
c
c        enddo
c        close(415)
c        close(416 + compteur_aimd)
c      endif

      end subroutine write_qm_inputs


      subroutine launch_qm_software()
      implicit none
      character*240  command

      inquire(file=qm_input_filename, exist=qm_input_file_found)
      if(.not. qm_input_file_found) then
        write(*,*) 'NO QM INPUT FILE FOUND'
      endif

      if (orca_qm .and. qm_input_file_found) then
        command='$(which orca) ' // qm_input_filename // '> '
     $            // qm_output_filename 
      endif

      call execute_command_line(command)
      
      end subroutine launch_qm_software


      subroutine get_gradient_from_qm()
      use atomsMirror
      use atmtyp, only: mass
      use units
      use domdec

      implicit none
      integer :: i,j,iostat
cc      integer :: line_length, pos, pos_start, pos_end
      integer :: number_atoms
      logical :: found_energy, found_gradient, found_atoms
      integer iqm, freeunit
      real(r_p), allocatable :: gradient_qm(:,:)
      character*40 energy_str 
      character*400 line

      write(compteur_aimd_str, '(I10)') compteur_aimd

      found_energy = .false.
      found_atoms = .false.
      found_gradient = .false.

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
              if(number_atoms .eq. n) then
                found_atoms=.true.  
              endif

            else if (index(line, 'current gradient in') /= 0 .and. 
     $            .not. found_gradient) then
              if(allocated(gradient_qm)) then
                deallocate(gradient_qm)
              endif
              allocate(gradient_qm(number_atoms,3))
              read(iqm,'(A)', iostat=iostat) line
              do i=1, number_atoms
                read(iqm,'(A)', iostat=iostat) line
                read(line,*) gradient_qm(i,1) 
                gradient_qm(i,1) = (gradient_qm(i,1) * hartree / bohr)
                read(iqm,'(A)', iostat=iostat) line
                read(line,*) gradient_qm(i,2)
                gradient_qm(i,2) = (gradient_qm(i,2) * hartree / bohr)
                read(iqm,'(A)', iostat=iostat) line
                read(line,*) gradient_qm(i,3) 
                gradient_qm(i,3) = (gradient_qm(i,3) * hartree / bohr)
              enddo
              found_gradient=.true.

            endif
          enddo
          close(iqm)

          if(found_gradient) then
            if(allocated(gradient_qm_t)) then
              deallocate(gradient_qm_t)
            endif
            allocate(gradient_qm_t(3, n))
            do i=1,nloc
              gradient_qm_t(1,i)=gradient_qm(i,1)
              gradient_qm_t(2,i)=gradient_qm(i,2)
              gradient_qm_t(3,i)=gradient_qm(i,3)
            enddo
          endif
        

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
        write(*,*) 'NO QM OUTPUT GRADIENT FILE FOUND'
        call fatal()
      endif

cc      if (found_atoms) then
cc        write(*,*) 'NUMBER OF ATOMS = ', nloc
cc      endif
      if (found_gradient) then
cc        write(*,*) 'Gradient values in engrad:'
cc        do i = 1, n
cc          write(*,'(A,I3,3F16.9)') 'Atom', i, gradient_qm(i, 1), 
cc     $                gradient_qm(i, 2), gradient_qm(i, 3)
cc        enddo
        write(*,*) 'Gradient values in transpose:'
        do i=1,nloc
          do j=1,3
            print*,j,i,gradient_qm_t(j,i) 
          enddo
        enddo
        write(*,*) ' '
      endif

cc      write(*,*) qm_input_filename, qm_output_filename
cc      write(*,*) qm_grad_filename
cc      WRITE(*,*) qm_init_filename

      compteur_aimd = compteur_aimd + 1

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

