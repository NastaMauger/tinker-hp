
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
      logical :: create_qm_file
      logical :: register_coord
      logical :: software_found
      logical qm_input_file_found
      logical qm_grad_file_found
      logical :: orca_qm, g16_qm, psi4_qm, pyscf_qm, qchem_qm
      real(r_p), allocatable :: gradient_qm_t(:,:)
      real :: energy_qm

      integer :: nprocs

      character*40 qm_input_filename
      character*40 qm_output_filename
      character*40 qm_grad_filename
      character*40 qm_coord_filename

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
          elseif (value .eq. 'GAUSSIAN') then
            g16_qm = .true.
          elseif (value .eq. 'PSI4') then
            psi4_qm = .true.
          elseif (value .eq. 'PYSCF') then
            pyscf_qm = .true.
          elseif (value .eq. 'QCHEM') then
            qchem_qm = .true.
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


      subroutine qm_inputs_name_gen() 
      use inform
      use atomsMirror
      use atmtyp
      use beads
      use domdec
      implicit none
      integer :: k
      integer iqm, freeunit
      integer :: ios, iostat
      character*40 qm_init_filename
      character*400 line
      character*3 :: k_beads

      qm_input_filename  = 'orca.inp'
      qm_output_filename = 'orca.out'
      qm_grad_filename   = 'orca.engrad'
      if(register_coord) then
        qm_coord_filename = 'coordinates.xyz'
      endif

      inquire(file=qm_input_filename, exist=qm_input_file_found)
      if(.not. qm_input_file_found) then
        write(*,*) 'NO QM INPUT FILE FOUND'
        call fatal()
      endif

      if(app_id .eq. pimd_a) then  
        if(create_qm_file) then
c       This part of the code allows to have as many QM inputs as beads
c       before going in write_qm_inputs. Hence these inputs are the same
c       than the one given by the user at the beginning of the simulation
          qm_init_filename   = 'orca.inp'
          open(unit=415, file=qm_init_filename, action='read'
     $              ,status='old')
          do k=1,nbeads
            write(k_beads,'(I3.3)') k
            qm_input_filename  = 'orca_beads' 
     $              //trim(adjustl(k_beads)) //'.inp'
            iqm = freeunit()
            open(unit=iqm, file=qm_input_filename, status='unknown',
     $            action='write')
            rewind(415)
            do while(.true.)
              read(415,'(A)', iostat=ios) line
              if (ios /= 0) exit
              write(iqm, '(A)') trim(line)
            enddo
            close(iqm)
            qm_output_filename = 'orca_beads' 
     $              //trim(adjustl(k_beads)) //'.out'
            qm_grad_filename   = 'orca_beads' 
     $              //trim(adjustl(k_beads)) //'.engrad'
            if(register_coord) then
              qm_coord_filename = 'coordinates_beads' 
     $              //trim(adjustl(k_beads))//'.xyz'
            endif
            close(415)
            call write_qm_inputs
            go to 14
          enddo
        endif
      endif

      call write_qm_inputs

   14 continue
      
      end subroutine qm_inputs_name_gen

      subroutine write_qm_inputs() 
      use atomsMirror
      use atmtyp
      use domdec
      implicit none
      integer :: i,k,iglob
      integer :: ios, iostat
      integer :: iqm, iqm_xyz, freeunit
      integer :: compteur_aimd = 0
      logical :: copy_file
      character*400 line
      character*240 filename_temp

      filename_temp = 'orca.inp_temp'
      copy_file=.true.
      iqm = freeunit()
      open(unit=415, file=qm_input_filename, action='read'
     $              ,status='old')
      open(iqm, file=filename_temp, status='unknown',
     $            action='write')
        if(register_coord) then
          iqm_xyz = freeunit()
          open(iqm_xyz, file=qm_coord_filename, status='unknown',
     $            action='write', position='append')
          write(iqm_xyz,'(I0)') n
          write(iqm_xyz,'(A,I0)') 'Coordinates for istep = ', 
     $                compteur_aimd
          do i=1,n
            iglob=glob(i)
            write(iqm_xyz, '(A,F18.12,1X, F18.12,1X,F18.12)') 
     $           name(iglob),x(iglob), y(iglob), z(iglob)
            enddo
          compteur_aimd = compteur_aimd + 1
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
     $                 name(iglob),x(iglob), y(iglob), z(iglob)
            enddo
            write(iqm,'(A)') '*'
            close(415,status='delete')
            close(iqm)
            call rename(filename_temp,qm_input_filename)
          endif

          if(register_coord) then
            close(iqm_xyz)
          endif
  
          if(copy_file) then
            write(iqm,'(A)') trim(line)
          endif
        enddo


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

      call qm_inputs_name_gen()

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
        write(*,*) 'NO QM GRADIENT FILE FOUND'
        call fatal()
      endif

cc      if (found_atoms) then
cc        write(*,*) 'NUMBER OF ATOMS = ', nloc
cc      endif
cc      if (found_gradient) then
cc        write(*,*) 'Gradient values in engrad:'
cc        do i = 1, n
cc          write(*,'(A,I3,3F16.9)') 'Atom', i, gradient_qm(i, 1), 
cc     $                gradient_qm(i, 2), gradient_qm(i, 3)
cc        enddo
cc        write(*,*) 'Gradient values in transpose:'
cc        do i=1,nloc
cc          do j=1,3
cc            print*,j,i,gradient_qm_t(j,i) 
cc          enddo
cc        enddo
cc        write(*,*) ' '
cc      endif

      end subroutine get_gradient_from_qm

      end module abinitio

