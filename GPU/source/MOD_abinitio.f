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
      logical :: software_found
      logical qm_input_file_found
      logical qm_grad_file_found
      logical :: orca_qm, g16_qm, psi4_qm, pyscf_qm, qchem_qm
      integer :: multiplicity
      integer :: charge

      integer :: nprocs

      integer :: compteur_aimd = 0

      character*40 filename_orca
      character*40 filename_g16 
      character*40 filename_psi4
      character*40 filename_pyscf
      character*40 filename_qchem
      character*40 qm_input_filename
      character*40 qm_grad_filename

      character*40 compteur_aimd_str


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


      nprocs        = 6         !Number of CPUs core that the QM software need

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


      subroutine write_qm_inputs() 
      implicit none
      integer :: ios, iostat
      logical :: copy_file
      character*400 line

      write(compteur_aimd_str, '(I10)') compteur_aimd


      if (compteur_aimd == 0) then
          if (orca_qm == .true.) then
            qm_input_filename = 'orca_' 
     $             // trim(adjustl(compteur_aimd_str)) // '.inp'
          
          else if(g16_qm == .true.) then
            qm_input_filename = 'g16_' 
     $             // trim(adjustl(compteur_aimd_str)) // '.inp'

          else if(psi4_qm == .true.) then
            qm_input_filename = 'psi4_' 
     $             // trim(adjustl(compteur_aimd_str)) // '.inp'

          else if(pyscf_qm == .true.) then
            qm_input_filename = 'pyscf_' 
     $             // trim(adjustl(compteur_aimd_str)) // '.inp'

          else if(qchem_qm == .true.) then
            qm_input_filename = 'qchem_' 
     $             // trim(adjustl(compteur_aimd_str)) // '.inp'
            
          endif
        inquire(file=qm_input_filename, exist=qm_input_file_found)
        if(.NOT. qm_input_file_found) then
          write(*,*) 'NEED A QM INPUT FILE'
          call fatal
        endif
      endif

      copy_file=.true.

      if (qm_input_file_found  .and.
     $        orca_qm) then
        open(unit=415, file=qm_input_filename, action='read')
        open(unit=416 + compteur_aimd, file='orca_1_beads.inp'
     $                , action='write', status='replace')
        do
          read(415,'(A)', iostat=ios) line
          if (ios /= 0) exit

          if(index(line, '*xyz') /=0 ) then
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
 
      call launch_qm_software
      call get_gradient_from_qm

      end subroutine write_qm_inputs


      subroutine launch_qm_software
      implicit none
      character*240  command


!!      if (compteur_aimd == 0) then
!!       if (orca_qm) then
!!        write(*,*) 'ORCA FILENAME = ', filename_orca 
!!               command = "$(which orca) filename_orca > orca_0.out"
!!               call system(command)
!!        endif
!!      endif
              

      end subroutine launch_qm_software


      subroutine get_gradient_from_qm
      use units
      implicit none
      integer :: i,j,k,iostat
      integer :: number_atoms
      integer :: dummy, count_lines, compteur, remainder,total
      integer :: pos, pos_start, pos_end
      logical :: start_count
      logical :: found_energy, found_gradient, found_atoms
      real :: energy
      integer :: line_length
      real(r_p), allocatable :: gradient_qm(:,:)
      real(r_p), allocatable :: grad(:,:)
      character*40 energy_str 
      character*400 line
      character*400 line_1, line_2, line_3
      character *1 dummy_char
      character*40 method

      found_energy = .false.
      found_atoms = .false.
      found_gradient = .false.


      if (orca_qm == .true.) then
        qm_grad_filename = 'orca_' 
     $         // trim(adjustl(compteur_aimd_str)) // '.engrad'
      
      else if(g16_qm == .true.) then
        qm_grad_filename = 'g16_' 
     $         // trim(adjustl(compteur_aimd_str)) // '.log'

      else if(psi4_qm == .true.) then
        qm_grad_filename = 'psi4_' 
     $         // trim(adjustl(compteur_aimd_str)) // '.inp.dat'

      else if(pyscf_qm == .true.) then
        qm_grad_filename = 'pyscf_' 
     $         // trim(adjustl(compteur_aimd_str)) // '.out'

      else if(qchem_qm == .true.) then
        qm_grad_filename = 'qchem_' 
     $         // trim(adjustl(compteur_aimd_str)) // '.out'
        
      endif



      inquire(file=qm_grad_filename, exist=qm_grad_file_found)
      if(qm_grad_file_found) then 
        write(*,*) 'READING GRADIENT FROM ', qm_grad_filename
        if(orca_qm) then
          open(425, file=qm_grad_filename, status='old')
          do 
            read(425, '(A)', iostat=iostat) line
            if (iostat /= 0) exit
            if (index(line, 'current total energy') /= 0 .and.
     $            .not. found_energy) then
              read(425,'(A)', iostat=iostat) line 
              read(425,'(A)', iostat=iostat) line 
              read(line,*) energy   
              found_energy=.true.

            else if (index(line, 'Number of atoms') /= 0 .and. 
     $            .not. found_atoms) then
              read(425,'(A)', iostat=iostat) line
              read(425,'(A)', iostat=iostat) line
              read(line,*) number_atoms
              found_atoms=.true.  

            else if (index(line, 'current gradient in') /= 0 .and. 
     $            .not. found_gradient) then
              allocate(gradient_qm(number_atoms,3))
              read(425,'(A)', iostat=iostat) line
              read(425,'(A)', iostat=iostat) line
              do i=1, number_atoms
                read(line,*) gradient_qm(i,1) 
                gradient_qm(i,1) = gradient_qm(i,1) * hartree / bohr
                read(425,'(A)', iostat=iostat) line
                read(line,*) gradient_qm(i,2)
                gradient_qm(i,2) = gradient_qm(i,2) * hartree / bohr
                read(425,'(A)', iostat=iostat) line
                read(line,*) gradient_qm(i,3) 
                gradient_qm(i,3) = gradient_qm(i,3) * hartree / bohr
              enddo
              found_gradient=.true.

            endif
          enddo
          close(425)

        else if(g16_qm) then
          open(426, file=qm_grad_filename, status='old')
          do 
            read(426, '(A)', iostat=iostat) line 
            if (iostat /= 0) exit
            if (index(line, '\MP2=') /= 0 .and. 
     $                  .not. found_energy) then
              method='MP2'
            endif
            if (index(line, '\MP2=') == 0 .and. 
     $         index(line, '\HF=') /=0 .and.
     $                  .not. found_energy) then
              method='HF'
            endif
            if (index(line, '\\' // trim(method) // '=') /= 0 .and.
     $                  .not. found_energy) then
              pos_start = index(line, trim(method) // '=')
              pos_start = pos_start + len(trim(method)) + 1 
              line_length=len_trim(line)
              do i = pos_start, line_length
                if (line(i:i) == '\\') then
                  pos_end = i-1 
                  write(*,*) pos_end
                  exit
                endif
              enddo
              energy_str = line(pos_start:pos_end)
              read(energy_str,*) energy
              found_energy=.true.

            else if (index(line, 'NAtoms=') /= 0 .and. 
     $            .not. found_atoms) then
              read(line(10:),*) number_atoms
              found_atoms=.true.

            else if (index(line, 'Forces (Hartrees/Bohr)') /= 0 .and. 
     $            .not. found_gradient) then
              allocate(gradient_qm(number_atoms,3))
              read(426,'(A)', iostat=iostat) line 
              read(426,'(A)', iostat=iostat) line 
              read(426,'(A)', iostat=iostat) line 
              do i=1, number_atoms
                read(line,*) dummy, dummy, gradient_qm(i,1) 
     $                  , gradient_qm(i,2), gradient_qm(i,3)
                read(426,'(A)', iostat=iostat) line
                gradient_qm(i,1) = gradient_qm(i,1) * hartree / bohr
                gradient_qm(i,2) = gradient_qm(i,2) * hartree / bohr
                gradient_qm(i,3) = gradient_qm(i,3) * hartree / bohr
              enddo
              found_gradient=.true. 

           endif
          enddo
        close(426)

        else if(psi4_qm) then
          open(427, file=qm_grad_filename, status='old')
          do 
            read(427, '(A)', iostat=iostat) line
            if (iostat /= 0) exit

            if ((index(line, '= energy(''MP2'')') /= 0 .or. 
     $           index(line, '= energy(''mp2'')') /= 0) .and. 
     $                  .not. found_energy) then
              method='MP2'
            endif
            if (index(line, 'energy(''scf'')') /= 0 .and. 
     $                  .not. found_energy) then
              method='HF'
            endif
 
            if (index(trim(line), 'MP2 Total Energy (a.u.)') /= 0 
     $            .and. .not. found_energy) then
                do i=1,6
                  read(427,'(A)', iostat=iostat) line 
                enddo
                read(line(40:), *) energy_str 
                read(energy_str,*) energy
                found_energy = .true.

            else if (index(line, 'Number of atoms') /= 0 .and. 
     $            .not. found_atoms) then
               read(line(40:),*) number_atoms
               found_atoms=.true.
   
            else if (index(line, '-Total gradient:') /= 0 .and. 
     $            .not. found_gradient) then
               allocate(gradient_qm(number_atoms,3))
               read(427,'(A)', iostat=iostat) line 
               read(427,'(A)', iostat=iostat) line 
               read(427,'(A)', iostat=iostat) line 
               do i=1, number_atoms
                 read(line,*) dummy, gradient_qm(i,1) 
     $                  , gradient_qm(i,2), gradient_qm(i,3)
                 read(427,'(A)', iostat=iostat) line
                 gradient_qm(i,1) = gradient_qm(i,1) * hartree / bohr
                 gradient_qm(i,2) = gradient_qm(i,2) * hartree / bohr
                 gradient_qm(i,3) = gradient_qm(i,3) * hartree / bohr
               enddo
               found_gradient=.true. 

             endif
           enddo
         close(427)
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

      if (found_energy) then
        write(*,*) 'ENERGY = ', energy
      endif
      if (found_atoms) then
        write(*,*) 'NUMBER OF ATOMS = ', number_atoms
      endif
      if (found_gradient) then
          write(*,*) 'Gradient values:'
          do i = 1, number_atoms 
              write(*,'(A,I3,3F16.9)') 'Atom', i, gradient_qm(i, 1), 
     $                gradient_qm(i, 2), gradient_qm(i, 3)
          enddo
      endif



      end subroutine get_gradient_from_qm

      end module abinitio

