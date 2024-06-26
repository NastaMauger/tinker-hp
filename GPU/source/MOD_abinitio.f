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
      logical :: orca_qm, g16_qm, psi4_qm, pyscf_qm, qchem_qm
      integer :: multiplicity
      integer :: charge
      character*11 :: method
      character*11 :: basis
      character*11 :: jobtype

      integer :: nprocs
      integer :: maxcore

      integer :: compteur_aimd = 0

      character*40 filename_orca
      character*40 filename_g16 
      character*40 filename_psi4
      character*40 filename_pyscf
      character*40 filename_qchem
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

      method        = 'B3LYP'
      basis         = '6-31g'
      multiplicity  = 1
      charge        = 0

      nprocs        = 6         !Number of CPUs core that the QM software need
      maxcore       = 100000    !100Gb en Mb

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
        case ('METHOD')
          read (string,*,iostat=ios) method 
        case('BASIS')
          read (string,*,iostat=ios) basis
        case('MULTIPLICITY')
          read (string,*,iostat=ios) multiplicity
        case('CHARGE')
          read (string,*,iostat=ios) charge
        case('NPROCS')
          read (string,*,iostat=ios) nprocs 
        case('MAXCORE')
          read (string,*,iostat=ios) maxcore
          maxcore=maxcore*1e3
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
            write(*, *) 'Unknown QM software', trim(value)
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

      write(compteur_aimd_str, '(I10)') compteur_aimd

      if (compteur_aimd == 0) then
          write(filename_orca, '(A,I0,A)') 'orca_', compteur_aimd, '.in'
          write(filename_g16, '(A,I0,A)') 'g16_', compteur_aimd, '.in'
          write(filename_psi4, '(A,I0,A)') 'psi4_', compteur_aimd, '.in'
          write(filename_pyscf, '(A,I0,A)') 'pyscf_',compteur_aimd,'.in'
          write(filename_qchem, '(A,I0,A)') 'qchem_',compteur_aimd,'.in'
!      else
!          write(filename_orca, '(A,I0,A)') 'orca_', compteur_aimd, 
!     $              '_beads', numberbeads, '.in'
!          write(filename_orca, '(A,I0,A)') 'orca_', compteur_aimd,
!     $              '_beads', numberbeads, '.in'
!          write(filename_g16, '(A,I0,A)') 'g16_', compteur_aimd, 
!     $              '_beads', numberbeads, '.in'
!          write(filename_psi4, '(A,I0,A)') 'psi4_', compteur_aimd,
!     $              '_beads', numberbeads, '.in'
!          write(filename_pyscf, '(A,I0,A)') 'pyscf_',compteur_aimd,
!     $              '_beads', numberbeads, '.in'
!          write(filename_qchem, '(A,I0,A)') 'qchem_',compteur_aimd,
!     $              '_beads', numberbeads, '.in'

        if (orca_qm == .true.) then
          jobtype='sp'
          open(415, file=trim(filename_orca))
          write(415,'(A, 1X, A, 1X, A, 1X, A, A, A)') '!',
     &     trim(method), trim(basis), trim(jobtype), 
     &     ' ENGRAD xyzfile NoFrozenCore'
          write(415,'(A,I0,A)') '%PAL NPROCS ', nprocs, ' END'
          write(415,'(A,I0)') '%MAXCORE ', maxcore
          write(415,'(A,I0,1X, I0)') '*xyz ', charge, multiplicity
          call getxyz
          write(415,'(A)') '*'
          close(415)

        elseif (g16_qm == .true.) then
          maxcore=maxcore/1e3
          jobtype='Force'
          open(416, file=trim(filename_g16))
          write(416,'(A)') '%chk=g16_0.chk'
          write(416,'(A,I0)') '%nprocshared=', nprocs
          write(416,'(A,I0,A,I0)') '%mem=', maxcore,'GB'
          write(416, '(A, A, A ,A, A, A )') '#P ', trim(method), '/', 
     &                         trim(basis), ' ', trim(jobtype)
          write(416,'(A)') ''
          write(416,'(A)') 'Initial Input'
          write(416,'(A)') ''
          write(416,'(I0,1X, I0)') charge, multiplicity
          call getxyz
          write(416,'(A)') ''
          close(416)

        elseif (psi4_qm == .true.) then
          maxcore=maxcore/1e3
          open(417, file=trim(filename_psi4))
          write(417,'(A,I0,A,I0)') 'memory ', maxcore,' GB'
          write(417,'(A)') ''
          write(417, '(A)') 'molecule {'
          write(417,'(I0,1X, I0)') charge, multiplicity
          call getxyz
          write(417, '(A)') '}'
          write(417,'(A)') ''
          write(417, '(A)') 'set {'
          write(417, '(A)') ' basis ' // trim(basis)
          write(417, '(A)') ' ' // trim(method) // '_type conv'
          write(417, '(A)') '}'
          write(417,'(A)') ''
          write(417, '(A)') '#Compute Gradient'
          write(417, '(A)') 'grad = gradient('''// trim(method) // ''')'
          write(417,'(A)') ''
          write(417, '(A)') '#Compute Energy'
          write(417, '(A)') 'ener = energy(''' // trim(method) // ''')'
          close(417)

        elseif (pyscf_qm == .true.) then
          open(418, file=trim(filename_pyscf))
          write(418,'(A)') 'from pyscf import gto, scf, mp, grad'
          write(418,'(A)') 'import numpy as np'
          write(418,'(A)') 'import sys' 
          write(418,'(A)') ''
          write(418,'(A)') 'sys.stdout = open("pyscf_'
     $            // trim(adjustl(compteur_aimd_str)) //'.out", "w")'
          write(418,'(A)') 'sys.stderr = sys.stdout'
          write(418,'(A)') ''
          write(418,'(A)') 'water = '''''''
          call getxyz
          write(418,'(A)') ' '''''''
          write(418,'(A)') 'mol = gto.Mole()'
          write(418,'(A)') 'mol.atom = water'
          write(418, '(A,A,A)') 'mol.basis = "', trim(basis), '"'
          write(418,'(A)') 'mol.build()'
          write(418, '(A)') ''
          write(418, '(A)') '#Perform HF calculation'
          write(418, '(A)') 'myhf = scf.RHF(mol)'
          write(418, '(A)') 'myhf.kernel()'
          write(418, '(A)') ''
          write(418, '(A)') '#Perform '// trim(method) // ' calculation'
          write(418, '(A)') 'mymp = mp.' // trim(method) // '(myhf)'
          write(418, '(A)') 'mymp.kernel()'
          write(418, '(A)') ''
          write(418, '(A)') '#Compute the '// trim(method)// ' gradient'
          write(418, '(A)')  trim(method) // 
     $                        '_grad = grad.' // trim(method) //
     $                        '.Gradients(mymp)'
          write(418, '(A)') 'gradient = ' // trim(method) //
     $                      '.grad.kernel()'
          write(418, '(A)') ''
          write(418, '(A)') 'sys.stdout.close()'
          write(418, '(A)') 'sys.stdout = sys.__stdout__'
          write(418, '(A)') 'sys.stderr = sys.__stderr__'
          close(418)
     
        elseif (qchem_qm == .true.) then
          jobtype='Force'
          open(419, file=trim(filename_qchem))
          write(419,'(A)') '$molecule'
          write(419,'(1X,I0,1X, I0)') charge, multiplicity
          call getxyz
          write(419,'(A)') '$end'
          write(419, '(A)') ''
          write(419, '(A)') '$rem'
          write(419, '(A,A)') ' jobtype ', trim(jobtype)
          write(419, '(A,A)') ' method ', trim(method)
          write(419, '(A,A)') ' basis ', trim(basis)
          write(419, '(A,A)') ' N_FROZEN_CORE 0' 
          write(419, '(A,A)') ' SYM_IGNORE = true' 
          write(419, '(A,I0)') ' MEM_TOTAL ', maxcore
          write(419,'(A)') '$end'
          close(419)
        endif
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
      integer :: i,iostat
      integer :: number_atoms
      integer :: dummy, count_lines
      integer :: pos, pos_start, pos_end
      logical :: start_count
      logical qm_grad_file_found
      logical :: found_energy, found_gradient, found_atoms
      real :: energy
      real(r_p), allocatable :: gradient_qm(:,:)
      character*40 energy_str 
      character*400 line
      character*1 dummy_char

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
     $         // trim(adjustl(compteur_aimd_str)) // '.in.dat'

      else if(pyscf_qm == .true.) then
        qm_grad_filename = 'pyscf_' 
     $         // trim(adjustl(compteur_aimd_str)) // '.out'
        
        
      endif



      inquire(file=qm_grad_filename, exist=qm_grad_file_found)
      if(qm_grad_file_found) then 
        write(*,*) 'READIND GRADIENT FROM ', qm_grad_filename
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
            if (index(line, trim(method) // '=') /= 0 .and.
     $                  .not. found_energy) then
!!!    NEED TO HANDLE THE CASE WHERE THE VALUE IS SEPARATED INTO 2 LINES  !!!
              pos_start = index(line, trim(method) // '=')
              pos_start = pos_start + len(trim(method)) + 1 
              pos_end = index(line(pos_start:), ' ')
              pos_end = pos_end + len(trim(method))
              energy_str = (line(pos_start:pos_end))
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
            if (index(line, trim(method) // ' Total Energy') /= 0 .and.
     $                  .not. found_energy) then
              read(line(40:),*) energy_str 
              found_energy=.true.

            else if (index(line, 'Number of atoms:') /= 0 .and. 
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

        else if(pyscf_qm) then
          open(428, file=qm_grad_filename, status='old')
          do 
            read(428, '(A)', iostat=iostat) line
            if (iostat /= 0) exit
            if (index(line, 'E(' // trim(method) // ') =') /= 0 .and.
     $                  .not. found_energy) then
              read(line(10:),*) energy
              found_energy=.true.

            elseif (index(line,trim(method) // ' gradients') /= 0 .and. 
     $                  .not. found_atoms) then
              start_count = .true.
              count_lines = 1
              do while (start_count)
                read(428, '(A)', iostat=iostat) line
                  if (iostat /= 0 .or. index(line, '----') /= 0) then
                      start_count = .false.
                  else
                      count_lines = count_lines + 1
                  endif
              end do
              number_atoms = count_lines - 2
              found_atoms = .true.
      
            rewind(428)
            elseif (index(line, 'gradients') /= 0 .and. 
     $            .not. found_gradient) then
              allocate(gradient_qm(number_atoms,3))
              read(428,'(A)', iostat=iostat) line
              read(428,'(A)', iostat=iostat) line
              do i=1, number_atoms
                read(line,*) dummy, dummy_char, gradient_qm(i,1) 
     $                  , gradient_qm(i,2), gradient_qm(i,3)
                read(428,'(A)', iostat=iostat) line
                gradient_qm(i,1) = gradient_qm(i,1) * hartree / bohr
                gradient_qm(i,2) = gradient_qm(i,2) * hartree / bohr
                gradient_qm(i,3) = gradient_qm(i,3) * hartree / bohr
              enddo
              found_gradient=.true. 
              
            endif
          enddo
        close(428)
        
        endif
      else
        write(*,*) 'NO QM GRADIENT FILE FOUND'
        call fatal()
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

