  subroutine edgestat_BdG
     use wmpi
     use para
     implicit none

     integer :: ierr, doslfile, dosrfile, dosbulkfile,nm

     ! general loop index
     integer :: i, j, io, ikp, nw_half, spindoslfile, spindosrfile

     real(dp) :: emin, emax,eta_s, w, k, time_start, time_end, s0(3), s1(3)

     real(dp), allocatable :: omega(:)

     real(dp), allocatable :: dos_l(:,:), dos_r(:,:), dos_l_only(:, :), dos_r_only(:,:)
     real(dp), allocatable :: dos_l_mpi(:,:), dos_r_mpi(:,:), dos_bulk(:,:), dos_bulk_mpi(:,:)

     complex(dp), allocatable :: GLL(:,:), GRR(:,:), GB (:,:), H00(:,:), H01(:,:), ones(:,:)
 
     ! Spin resolved component
     REAL(DP),    ALLOCATABLE  :: sx_l(:, :), sy_l(:, :), sz_l(:, :)
     REAL(DP),    ALLOCATABLE  :: sx_r(:, :), sy_r(:, :), sz_r(:, :)
     REAL(DP),    ALLOCATABLE  :: sx_l_mpi(:, :), sy_l_mpi(:, :), sz_l_mpi(:, :)
     REAL(DP),    ALLOCATABLE  :: sx_r_mpi(:, :), sy_r_mpi(:, :), sz_r_mpi(:, :)
     COMPLEX(DP), ALLOCATABLE  :: sigma_x(:,:), sigma_y(:,:), sigma_z(:,:)
     COMPLEX(DP), ALLOCATABLE  :: ctemp(:,:)
     
     nm=2*Np*nslab*num_wann
     allocate( omega(omeganum), dos_l(knv2, omeganum), dos_r(knv2, omeganum))
     allocate( dos_l_only(knv2, omeganum), dos_r_only(knv2, omeganum))
     allocate( dos_l_mpi(knv2, omeganum), dos_r_mpi(knv2, omeganum))
     allocate( dos_bulk(knv2, omeganum), dos_bulk_mpi(knv2, omeganum))
     omega=0d0;  dos_l=0d0;  dos_r=0d0;  dos_l_only=0d0;  dos_r_only=0d0
     dos_l_mpi=0d0;  dos_r_mpi=0d0;  dos_bulk=0d0;  dos_bulk_mpi=0d0

     ALLOCATE( sx_l(knv2, omeganum), sy_l(knv2, omeganum), sz_l(knv2, omeganum))
     ALLOCATE( sx_l_mpi(knv2, omeganum), sy_l_mpi(knv2, omeganum), sz_l_mpi(knv2, omeganum))
     ALLOCATE( sx_r(knv2, omeganum), sy_r(knv2, omeganum), sz_r(knv2, omeganum))
     ALLOCATE( sx_r_mpi(knv2, omeganum), sy_r_mpi(knv2, omeganum), sz_r_mpi(knv2, omeganum))
     ALLOCATE( sigma_x(nm,nm), sigma_y(nm,nm),sigma_z(nm,nm),ctemp(nm,nm))
     sigma_x      = 0d0;      sigma_y      = 0d0;      sigma_z      = 0d0
     sx_l         = 0d0;      sy_l         = 0d0;      sz_l         = 0d0
     sx_r         = 0d0;      sy_r         = 0d0;      sz_r         = 0d0
     sx_l_mpi     = 0d0;      sy_l_mpi     = 0d0;      sz_l_mpi     = 0d0
     sx_r_mpi     = 0d0;      sy_r_mpi     = 0d0;      sz_r_mpi     = 0d0


    ! eta_s=(omegamax- omegamin)/dble(omeganum)*1.0d0
     eta=Eta_Arc

     do i= 1, omeganum
        omega(i)=omegamin+(i-1)*(omegamax-omegamin)/dble(omeganum)
     enddo

     if (index(Particle,'phonon')/=0) then
        do i= 1, omeganum
           omega(i)= omega(i)*omega(i)
        enddo
     endif

     allocate(GLL(nm,nm), GRR(nm,nm),GB (nm,nm))
     allocate(H00(nm,nm), H01(nm,nm),ones(nm,nm))
     GLL= 0d0; GRR= 0d0; GB = 0d0; H00= 0d0; H01= 0d0; ones= 0d0

     nw_half = nslab*Num_wann/2
     do i=1, Np
        do j=1, nw_half
           sigma_x( nslab*Num_wann*(i-1)+j        , nslab*Num_wann*(i-1)+j+nw_half ) =  1.0d0
           sigma_x( nslab*Num_wann*(i-1)+j+nw_half, nslab*Num_wann*(i-1)+j         ) =  1.0d0
           sigma_y( nslab*Num_wann*(i-1)+j        , nslab*Num_wann*(i-1)+j+nw_half ) = -zi
           sigma_y( nslab*Num_wann*(i-1)+j+nw_half, nslab*Num_wann*(i-1)+j         ) =  zi
           sigma_z( nslab*Num_wann*(i-1)+j        , nslab*Num_wann*(i-1)+j         ) =  1.0d0
           sigma_z( nslab*Num_wann*(i-1)+j+nw_half, nslab*Num_wann*(i-1)+j+nw_half ) = -1.0d0
           sigma_x( nm/2+nslab*Num_wann*(i-1)+j,nm/2+nslab*Num_wann*(i-1)+j+nw_half ) =  1.0d0
           sigma_x( nm/2+nslab*Num_wann*(i-1)+j+nw_half, nm/2+nslab*Num_wann*(i-1)+j) =  1.0d0
           sigma_y( nm/2+nslab*Num_wann*(i-1)+j,nm/2+nslab*Num_wann*(i-1)+j+nw_half ) = -zi
           sigma_y( nm/2+nslab*Num_wann*(i-1)+j+nw_half, nm/2+nslab*Num_wann*(i-1)+j) =  zi
           sigma_z( nm/2+nslab*Num_wann*(i-1)+j, nm/2+nslab*Num_wann*(i-1)+j) =  1.0d0
           sigma_z( nm/2+nslab*Num_wann*(i-1)+j+nw_half,nm/2+nslab*Num_wann*(i-1)+j+nw_half ) = -1.0d0
        enddo 
     enddo


     do i=1,nm
        ones(i,i)=1.0d0
     enddo
     time_start= 0d0
     time_end= 0d0
     if ( cpuid==0 )then
         write(stdout, *) '  EdgeSS, ik',  '        Nk', '   time left'
        endif
     do ikp= 1+cpuid, knv2, num_cpu
        k= k2_path(ikp,1)

        call now(time_start)
        call ham_qlayer2qlayerribbon_BdG(k,H00,H01)

        do j = 1, omeganum
           w=omega(j)

           call surfgreen_edge(w,GLL,GRR,GB,H00,H01,ones,eta_s)
           if(cpuid==0) print *,j
           ! calculate spectral function
           do i= 1, NslabTopOrbitals
              io= SlabTopOrbitals(i)
              dos_l(ikp, j)=dos_l(ikp,j)- aimag(GLL(io,io))
           enddo ! i
           do i= 1,NslabTopOrbitals
              io= nm/2+SlabTopOrbitals(i)
              dos_l(ikp, j)=dos_l(ikp,j)- aimag(GLL(io,io))
           enddo ! i
           do i= 1, NslabBottomOrbitals
              io= nm/2-nslab*num_wann+SlabBottomOrbitals(i)
              dos_r(ikp, j)=dos_r(ikp,j)- aimag(GRR(io,io))
           enddo ! i

           do i= 1,NslabBottomOrbitals
              io= nm-nslab*num_wann+SlabBottomOrbitals(i)
              dos_r(ikp, j)=dos_r(ikp,j)- aimag(GRR(io,io))
           enddo ! i

           do i= 1, nm
              dos_bulk(ikp, j)=dos_bulk(ikp,j)- aimag(GB(i,i))
           enddo ! i

            ! Spin resolved sprectrum
    !        call mat_mul(nm,gll,sigma_x,ctemp)
    !        do i = 1, NSlabtopOrbitals
    !            io = SlabTopOrbitals(i)
    !            sx_l_mpi(ikp, j) = sx_l_mpi(ikp, j)- AIMAG(ctemp(io,io))
    !        enddo ! i
    !        call mat_mul(nm,gll,sigma_y,ctemp)
    !        do i = 1, NSlabtopOrbitals
    !            io = SlabTopOrbitals(i)
    !            sy_l_mpi(ikp, j) = sy_l_mpi(ikp, j)- AIMAG(ctemp(io,io))
    !        enddo !
    !        call mat_mul(nm,gll,sigma_z,ctemp)
    !        do i = 1, NSlabtopOrbitals
    !            io = SLabTopOrbitals(i)
    !            sz_l_mpi(ikp, j) = sz_l_mpi(ikp, j)- AIMAG(ctemp(io,io))
    !        enddo ! i
    !        call mat_mul(nm,grr,sigma_x,ctemp)
    !        do i = 1, NSlabtopOrbitals
    !            io = SlabTopOrbitals(i)
    !            sx_r_mpi(ikp, j) = sx_r_mpi(ikp, j)- AIMAG(ctemp(io,io))
    !        enddo ! i
    !        call mat_mul(nm,grr,sigma_y,ctemp)
    !        do i = 1, NSlabtopOrbitals
    !            io = SlabTopOrbitals(i)
    !            sy_r_mpi(ikp, j) = sy_r_mpi(ikp, j)- AIMAG(ctemp(io,io))
    !        enddo !
    !        call mat_mul(nm,grr,sigma_z,ctemp)
    !        do i = 1, NSLabtopOrbitals
    !            io = SlabTopOrbitals(i)
    !            sz_r_mpi(ikp, j) = sz_r_mpi(ikp, j)- AIMAG(ctemp(io,io))
    !        enddo ! i


        enddo ! j
        call now(time_end)

        if ( cpuid==0 )then 
         write(stdout, *)  ikp, knv2,(knv2-ikp)*(time_end- time_start)/num_cpu, ' s'
        endif
     enddo ! ikp

!> we do have to do allreduce operation
!#ifdef MPI
!    call mpi_allreduce(sx_l_mpi, sx_l, SIZE(sx_l), mpi_double_precision,&
!    mpi_sum, mpi_comm_world, ierr)
!    call mpi_allreduce(sy_l_mpi, sy_l, SIZE(sy_l), mpi_double_precision,&
!    mpi_sum, mpi_comm_world, ierr)
!    call mpi_allreduce(sz_l_mpi, sz_l, SIZE(sz_l), mpi_double_precision,&
!    mpi_sum, mpi_comm_world, ierr)
!    call mpi_allreduce(sx_r_mpi, sx_r, SIZE(sx_r), mpi_double_precision,&
!    mpi_sum, mpi_comm_world, ierr)
!    call mpi_allreduce(sy_r_mpi, sy_r, SIZE(sy_r), mpi_double_precision,&
!    mpi_sum, mpi_comm_world, ierr)
!    call mpi_allreduce(sz_r_mpi, sz_r, SIZE(sz_r), mpi_double_precision,&
!    mpi_sum, mpi_comm_world, ierr)
!#else
!    sx_l     = sx_l_mpi;      sy_l     = sy_l_mpi;      sz_l     = sz_l_mpi
!    sx_r     = sx_r_mpi;      sy_r     = sy_r_mpi;      sz_r     = sz_r_mpi
!#endif

    deallocate( sx_l_mpi, sy_l_mpi, sz_l_mpi )
    deallocate( sx_r_mpi, sy_r_mpi, sz_r_mpi )
    deallocate( sigma_x, sigma_y, sigma_z, ctemp )


     !> we don't have to do allreduce operation
#if defined (MPI)
     call mpi_reduce(dos_l, dos_l_mpi, size(dos_l), mpi_double_precision,&
                     mpi_sum, 0, mpi_comm_world, ierr)
     call mpi_reduce(dos_r, dos_r_mpi, size(dos_r), mpi_double_precision,&
                     mpi_sum, 0, mpi_comm_world, ierr)
     call mpi_reduce(dos_bulk, dos_bulk_mpi, size(dos_bulk), mpi_double_precision,&
                     mpi_sum, 0, mpi_comm_world, ierr)
#else
     dos_l_mpi= dos_l
     dos_r_mpi= dos_r
     dos_bulk_mpi= dos_bulk
#endif

     dos_l=log(abs(dos_l_mpi))
     dos_r=log(abs(dos_r_mpi))
     dos_bulk=log(abs(dos_bulk_mpi)+eps9)
     do ikp=1, knv2
        do j=1, omeganum
           dos_l_only(ikp, j)= dos_l_mpi(ikp, j)- dos_bulk_mpi(ikp, j)
           if (dos_l_only(ikp, j)<0) dos_l_only(ikp, j)=eps9
           dos_r_only(ikp, j)= dos_r_mpi(ikp, j)- dos_bulk_mpi(ikp, j)
           if (dos_r_only(ikp, j)<0) dos_r_only(ikp, j)=eps9
        enddo
     enddo

     outfileindex= outfileindex+ 1
     doslfile= outfileindex
     outfileindex= outfileindex+ 1
     dosrfile= outfileindex
     outfileindex= outfileindex+ 1
     dosbulkfile= outfileindex
     outfileindex = outfileindex+ 1
!     spindoslfile = outfileindex
     outfileindex = outfileindex+ 1
!     spindosrfile = outfileindex

     !> deal with phonon system
     !> for phonon system, omega should be changed to omega^2
     if (index(Particle,'phonon')/=0) then
        do i= 1, omeganum
           omega(i)= sqrt(omega(i))
        enddo
     endif

 
     ! Write surface state to files
     IF (cpuid .eq. 0) THEN
         open(unit=doslfile    , file='dos.dat_l')
         open(unit=dosrfile    , file='dos.dat_r')
         open(unit=dosbulkfile , file='dos.dat_bulk')
!         open(unit=spindoslfile, file='spindos.dat_l')
!         open(unit=spindosrfile, file='spindos.dat_r')
         write(doslfile, '("#", a12, 6a17)') ' k(1/Ang)', ' E(eV)', 'dos_l', 'dos_l_only'
         write(dosrfile, '("#", a12, 6a17)') ' k(1/Ang)', ' E(eV)', 'dos_r', 'dos_r_only'
         write(dosbulkfile, '("#", a12, 6a17)') ' k(1/Ang)', ' E(eV)', 'dos_bulk'
   !      write(spindoslfile, '("#", a)') ' spin dos_l, the axis is rotated as '
   !      write(spindoslfile, '("#", a)') " x is along R1', z is along R1'xR2', y is along z x y"
   !      write(spindoslfile, '("#", a12, 6a17)') ' k(1/Ang)', ' E(eV)', 'sx', 'sy', 'sz'
   !      write(spindosrfile, '("#", a)') ' spin dos_r, the axis is rotated as '
   !      write(spindosrfile, '("#", a)') " x is along R1', z is along R1'xR2', y is along z x y"
   !      write(spindosrfile, '("#", a12, 6a17)') ' k(1/Ang)', ' E(eV)', 'sx', 'sy', 'sz'
         do ikp = 1, knv2
            do j = 1, omeganum
                write(doslfile,    2002) k2len(ikp), omega(j)/eV2Hartree, dos_l(ikp, j), dos_l_only(ikp, j)
                write(dosrfile,    2002) k2len(ikp), omega(j)/eV2Hartree, dos_r(ikp, j), dos_r_only(ikp, j)
                write(dosbulkfile, 2003) k2len(ikp), omega(j)/eV2Hartree, dos_bulk(ikp, j)
         !       s0(1)=sx_l(ikp, j); s0(2)=sy_l(ikp, j); s0(3)=sz_l(ikp, j); 
         !       call rotate(s0, s1)
         !       write(spindoslfile,2001) k2len(ikp), omega(j)/eV2Hartree, s1
         !       s0(1)=sx_r(ikp, j); s0(2)=sy_r(ikp, j); s0(3)=sz_r(ikp, j); 
         !       call rotate(s0, s1)
         !       write(spindosrfile,2001) k2len(ikp), omega(j)/eV2Hartree, s1
             ENDDO
             write(doslfile    , *)
             write(dosrfile    , *)
             write(dosbulkfile , *)
          !   write(spindoslfile, *)
          !   write(spindosrfile, *)
         ENDDO
         close(doslfile)
         close(dosrfile)
         close(dosbulkfile)
       !  close(spindoslfile)
       !  close(spindosrfile)
         write(stdout,*)'ndim',ndim
         write(stdout,*) 'knv2,omeganum,eta',knv2, omeganum, eta/eV2Hartree
         write(stdout,*)'calculate density of state successfully'
     ENDIF



2001 FORMAT(5(1X,E16.8))
2002 FORMAT(4(1X,E16.8))
2003 FORMAT(3(1X,E16.8))

     emin= minval(omega)/eV2Hartree
     emax= maxval(omega)/eV2Hartree
     !> write script for gnuplot
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='surfdos_l.gnu')
        write(outfileindex, '(a)')"set encoding iso_8859_1"
        write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
        write(outfileindex, '(a)')"#set output 'surfdos_l.eps'"
        write(outfileindex, '(3a)')'#set terminal  pngcairo truecolor enhanced', &
           '  font ", 60" size 1920, 1680'
        write(outfileindex, '(3a)')'set terminal  png truecolor enhanced', &
           ' font ", 60" size 1920, 1680'
        write(outfileindex, '(a)')"set output 'surfdos_l.png'"
        write(outfileindex,'(2a)') 'set palette defined (-10 "#194eff", ', &
           '0 "white", 10 "red" )'
        write(outfileindex, '(a)')'#set palette rgbformulae 33,13,10'
        write(outfileindex, '(a)')'set style data linespoints'
        write(outfileindex, '(a)')'set size 0.8, 1'
        write(outfileindex, '(a)')'set origin 0.1, 0'
        write(outfileindex, '(a)')'unset ztics'
        write(outfileindex, '(a)')'unset key'
        write(outfileindex, '(a)')'set pointsize 0.8'
        write(outfileindex, '(a)')'set pm3d'
        write(outfileindex, '(a)')'#set view equal xyz'
        write(outfileindex, '(a)')'set view map'
        write(outfileindex, '(a)')'set border lw 3'
        write(outfileindex, '(a)')'#set cbtics font ",48"'
        write(outfileindex, '(a)')'#set xtics font ",48"'
        write(outfileindex, '(a)')'#set ytics font ",48"'
        write(outfileindex, '(a)')'#set ylabel font ",48"'
        write(outfileindex, '(a)')'set ylabel "Energy (eV)"'
        write(outfileindex, '(a)')'#set xtics offset 0, -1'
        write(outfileindex, '(a)')'#set ylabel offset -6, 0 '
        write(outfileindex, '(a, f18.5, a)')'set xrange [0: ', maxval(k2len), ']'
        write(outfileindex, '(a, f18.5, a, f18.5, a)')'set yrange [', emin, ':', emax, ']'
        write(outfileindex, 202, advance="no") (k2line_name(i), k2line_stop(i), i=1, nk2lines)
        write(outfileindex, 203)k2line_name(nk2lines+1), k2line_stop(nk2lines+1)

        do i=1, nk2lines-1
           write(outfileindex, 204)k2line_stop(i+1), emin, k2line_stop(i+1), emax
        enddo
        write(outfileindex, '(a)')'set pm3d interpolate 2,2'
        write(outfileindex, '(2a)')"splot 'dos.dat_l' u 1:2:3 w pm3d"

      !  write(outfileindex, '(3a)')'set terminal png truecolor enhanced', &
      !     ' font ", 30" size 1920, 1680'
      !  write(outfileindex, '(a)')"set output 'spindos_l.png'"
      !  write(outfileindex, '(a)')"set multiplot layout 3, 1"
      !  write(outfileindex, '(a)')"set title 'sx'"
      !  write(outfileindex, '(2a)')"splot 'spindos.dat_l' u 1:2:3 w pm3d "
      !  write(outfileindex, '(a)')"set title 'sy'"
      !  write(outfileindex, '(2a)')"splot 'spindos.dat_l' u 1:2:4 w pm3d"
      !  write(outfileindex, '(a)')"set title 'sz'"
      !  write(outfileindex, '(2a)')"splot 'spindos.dat_l' u 1:2:5 w pm3d"
 

        close(outfileindex)

     endif

     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='surfdos_l_only.gnu')
        write(outfileindex, '(a)')"set encoding iso_8859_1"
        write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
        write(outfileindex, '(a)')"#set output 'surfdos_l.eps'"
        write(outfileindex, '(3a)')'#set terminal  pngcairo truecolor enhanced', &
           '  font ", 60" size 1920, 1680'
        write(outfileindex, '(3a)')'set terminal  png truecolor enhanced', &
           ' font ", 60" size 1920, 1680'
        write(outfileindex, '(a)')"set output 'surfdos_l_only.png'"
        write(outfileindex,'(2a)') 'set palette defined (0  "white", ', &
           '6 "red", 20 "black" )'
        write(outfileindex, '(a)')'#set palette rgbformulae 33,13,10'
        write(outfileindex, '(a)')'set style data linespoints'
        write(outfileindex, '(a)')'set size 0.8, 1'
        write(outfileindex, '(a)')'set origin 0.1, 0'
        write(outfileindex, '(a)')'unset ztics'
        write(outfileindex, '(a)')'unset key'
        write(outfileindex, '(a)')'set pointsize 0.8'
        write(outfileindex, '(a)')'set pm3d'
        write(outfileindex, '(a)')'#set view equal xyz'
        write(outfileindex, '(a)')'set view map'
        write(outfileindex, '(a)')'set border lw 3'
        write(outfileindex, '(a)')'#set cbtics font ",48"'
        write(outfileindex, '(a)')'#set xtics font ",48"'
        write(outfileindex, '(a)')'#set ytics font ",48"'
        write(outfileindex, '(a)')'unset cbtics'
        write(outfileindex, '(a)')'#set ylabel font ",48"'
        write(outfileindex, '(a)')'set ylabel "Energy (eV)"'
        write(outfileindex, '(a)')'#set xtics offset 0, -1'
        write(outfileindex, '(a)')'#set ylabel offset -6, 0 '
        write(outfileindex, '(a, f18.5, a)')'set xrange [0: ', maxval(k2len), ']'
        write(outfileindex, '(a, f18.5, a, f18.5, a)')'set yrange [', emin, ':', emax, ']'
        write(outfileindex, 202, advance="no") (k2line_name(i), k2line_stop(i), i=1, nk2lines)
        write(outfileindex, 203)k2line_name(nk2lines+1), k2line_stop(nk2lines+1)

        do i=1, nk2lines-1
           write(outfileindex, 204)k2line_stop(i+1), emin, k2line_stop(i+1), emax
        enddo
        write(outfileindex, '(a)')'set pm3d interpolate 2,2'
        write(outfileindex, '(2a)')"splot 'dos.dat_l' u 1:2:4 w pm3d"
        CLOSE(outfileindex)

     endif

     !> write script for gnuplot
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='surfdos_r_only.gnu')
        write(outfileindex, '(a)')"set encoding iso_8859_1"
        write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
        write(outfileindex, '(a)')"#set output 'surfdos_r.eps'"
        write(outfileindex, '(3a)')'#set terminal pngcairo truecolor enhanced', &
           '  font ", 60" size 1920, 1680'
        write(outfileindex, '(3a)')'set terminal png truecolor enhanced', &
           ' font ", 60" size 1920, 1680'
        write(outfileindex, '(a)')"set output 'surfdos_r_only.png'"
        write(outfileindex,'(2a)') 'set palette defined (0  "white", ', &
           '6 "red", 20 "black" )'
        write(outfileindex, '(a)')'#set palette rgbformulae 33,13,10'
        write(outfileindex, '(a)')'set style data linespoints'
        write(outfileindex, '(a)')'unset ztics'
        write(outfileindex, '(a)')'unset key'
        write(outfileindex, '(a)')'set pointsize 0.8'
        write(outfileindex, '(a)')'set pm3d'
        write(outfileindex, '(a)')'set border lw 3'
        write(outfileindex, '(a)')'set size 0.8, 1'
        write(outfileindex, '(a)')'set origin 0.1, 0'
        write(outfileindex, '(a)')'#set size ratio -1'
        write(outfileindex, '(a)')'#set view equal xyz'
        write(outfileindex, '(a)')'set view map'
        write(outfileindex, '(a)')'#set cbtics font ",48"'
        write(outfileindex, '(a)')'#set xtics font ",48"'
        write(outfileindex, '(a)')'#set ytics font ",48"'
        write(outfileindex, '(a)')'#set ylabel font ",48"'
        write(outfileindex, '(a)')'set ylabel "Energy (eV)"'
        write(outfileindex, '(a)')'unset cbtics'
        write(outfileindex, '(a)')'#set xtics offset 0, -1'
        write(outfileindex, '(a)')'#set ylabel offset -6, 0 '
        write(outfileindex, '(a, f18.5, a)')'set xrange [0: ', maxval(k2len), ']'
        write(outfileindex, '(a, f18.5, a, f18.5, a)')'set yrange [', emin, ':', emax, ']'
        write(outfileindex, 202, advance="no") (k2line_name(i), k2line_stop(i), i=1, nk2lines)
        write(outfileindex, 203)k2line_name(nk2lines+1), k2line_stop(nk2lines+1)

        do i=1, nk2lines-1
           write(outfileindex, 204)k2line_stop(i+1), emin, k2line_stop(i+1), emax
        enddo
        write(outfileindex, '(a)')'set pm3d interpolate 2,2'
        write(outfileindex, '(2a)')"splot 'dos.dat_r' u 1:2:4 w pm3d"

        CLOSE(outfileindex)

     endif


     !> write script for gnuplot
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='surfdos_r.gnu')
        write(outfileindex, '(a)')"set encoding iso_8859_1"
        write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
        write(outfileindex, '(a)')"#set output 'surfdos_r.eps'"
        write(outfileindex, '(3a)')'#set terminal pngcairo truecolor enhanced', &
           '  font ", 60" size 1920, 1680'
        write(outfileindex, '(3a)')'set terminal png truecolor enhanced', &
           ' font ", 60" size 1920, 1680'
        write(outfileindex, '(a)')"set output 'surfdos_r.png'"
        write(outfileindex,'(2a)') 'set palette defined (-10 "#194eff", ', &
           '0 "white", 10 "red" )'
        write(outfileindex, '(a)')'#set palette rgbformulae 33,13,10'
        write(outfileindex, '(a)')'set style data linespoints'
        write(outfileindex, '(a)')'unset ztics'
        write(outfileindex, '(a)')'unset key'
        write(outfileindex, '(a)')'set pointsize 0.8'
        write(outfileindex, '(a)')'set pm3d'
        write(outfileindex, '(a)')'set border lw 3'
        write(outfileindex, '(a)')'set size 0.8, 1'
        write(outfileindex, '(a)')'set origin 0.1, 0'
        write(outfileindex, '(a)')'#set size ratio -1'
        write(outfileindex, '(a)')'#set view equal xyz'
        write(outfileindex, '(a)')'set view map'
        write(outfileindex, '(a)')'#set cbtics font ",48"'
        write(outfileindex, '(a)')'#set xtics font ",48"'
        write(outfileindex, '(a)')'#set ytics font ",48"'
        write(outfileindex, '(a)')'#set ylabel font ",48"'
        write(outfileindex, '(a)')'set ylabel "Energy (eV)"'
        write(outfileindex, '(a)')'#set xtics offset 0, -1'
        write(outfileindex, '(a)')'#set ylabel offset -6, 0 '
        write(outfileindex, '(a, f18.5, a)')'set xrange [0: ', maxval(k2len), ']'
        write(outfileindex, '(a, f18.5, a, f18.5, a)')'set yrange [', emin, ':', emax, ']'
        write(outfileindex, 202, advance="no") (k2line_name(i), k2line_stop(i), i=1, nk2lines)
        write(outfileindex, 203)k2line_name(nk2lines+1), k2line_stop(nk2lines+1)

        do i=1, nk2lines-1
           write(outfileindex, 204)k2line_stop(i+1), emin, k2line_stop(i+1), emax
        enddo
        write(outfileindex, '(a)')'set pm3d interpolate 2,2'
        write(outfileindex, '(2a)')"splot 'dos.dat_r' u 1:2:3 w pm3d"

      !  write(outfileindex, '(3a)')'set terminal png truecolor enhanced', &
      !     ' font ", 30" size 1920, 1680'
      !  write(outfileindex, '(a)')"set output 'spindos_r.png'"
      !  write(outfileindex, '(a)')"set multiplot layout 3, 1"
      !  write(outfileindex, '(a)')"set title 'sx'"
      !  write(outfileindex, '(2a)')"splot 'spindos.dat_r' u 1:2:3 w pm3d "
      !  write(outfileindex, '(a)')"set title 'sy'"
      !  write(outfileindex, '(2a)')"splot 'spindos.dat_r' u 1:2:4 w pm3d"
      !  write(outfileindex, '(a)')"set title 'sz'"
      !  write(outfileindex, '(2a)')"splot 'spindos.dat_r' u 1:2:5 w pm3d"
        close(outfileindex)

     endif

     !> write script for gnuplot
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='surfdos_bulk.gnu')
        write(outfileindex, '(a)')"set encoding iso_8859_1"
        write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
        write(outfileindex, '(a)')"#set output 'surfdos_bulk.eps'"
        write(outfileindex, '(3a)')'#set terminal  pngcairo truecolor enhanced', &
           '  font ", 60" size 1920, 1680'
        write(outfileindex, '(3a)')'set terminal  png truecolor enhanced', &
           ' font ", 60" size 1920, 1680'
        write(outfileindex, '(a)')"set output 'surfdos_bulk.png'"
        write(outfileindex,'(2a)') 'set palette defined (-10 "#194eff", ', &
           '0 "white", 10 "red" )'
        write(outfileindex, '(a)')'#set palette rgbformulae 33,13,10'
        write(outfileindex, '(a)')'set style data linespoints'
        write(outfileindex, '(a)')'set size 0.8, 1'
        write(outfileindex, '(a)')'set origin 0.1, 0'
        write(outfileindex, '(a)')'unset ztics'
        write(outfileindex, '(a)')'unset key'
        write(outfileindex, '(a)')'set pointsize 0.8'
        write(outfileindex, '(a)')'set pm3d'
        write(outfileindex, '(a)')'#set view equal xyz'
        write(outfileindex, '(a)')'set view map'
        write(outfileindex, '(a)')'set border lw 3'
        write(outfileindex, '(a)')'#set cbtics font ",48"'
        write(outfileindex, '(a)')'#set xtics font ",48"'
        write(outfileindex, '(a)')'#set ytics font ",48"'
        write(outfileindex, '(a)')'#set ylabel font ",48"'
        write(outfileindex, '(a)')'unset cbtics'
        write(outfileindex, '(a)')'set ylabel "Energy (eV)"'
        write(outfileindex, '(a)')'#set xtics offset 0, -1'
        write(outfileindex, '(a)')'#set ylabel offset -6, 0 '
        write(outfileindex, '(a, f18.5, a)')'set xrange [0: ', maxval(k2len), ']'
        write(outfileindex, '(a, f18.5, a, f18.5, a)')'set yrange [', emin, ':', emax, ']'
        write(outfileindex, 202, advance="no") (k2line_name(i), k2line_stop(i), i=1, nk2lines)
        write(outfileindex, 203)k2line_name(nk2lines+1), k2line_stop(nk2lines+1)

        do i=1, nk2lines-1
           write(outfileindex, 204)k2line_stop(i+1), emin, k2line_stop(i+1), emax
        enddo
        write(outfileindex, '(a)')'set pm3d interpolate 2,2'
        write(outfileindex, '(2a)')"splot 'dos.dat_bulk' u 1:2:(exp($3)) w pm3d"
        CLOSE(outfileindex)

     endif

     202 format('set xtics (',:20('"',A1,'" ',F8.5,','))
     203 format(A1,'" ',F8.5,')')
     204 format('set arrow from ',F8.5,',',F10.5,' to ',F8.5,',',F10.5, ' nohead front lw 3')

#if defined (MPI)
     call mpi_barrier(mpi_cmw, ierr)
#endif

     deallocate( omega)
     deallocate( dos_l)
     deallocate( dos_r)
     deallocate( dos_l_only)
     deallocate( dos_r_only)
     deallocate( dos_l_mpi)
     deallocate( dos_r_mpi)
     deallocate( dos_bulk)
     deallocate( dos_bulk_mpi)
     deallocate(GLL)
     deallocate(GRR)
     deallocate(GB )
     deallocate(H00)
     deallocate(H01)
     deallocate(ones)

  return
  end subroutine edgestat_BdG
