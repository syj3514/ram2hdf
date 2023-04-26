module readramses
    use omp_lib
    implicit none
    real(kind=8),    dimension(:,:),   allocatable :: real_table
    integer(kind=4), dimension(:,:),   allocatable :: integer_table
    integer(kind=8), dimension(:,:),   allocatable :: long_table
    integer(kind=1), dimension(:,:),   allocatable :: byte_table
    real(kind=8),    dimension(:),     allocatable :: sfrs
    integer :: nhvar
    integer :: ncell_tot, ndm, nsink, nstar, ncloud
    integer(kind=8) :: npart_tot
    real(kind=8) :: aexp

contains

!#####################################################################
    subroutine read_part(repo, iout, cpu_list, verbose, nthread)
!#####################################################################
        implicit none

        integer :: i, icpu
        integer :: ncpu, ndim
        integer :: npart, nthread

        integer :: part_n
        logical :: ok

        character(len=128) :: file_path

        character(len=128),    intent(in) :: repo
        integer,               intent(in) :: iout, verbose
        integer, dimension(:), intent(in) :: cpu_list

        part_n = 30

        file_path = part_filename(repo, iout, cpu_list(1))

        ! Step 1: Verify there is file
        if(verbose>1)write(*,'(a)') '> Check if there is file: '//trim(file_path)
        inquire(file=file_path, exist=ok)
        if ( .not. ok ) then
            write(*,'(a)') file_path, ' not found.'
            stop
        endif

        ! Step 2: Count the total number of particles.
        if(verbose>1)write(*,'(a)') '> Check the total number of stars'
        open(unit=part_n, file=file_path, status='old', form='unformatted')
        read(part_n) ncpu
        read(part_n) ndim
        call skip_read(part_n, 2)
        read(part_n) nstar
        call skip_read(part_n, 2)
        read(part_n) nsink
        close(part_n)
        ncloud = 2109*nsink

        if(verbose>1)write(*,'(a)') '> Check the total number of all particles...'
        npart_tot = 0
        if(verbose>1)write(6, '(a)', advance='no') '>> Progress: '
        !$OMP PARALLEL DO SHARED(npart_tot) PRIVATE(i, icpu, npart, part_n) NUM_THREADS(nthread)
        do i = 1, SIZE(cpu_list)
            !$OMP CRITICAL
            if(verbose>1)call progress_bar(i, SIZE(cpu_list))
            !$OMP END CRITICAL
            icpu = cpu_list(i)
            part_n = 30 + omp_get_thread_num()
            call get_npart(part_n, part_filename(repo, iout, icpu), npart)
            !$OMP ATOMIC
            npart_tot = npart_tot + npart
        end do
        !$OMP END PARALLEL DO
        ndm = npart_tot - nstar - ncloud
        if(verbose>1)write(*,'(a)') '> Done'
    end subroutine read_part


!#####################################################################
    subroutine read_cell(repo, iout, cpu_list, verbose)
!#####################################################################
        implicit none

        integer :: i, icpu, jcpu, ilevel, jdim, igrid
        integer :: ncpu, ndim, nlevelmax, nboundary, ngrid_a, ngridmax

        integer :: amr_n, twotondim, skip_amr
        logical :: ok
        character(LEN=10) :: ordering
        character(len=128) :: file_path

        ! temporary arrays
        integer(kind=4), dimension(:,:),   allocatable :: ngridfile
        integer(kind=4), dimension(:),     allocatable :: son

        character(len=128),    intent(in) :: repo
        integer,               intent(in) :: iout, verbose
        integer, dimension(:), intent(in) :: cpu_list

        amr_n = 10

        ordering='hilbert'

        file_path = amr_filename(repo, iout, cpu_list(1))

        ! Step 1: Check if there is file
        inquire(file=file_path, exist=ok)
        if(verbose>0)write(*,'(a)') '> Check if there is file: ',file_path
        if ( .not. ok ) then
            write(*,'(a)')'File not found in repo: '//trim(file_path)
            stop
        endif

        ! Step 2: Read first file for header
        if(verbose>0)write(*,'(a)') '> Read first file for header'
        open(unit=amr_n, file=amr_filename(repo, iout, cpu_list(1)), status='old', form='unformatted')
        read(amr_n) ncpu
        read(amr_n) ndim
        read(amr_n)
        read(amr_n) nlevelmax
        read(amr_n)
        read(amr_n) nboundary

        allocate(ngridfile(1:ncpu+nboundary, 1:nlevelmax))
        close(amr_n)

        twotondim = 2**ndim
        skip_amr = 3 * (2**ndim + ndim) + 1

        ! Check total number of grids
        if(verbose>0)write(*,'(a)') '> Check total number of grids...'
        ncell_tot = 0
        ngridmax = 0
        if(verbose>1)write(6, '(a)', advance='no') '>> Progress: '
        do i = 1, SIZE(cpu_list)
            if(verbose>1)call progress_bar(i, SIZE(cpu_list))
            icpu = cpu_list(i)
            open(unit=amr_n, file=amr_filename(repo, iout, icpu), status='old', form='unformatted')
            call skip_read(amr_n, 21)
            ! Read grid numbers
            read(amr_n) ngridfile(1:ncpu,1:nlevelmax)
            ngridmax=maxval(ngridfile)!!!
            call skip_read(amr_n, 7)
            ! For non-periodic boundary conditions (not tested!)
            if(nboundary>0) then
                call skip_read(amr_n, 3)
            endif
            allocate(son(1:ngridmax))
            do ilevel = 1, nlevelmax
                ngrid_a = ngridfile(icpu, ilevel)
                ! Loop over domains
                do jcpu = 1, nboundary + ncpu
                    if(ngridfile(jcpu, ilevel) > 0) then
                        call skip_read(amr_n, 3)
                        ! Read grid center
                        if(jcpu == icpu) then
                            call skip_read(amr_n, ndim) ! Skip position
                            read(amr_n) ! Skip father index
                            call skip_read(amr_n, 2*ndim) ! Skip nbor index
                            ! Read son index to check refinement
                            do jdim = 1, twotondim
                                read(amr_n) son(1:ngrid_a)
                                do igrid=1, ngrid_a
                                    if(son(igrid) == 0) ncell_tot = ncell_tot+1
                                end do
                            end do
                            call skip_read(amr_n, 2*twotondim) ! Skip cpu, refinement map
                        else
                            call skip_read(amr_n, skip_amr)
                        end if
                    end if
                end do
            end do
            deallocate(son)
            close(amr_n)
        end do
        if(verbose>0)write(*,'(a)') '> Done'
    end subroutine read_cell

! !#####################################################################
!     subroutine sfrcell(cx, cy, cz, cdx, sage, sx, sy, sz, sm, timewindow)
! !#####################################################################
!         implicit none

!         real(kind=8), dimension(:), intent(in) :: cx, cy, cz, cdx
!         real(kind=8), dimension(:), intent(in) :: sx, sy, sz, sm, sage
!         real(kind=8), intent(in) :: timewindow

!         integer(kind=4) :: i, j, nmatches
!         real(kind=8) :: cxmin, cxmax, cymin, cymax, czmin, czmax, ysm_sum
!         logical :: condx, condy, condz, cond
!         integer(kind=4), dimension(size(sx)) :: cell_indices, ystar_indices
!         real(kind=8), dimension(size(sx)) :: ysm
!         integer(kind=4), dimension(size(sx)) :: match_indices

!         ! Loop over cells
!         if(allocated(sfrs)) deallocate(sfrs)
!         allocate(sfrs(size(cx)))
!         sfrs = 0.0d0
!         do i = 1, size(cx)
!             cxmin = cx(i) - cdx(i)/2
!             cxmax = cx(i) + cdx(i)/2
!             cymin = cy(i) - cdx(i)/2
!             cymax = cy(i) + cdx(i)/2
!             czmin = cz(i) - cdx(i)/2
!             czmax = cz(i) + cdx(i)/2

!             ! Loop over young stars within time window
!             nmatches = 0
!             ysm_sum = 0.0d0
!             do j = 1, size(sx)
!                 if (sage(j) < timewindow) then
!                     condx = (cxmin <= sx(j)) .and. (cxmax > sx(j))
!                     condy = (cymin <= sy(j)) .and. (cymax > sy(j))
!                     condz = (czmin <= sz(j)) .and. (czmax > sz(j))
!                     cond = condx .and. condy .and. condz
!                     if (cond) then
!                         nmatches = nmatches + 1
!                         ystar_indices(nmatches) = j
!                         ysm(nmatches) = sm(j) / 1d8
!                     endif
!                 endif
!             end do
            
!             ! Compute the SFR for this cell
!             if (nmatches > 0) then
!                 cell_indices(1:nmatches) = i
!                 match_indices(1:nmatches) = ystar_indices(1:nmatches)
!                 sfrs(i) = sum(ysm(1:nmatches))
!             endif
!         end do

!     end subroutine sfrcell

!#####################################################################
    subroutine sfrincell(cx, cy, cz, cdx, sage, sx, sy, sz, sm, timewindow, nthread)
!#####################################################################
        use omp_lib
        implicit none

        real(kind=8), dimension(:), intent(in) :: cx, cy, cz, cdx
        real(kind=8), dimension(:), intent(in) :: sx, sy, sz, sm, sage
        real(kind=8), intent(in) :: timewindow

        integer(kind=4) :: i, j, nmatches, nthread
        real(kind=8) :: cxmin, cxmax, cymin, cymax, czmin, czmax, ysm_sum
        logical :: condx, condy, condz, cond
        integer(kind=4), dimension(size(sx)) :: cell_indices, ystar_indices
        real(kind=8), dimension(size(sx)) :: ysm
        integer(kind=4), dimension(size(sx)) :: match_indices

        ! Loop over cells
        if(allocated(sfrs)) deallocate(sfrs)
        allocate(sfrs(size(cx)))
        sfrs = 0.0d0
        !$OMP PARALLEL DO SHARED(cx,cy,cz,cdx,sage,sx,sy,sz,sm,timewindow,sfrs) &
        !$OMP PRIVATE(i,j,cxmin,cxmax,cymin,cymax,czmin,czmax,condx,condy,condz,cond) &
        !$OMP PRIVATE(nmatches,ystar_indices,ysm_sum,cell_indices,match_indices) &
        !$OMP SCHEDULE(dynamic) &
        !$OMP NUM_THREADS(nthread)
        do i = 1, size(cx)
            cxmin = cx(i) - cdx(i)/2
            cxmax = cx(i) + cdx(i)/2
            cymin = cy(i) - cdx(i)/2
            cymax = cy(i) + cdx(i)/2
            czmin = cz(i) - cdx(i)/2
            czmax = cz(i) + cdx(i)/2

            ! Loop over young stars within time window
            nmatches = 0
            ysm_sum = 0.0d0
            do j = 1, size(sx)
                if (sage(j) < timewindow) then
                    condx = (cxmin <= sx(j)) .and. (cxmax > sx(j))
                    condy = (cymin <= sy(j)) .and. (cymax > sy(j))
                    condz = (czmin <= sz(j)) .and. (czmax > sz(j))
                    cond = condx .and. condy .and. condz
                    if (cond) then
                        nmatches = nmatches + 1
                        ystar_indices(nmatches) = j
                        ysm(nmatches) = sm(j) / 1d8
                    endif
                endif
            end do
            
            ! Compute the SFR for this cell
            if (nmatches > 0) then
                cell_indices(1:nmatches) = i
                match_indices(1:nmatches) = ystar_indices(1:nmatches)
                sfrs(i) = sum(ysm(1:nmatches))
            endif
        end do
        !$OMP END PARALLEL DO

    end subroutine sfrincell

! !#####################################################################
!     subroutine sfrcell_omp(cx, cy, cz, cdx, sage, sx, sy, sz, sm, timewindow)
! !#####################################################################
!         use omp_lib
!         implicit none

!         real(kind=8), dimension(:), intent(in) :: cx, cy, cz, cdx
!         real(kind=8), dimension(:), intent(in) :: sx, sy, sz, sm, sage
!         real(kind=8), intent(in) :: timewindow

!         integer(kind=4)      :: i, j, nmatches
!         real(kind=8) :: cxmin, cxmax, cymin, cymax, czmin, czmax
!         logical :: condx, condy, condz, cond
!         integer(kind=4), dimension(:),   allocatable :: cell_indices, ystar_indices
!         real(kind=8), dimension(:),   allocatable    :: ysm

!         ! Allocate temporary array for indices of matching young stars
!         integer(kind=4), dimension(:), allocatable :: match_indices
!         integer(kind=4), dimension(:), allocatable :: star_occupied

!         if(allocated(sfrs)) deallocate(sfrs)
!         allocate(sfrs(1:SIZE(cx)))
!         allocate(ysm(1:SIZE(sx)))
!         allocate(cell_indices(1:SIZE(sx)))
!         allocate(ystar_indices(1:SIZE(sx)))
!         allocate(star_occupied(1:SIZE(sx)))
!         do j = 1, SIZE(sx)
!             star_occupied(j) = 0
!             ystar_indices(j) = 0
!             cell_indices(j) = 0
!             ysm(j) = 0.d0
!         end do

!         ! Loop over cells
!         CALL OMP_SET_NUM_THREADS(32)
!         !$OMP PARALLEL DEFAULT(NONE) &
!         !$OMP SHARED(cx,cy,cz,cdx,sage,sx,sy,sz,sm,timewindow,sfrs) &
!         !$OMP PRIVATE(i,j,cxmin,cxmax,cymin,cymax,czmin,czmax,condx,condy,condz,cond,match_indices) &
!         !$OMP DO SCHEDULE(static)
!         do i = 1, SIZE(cx)
!             cxmin = cx(i) - cdx(i)/2
!             cxmax = cx(i) + cdx(i)/2
!             cymin = cy(i) - cdx(i)/2
!             cymax = cy(i) + cdx(i)/2
!             czmin = cz(i) - cdx(i)/2
!             czmax = cz(i) + cdx(i)/2

!             ! Loop over young stars
!             sfrs(i) = 0d0
!             nmatches = 0
!             do j = 1, SIZE(sx)
!                 if (sage(j) < timewindow) then
!                     if (star_occupied(j) < 1) then
!                         condx = (cxmin <= sx(j)) .and. (cxmax > sx(j))
!                         condy = (cymin <= sy(j)) .and. (cymax > sy(j))
!                         condz = (czmin <= sz(j)) .and. (czmax > sz(j))
!                         cond = condx .and. condy .and. condz
!                         if (cond) then
!                             nmatches = nmatches + 1
!                             cell_indices(j) = i
!                             ystar_indices(j) = j
!                             ysm(j) = sm(j)/1d8
!                             star_occupied(j) = 1
!                         endif
!                     endif
!                 endif
!             end do
!             ! Compute the SFR for this cell
!                 if (nmatches > 0) then
!                     cell_indices(1:nmatches) = i
!                     match_indices(1:nmatches) = ystar_indices(1:nmatches)
!                     sfrs(i) = sum(ysm(1:nmatches))
!                 endif
!         end do
!         !$OMP END PARALLEL DO
!     end subroutine sfrcell_omp

! USE omp_lib
! call omp_get_num_threads()
!#####################################################################
    subroutine close()
!#####################################################################
        ! Clean old data table is used more then once
        implicit none
        if(allocated(real_table)) deallocate(real_table)
        if(allocated(integer_table)) deallocate(integer_table)
        if(allocated(byte_table)) deallocate(byte_table)
        if(allocated(long_table)) deallocate(long_table)
        if(allocated(sfrs)) deallocate(sfrs)
    end subroutine close

!#####################################################################
    subroutine skip_read(unit,nskip)
!#####################################################################
        ! skip the given number of reads

        implicit none
        integer,intent(in) :: unit, nskip
        integer :: i
        do i=1,nskip
            read(unit)
        end do
    end subroutine skip_read

!#####################################################################
    character(len=5) function charind(iout)
!#####################################################################
        implicit none
        integer, intent(in) :: iout

        write(charind, '(I0.5)') iout

    end function charind

!#####################################################################
    character(len=128) function part_filename(repo, iout, icpu)
!#####################################################################
        implicit none
        character(len=128), intent(in)  :: repo
        integer,            intent(in)  :: iout, icpu

        part_filename = TRIM(repo)//'/output_'//charind(iout)//'/part_'//charind(iout)//'.out'//charind(icpu)
    end function part_filename

!#####################################################################
    character(len=128) function amr_filename(repo, iout, icpu)
!#####################################################################
        implicit none
        character(len=128), intent(in)  :: repo
        integer,            intent(in)  :: iout, icpu

        amr_filename = TRIM(repo)//'/output_'//charind(iout)//'/amr_'//charind(iout)//'.out'//charind(icpu)
    end function amr_filename


!#####################################################################
    subroutine progress_bar(iteration, maximum)
!#####################################################################
        implicit none
        integer :: iteration, maximum
        if(iteration == maximum) then
            write(6, '(a)') '# - Done!'
        elseif(MOD(iteration, MAX(1, maximum/50)) == 0) then
            write(6, '(a)', advance='no') '#'
        end if
    end subroutine progress_bar

!#####################################################################
    subroutine get_npart(part_n, file_path, npart)
!#####################################################################
        ! Suggested by ChatGPT
        implicit none
        integer, intent(in) :: part_n
        character(len=*), intent(in) :: file_path
        integer, intent(out) :: npart
        open(unit=part_n, file=file_path, status='old', form='unformatted')
        call skip_read(part_n, 2)
        read(part_n) npart
        close(part_n)
        return
    end subroutine get_npart

end module readramses