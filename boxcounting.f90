! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! RAY-CHAOS
!
! Ray-chaos is a ray-tracer tailor-made for investigating classical dynamics in
! scattering systems, i.e. open billiards.
!
! boxcounting.f90 calculates the log(scale) and log(number-of-boxes) for calculating the fractal
! dimension of the scattering fractal
!
! Compile with: gfortran parameters.f90 shapes.f90 functions.f90 boxcounting.f90 -O3 -fopenmp -o rays.out
! Run with: ./rays.out
!
! Needs link to OpenMP
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
program bucketmodel
  use parameters
  use shapes
  use functions
  !$ use omp_lib
  implicit none

  integer(wp)               :: j, k, l, hours, minutes, seconds, PSOScount, survived, chunk
  integer(wp), allocatable  :: survivers(:), limits(:)
  real(wp)                  :: timerstart, timerstop, time, prob, cpl
  real(wp), dimension(2)    :: ray, cray
  real(sp), allocatable     :: matrix(:,:)
  character(len=90)         :: filename
  logical                   :: tir, abson
  
  write(*,FMT='(A12,I0,A1,I0)') " Matrix is: ", grid, "x", grid
  write(*,FMT='(A12,I0)')       " #rays    : ", grid*grid
  print *, "____________________________________________________________________________"
  print *, ""

  !$ timerstart = omp_get_wtime()

  ! File for saving survival rate
  open(unit = 2, file = "tsurvival.txt", form='formatted')
  write(2, FMT='(2A15)') "psoscount", "srate"

  ! Check grid for BCM. How many times is it divisble?
  call checkGrid(grid)
  
  ! Define start for ray
  if (initloc .eq. "B ") then
     cray(2) = 0.0001_wp
     ray = [tan(inang*PI/180._wp),1.0_wp]
     abson = .true.
  elseif (initloc .eq. "C ") then
     print *, "Rays are starting at the top. Check the parameters."
     call exit(1)
  end if
  prob = 1._wp
  tir  = .false.

  ! Set potential
  skew = oldskew ! For varying shape
  skew = n1
  call setPotential(potentialtype)

  ! Allocate system
  allocate(matrix(grid,grid))
  matrix = 0._wp

  print *, "Generating the system..."
  survived = 0

  ! This is only if you want to set chunk size for OMP dynamic schedule
  chunk = grid/16

  ! Setup survival counting system
  allocate(limits(boxcntlim+1))
  allocate(survivers(boxcntlim+1))
  do l = 1, boxcntlim+1
     limits(l) = l-1
  end do
  survivers = 0
  survived = 0

  !$OMP parallel ! num_threads(1)
  !$ write(*,FMT='(A4,I2,A13)') "Red ", omp_get_thread_num(), " standing by"
  !$OMP do private (PSOScount,PSOSrays,erays), firstprivate(ray,cray,tir), reduction(+:nextid,survived,deserter), &
  !$OMP & schedule(dynamic)
  do k = 1, grid
     ! Print progress
     if ((mod(k,grid/prog) .eq. 0) .and. (grid .ge. prog)) then
        print *, "Progress: ", (k*100)/grid, "%"
     end if
     ray(1) = tan((iang+(k-1)*(fang-iang)/(grid-1))*PI/180._wp) ! Set angle

     do j = 1, grid
        cray(1) = xi + (j-1)*(xf-xi)/(grid-1) ! Set position
        
        if (n1 .lt. n2) then
           ! skip calculation of bouncing ball orbit (TIR condition at n1/n2 interface)
           if (abs(iang+(k-1)*(fang-iang)/(grid-1)) .gt. asin(n1/n2)*180._wp/PI) then
              PSOScount = boxcntlim
              call checkForSurvivors(PSOScount,limits,survivers)
           else
              call createInsideRay(ray,cray,prob,nextid,initn,initloc,PSOScount,tir,abson,cpl)
              call checkForSurvivors(PSOScount,limits,survivers)
           end if
        else
           call createInsideRay(ray,cray,prob,nextid,initn,initloc,PSOScount,tir,abson,cpl)
           call checkForSurvivors(PSOScount,limits,survivers)
        end if
        
        ! Only fractal part or whole bouncespace?
        if (PSOScount .ge. boxcntlim) then
           matrix(j,k) = boxcntlim!PSOScount
           survived = survived+1
        end if
        matrix(j,k) = PSOScount
        nextid = nextid+1
     end do
  end do
  !$OMP end do
  !$OMP end parallel

  ! Storing the survival data
  do l = 1, boxcntlim+1
     write(2,*) limits(l), dble(survivers(l))/(grid**2)
  end do

  print *, "Doing boxcounting..."
  ! Optional save the fractal
  if (savefractal .eqv. .true.) then
     ! Open binary file to store shape
     open(unit = 2, file = "tfractal.dat", form='unformatted', access='stream')
     print *, "Write .dat..."
     matrix = transpose(matrix)
     do j = 1, grid
        write(2) real(matrix(j,:))
     end do
     close(2)

     ! Store as .txt
     open(unit = 2, file = "tfractal.txt", form='formatted')
     print *, "Write .txt..."
     do k = 1, grid
        ray(1) = (iang+(k-1)*(fang-iang)/(grid-1))
        do j = 1, grid
           cray(1) = xi + (j-1)*(xf-xi)/(grid-1)
           write(2,*) k, j, matrix(k,j), cray(1), ray(1)
        end do
     end do
     close(2)
  end if

  ! BCM
  call boxCountingMethod(matrix,grid)

  call deallocatePotential(potential,point,edgecoords,edgeid,potsid,rileft,riright)
  deallocate(matrix)

  !$ timerstop = omp_get_wtime()

  ! Write useful info to screen
  time    = timerstop-timerstart
  hours   = int8(time/3600)
  minutes = int8((time-hours*3600)/60)
  seconds = int8(time-hours*3600-minutes*60)
  ! Write log entry
  call logparameters(hours,minutes,seconds,"box ")
  print *, "==============================================================="
  print *, "Run time:               ", timerstop-timerstart
  print *, "==============================================================="
  write(*,FMT='(T2,A29,I0)') "Number of rays sent in:      ", nextid-100
  write(*,FMT='(T2,A29,I0)') "Number of rays exiting -> C: ", erays
  write(*,FMT='(T2,A29,I0)') "Number of bounces on PSOS:   ", PSOSrays
  write(*,FMT='(T2,A29,I0)') "Box counting limit:          ", boxcntlim
  write(*,FMT='(T2,A29,I0)') "Number of surviviors:        ", survived
  write(*,FMT='(T2,A29,I0)') "Number of deserters:         ", deserter
  write(*,FMT='(T2,A6,A2)') "PSOS: ", PSOS
  write(*,FMT='(T2,A14,L1)') "Raysplitting: ", raysplitting
  write(*,FMT='(T2,A18,I0)') "Min PSOS bounces: ", PSOSlim
  write(*,FMT='(T2,A10)', advance="no") "Run time: "
  write(*,FMT='(I2,A6)') hours, " hours"
  write(*,FMT='(T12,I2,A8)') minutes, " minutes"
  write(*,FMT='(T12,I2,A8)') seconds, " seconds"
  print *, "==============================================================="
end program bucketmodel
