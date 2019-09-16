! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! RAY-CHAOS
!
! Ray-chaos is a ray-tracer tailor-made for investigating classical dynamics in
! scattering systems, i.e. open billiards.
!
! buckettrace.f90 calculates the average path length of rays in a scatterer
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
  ! "!$" is the notation used by OpenMP
  !$ use omp_lib
  implicit none

  integer(wp)               :: j, hours, minutes, seconds, rayinitialid, k, l
  real(wp)                  :: timerstart, timerstop, time, prob, cpl
  real(wp), dimension(2)    :: ray, cray
  real(wp), allocatable     :: veff(:,:), vcpl(:,:)
  logical                   :: tir, abson
  character(len=90)         :: filename
  
  ! Seed RNG
  call random_seed

  ! Define start for ray
  if (initloc .eq. "B ") then
     cray(2) = 0.0001_wp
     ray = [tan(inang*PI/180._wp),1.0_wp]
     abson = .true.
  elseif (initloc .eq. "C ") then
     cray(2) = height-0.0001_wp
     ray = [tan(inang*PI/180._wp),-1.0_wp]
     abson = .false.
  end if
  prob = 1._wp
  tir  = .false.
  l = 1 ! A counter for the ray start grid
  k = 1

  ! Open a file to store scattering coordinates
  open(unit = pid, file = "tphasespace.txt", form='formatted')
  write(pid, FMT='(2A15)') "x", "sin(theta)"

  ! Open a file to store optical path length
  open(unit = oplid, file = "topl.txt", form='formatted')
  write(oplid, FMT='(2A15)') "x", "opl"

  ! Open a file to store conversion path length
  open(unit = cplid, file = "tcpl.txt", form='formatted')
  write(cplid, FMT='(2A15)') "x", "cpl"

  ! Open a file to store efficiency
  open(unit = effid, file = "tefficiency.txt", form='formatted')
  write(effid, FMT='(2A15)') "x", "eff"

  print *, "Initial coord: ", cray
  print *, "____________________________________________________________________________"
  print *, ""
  rayinitialid = 1

  skew = n1
  ! Initialize potential
  call setPotential(potentialtype)

  print *, "Tracing rays..."
  !$ timerstart = omp_get_wtime()
  do j = 1, nrays
     probloss = 0._wp
     cpl = 0._wp
     truncray = 0._wp

     if ((nrays .ne. 1) .and. (startagrid .eqv. .false.)) then
        cray(1) = xi + (j-1)*(xf-xi)/(nrays-1)

        if ((mod(j,nrays/prog) .eq. 0) .and. (nrays .ge. prog)) then
           print *, "Progress: ", (j*100)/nrays, "%"
        end if
     elseif ((nrays .eq. 1) .and. (startagrid .eqv. .false.)) then
        ! Set manually
        cray(1) = xi + (1600-1)*(xf-xi)/(grid-1) ! Set s coordinate

        print *, "Initial direction ray: ", ray
        print *, "Initial coordinate:    ", cray
        print *, "inang: ", inang
     elseif ((nrays .ne. 1) .and. (startagrid .eqv. .true.)) then
        ! Start rays on a set grid. nrays must be square
        if (abs(sqrt(dble(nrays))-nint(sqrt(dble(nrays)))) .ge. eps) then
           print *, "Please choose a square number"
           call exit(1)
        end if

        if ((mod(j-1,nint(sqrt(dble(nrays)))) .eq. 0) .and. (j .ne. 1)) then
           k = 1
           l = l+1
           print *, "j switch: ", j
        end if
        cray(1) = xi + (k-1)*(xf-xi)/(sqrt(dble(nrays))-1)
        ray(1)  = tan((iang+(l-1)*(fang-iang)/(sqrt(dble(nrays))-1))*PI/180._wp)
        k = k+1

        ! Open a file to store scattering coordinates
        close(pid) ! close old file
        write(filename, '("tphasespace-", I0, ".txt")') nextid
        open(unit = pid, file = filename, form='formatted')
        write(pid, FMT='(2A15)') "x", "sin(theta)"
     end if

     call createRay(ray,cray,prob,nextid,initn,initloc,tir,abson,cpl)
     nextid = nextid+1

     ! write efficiency
     write(effid,*) cray(1), probloss, truncray
     ! write CPL
     write(cplid,*) cray(1), cpl
     
  end do
  !$ timerstop = omp_get_wtime()

  ! Deallocate
  call deallocatePotential(potential,point,edgecoords,edgeid,potsid,rileft,riright)

  ! Plot OPL
  if ((raysplitting .eqv. .false.) .and. (nrays .gt. 1)) then
     call system('gnuplot gopl.gp')
     call system('gnuplot gcpl.gp')
     print *, "OPL and CPL is plotted"
  end if

  ! Read efficiency data from file
  rewind(effid)
  read(effid,*) ! Reads descriptors/file head
  ! Read cpl data from file
  rewind(cplid)
  read(cplid,*) ! Reads descriptors/file head
  
  allocate(veff(nrays,2))
  allocate(vcpl(nrays,2))
  do j = 1, nrays
     read(effid, *) veff(j,:)
     read(cplid, *) vcpl(j,:)
  end do

  close(effid)
  close(pid)
  close(oplid)
  close(cplid)

  ! Write useful info to screen
  time    = timerstop-timerstart
  hours   = int8(time/3600)
  minutes = int8((time-hours*3600)/60)
  seconds = int8(time-hours*3600-minutes*60)
  print *, "==============================================================="
  print *, "Run time:               ", timerstop-timerstart
  print *, "==============================================================="
  write(*,FMT='(T2,A29,F6.4)') "Total efficiency (avg):      ", sum(veff(:,2))/nrays
  write(*,FMT='(T2,A,F8.6)') "Loss due to truncation:      ", truncloss/nrays
  write(*,FMT='(T2,A,F10.4)') "Mean pathlength (microns):   ", sum(vcpl(:,2))/nrays
  write(*,FMT='(T2,A29,I0)') "Number of rays sent in:      ", nextid-100
  write(*,FMT='(T2,A29,I0)') "Number of rays exiting -> C: ", erays
  write(*,FMT='(T2,A29,I0)') "Number of bounces on PSOS:   ", PSOSrays
  write(*,FMT='(T2,A6,A2)') "PSOS: ", PSOS
  write(*,FMT='(T2,A14,L1)') "Raysplitting: ", raysplitting
  write(*,FMT='(T2,A18,I0)') "Min PSOS bounces: ", PSOSlim
  write(*,FMT='(T2,A10)', advance="no") "Run time: "
  write(*,FMT='(I2,A6)') hours, " hours"
  write(*,FMT='(T12,I2,A8)') minutes, " minutes"
  write(*,FMT='(T12,I2,A8)') seconds, " seconds"
  print *, "==============================================================="
end program bucketmodel
