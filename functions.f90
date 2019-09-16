module functions
  use parameters
  use shapes
  !$ use omp_lib
  implicit none

contains

  ! Create a ray
  recursive subroutine createRay(iray,icray,iprob,nextid,ic,ibound,tir,iabson,cpl)
    real(wp), dimension(2), intent(in)       :: iray, icray
    real(wp), dimension(2)                   :: isect, ray, cray, temp
    real(wp), intent(in)                     :: iprob, ic
    real(wp), intent(inout)                  :: cpl
    real(wp)                                 :: prob, ncurr, nold, opl, seglen, segprobloss, nnext
    integer(wp), dimension(1)                :: listel
    integer(wp), intent(in)                  :: nextid
    integer(wp)                              :: rayid, bnum, PSOScount, j
    character(len=90)                        :: filename
    character(len=4), intent(in)             :: ibound
    character(len=4)                         :: boundary
    logical                                  :: tir, abson, iabson, nodetect

    ! print *, "hello", nextid
    rayid = nextid
    bnum = 0
    nodetect = .false.
    ray  = iray
    cray = icray
    abson = iabson
    opl  = 0
    prob = iprob
    ncurr = ic
    PSOScount = 0
    boundary = ibound

    if ((trackrays .eqv. .true.) .and. (nextid .lt. savepaths) .or. (nextid .eq. specray)) then
       ! Open a file to store scattering coordinates
       write(filename, '("tscatteringcoords-no-", I0, ".txt")') rayid
       open(unit = rayid, file = filename, form='formatted')
       write(rayid, FMT='(5A15)') "x", "y", "z", "prob"!,  "rayid" , "coord"
       write(rayid,*) cray(1), cray(2), 100, prob!, rayid, int8(1)
    end if

    do while (prob .gt. problim)
       nold = ncurr ! Update index of refraction

       call detectHit(ray,cray,boundary,isect,ncurr,tir,nodetect)

       if (nodetect .eqv. .true.) then
          exit
       end if

       ! Find the index of refraction which the ray is heading into
       if (any(bucketid .eq. boundary) .eqv. .false.) then
          listel = pack([(j,j=1,pots)],potsid .eq. boundary) ! Go from boundary (char) to integer
          call findRefInd(ncurr,nnext,listel(1))
       end if
       temp = cray
       cray = isect

       ! Calculate segment length
       seglen = sqrt((temp(1)-cray(1))**2+(temp(2)-cray(2))**2)

       ! Calculate probability loss for segment
       if ((abson .eqv. .true.) .and. (absorption .eqv. .true.)) then
          ! cpl = cpl+seglen*nold*prob ! Use this if you want CPL weighted
          segprobloss = exp(-seglen*4._wp*PI*ni/lambda) ! Fraction lost travelling the segment
          probloss = probloss+prob*(1-segprobloss) ! Accumulated loss along all absorbing segments
          prob = prob*segprobloss ! Reduce probability accordingly to segment loss
       end if

       ! Calculate OPL and GPL (cpl is the old variable name, GPL = Geometric Path Length)
       if ((nold .ne. 1._wp) .and. (optical .eqv. .true.) .and. (weighted .eqv. .false.)) then
          opl = opl+seglen*nold
          cpl = cpl+seglen*nold ! Use this if you want OPL
       elseif ((nold .ne. 1._wp) .and. (optical .eqv. .false.) .and. (weighted .eqv. .false.)) then
          opl = opl+seglen
          cpl = cpl+seglen
       elseif ((nold .ne. 1._wp) .and. (optical .eqv. .true.) .and. (weighted .eqv. .true.)) then
          opl = opl+seglen*nold
          cpl = cpl+seglen*nold*prob
       elseif ((nold .ne. 1._wp) .and. (optical .eqv. .false.) .and. (weighted .eqv. .true.)) then
          opl = opl+seglen
          cpl = cpl+seglen*prob
       end if

       ! Scatter the ray
       call scatterRay(ray,boundary,cray,prob,ncurr,tir,abson,cpl,nnext,listel)

       if (boundary .eq. PSOS) then
          PSOScount = PSOScount+1
       end if

       ! Add to phasespace
       if ((boundary .eq. PSOS) .and. (PSOScount .gt. PSOSlim) .and. (storephasspc .eqv. .true.)) then
          call addPhase(ray,isect)
       end if

       ! Turn absorption on/off
       if (any(abslist .eq. boundary) .and. (tir .eqv. .false.)) then
          if (abson .eqv. .true.) then
             abson = .false.
          else
             abson = .true.
          end if
       end if
       
       ! Write trajectory and probability to file
       if ((trackrays .eqv. .true.) .and. (rayid .lt. savepaths) .or. (nextid .eq. specray)) then
          write(rayid,*) cray(1), cray(2), 100, prob!, rayid, bnum
       end if
       bnum = bnum+1

       ! Exit if at ceiling
       if (boundary .eq. "C") then
          erays = erays+1
          exit
       elseif (((boundary .eq. "L") .or. (boundary .eq. "R")) .and. (periodicbcs .eqv. .true.)) then
          call periodicShift(boundary,cray)
          if ((trackrays .eqv. .true.) .and. (rayid .lt. savepaths) .or. (nextid .eq. specray)) then
             write(rayid,*) cray(1), cray(2), 100, prob!, rayid, bnum
          end if
          bnum = bnum+1
       end if

       ! Exit if too many bounces
       if (bnum .gt. bouncelim) then
          print *, "Ray exeeded bounce limit. Too many bounces."
          exit
       end if

    end do

    if (boundary .ne. "C") then
       truncloss = truncloss+prob
       truncray = truncray+prob
    end if

    ! write OPL
    if ((raysplitting .eqv. .false.)) then
       write(oplid,*) icray(1), opl
    end if

    if ((trackrays .eqv. .true.) .and. (nextid .lt. savepaths) .or. (nextid .eq. specray)) then
       close(rayid)
    end if
  end subroutine createRay

  ! Create an inside ray
  recursive subroutine createInsideRay(iray,icray,iprob,nextid,ic,ibound,PSOScount,tir,iabson,cpl)
    real(wp), dimension(:), intent(in)       :: iray, icray
    real(wp), dimension(2)                   :: isect, ray, cray, temp
    real(wp), intent(in)                     :: iprob, ic
    real(wp), intent(inout)                  :: cpl
    real(wp)                                 :: prob, ncurr, nold, opl, nnext
    integer(wp), intent(in)                  :: nextid
    integer(wp), intent(inout)               :: PSOScount
    integer(wp), dimension(1)                :: listel
    integer(wp)                              :: rayid, bnum, j
    character(len=90)                        :: filename
    character(len=4), intent(in)             :: ibound
    character(len=4)                         :: boundary
    logical                                  :: tir, abson, iabson, nodetect

    ! print *, "hello", nextid
    rayid = nextid
    bnum = 0
    nodetect = .false.
    ray  = iray
    cray = icray
    abson = iabson
    opl  = 0
    prob = iprob
    ncurr = ic
    PSOScount = 0
    boundary = ibound

    do while ((PSOScount .le. boxcntlim))
       prob = 1._wp ! To keep other functions happy
       nold = ncurr
       call detectHit(ray,cray,boundary,isect,ncurr,tir,nodetect)

       if (nodetect .eqv. .true.) then
          exit
       end if

       ! Find the index of refraction which the ray is heading into
       if (any(bucketid .eq. boundary) .eqv. .false.) then
          listel = pack([(j,j=1,pots)],potsid .eq. boundary) ! Go from boundary (char) to integer
          call findRefInd(ncurr,nnext,listel(1))
       end if
       
       temp = cray
       cray = isect
       
       call scatterRay(ray,boundary,cray,prob,ncurr,tir,abson,cpl,nnext,listel)

       bnum = bnum+1
       
       if (boundary .eq. PSOS) then
          PSOScount = PSOScount+1 ! Per ray
          !$OMP atomic
          PSOSrays = PSOSrays+1 ! Global
          !$OMP end atomic
       end if

       ! Exit if at ceiling
       if (boundary .eq. "C") then
          !$OMP atomic
          erays = erays+1 ! Global
          !$OMP end atomic
          if (ncurr .ne. n0) then
             !$OMP atomic
             deserter = deserter+1 ! Global
             !$OMP end atomic
          end if
          exit
       end if
    end do

  end subroutine createInsideRay

  ! Scatter ray depending on the boundary
  subroutine scatterRay(ray,boundary,isect,prob,ncurr,tir,abson,cpl,nnext,listel)
    real(wp), dimension(:), intent(inout)    :: ray
    real(wp), dimension(:), intent(inout)    :: isect
    real(wp), intent(inout)                  :: prob, ncurr, cpl
    real(wp), intent(in)                     :: nnext
    integer(wp), dimension(1), intent(in)    :: listel
    character(len=4), intent(inout)          :: boundary
    logical                                  :: tir, abson

    if (boundary .eq. "R") then
       if (periodicbcs .eqv. .false.) then
          ray  = ray - 2._wp*dot_product(ray,rwall)*rwall
          tir = .false.
       end if
    elseif (boundary .eq. "L") then
       if (periodicbcs .eqv. .false.) then
          ray  = ray - 2._wp*dot_product(ray,lwall)*lwall
          tir = .false.
       end if
    elseif (boundary .eq. "B") then
       ray  = ray - 2._wp*dot_product(ray,bwall)*bwall
       tir = .false.
       ! call randomMirror(ray,bwall)
    elseif (boundary .eq. "C") then
       ray  = ray - 2._wp*dot_product(ray,cwall)*cwall
       tir = .false.
    else
       call refractray(ray,boundary,isect,prob,ncurr,tir,abson,cpl,nnext,listel)
    end if
  end subroutine scatterRay

  ! Add coordinate to phase space
  subroutine addPhase(ray,isect)
    real(wp), dimension(:), intent(in)    :: ray, isect
    real(wp), dimension(2)                :: normal
    real(wp)                              :: theta

    !! normal must correspond to the normal of PSOS
    ! PSOSnorm is a global variable defined in parameters.f90
    normal = PSOSnorm

    ! Angle between ray and norm of PSOS
    if (dot_product(ray,normal)/(norm2(ray)*norm2(normal)) .gt. 1._wp) then
       theta = acos(1._wp)
    elseif (dot_product(ray,normal)/(norm2(ray)*norm2(normal)) .lt. -1._wp) then
       theta = acos(-1._wp)
    else
       theta = acos(dot_product(ray,normal)/(norm2(ray)*norm2(normal)))
    end if

    write(pid,*) isect(1), sin(theta)
    PSOSrays = PSOSrays+1
  end subroutine addPhase

  ! Shift ray if periodic boundary conditions are needed
  subroutine periodicShift(boundary,cray)
    real(wp), dimension(:), intent(inout)    :: cray
    character(len=4), intent(inout)          :: boundary

    if (boundary .eq. "L") then
       cray(1) = width
       boundary = "R"
    elseif (boundary .eq. "R") then
       cray(1) = 0._wp
       boundary = "L"
    end if
  end subroutine periodicShift

  ! Function to refract ray with Snell's law or TIR
  subroutine refractRay(ray,boundary,isect,prob,ncurr,tir,abson,cpl,nnext,listel)
    real(wp), dimension(:), intent(inout)    :: ray
    real(wp), dimension(:), intent(in)       :: isect
    real(wp), dimension(2)                   :: normal, reflray
    integer(wp), dimension(1), intent(in)    :: listel
    real(wp), dimension(2,2)                 :: rotmat
    real(wp), intent(inout)                  :: prob, ncurr, cpl
    real(wp), intent(in)                     :: nnext
    real(wp)                                 :: theta, phi, rot, reflprob
    character(len=4)                         :: boundary
    logical                                  :: tir, abson

    ! Make a copy of ray which will be reflected for raysplitting
    reflray = ray

    ! Find normal
    if (scan(boundary, "X", kind=8) .eq. 0) then
       normal = [-potential(listel(1),2),potential(listel(1),1)]/&
            (sqrt(potential(listel(1),1)**2+potential(listel(1),2)**2))
    elseif (scan(boundary, "X", kind=8) .eq. 1) then
       call findConicNormal(normal,isect(1),isect(2),boundary)
    end if

    ! Make sure that argument for Snells is in the intervall [-1,1]
    if (dot_product(ray,normal)/(norm2(ray)*norm2(normal)) .gt. 1._wp) then
       theta = acos(dot_product(ray,normal)/(norm2(ray)*norm2(normal))-eps)
    elseif (dot_product(ray,normal)/(norm2(ray)*norm2(normal)) .lt. -1._wp) then
       theta = acos(dot_product(ray,normal)/(norm2(ray)*norm2(normal))+eps)
    else
       theta = acos(dot_product(ray,normal)/(norm2(ray)*norm2(normal)))
    end if

    ! Make sure that the incoming angle theta is correct
    ! and that the normal is pointing in the same general direction
    if (theta .gt. PI/2._wp) then
       theta = PI-theta
       normal = -normal
    end if

    ! This structure handles refraction (with or without ray splitting) and total internal reflection
    if ((ncurr/nnext)*sin(theta) .le. 1._wp+eps) then
       phi = asin((ncurr/nnext)*sin(theta))
       rot = abs(theta-phi) ! "rotation" needed for refraction

       ! Calculate 2D "cross product" to determine rotation direction
       if (((ray(1)*normal(2)-ray(2)*normal(1) .gt. 0._wp) .and. (ncurr .lt. nnext)) .or. &
            & ((ray(1)*normal(2)-ray(2)*normal(1) .lt. 0._wp) .and. (ncurr .gt. nnext))) then
          ! Rotate counterclockwise if CP is positive
          rotmat = rotate(rot)
          ray = matmul(rotmat,ray)
       elseif (((ray(1)*normal(2)-ray(2)*normal(1) .gt. 0._wp) .and. (ncurr .gt. nnext)) .or. &
            & ((ray(1)*normal(2)-ray(2)*normal(1) .lt. 0._wp) .and. (ncurr .lt. nnext))) then
          ! Rotate clockwise if CP is positive
          rotmat = rotate(-rot)
          ray = matmul(rotmat,ray)
       end if

       ! Call createRay to calculate a reflected ray
       reflprob = prob*rte(theta,nnext/ncurr) ! includes TE polarization only
       if ((reflprob .gt. problim) .and. (raysplitting .eqv. .true.)) then
          reflray  = reflray - 2._wp*dot_product(reflray,normal)*normal
          nextid = nextid+1
          call createRay(reflray,isect,reflprob,nextid,ncurr,boundary,tir,abson,cpl)
       end if

       ! Update prob
       prob = prob*tte(theta,nnext/ncurr)

       ! Switch index of refraction
       ncurr = nnext
       tir = .false.
       
    else ! Total internal reflection
       ray = ray - 2._wp*dot_product(ray,normal)*normal
       tir = .true.
    end if

  end subroutine refractRay

  ! Find new index of refraction
  subroutine findRefInd(ncurr,nnext,idnum)
    real(wp), intent(in)       :: ncurr
    real(wp), intent(out)      :: nnext
    integer(wp), intent(in)    :: idnum

    ! swap refractive index
    if (ncurr .eq. rileft(idnum)) then
       nnext = riright(idnum)
    else
       nnext = rileft(idnum)
    end if
  end subroutine findRefInd

  ! Box Counting Method
  subroutine boxCountingMethod(matrix,grid)
    real(sp), intent(in)      :: matrix(:,:)
    integer(wp), intent(in)   :: grid
    integer(wp)               :: boxcount, boxmesh, boxgrid, cc, cr, div, k, j
    real(wp)                  :: scale
    character(len=90)         :: filename

    ! Open a file to store ln(N) and ln(M)
    open(unit = 3, file = "tboxcount.txt", form='formatted')
    write(3, FMT='(3A15)') "M", "N"

    ! Open a file to store ln(N) and ln(M) 
    write(filename, '("tboxcount-gen-", I0, ".txt")') grid
    open(unit = 4, file = filename, form='formatted')
    write(4, FMT='(3A15)') "M", "N"

    boxgrid = grid
    boxmesh = 1
    scale   = 1
    div     = 2
    
    print *, "Count boxes"
    do while (boxmesh .lt. grid)
       ! Set the initial values
       boxcount = 0
       do while (mod(grid,div) .ne. 0)
          div = div+1
       end do
       
       scale = div
       boxmesh = div  ! Number of boxes in each dir
       boxgrid = grid/boxmesh ! Grid of 1 box (boxgrid*boxgrid)\
       
       if (boxgrid*boxmesh .ne. grid) then
          print *, "There is a mismatch between the grid and boxgrid"
          call exit(1)
       end if
       
       print *, "---------------------------------"
       print *, "boxmesh: ", boxmesh
       print *, "boxgrid: ", boxgrid
       print *, "scale:   ", scale

       ! BCM
       do k = 1, boxmesh
          cc = (k-1)*boxgrid
          do j = 1, boxmesh
             cr = (j-1)*boxgrid
             call loopSection(matrix,boxgrid,boxcount,cc,cr)
          end do
       end do

       ! Write number of boxes
       write(*,FMT='(A8,I0,A8)') "I found ", boxcount, " box(es)"
       if (boxcount .eq. 0) boxcount=1 ! For the log
       write(3,*) log(dble(scale)), log(dble(boxcount))
       write(4,*) log(dble(scale)), log(dble(boxcount))

       div = div+1
    end do
    close(3)
    close(4)
  end subroutine boxCountingMethod

  ! Check how many times it is possible to section a grid into smaller grids
  subroutine checkGrid(grid)
    integer(wp), intent(in)  :: grid
    integer(wp)              :: div, boxmesh, boxgrid, iter
    ! real(wp)                 :: scale

    boxgrid = grid
    boxmesh = 1
    ! scale   = 1
    div     = 2
    iter    = 0
    do while (boxmesh .lt. grid)
       do while (mod(grid,div) .ne. 0)
          div = div+1
       end do

       boxmesh = div  ! Number of boxes in each dir
       boxgrid = grid/boxmesh ! Grid of 1 box (boxgrid*boxgrid)\
       if (boxgrid*boxmesh .ne. grid) then
          print *, "There is a mismatch between the grid and boxgrid"
          print *, "div: ", div
          print *, "boxmesh: ", boxmesh
          print *, "boxgrid: ", boxgrid
          print *, "prod: ", boxmesh*boxgrid
          call exit(1)
       end if
       
       iter = iter+1
       div  = div+1
    end do
    
    print *, "--------------------------------------------"
    write(*,FMT='(A11,I0,A14,I0,A8)') "The number ", grid, " is divisible ", iter, " time(s)"
    print *, "--------------------------------------------"
  end subroutine checkGrid

  ! Loop over a section of a matrix to look for elements above boxcntlim
  subroutine loopSection(matrix,boxsize,boxcount,cc,cr)
    real(sp), intent(in)       :: matrix(:,:)
    integer(wp), intent(in)    :: boxsize, cc, cr
    integer(wp)                :: j, k
    integer(wp), intent(inout) :: boxcount
    logical                    :: foundone

    foundone = .false.
    
    ! Loop over section
    do k = 1+cc, boxsize+cc
       do j = 1+cr, boxsize+cr
          ! Check elements inside the box if there a part of the fractal there or not
          if (matrix(j,k) .ge. boxcntlim) then
             boxcount = boxcount+1
             foundone = .true.
             exit
          end if
       end do
       
       if (foundone .eqv. .true.) exit
    end do
  end subroutine loopSection

  ! Find coordinates of intersection (slower)
  pure subroutine intersection(ray,cray,wall,cwall,isect)
    real(wp), dimension(:), intent(in)      :: ray, cray, wall, cwall
    real(wp), dimension(:), intent(inout)   :: isect
    real(wp)                                :: det, t
    det = ray(2)*wall(1)-ray(1)*wall(2)

    if (det .ne. 0._wp) then
       t = (1._wp/det)*(-ray(2)*(cwall(1)-cray(1))+ray(1)*(cwall(2)-cray(2)))
       isect(1) = cwall(1)+t*wall(1)
       isect(2) = cwall(2)+t*wall(2)
    else
       isect = [2._wp*width,2._wp*height] ! Out of bounds
    end if
    
  end subroutine intersection

  ! Find coordinates of intersection (faster). Using Cramer's rule
  pure subroutine intersection2(rayx,rayy,crayx,crayy,wallx,wally,cwallx,cwally,isectx,isecty)
    real(wp), intent(in)      :: rayx, rayy, crayx, crayy, wallx, wally, cwallx, cwally
    real(wp), intent(inout)   :: isectx, isecty
    real(wp)                  :: det, t
    det = rayy*wallx-rayx*wally

    if (det .ne. 0._wp) then
       t = (1._wp/det)*(-rayy*(cwallx-crayx)+rayx*(cwally-crayy))
       isectx = cwallx+t*wallx
       isecty = cwally+t*wally
    else
       isectx = 2._wp*width
       isecty = 2._wp*height ! Out of bounds
    end if
  end subroutine intersection2

  ! Determines the distances to boundaries
  subroutine detectHit(ray,cray,lastboundary,isect,ncurr,tir,nodetect)
    integer(wp)                           :: k
    real(wp), dimension(2), intent(in)    :: ray, cray
    real(wp), dimension(2), intent(inout) :: isect
    real(wp), intent(in)                  :: ncurr
    real(wp), dimension(bndrs)            :: dist
    real(wp), dimension(bndrs,2)          :: isectList
    real(wp), dimension(bndrs,4)          :: distlist
    character(len=4), intent(inout)       :: lastboundary
    logical, intent(in)                   :: tir
    logical, intent(inout)                :: nodetect

    ! Find intersection with every boundary of bucket
    do k = 1, 4
       call intersection2(ray(1),ray(2),cray(1),cray(2),bucket(k,1),bucket(k,2),&
            edgecoords(k,1),edgecoords(k,2),isectList(k,1),isectList(k,2))
       dist(k) = distance(cray(1),cray(2),isectList(k,1),isectList(k,2))
    end do

    ! Find intersection with every boundary of potential
    do k = 5, bndrs
       ! Decide between straight and curved pieces
       if (scan(edgeid(k), "X", kind=8) .eq. 0) then
          call intersection2(ray(1),ray(2),cray(1),cray(2),potential(k-4,1),potential(k-4,2),&
               point(k-4,1),point(k-4,2),isectList(k,1),isectList(k,2))
          dist(k) = distance(cray(1),cray(2),isectList(k,1),isectList(k,2))
       elseif (scan(edgeid(k), "X", kind=8) .eq. 1) then
          call conicIntersection(ray(1),ray(2),cray(1),cray(2),isectList(k,1),isectList(k,2),edgeid(k))
          dist(k) = distance(cray(1),cray(2),isectList(k,1),isectList(k,2))
       end if
    end do

    ! Make distlist (could benefit from OOP)
    do k = 1, bndrs
       distlist(k,1) = dist(k)        ! Distance to intersection
       distlist(k,2) = isectList(k,1) ! x-coordinate of intersection
       distlist(k,3) = isectList(k,2) ! y-coordinate of intersection
       distlist(k,4) = real(k)        ! Identity (number corresponding to edgeid)
    end do

    ! Sort the system from shortest distance to intersection
    call insertionSort4cols(distlist,1)

    if (tir .eqv. .true.) then
       call pickIntersectionTIR(ray,cray,distlist,isect,ncurr,lastboundary)
    elseif (tir .eqv. .false.) then
       call pickIntersectionREF(ray,cray,distlist,isect,lastboundary,nodetect)
    end if

  end subroutine detectHit

  subroutine pickIntersectionREF(ray,cray,distlist,isect,lastboundary,nodetect)
    real(wp), dimension(2), intent(in)        :: ray, cray
    real(wp), dimension(2), intent(inout)     :: isect
    real(wp), dimension(bndrs,4), intent(in)  :: distlist
    real(wp)                                  :: sgnx, sgny
    integer(wp)                               :: j, k, numofpots
    character(len=4), intent(inout)           :: lastboundary
    character(len=4)                          :: el
    logical                                   :: detection
    logical, intent(inout)                    :: nodetect

    detection = .false.

    do k = 1, bndrs
       j = int(distlist(k,4))

       ! Skip missed boundaries
       if (scan(edgeid(j), "X", kind=8) .ne. 1) then
          if (isOnLine(edgecoords(j,1),edgecoords(j,2),edgecoords(j,3),edgecoords(j,4),&
               distlist(k,2),distlist(k,3)) .eqv. .false.) then
             cycle
          end if
       end if

       el = edgeid(j) ! Identity of the closest intersection

       ! Handle direction condition
       sgnx = sign(1._wp,distlist(k,2)-cray(1))
       sgny = sign(1._wp,distlist(k,3)-cray(2))
       ! Ray precision
       if ((ray(1) .lt. 0._wp+eps15) .and. (ray(1) .gt. 0._wp-eps15)) then
          sgnx = sign(1._wp,ray(1))
       end if
       if ((ray(2) .lt. 0._wp+eps15) .and. (ray(2) .gt. 0._wp-eps15)) then
          sgny = sign(1._wp,ray(2))
       end if

       ! The first value is very finicky
       ! It defines the minimum distance for multiple boundary hit
       if ((distlist(k,1) .gt. 4e-16_wp) .and. &
            (sign(1._wp,ray(1)) .eq. sgnx) .and. &
            (sign(1._wp,ray(2)) .eq. sgny) .and. &
            (distlist(k,2) .ge. 0._wp) .and. (distlist(k,3) .ge. 0._wp) .and. &
            (distlist(k,2) .le. width) .and. (distlist(k,3) .le. height)) then
          if (el .eq. lastboundary) then
             if ((scan(edgeid(j), "X", kind=8) .eq. 1) .and. (distlist(k,1) .lt. 2e-14_wp)) cycle
             if (scan(edgeid(j), "X", kind=8) .ne. 1) cycle
          end if
          lastboundary = el
          isect(1) = distlist(k,2)
          isect(2) = distlist(k,3)
          detection = .true.
          exit
       end if
    end do

    ! In case something went wrong
    if (detection .eqv. .false.) then
       print *, "No detection in pickIntersectionREF"
       print *, "Current ray is:   ", nextid
       print *, "Direction of ray: ", ray
       print *, "Current point     ", cray
       print *, "Current boundary: ", lastboundary
       
       do k = 1, bndrs
          print *, distlist(k,1), distlist(k,2), distlist(k,3), edgeid(int(distlist(k,4)))
       end do
       numofpots = size(potsid)
       if ((lastboundary .ne. "L") .and. (lastboundary .ne. "R") .and. &
            (lastboundary .ne. potsid(1)) .and. (lastboundary .ne. potsid(numofpots))) then
          print *, "???"
          call exit(1)
       end if
       
       nodetect = .true. ! This shall terminate the ray in the create functions
    end if
    
  end subroutine pickIntersectionREF

  subroutine pickIntersectionTIR(ray,cray,distlist,isect,ncurr,lastboundary)
    real(wp), dimension(2), intent(in)        :: ray, cray
    real(wp), dimension(2), intent(inout)     :: isect
    real(wp), dimension(bndrs,4), intent(in)  :: distlist
    real(wp), intent(in)                      :: ncurr
    real(wp)                                  :: sgnx, sgny, nalt
    integer(wp), dimension(1)                 :: listel
    integer(wp)                               :: j, k
    character(len=4), intent(inout)           :: lastboundary
    character(len=4)                          :: el
    logical                                   :: detection

    detection = .false.

    ! Skip missed boundaries
    do k = 1, bndrs
       j = int(distlist(k,4))
       if (scan(edgeid(j), "X", kind=8) .ne. 1) then
          if (isOnLine(edgecoords(j,1),edgecoords(j,2),edgecoords(j,3),edgecoords(j,4),&
               distlist(k,2),distlist(k,3)) .eqv. .false.) then
             cycle
          end if
       end if

       el = edgeid(j) ! Identity of the closest intersection

       ! Handle direction condition
       sgnx = sign(1._wp,distlist(k,2)-cray(1))
       sgny = sign(1._wp,distlist(k,3)-cray(2))
       ! Ray precision
       if ((ray(1) .lt. 0._wp+eps15) .and. (ray(1) .gt. 0._wp-eps15)) then
          sgnx = sign(1._wp,ray(1))
       end if
       if ((ray(2) .lt. 0._wp+eps15) .and. (ray(2) .gt. 0._wp-eps15)) then
          sgny = sign(1._wp,ray(2))
       end if

       if ((distlist(k,1) .gt. 0._wp) .and. &
            (sign(1._wp,ray(1)) .eq. sgnx) .and. &
            (sign(1._wp,ray(2)) .eq. sgny) .and. &
            (distlist(k,2) .ge. 0._wp) .and. (distlist(k,3) .ge. 0._wp) .and. &
            (distlist(k,2) .le. width) .and. (distlist(k,3) .le. height)) then
          if (el .eq. lastboundary) then
             if ((scan(el, "L", kind=8) .eq. 1) .or. (any(bucketid .eq. lastboundary))) then
                cycle
             end if
             listel = pack([(j,j=1,pots)],potsid .eq. el)
             call findRefInd(ncurr,nalt,listel(1))
             if (ncurr .le. nalt) then
                cycle
             end if
          end if
          lastboundary = el
          isect(1) = distlist(k,2)
          isect(2) = distlist(k,3)
          detection = .true.
          exit
       end if
    end do

    ! In case something went wrong
    if (detection .eqv. .false.) then
       print *, "No detection in pickIntersectionTIR"
       print *, "Current ray is:   ", nextid
       print *, "Direction of ray: ", ray
       print *, "Current point     ", cray
       print *, "Current boundary: ", lastboundary
       do k = 1, bndrs
          print *, distlist(k,1), distlist(k,4), edgeid(int(distlist(k,4)))
       end do
       call exit(1)
    end if
  end subroutine pickIntersectionTIR

  ! Return a reflected ray in the opposite direction with a random angle [0+eps,180-eps]
  subroutine randomMirror(ray,wall)
    real(wp), dimension(2), intent(inout)    :: ray
    real(wp), dimension(2), intent(in)       :: wall
    real(wp), dimension(2,2)                 :: rotmat
    real(wp)                                 :: rng, uni

    ! Maxwell Demon Mirror
    ! call random_number(rng)
    ! rng = (PI-2._wp*eps)*rng-PI/2._wp+eps ! rng is uniformly distributed in [-pi+eps,pi-eps]

    ! Cosine mirror
    call random_number(rng)
    call random_number(uni)
    uni = sign(1._wp,uni-0.5_wp)
    rng = uni*asin(rng)
    
    rotmat = rotate(rng)
    ray = matmul(rotmat,wall) ! ray is now the randomly rotated wall normal
    
  end subroutine randomMirror
    
  ! Check if a point (xsect,ysect) is on the line defined by edgecoords
  pure function isOnLine(ec1,ec2,ec3,ec4,xsect,ysect)
    real(wp), intent(in)                  :: xsect, ysect, ec1, ec2, ec3, ec4
    real(wp)                              :: A, B, C
    logical                               :: isOnLine

    A = distance(ec1,ec2,ec3,ec4)
    B = distance(ec1,ec2,xsect,ysect)
    C = distance(ec3,ec4,xsect,ysect)
    if ((A .lt. B+C+eps) .and. (A .gt. B+C-eps)) then
       isOnLine = .true.
    else
       isOnLine = .false.
    end if
  end function isOnLine

  ! Calc the distance between two points p and q
  pure function distance(p1,p2,q1,q2)
    real(wp), intent(in)    :: p1, p2, q1, q2
    real(wp)                :: distance
    distance = sqrt((p1-q1)**2+(p2-q2)**2)
  end function distance

  ! Call the functions to build the chosen potential
  subroutine setPotential(type,k)
    character(len=4), intent(in)    :: type
    integer(wp), optional           :: k

    if ((type .eq. "1p1l") .and. (present(k) .eqv. .true.)) then
       call initialize1p1l(k)
    elseif ((type .eq. "1p1l") .and. (present(k) .eqv. .false.)) then
       call initialize1p1l()
       
    elseif ((type .eq. "1p2l") .and. (present(k) .eqv. .true.)) then
       call initialize1p2l(k)
    elseif ((type .eq. "1p2l") .and. (present(k) .eqv. .false.)) then
       call initialize1p2l()
       
    elseif ((type .eq. "1p4l") .and. (present(k) .eqv. .true.)) then
       call initialize1p4l(k)
    elseif ((type .eq. "1p4l") .and. (present(k) .eqv. .false.)) then
       call initialize1p4l()
       
    elseif ((type .eq. "1pnl") .and. (present(k) .eqv. .true.)) then
       call initialize1pnl(k)
    elseif ((type .eq. "1pnl") .and. (present(k) .eqv. .false.)) then
       call initialize1pnl()
       
    elseif ((type .eq. "2p3l") .and. (present(k) .eqv. .true.)) then
       call initialize2p3l(k)
    elseif ((type .eq. "2p3l") .and. (present(k) .eqv. .false.)) then
       call initialize2p3l()
       
    elseif ((type .eq. "2p2l") .and. (present(k) .eqv. .true.)) then
       call initialize2p2l(k)
    elseif ((type .eq. "2p2l") .and. (present(k) .eqv. .false.)) then
       call initialize2p2l()
       
    elseif ((type .eq. "2p6l") .and. (present(k) .eqv. .true.)) then
       call initialize2p6l(k)
    elseif ((type .eq. "2p6l") .and. (present(k) .eqv. .false.)) then
       call initialize2p6l()
       
    elseif ((type .eq. "1p1c") .and. (present(k) .eqv. .true.)) then
       call initialize1p1c(k)
    elseif ((type .eq. "1p1c") .and. (present(k) .eqv. .false.)) then
       call initialize1p1c()

    elseif ((type .eq. "1p2c") .and. (present(k) .eqv. .true.)) then
       call initialize1p2c(k)
    elseif ((type .eq. "1p2c") .and. (present(k) .eqv. .false.)) then
       call initialize1p2c()
    elseif ((type .eq. "1p2cdisk") .and. (present(k) .eqv. .true.)) then
       call initialize1p2cdisk(k)
    elseif ((type .eq. "1p2cdisk") .and. (present(k) .eqv. .false.)) then
       call initialize1p2cdisk()
    else
       print *, type, " is not a valid potential"
       call exit(1)
    end if
  end subroutine setPotential

  ! Sorts Nx4 matrix, quicksort implementation
  recursive pure subroutine quicksort4cols(a, first, last, sortbycol)
    implicit none
    real(wp), dimension(:,:), intent(inout)  :: a
    real(wp)                                 :: x
    real(wp), dimension(4)                   :: t
    integer(wp), intent(in)                  :: first, last, sortbycol
    integer(wp)                              :: i, j

    x = a((first+last)/2,sortbycol)
    i = first
    j = last
    do
       do while (a(i,sortbycol) < x)
          i=i+1
       end do
       do while (x < a(j,sortbycol))
          j=j-1
       end do
       if (i >= j) exit
       t = a(i,:);  a(i,:) = a(j,:);  a(j,:) = t
       i=i+1
       j=j-1
    end do
    if (first < i-1) call quicksort4cols(a, first, i-1, sortbycol)
    if (j+1 < last)  call quicksort4cols(a, j+1, last, sortbycol)
  end subroutine quicksort4cols

  ! Sorts array, insertionsort implementation
  subroutine insertionSort4cols(a, sortbycol)
    real(wp), dimension(:,:), intent(inout)   :: a
    real(wp), dimension(4)                    :: temp
    integer(wp)                               :: i, j
    integer(sp), intent(in)                   :: sortbycol

    do i = 2, size(a,1)
       j = i - 1
       temp = a(i,:)
       do while (j .ge. 1 .and. a(j,sortbycol) .gt. temp(sortbycol))
          a(j+1,:) = a(j,:)
          j = j - 1
          if (j .eq. 0) exit
       end do
       a(j+1,:) = temp
    end do
  end subroutine insertionSort4cols
  
  ! Defines a rotation matrix
  function rotate(theta)
    real(wp), dimension(2,2)    :: rotate
    real(wp), intent(in)        :: theta
    
    rotate(1,1) = cos(theta)
    rotate(2,1) = sin(theta)
    rotate(1,2) = -sin(theta)
    rotate(2,2) = cos(theta)
    
  end function rotate

  ! Transverse Electric reflection coefficient from Pedrotti p.495-496
  function rte(theta,n)
    real(wp)                :: rte, tcos
    real(wp), intent(in)    :: theta, n
    real(wp)                :: temp
    
    temp = sqrt(n**2-sin(theta)**2)
    tcos = cos(theta)
    rte =  abs((tcos-temp)/(tcos+temp))**2
    
    if (rte .gt. 1._wp+eps) then
       print *, "rte: ", rte
       print *, "rte > 1, something is wrong"
       call exit(1)
    end if
    
  end function rte

  ! Transverse Electric transmission coefficient from Pedrotti p.495-496
  function tte(theta,n)
    real(wp)                :: tte, tcos
    real(wp), intent(in)    :: theta, n
    real(wp)                :: temp

    temp = sqrt(n**2-sin(theta)**2)
    tcos = cos(theta)
    tte  = (temp/tcos)*abs((2._wp*tcos)/(tcos+temp))**2

    if (tte .gt. 1._wp+eps) then
       print *, "tte: ", tte
       print *, "Nextid: ", nextid
       print *, "tte > 1, something is wrong"
       call exit(1)
    end if

  end function tte

  ! Conic-Line intersection
  subroutine conicIntersection(rayx,rayy,crayx,crayy,isectx,isecty,boundary)
    real(wp), intent(in)          :: rayx, rayy, crayy, crayx
    real(wp), intent(inout)       :: isectx, isecty
    real(wp)                      :: discrim, a, b, c, t1, t2, ix1, ix2, iy1, iy2, d1, d2, sax, say, xoff, yoff, epsi
    character(len=4), intent(in)  :: boundary

    ! Local eps
    epsi = 1.5e-14_wp
    
    if (boundary .eq. "X1") then
       xoff = ellipse(1,1)
       yoff = ellipse(1,2)
       sax = ellipse(1,3)**2
       say = ellipse(1,4)**2
    elseif (boundary .eq. "X2") then
       xoff = ellipse(2,1)
       yoff = ellipse(2,2)
       sax = ellipse(2,3)**2
       say = ellipse(2,4)**2
    end if
    
    ! Coefficients of a line intersection an ellipse quadric
    a = say*rayx*rayx+sax*rayy*rayy
    b = 2._wp*say*rayx*(crayx-xoff)+2._wp*sax*rayy*(crayy-yoff)
    c = say*(crayx-xoff)**2+sax*(crayy-yoff)**2-sax*say

    discrim = b**2-4*a*c

    if (discrim .gt. 0._wp) then
       t1 = (-b+sqrt(discrim))/(2._wp*a)
       t2 = (-b-sqrt(discrim))/(2._wp*a)

       ix1 = crayx+t1*rayx
       iy1 = crayy+t1*rayy
       ix2 = crayx+t2*rayx
       iy2 = crayy+t2*rayy

       d1 = distance(crayx,crayy,ix1,iy1)
       d2 = distance(crayx,crayy,ix2,iy2)
       
       ! The case where the ray scatters on the curved surface
       if      ((d1 .lt. epsi) .and. (t2 .gt. 0._wp+epsi) .and. (iy2 .ge. yoff)) then
          isectx = ix2
          isecty = iy2
          return
       else if ((d2 .lt. epsi) .and. (t1 .gt. 0._wp+epsi) .and. (iy1 .ge. yoff)) then
          isectx = ix1
          isecty = iy1
          return
       end if

       ! The case where the ray comes from elsewhere
       if (d1 .lt. d2) then
          if ((t1 .gt. epsi) .and. (iy1 .ge. yoff)) then
             isectx = ix1
             isecty = iy1
          elseif ((t2 .gt. epsi) .and. (iy2 .ge. yoff)) then
             isectx = ix2
             isecty = iy2
          else
             isectx = 2._wp*width
             isecty = 2._wp*height
          end if
       elseif (d2 .lt. d1) then
          if ((t2 .gt. epsi) .and. (iy2 .ge. yoff)) then
             isectx = ix2
             isecty = iy2
          elseif ((t1 .gt. epsi) .and. (iy1 .ge. yoff)) then
             isectx = ix1
             isecty = iy1
          else
             isectx = 2._wp*width
             isecty = 2._wp*height
          end if
       elseif (d1 .eq. d2) then
          if (iy1 .ge. yoff) then
             isectx = ix1
             isecty = iy1
          elseif (iy2 .ge. yoff) then
             isectx = ix2
             isecty = iy2
          end if
       else
          print *, "Failure in conicSection!"
          print *, d1, d2
          print *, crayx, crayy
          print *, ix1, iy1
          print *, ix2, iy2
          call exit(1)
       end if

    else ! Out of bounds
       isectx = 2._wp*width
       isecty = 2._wp*height
    end if
  end subroutine conicIntersection

  ! Count surviving rays
  subroutine checkForSurvivors(PSOScount,limits,survived)
    integer(wp)                   :: l
    integer(wp), intent(inout)    :: survived(:)
    integer(wp), intent(in)       :: limits(:), PSOScount

    l = 1
    do while ((PSOScount .ge. limits(l)) .and. (l .le. boxcntlim+1))
       !$OMP atomic
       survived(l) = survived(l)+1
       !$OMP end atomic
       if (l .le. boxcntlim) then
          l = l+1
       else
          exit
       end if
    end do
  end subroutine checkForSurvivors
  
  ! Find the normal to a conic shape
  subroutine findConicNormal(normal, crayx, crayy, boundary)
    real(wp), dimension(2), intent(inout)  :: normal
    real(wp), intent(in)                   :: crayx, crayy
    real(wp)                               :: t, sax, say, xoff, yoff
    character(len=4)                       :: boundary

    if (boundary .eq. "X1") then
       xoff = ellipse(1,1)
       yoff = ellipse(1,2)
       sax = ellipse(1,3)
       say = ellipse(1,4)
    elseif (boundary .eq. "X2") then
       xoff = ellipse(2,1)
       yoff = ellipse(2,2)
       sax = ellipse(2,3)
       say = ellipse(2,4)
    end if

    if ((crayx-xoff)/sax .gt. 1._wp) then
       t = acos(1._wp)
    elseif ((crayx-xoff)/sax .lt. -1._wp) then
       t = acos(-1._wp)
    else
       t = acos((crayx-xoff)/sax)
    end if

    normal(1) = say*cos(t)
    normal(2) = sax*sin(t)
    normal = normal/sqrt(normal(1)**2+normal(2)**2)
  end subroutine findConicNormal

  subroutine readFileLength(fileid,rows)
    integer(wp), intent(inout)    :: rows
    integer(wp), intent(in)       :: fileid
    integer(wp)                   :: io

    ! Count the number of lines in a file
    rows = 0
    rewind(fileid)
    do
       read(fileid,*,iostat=io)
       if (io .ne. 0) then
          ! print *, "eof reached, all good"
          exit
       end if
       rows = rows + 1
    end do
    rewind(fileid)

    ! Now rows is the total number of lines in the file
    ! not interested in the descriptors (1 line header), so:
    rows = rows-1
    
  end subroutine readFileLength

  ! Create logbook entry
  subroutine logParameters(hours,minutes,seconds,programtype)
    implicit none
    integer(wp), intent(in)          :: hours, minutes, seconds
    real(wp)                         :: tsec
    character(len=4), intent(in)     :: programtype
    integer(wp), dimension(8)        :: dt

    tsec = hours*3600+minutes*60+seconds
    call date_and_time(values=dt)

    ! Create logbook
    open(unit = 25, file = "logbook.txt", form='formatted')
    write(25, FMT='(A62)') "______________________________________________________________"
    write(25, FMT='(A62)') ""
    write(25, FMT='(A6,I0,A1,I0,A1,I0)') "Date: ", dt(1), "-", dt(2), "-", dt(3)
    write(25, FMT='(A6,I0,A1,I0)')       "Time: ", dt(5), ":", dt(6)
    write(25, FMT='(A62)') "______________________________________________________________"
    write(25, FMT='(A62)') ""
    write(25, FMT='(A26,F4.2)')     "System width:             ", width
    write(25, FMT='(A26,F4.2)')     "System height:            ", height
    write(25, FMT='(A26,ES8.2E2,A2)')"Wavelength:               ", lambda/1e6_wp, " m"
    write(25, FMT='(A26,I0)')       "Variation chosen:         ", specshap
    write(25, FMT='(A26,L)')        "Raysplitting:             ", raysplitting
    write(25, FMT='(A26,I0)')       "Max bounce limit:         ", bouncelim
    write(25, FMT='(A26,F4.2)')     "Index of refraction 0:    ", n0
    write(25, FMT='(A26,F4.2)')     "Index of refraction 1:    ", n1
    write(25, FMT='(A26,F4.2)')     "Index of refraction 2:    ", n2

    select case (programtype)
    case ("trac")
       write(25, FMT='(A1)')           " " !spacer
       write(25, FMT='(A8)')           "Tracing:"
       write(25, FMT='(A26,I0)')       "Number of rays:           ", nrays
       write(25, FMT='(A26,I0)')       "Number of variations:     ", vars
       write(25, FMT='(A26,I0,A1,I0)') "Total number of rays:     ", nrays*vars
       write(25, FMT='(A26,I0)')       "Number of rays escaping:  ", erays
       write(25, FMT='(A26,I0)')       "Number of PSOS bounces:   ", PSOSrays
    case ("box ")
       write(25, FMT='(A1)')           " " !spacer
       write(25, FMT='(A4)')           "BCM:"
       write(25, FMT='(A26,I0,A1,I0)') "Poincare grid:            ", grid, "x", grid
       write(25, FMT='(A26,I0)')       "Box counting limit:       ", boxcntlim
       write(25, FMT='(A26,I0,A1,I0)') "Total number of rays:     ", grid**2
       write(25, FMT='(A26,I0)')       "Number of rays escaping:  ", erays
       write(25, FMT='(A26,I0)')       "Number of PSOS bounces:   ", PSOSrays
    end select

    write(25, FMT='(A1)')              " " !spacer
    write(25, FMT='(A26,A4)')          "Potential type:           ", potentialtype

    write(25, FMT='(A1)')           " " !spacer
    write(25, FMT='(A10)', advance="no") "Run time: "
    write(25, FMT='(I2,A6)') hours, " hours"
    write(25, FMT='(T11,I2,A8)') minutes, " minutes"
    write(25, FMT='(T11,I2,A8)') seconds, " seconds"
    write(25, FMT='(A1)')           " " !spacer
    close(unit = 25)
  end subroutine logParameters

end module functions
