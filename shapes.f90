module shapes
  use parameters
  implicit none
  !_______________________________________________________________________________________
  !
  ! Shapes
  !_______________________________________________________________________________________
  
  integer(sp)    :: pline1, pline2, pots, bndrs, curves
  
  ! Define bucket
  real(wp), dimension(2), parameter            :: vert = [0._wp,1._wp] ! Left and right wall
  real(wp), dimension(2), parameter            :: hor  = [1._wp,0._wp] ! Ceiling and bottom wall
  real(wp), dimension(4,2), parameter          :: bucket   = transpose(reshape([vert,vert,hor,hor], [2,4]))
  character(len=1), dimension(4), parameter    :: bucketid = [ "R", "L", "B", "C"]

  ! Define bucket, never change
  real(wp), dimension(2), parameter       :: sp1 = [width,0._wp]  ! Lower right
  real(wp), dimension(2), parameter       :: sp2 = [width,height] ! Top right
  real(wp), dimension(2), parameter       :: sp3 = [0._wp,height] ! Top left
  real(wp), dimension(2), parameter       :: sp4 = [0._wp,0._wp]  ! Lower left

  real(wp), dimension(1,4), parameter     :: sR = transpose(reshape([ sp1,sp2 ], [4,1]))
  real(wp), dimension(1,4), parameter     :: sL = transpose(reshape([ sp4,sp3 ], [4,1]))
  real(wp), dimension(1,4), parameter     :: sB = transpose(reshape([ sp4,sp1 ], [4,1]))
  real(wp), dimension(1,4), parameter     :: sC = transpose(reshape([ sp3,sp2 ], [4,1]))

  ! Define unit normals for bucket depending on width and height, never change
  real(wp), dimension(2), parameter     :: rwall = [-1,0] ! Right wall
  real(wp), dimension(2), parameter     :: lwall = [1,0]  ! Left wall
  real(wp), dimension(2), parameter     :: bwall = [0,1]  ! Bottom wall
  real(wp), dimension(2), parameter     :: cwall = [0,-1] ! Ceiling

  ! Define potential in terms of points from left to right. Minimum 2 points => 1 line
  real(wp), parameter    :: oldskew = -(iskew+(specshap-1)*(fskew-iskew)/(nshapes-1))
  real(wp)               :: skew

  real(wp), allocatable             :: potential(:,:), point(:,:), edgecoords(:,:), ellipse(:,:)
  real(wp), allocatable             :: rileft(:), riright(:)
  character(len=4), allocatable     :: edgeid(:), potsid(:), abslist(:)

contains
  subroutine initialize1p1c(k)
    real(wp), dimension(2)    :: p1,p2,p3,p4,l1,l2,l3
    integer(wp)               :: j
    integer(wp), optional     :: k

    ! Set numbers
    pline1 = 3 ! Total number of potential lines in line 1
    curves = 1 ! Total number of curved potentials
    pots   = pline1+curves
    bndrs  = pots+4 ! Total number of lines
    
    call p1c1(p1,p2,p3,p4,l1,l2,l3,k)
    
    ! Allocate
    allocate(potential(pots-1,2))
    allocate(point(pots-1,2))
    allocate(edgecoords(bndrs-1,4))
    allocate(edgeid(bndrs))
    allocate(potsid(pots))
    allocate(rileft(pots))
    allocate(riright(pots))
    allocate(abslist(pline1))

    ! If the varying part is index of refraction in the dome. Otherwise comment line below
    n1 = skew
    
    ! Set
    potential = transpose(reshape([ l1,l2,l3 ], [2,pots-1])) ! List of direction vectors
    point = transpose(reshape([ p2,p3,p4 ], [2,pots-1]))    ! List of point on lines
    edgecoords = transpose(reshape([ sp1,sp2,sp4,sp3,sp4,sp1,sp3,sp2,p1,p2,p2,p3,p3,p4 ], [4,bndrs-1])) ! -1 because curved potential works differently. There are no edgecoords needed, but it still counts as a "potential line"
    edgeid = [ "R ","L ","B ","C ","L1","L2","L3","X1" ]
    potsid = [ "L1","L2","L3","X1" ]
    rileft  = [ n2, n2, n2, n1 ]
    riright = [ n0, n1, n0, n0 ]
    
    ! Set which lines should turn absorption on/off
    do j = 1, pline1
       abslist(j) = potsid(j)
    end do
  end subroutine initialize1p1c

  subroutine p1c1(p1,p2,p3,p4,l1,l2,l3,k)
    real(wp)       :: dx, x1, x2, base, tinybit, arg, xoff, yoff, semiaxisx, semiaxisy, radius
    integer(wp)    :: j, steps
    real(wp), dimension(2), intent(inout)    :: p1,p2,p3,p4,l1,l2,l3
    character(len=90)                        :: filename
    integer(wp), optional                    :: k

    ! Set globals
    radius = 0.4999_wp*width
    xoff = 0.5_wp*width
    yoff = 2.0_wp
    tinybit = xoff-radius
    semiaxisx = 1.0_wp*radius
    semiaxisy = 5*0.69216923076923076_wp ! best at n=2, thick=400nm

    allocate(ellipse(curves,4))
    ellipse = transpose(reshape([ xoff,yoff,semiaxisx,semiaxisy ], [4,1]))
    
    if (present(k) .eqv. .false.) then
       print *, "Ellipse semiaxisy: ", semiaxisy
    end if

    if (semiaxisy+yoff .ge. height) then
       print *, "The dome tip is too far up! :", semiaxisy+yoff
       print *, "C is located at:       ", height
       call exit(1)
    end if
    
    base = yoff
    p1 = [0._wp,base]
    p2 = [tinybit,base]
    p3 = [xoff+radius,base]
    p4 = [width,base]

    ! Lines
    l1 = [p2(1)-p1(1),p2(2)-p1(2)]
    l2 = [p3(1)-p2(1),p3(2)-p2(2)]
    l3 = [p4(1)-p3(1),p4(2)-p3(2)]

    steps = 500
    
    ! Print potential
    x1 = ellipse(1,1)-ellipse(1,3)
    x2 = ellipse(1,1)+ellipse(1,3)
    dx = (x2-x1)/(steps-1)
    open(unit = 1, file = "tpotential1.txt", form='formatted')
    write(1,FMT='(2A15)') "x", "y"
    do j = 1, steps
       arg = abs(1._wp-(((x1+dx*(j-1))-ellipse(1,1))/ellipse(1,3))**2)
       write(1,*) x1+dx*(j-1), ellipse(1,4)*sqrt(arg)+ellipse(1,2)
    end do
    close(1)

    open(unit = 1, file = "tpotential2.txt", form='formatted')
    write(1,FMT='(2A15)') "x", "y"
    write(1,*) p1(1), p1(2)
    write(1,*) p2(1), p2(2)
    write(1,*) p3(1), p3(2)
    write(1,*) p4(1), p4(2)
    close(1)

    ! If multi is ran
    if (present(k)) then
       write(filename, '("tpotential1-", I0, ".txt")') k
       open(unit = 1, file = filename, form='formatted')
       write(1,FMT='(2A15)') "x", "y"
       do j = 1, steps
          arg = abs(1._wp-(((x1+dx*(j-1))-xoff)/semiaxisx)**2)
          write(1,*) x1+dx*(j-1), semiaxisy*sqrt(arg)+yoff
       end do
       close(1)

       write(filename, '("tpotential2-", I0, ".txt")') k
       open(unit = 1, file = filename, form='formatted')
       write(1,FMT='(2A15)') "x", "y"
       write(1,*) p1(1), p1(2)
       write(1,*) p2(1), p2(2)
       write(1,*) p3(1), p3(2)
       write(1,*) p4(1), p4(2)
       close(1)
    end if
  end subroutine p1c1

  subroutine initialize1p2c(k)
    real(wp), dimension(2)    :: p1,p2,p3,p4,p5,p6,l1,l2,l3,l4,l5
    integer(wp)               :: j
    integer(wp), optional     :: k

    ! Set numbers
    pline1 = 5 ! Total number of potential lines in line 1
    curves = 2 ! Total number of curved potentials
    pots   = pline1+curves
    bndrs  = pots+4 ! Total number of lines
    
    call p1c2(p1,p2,p3,p4,p5,p6,l1,l2,l3,l4,l5,k)
    
    ! Allocate
    allocate(potential(pline1,2))
    allocate(point(pline1,2))
    allocate(edgecoords(bndrs-curves,4))
    allocate(edgeid(bndrs))
    allocate(potsid(pots))
    allocate(rileft(pots))
    allocate(riright(pots))
    allocate(abslist(pline1))

    ! Set
    potential = transpose(reshape([ l1,l2,l3,l4,l5 ], [2,pline1])) ! List of direction vectors
    point = transpose(reshape([ p2,p3,p4,p5,p6 ], [2,pline1]))    ! List of point on lines
    edgecoords = transpose(reshape([ sp1,sp2,sp4,sp3,sp4,sp1,sp3,sp2,p1,p2,p2,p3,p3,p4,p4,p5,p5,p6 ], [4,bndrs-curves])) ! -1 because curved potential works differently. There are no edgecoords needed, but it still counts as a "potential line"
    edgeid = [ "R ","L ","B ","C ","L1","L2","L3","L4","L5","X1","X2" ]
    potsid = [ "L1","L2","L3","L4","L5","X1","X2" ]
    rileft  = [ n1, n2, n1, n2, n1, n2, n2 ]
    riright = [ n0, n1, n0, n1, n0, n0, n0 ]
    
    ! Set which lines should turn absorption on/off
    do j = 1, pline1
       abslist(j) = potsid(j)
    end do
  end subroutine initialize1p2c

  subroutine p1c2(p1,p2,p3,p4,p5,p6,l1,l2,l3,l4,l5,k)
    real(wp)       :: dx, x1, x2, base, tinybit1, tinybit2, arg, axx, axy, semiaxisx, semiaxisy, radius, xoff, yoff
    integer(wp)    :: j, steps
    real(wp), dimension(2), intent(inout)    :: p1,p2,p3,p4,p5,p6,l1,l2,l3,l4,l5
    character(len=90)                        :: filename
    integer(wp), optional                    :: k

    ! Set globals
    radius = 0.4999_wp*width
    xoff = 0.5_wp*width
    yoff = normalizedthickness!*height
    axx = 0.249_wp*width
    axy = -axx*skew!2.5_wp
    axy = 0.7_wp*radius!2.5_wp
    print *, axy/axx
    allocate(ellipse(curves,4))
    ! ellipse(x,:) = [ xoff,yoff,axx,axy ] How to input
    ellipse(1,:) = [ 0.25_wp*width,yoff,axx,axy*1.0_wp ]
    ellipse(2,:) = [ 0.75_wp*width,yoff,axx,axy*1.0_wp ]
    tinybit1 = ellipse(1,1)-ellipse(1,3)
    tinybit2 = width-ellipse(2,1)-ellipse(2,3)!ellipse(2,1)-ellipse(2,3)

    ! Points
    base = yoff
    p1 = [0._wp,base]
    p2 = [tinybit1,base]
    p3 = [ellipse(1,1)+ellipse(1,3),base]
    p4 = [p3(1)+tinybit1+tinybit2,base]
    p5 = [ellipse(2,1)+ellipse(2,3),base]
    p6 = [width,base]

    ! Lines
    l1 = [ p2(1)-p1(1),p2(2)-p1(2) ]
    l2 = [ p3(1)-p2(1),p3(2)-p2(2) ]
    l3 = [ p4(1)-p3(1),p4(2)-p3(2) ]
    l4 = [ p5(1)-p4(1),p5(2)-p4(2) ]
    l5 = [ p6(1)-p5(1),p6(2)-p5(2) ]

    steps = 500
    x1 = xoff-radius
    x2 = xoff+radius
    dx = (x2-x1)/(steps-1)

    ! Print potential
    x1 = ellipse(1,1)-ellipse(1,3)
    x2 = ellipse(1,1)+ellipse(1,3)
    dx = (x2-x1)/(steps-1)
    open(unit = 1, file = "tpotential1.txt", form='formatted')
    write(1,FMT='(2A15)') "x", "y"
    do j = 1, steps
       arg = abs(1._wp-(((x1+dx*(j-1))-ellipse(1,1))/ellipse(1,3))**2)
       write(1,*) x1+dx*(j-1), ellipse(1,4)*sqrt(arg)+ellipse(1,2)
    end do
    close(1)

    x1 = ellipse(2,1)-ellipse(2,3)
    x2 = ellipse(2,1)+ellipse(2,3)
    dx = (x2-x1)/(steps-1)
    open(unit = 1, file = "tpotential2.txt", form='formatted')
    write(1,FMT='(2A15)') "x", "y"
    do j = 1, steps
       arg = abs(1._wp-(((x1+dx*(j-1))-ellipse(2,1))/ellipse(2,3))**2)
       write(1,*) x1+dx*(j-1), ellipse(2,4)*sqrt(arg)+ellipse(2,2)
    end do
    close(1)

    open(unit = 1, file = "tpotential3.txt", form='formatted')
    write(1,FMT='(2A15)') "x", "y"
    write(1,*) p1(1), p1(2)
    write(1,*) p2(1), p2(2)
    write(1,*) p3(1), p3(2)
    write(1,*) p4(1), p4(2)
    write(1,*) p5(1), p5(2)
    write(1,*) p6(1), p6(2)
    close(1)

    ! If multi is ran
    if (present(k)) then
       write(filename, '("tpotential1-", I0, ".txt")') k
       open(unit = 1, file = filename, form='formatted')
       write(1,FMT='(2A15)') "x", "y"
       do j = 1, steps
          arg = abs(1._wp-(((x1+dx*(j-1))-xoff)/axx)**2)
          write(1,*) x1+dx*(j-1), axy*sqrt(arg)+yoff
       end do
       close(1)

       write(filename, '("tpotential2-", I0, ".txt")') k
       open(unit = 1, file = filename, form='formatted')
       write(1,FMT='(2A15)') "x", "y"
       write(1,*) p1(1), p1(2)
       write(1,*) p2(1), p2(2)
       write(1,*) p3(1), p3(2)
       write(1,*) p4(1), p4(2)
       write(1,*) p5(1), p5(2)
       write(1,*) p6(1), p6(2)
       close(1)
    end if
  end subroutine p1c2

  subroutine initialize1p2cdisk(k)
    real(wp), dimension(2)    :: p1,p2,p3,p4,p5,p6,l1,l2,l3,l4,l5
    integer(wp)               :: j
    integer(wp), optional     :: k

    ! Set numbers
    pline1 = 5 ! Total number of potential lines in line 1
    curves = 2 ! Total number of curved potentials
    pots   = pline1+curves
    bndrs  = pots+4 ! Total number of lines
    
    call p1c2disk(p1,p2,p3,p4,p5,p6,l1,l2,l3,l4,l5,k)
    
    ! Allocate
    allocate(potential(pline1,2))
    allocate(point(pline1,2))
    allocate(edgecoords(bndrs-curves,4))
    allocate(edgeid(bndrs))
    allocate(potsid(pots))
    allocate(rileft(pots))
    allocate(riright(pots))
    allocate(abslist(pline1))

    ! Set
    potential = transpose(reshape([ l1,l2,l3,l4,l5 ], [2,pline1])) ! List of direction vectors
    point = transpose(reshape([ p2,p3,p4,p5,p6 ], [2,pline1]))    ! List of point on lines
    edgecoords = transpose(reshape([ sp1,sp2,sp4,sp3,sp4,sp1,sp3,sp2,p1,p2,p2,p3,p3,p4,p4,p5,p5,p6 ], [4,bndrs-curves])) ! -1 because curved potential works differently. There are no edgecoords needed, but it still counts as a "potential line"
    edgeid = [ "R ","L ","B ","C ","L1","L2","L3","L4","L5","X1","X2" ]
    potsid = [ "L1","L2","L3","L4","L5","X1","X2" ]
    rileft  = [ n1, n2, n1, n2, n1, n2, n2 ]
    riright = [ n0, n1, n0, n1, n0, n0, n0 ]
    
    ! Set which lines should turn absorption on/off
    do j = 1, pline1
       abslist(j) = potsid(j)
    end do
  end subroutine initialize1p2cdisk

  subroutine p1c2disk(p1,p2,p3,p4,p5,p6,l1,l2,l3,l4,l5,k)
    real(wp)       :: dx, x1, x2, base, tinybit1, tinybit2, arg, axx, axy, semiaxisx, semiaxisy, radius, xoff, yoff
    integer(wp)    :: j, steps
    real(wp), dimension(2), intent(inout)    :: p1,p2,p3,p4,p5,p6,l1,l2,l3,l4,l5
    character(len=90)                        :: filename
    integer(wp), optional                    :: k

    ! Set globals
    radius = 0.4999_wp*width
    xoff = 0.5_wp*width
    yoff = normalizedthickness!*height
    axx = 0.249_wp*width
    axy = -axx*skew!2.5_wp
    axy = 0.7_wp*radius!2.5_wp
    print *, axy/axx
    allocate(ellipse(curves,4))
    ellipse(1,:) = [ 0.25_wp*width,yoff,axx,axy*1.0_wp ]
    ellipse(2,:) = [ 0.75_wp*width,yoff,axx,axy*1.0_wp ]
    tinybit1 = ellipse(1,1)-ellipse(1,3)
    tinybit2 = width-ellipse(2,1)-ellipse(2,3)!ellipse(2,1)-ellipse(2,3)

    ! Points
    base = yoff
    p1 = [0._wp,base]
    p2 = [tinybit1,base]
    p3 = [ellipse(1,1)+ellipse(1,3),base]
    p4 = [p3(1)+tinybit1+tinybit2,base]
    p5 = [ellipse(2,1)+ellipse(2,3),base]
    p6 = [width,base]

    ! Lines
    l1 = [ p2(1)-p1(1),p2(2)-p1(2) ]
    l2 = [ p3(1)-p2(1),p3(2)-p2(2) ]
    l3 = [ p4(1)-p3(1),p4(2)-p3(2) ]
    l4 = [ p5(1)-p4(1),p5(2)-p4(2) ]
    l5 = [ p6(1)-p5(1),p6(2)-p5(2) ]

    steps = 500
    x1 = xoff-radius
    x2 = xoff+radius
    dx = (x2-x1)/(steps-1)

    ! Print potential
    x1 = ellipse(1,1)-ellipse(1,3)
    x2 = ellipse(1,1)+ellipse(1,3)
    dx = (x2-x1)/(steps-1)
    open(unit = 1, file = "tpotential1.txt", form='formatted')
    write(1,FMT='(2A15)') "x", "y"
    do j = 1, steps
       arg = abs(1._wp-(((x1+dx*(j-1))-ellipse(1,1))/ellipse(1,3))**2)
       write(1,*) x1+dx*(j-1), ellipse(1,4)*sqrt(arg)+ellipse(1,2)
    end do
    close(1)

    x1 = ellipse(2,1)-ellipse(2,3)
    x2 = ellipse(2,1)+ellipse(2,3)
    dx = (x2-x1)/(steps-1)
    open(unit = 1, file = "tpotential2.txt", form='formatted')
    write(1,FMT='(2A15)') "x", "y"
    do j = 1, steps
       arg = abs(1._wp-(((x1+dx*(j-1))-ellipse(2,1))/ellipse(2,3))**2)
       write(1,*) x1+dx*(j-1), ellipse(2,4)*sqrt(arg)+ellipse(2,2)
    end do
    close(1)

    open(unit = 1, file = "tpotential3.txt", form='formatted')
    write(1,FMT='(2A15)') "x", "y"
    write(1,*) p1(1), p1(2)
    write(1,*) p2(1), p2(2)
    write(1,*) p3(1), p3(2)
    write(1,*) p4(1), p4(2)
    write(1,*) p5(1), p5(2)
    write(1,*) p6(1), p6(2)
    close(1)

    ! If multi is ran
    if (present(k)) then
       write(filename, '("tpotential1-", I0, ".txt")') k
       open(unit = 1, file = filename, form='formatted')
       write(1,FMT='(2A15)') "x", "y"
       do j = 1, steps
          arg = abs(1._wp-(((x1+dx*(j-1))-xoff)/axx)**2)
          write(1,*) x1+dx*(j-1), axy*sqrt(arg)+yoff
       end do
       close(1)

       write(filename, '("tpotential2-", I0, ".txt")') k
       open(unit = 1, file = filename, form='formatted')
       write(1,FMT='(2A15)') "x", "y"
       write(1,*) p1(1), p1(2)
       write(1,*) p2(1), p2(2)
       write(1,*) p3(1), p3(2)
       write(1,*) p4(1), p4(2)
       write(1,*) p5(1), p5(2)
       write(1,*) p6(1), p6(2)
       close(1)
    end if
  end subroutine p1c2disk
  
  subroutine initialize1p1l(k)
    real(wp), dimension(2)    :: p1,p2,l1
    integer(wp), optional     :: k

    ! Set numbers
    pline1 = 1 ! Total number of potential lines in line 1
    pots   = pline1
    bndrs  = pots+4 ! Total number of lines
    
    call p1l1(p1,p2,l1)
    
    ! Allocate
    allocate(potential(pots,2))
    allocate(point(pots,2))
    allocate(edgecoords(bndrs,4))
    allocate(edgeid(bndrs))
    allocate(potsid(pots))
    allocate(rileft(pots))
    allocate(riright(pots))
    allocate(abslist(pots))

    ! Set
    potential = transpose(reshape([ l1 ], [2,pots])) ! List of direction vectors
    point = transpose(reshape([ p2 ], [2,pots]))    ! List of point on lines
    edgecoords = transpose(reshape([ sp1,sp2,sp4,sp3,sp4,sp1,sp3,sp2,p1,p2 ], [4,bndrs]))
    edgeid = [ "R   ","L   ","B   ","C   ","L1  " ]
    potsid = [ "L1  " ]
    rileft  = [ n2 ]
    riright = [ n0 ]
    abslist = potsid
  end subroutine initialize1p1l

  subroutine p1l1(p1,p2,l1)
    real(wp), dimension(2), intent(inout)    :: p1,p2,l1
    real(wp)                                 :: base

    base = 3.5_wp!0.3_wp*height
    p1 = [0._wp,base]
    p2 = [width,base]

    ! Lines
    l1   = [p2(1)-p1(1),p2(2)-p1(2)]

    ! Print potential
    open(unit = 1, file = "tpotential1.txt", form='formatted')
    write(1,FMT='(2A15)') "x", "y"
    write(1,*) p1(1), p1(2)
    write(1,*) p2(1), p2(2)
    close(1)
  end subroutine p1l1

  subroutine initialize1p2l(k)
    real(wp), dimension(2)    :: p1,p2,p3,l1,l2
    integer(wp), optional     :: k

    ! Set numbers
    pline1 = 2 ! Total number of potential lines in line 1
    pots   = pline1
    bndrs  = pots+4 ! Total number of lines
    
    call p1l2(p1,p2,p3,l1,l2,k)
    
    ! Allocate
    allocate(potential(pots,2))
    allocate(point(pots,2))
    allocate(edgecoords(bndrs,4))
    allocate(edgeid(bndrs))
    allocate(potsid(pots))
    allocate(rileft(pots))
    allocate(riright(pots))
    allocate(abslist(pline1))

    ! Set
    potential = transpose(reshape([ l1,l2 ], [2,pots])) ! List of direction vectors
    point = transpose(reshape([ p2,p3 ], [2,pots]))    ! List of point on lines
    edgecoords = transpose(reshape([ sp1,sp2,sp4,sp3,sp4,sp1,sp3,sp2,p1,p2,p2,p3 ], [4,bndrs]))
    edgeid = [ "R ","L ","B ","C ","L1","L2" ]
    potsid = [ "L1","L2" ]
    abslist = potsid
    rileft  = [ n1, n1 ]
    riright = [ n0, n0 ]
  end subroutine initialize1p2l

  subroutine p1l2(p1,p2,p3,l1,l2,k)
    real(wp), dimension(2), intent(inout)  :: p1,p2,p3,l1,l2
    real(wp)                               :: base
    integer(wp), optional                  :: k
    character(len=90)                      :: filename

    base = 0.2_wp
    ! Top layer
    p1 = [0._wp,0.7_wp]
    p2 = [0.01_wp-skew,base]
    p3 = [width,0.7_wp]

    print *, "Low point is at: ", p2

    ! Lines
    l1 = [p2(1)-p1(1),p2(2)-p1(2)]
    l2 = [p3(1)-p2(1),p3(2)-p2(2)]

    if (present(k)) then
       ! Print potentials
       write(filename, '("tpotential1-", I0, ".txt")') k
       ! print *, "filename: ", filename
       open(unit = 1, file = filename, form='formatted')
       write(1,FMT='(2A15)') "x", "y"
       write(1,*) p1(1), p1(2)
       write(1,*) p2(1), p2(2)
       write(1,*) p3(1), p3(2)
       close(1)
    end if

    ! Print potential
    open(unit = 1, file = "tpotential1.txt", form='formatted')
    write(1,FMT='(2A15)') "x", "y"
    write(1,*) p1(1), p1(2)
    write(1,*) p2(1), p2(2)
    write(1,*) p3(1), p3(2)
    close(1)
  end subroutine p1l2

  subroutine initialize1p4l(k)
    real(wp), dimension(2)    :: p1,p2,p3,p4,p5,l1,l2,l3,l4
    integer(wp), optional     :: k

    ! Set numbers
    pline1 = 4 ! Total number of potential lines in line 1
    pots   = pline1
    bndrs  = pots+4 ! Total number of lines
    
    call p1l4(p1,p2,p3,p4,p5,l1,l2,l3,l4,k)
    
    ! Allocate
    allocate(potential(pots,2))
    allocate(point(pots,2))
    allocate(edgecoords(bndrs,4))
    allocate(edgeid(bndrs))
    allocate(potsid(pots))
    allocate(rileft(pots))
    allocate(riright(pots))

    ! Set
    potential = transpose(reshape([ l1,l2,l3,l4 ], [2,pots])) ! List of direction vectors
    point = transpose(reshape([ p2,p3,p4,p5 ], [2,pots]))    ! List of point on lines
    edgecoords = transpose(reshape([ sp1,sp2,sp4,sp3,sp4,sp1,sp3,sp2,p1,p2,p2,p3,p3,p4,p4,p5 ], [4,bndrs]))
    edgeid = [ "R ","L ","B ","C ","L1","L2","L3","L4" ]
    potsid = [ "L1","L2","L3","L4" ]
    rileft  = [ n1, n1, n1, n1 ]
    riright = [ n0, n0, n0, n0 ]
  end subroutine initialize1p4l

  subroutine p1l4(p1,p2,p3,p4,p5,l1,l2,l3,l4,k)
    real(wp), dimension(2), intent(inout)  :: p1,p2,p3,p4,p5,l1,l2,l3,l4
    real(wp)                               :: base
    integer(wp), optional                  :: k
    character(len=90)                      :: filename

    base = 0.217_wp
    p1 = [0._wp,0.7_wp]
    p2 = [0.01_wp-skew,base]
    p3 = [width/2._wp,0.7_wp]
    p4 = [0.01_wp-skew+width/2._wp,base]
    p5 = [width,0.7_wp]

    print *, "Low point 1 is at: ", p2
    print *, "Low point 2 is at: ", p4

    ! Lines
    l1 = [p2(1)-p1(1),p2(2)-p1(2)]
    l2 = [p3(1)-p2(1),p3(2)-p2(2)]
    l3 = [p4(1)-p3(1),p4(2)-p3(2)]
    l4 = [p5(1)-p4(1),p5(2)-p4(2)]

    if (present(k)) then
       ! Print potentials
       write(filename, '("tpotential1-", I0, ".txt")') k
       ! print *, "filename: ", filename
       open(unit = 1, file = filename, form='formatted')
       write(1,FMT='(2A15)') "x", "y"
       write(1,*) p1(1), p1(2)
       write(1,*) p2(1), p2(2)
       write(1,*) p3(1), p3(2)
       write(1,*) p4(1), p4(2)
       write(1,*) p5(1), p5(2)
       close(1)
    end if

    ! Print potential
    open(unit = 1, file = "tpotential1.txt", form='formatted')
    write(1,FMT='(2A15)') "x", "y"
    write(1,*) p1(1), p1(2)
    write(1,*) p2(1), p2(2)
    write(1,*) p3(1), p3(2)
    write(1,*) p4(1), p4(2)
    write(1,*) p5(1), p5(2)
    close(1)
  end subroutine p1l4

  subroutine initialize1pnl(k)
    real(wp)                  :: amp, ran, baseline, area, var, safety, bin, range
    real(wp), allocatable     :: xdir(:), ydir(:)
    integer(wp), optional     :: k
    integer(wp)               :: j, l, randompulls
    integer(sp), allocatable  :: seed(:)
    integer(sp)               :: n
    character(len=4)          :: pname
    character(len=90)         :: filename

    ! Set numbers
    amp = roughamp
    !pline1 = 2 ! Total number of potential lines in line 1
    if (present(k)) then
       pots = k!*randomlines
    else
       pots = randomlines
    end if
    bndrs  = pots+4 ! Total number of lines
    allocate(xdir(pots+1))
    allocate(ydir(pots+1))

    ! call random_seed(size=n)
    ! allocate(seed(n))
    ! call random_seed(get=seed)
    ! print *, seed

    ! Optional seeding of RNG (probably not needed)
    ! Skip a few numbers of RNG
    if (skiprandom .ne. 0) then ! 841
       randompulls = skiprandom*2*14
       do j = 1, randompulls
          call random_number(ran)
       end do
    end if

    baseline = normalizedthickness
    var = 0.045_wp ! variance?

    ! Gaussian y, approx normal x
    xdir(1) = 0._wp
    xdir(pots+1) = width
    ! call random_number(ran)
    ! ran = baseline+(ran-0.5_wp)*amp*2._wp
    ! call pullGauss(ran)
    ! ran = ran*var+baseline
    ydir(1) = baseline
    do j = 2, pots
       call random_number(ran)
       ! xdir(j) = ran+xdir(j-1)
       xdir(j) = (width/pots)*(j-1)+(ran-0.5_wp)*(width/(2._wp*pots))
       call random_number(ran)
       ran = baseline+(ran-0.5_wp)*amp*2._wp ! Random number [-amp,amp]+baseline
       ! ran = baseline+ran*amp*2._wp*(-1)**j
       ! call pullGauss(ran)
       ! ran = ran*var+baseline
       ydir(j) = ran
    end do
    ydir(pots+1) = baseline
    xdir = (width*xdir)/xdir(pots+1)
    ! do j = 1, pots+1
    !    print *, xdir(j), (width/pots)*(j-1)
    ! end do

    ! Calculate area of potential
    area = 0._wp
    do j = 1, pots
       area = area+(xdir(j+1)-xdir(j))*(ydir(j)+ydir(j+1))/2._wp ! area of a trapezoid
    end do
    ! print *, "Area under random potential: ", area

    ! Rescale potential
    ydir = (baseline*width/area)*ydir

    ! Check the new are of the potential
    area = 0._wp
    do j = 1, pots
       area = area+(xdir(j+1)-xdir(j))*(ydir(j)+ydir(j+1))/2._wp ! area of a trapezoid
    end do
    ! print *, "Corrected area is now: ", area

    ! Gaussian y, stretched random with safety x
    ! bin = width/pots
    ! safety = 0.1_wp*bin
    ! range = bin-safety
    ! xdir(1) = 0._wp
    ! call pullGauss(ran)
    ! ran = ran*var+baseline
    ! ydir(1) = ran
    ! do j = 1, pots-1
    !    call random_number(ran)
    !    ran = ran*range
    !    xdir(j+1) = (j-1)*bin+safety+ran
    !    call pullGauss(ran)
    !    ran = ran*var+baseline
    !    ydir(j+1) = ran
    ! end do
    ! xdir(pots+1) = width
    ! call pullGauss(ran)
    ! ran = ran*var+baseline
    ! ydir(pots+1) = ran
    
    ! Allocate
    allocate(potential(pots,2))
    allocate(point(pots,2))
    allocate(edgecoords(bndrs,4))
    allocate(edgeid(bndrs))
    allocate(potsid(pots))
    allocate(rileft(pots))
    allocate(riright(pots))
    allocate(abslist(pots))

    ! Set edgecoordinates of each line
    edgecoords(1,:) = sR(1,:)
    edgecoords(2,:) = sL(1,:)
    edgecoords(3,:) = sB(1,:)
    edgecoords(4,:) = sC(1,:)
    do j = 1, pots
       edgecoords(j+4,:) = [ xdir(j),ydir(j),xdir(j+1),ydir(j+1) ]
    end do

    ! Set: identification, vector, a point and id of each line
    edgeid(1) = "R   "
    edgeid(2) = "L   "
    edgeid(3) = "B   "
    edgeid(4) = "C   "
    do j = 1, pots
       potential(j,1) = xdir(j+1)-xdir(j)
       potential(j,2) = ydir(j+1)-ydir(j)
       point(j,1) = xdir(j)
       point(j,2) = ydir(j)
       write(pname, '("L", I0)') j
       edgeid(j+4) = pname
       potsid(j) = pname
       ! potsid = [ "L1","L2" ]
    end do

    ! Set index of refraction of each line
    rileft  = n1
    riright = n0

    ! Set which lines should turn absorption on/off
    abslist = potsid

    ! Store potential
    open(unit = 1, file = "tpotential1.txt", form='formatted')
    write(1,FMT='(2A15)') "x", "y"
    do j = 1, pots+1
       write(1,*) xdir(j), ydir(j)
    end do
    close(1)
    call system('gnuplot gsneakpeak.gp')

    if (present(k)) then
       write(filename, '("tpotential1-", I0, ".txt")') k
       open(unit = 1, file = filename, form='formatted')
       write(1,FMT='(2A15)') "x", "y"
       do j = 1, pots+1
          write(1,*) xdir(j), ydir(j)
       end do
       close(1)
    end if

    deallocate(xdir,ydir)
  end subroutine initialize1pnl

  ! Pull a gaussian distributed number
  subroutine pullGauss(gauss)
    real(wp), intent(inout)    :: gauss
    real(wp)                   :: u, v

    ! Pull a random number using the Box-Muller method
    call random_number(u)
    call random_number(v)
    gauss = sqrt(-2._wp*log(u))*cos(PI*2._wp*v)
  end subroutine pullGauss

  subroutine initialize2p3l(k)
    real(wp), dimension(2)    :: p1,p2,p3,p4,p5,l1,l2,l3
    integer(wp), optional                   :: k

    ! Set numbers
    pline1 = 2 ! Total number of potential lines in line 1
    pline2 = 1 ! Total number of potential lines in line 2
    pots   = pline1+pline2
    bndrs  = pots+4 ! Total number of lines
    
    call p2l3(p1,p2,p3,p4,p5,l1,l2,l3,k)
    
    ! Allocate
    allocate(potential(pots,2))
    allocate(point(pots,2))
    allocate(edgecoords(bndrs,4))
    allocate(edgeid(bndrs))
    allocate(potsid(pots))
    allocate(rileft(pots))
    allocate(riright(pots))

    ! Set
    potential = transpose(reshape([ l1,l2,l3 ], [2,pots])) ! List of direction vectors
    point = transpose(reshape([ p2,p3,p5 ], [2,pots]))    ! List of point on lines
    edgecoords = transpose(reshape([ sp1,sp2,sp4,sp3,sp4,sp1,sp3,sp2,p1,p2,p2,p3,p4,p5 ], [4,bndrs]))
    edgeid = [ "R ","L ","B ","C ","L1","L2","S1" ]
    potsid = [ "L1","L2","S1" ]
    rileft  = [ n1, n2, n2 ]
    riright = [ n0, n0, n1 ]
  end subroutine initialize2p3l

  subroutine p2l3(p1,p2,p3,p4,p5,l1,l2,l3,k)
    real(wp), dimension(2), intent(inout)  :: p1,p2,p3,p4,p5,l1,l2,l3
    real(wp)                               :: base
    integer(wp), optional                  :: k
    character(len=90)                      :: filename

    base = 0.7_wp
    ! Top layer
    p1 = [0._wp,base]
    p2 = [0.00001_wp-skew,0.3_wp]
    p3 = [width,base]
    print *, "Low point at: ", p2(1)

    ! Divide
    p4 = p2
    p5 = [width/2._wp,0.0_wp]

    ! Lines
    l1 = [p2(1)-p1(1),p2(2)-p1(2)]
    l2 = [p3(1)-p2(1),p3(2)-p2(2)]
    l3 = [p5(1)-p4(1),p5(2)-p4(2)]

    if (present(k)) then
       ! Print potentials
       write(filename, '("tpotential1-", I0, ".txt")') k
       ! print *, "filename: ", filename
       open(unit = 1, file = filename, form='formatted')
       write(1,FMT='(2A15)') "x", "y"
       write(1,*) p1(1), p1(2)
       write(1,*) p2(1), p2(2)
       write(1,*) p3(1), p3(2)
       close(1)

       write(filename, '("tpotential2-", I0, ".txt")') k
       ! print *, "filename: ", filename
       open(unit = 2, file = filename, form='formatted')
       write(2,FMT='(2A15)') "x", "y"
       write(2,*) p4(1), p4(2)
       write(2,*) p5(1), p5(2)
       close(2)
    end if

    ! Print potential
    open(unit = 1, file = "tpotential1.txt", form='formatted')
    write(1,FMT='(2A15)') "x", "y"
    write(1,*) p1(1), p1(2)
    write(1,*) p2(1), p2(2)
    write(1,*) p3(1), p3(2)
    close(1)

    open(unit = 2, file = "tpotential2.txt", form='formatted')
    write(2,FMT='(2A15)') "x", "y"
    write(2,*) p4(1), p4(2)
    write(2,*) p5(1), p5(2)
    close(2)
  end subroutine p2l3
  
  subroutine initialize2p2l(k)
    real(wp), dimension(2)    :: p1,p2,p3,p4,l1,l2
    integer(wp), optional     :: k

    ! Set numbers
    pline1 = 1 ! Total number of potential lines in line 1
    pline2 = 1 ! Total number of potential lines in line 2
    pots   = pline1+pline2
    bndrs  = pots+4 ! Total number of lines
    
    call p2l2(p1,p2,p3,p4,l1,l2,k)
    
    ! Allocate
    allocate(potential(pots,2))
    allocate(point(pots,2))
    allocate(edgecoords(bndrs,4))
    allocate(edgeid(bndrs))
    allocate(potsid(pots))
    allocate(rileft(pots))
    allocate(riright(pots))

    ! Set
    potential = transpose(reshape([ l1,l2 ], [2,pots])) ! List of direction vectors
    point = transpose(reshape([ p2,p4 ], [2,pots]))    ! List of point on lines
    edgecoords = transpose(reshape([ sp1,sp2,sp4,sp3,sp4,sp1,sp3,sp2,p1,p2,p3,p4 ], [4,bndrs]))
    edgeid = [ "R ","L ","B ","C ","L1","b1" ]
    potsid = [ "L1","b1" ]
    rileft  = [ n1, n2 ]
    riright = [ n0, n1 ]
  end subroutine initialize2p2l

  subroutine p2l2(p1,p2,p3,p4,l1,l2,k)
    real(wp), dimension(2), intent(inout)    :: p1,p2,p3,p4,l1,l2
    integer(wp), optional                    :: k
    character(len=90)                        :: filename

    ! Top layer
    p1  = [0.0_wp,0.3_wp]
    p2  = [0.5_wp,0.0_wp]
    
    !Bottom layer
    p3 = [0.0_wp,0.3_wp]
    p4 = [width,0.7_wp]

    ! Lines
    l1   = [p2(1)-p1(1),p2(2)-p1(2)]
    l2   = [p3(1)-p4(1),p3(2)-p4(2)]    

    if (present(k)) then
       ! Print potentials
       write(filename, '("tpotential1-", I0, ".txt")') k
       print *, "filename: ", filename
       open(unit = 1, file = filename, form='formatted')
       write(1,FMT='(2A15)') "x", "y"
       write(1,*) p1(1), p1(2)
       write(1,*) p2(1), p2(2)
       close(1)

       write(filename, '("tpotential2-", I0, ".txt")') k
       print *, "filename: ", filename
       open(unit = 2, file = filename, form='formatted')
       write(2,FMT='(2A15)') "x", "y"
       write(2,*) p3(1),  p3(2)
       write(2,*) p4(1),  p4(2)
       close(2)
    end if

    ! Print potential
    open(unit = 1, file = "tpotential1.txt", form='formatted')
    write(1,FMT='(2A15)') "x", "y"
    write(1,*) p1(1), p1(2)
    write(1,*) p2(1), p2(2)
    close(1)

    open(unit = 2, file = "tpotential2.txt", form='formatted')
    write(2,FMT='(2A15)') "x", "y"
    write(2,*) p3(1), p3(2)
    write(2,*) p4(1), p4(2)
    close(2)
  end subroutine p2l2

  subroutine initialize2p6l(k)
    real(wp), dimension(2)    :: p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14, &
         l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12
    integer(wp), optional     :: k

    ! Set numbers
    pline1 = 6 ! Total number of potential lines in line 1
    pline2 = 6 ! Total number of potential lines in line 2
    pots   = pline1+pline2
    bndrs  = pots+4 ! Total number of lines
    

    call p2l6(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14, &
         l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,k)
    
    ! Allocate
    allocate(potential(pots,2))
    allocate(point(pots,2))
    allocate(edgecoords(bndrs,4))
    allocate(edgeid(bndrs))
    allocate(potsid(pots))
    allocate(rileft(pots))
    allocate(riright(pots))

    ! Set
    potential = transpose(reshape([ l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12 ], [2,pots])) ! List of direction vectors
    point = transpose(reshape([ p1,p2,p3,p4,p5,p6,p8,p9,p10,p11,p12,p13 ], [2,pots]))    ! List of point on lines
    edgecoords = transpose(reshape([ sp1,sp2,sp4,sp3,sp4,sp1,sp3,sp2,p1,p2,p2,p3,p3,p4,p4,p5,p5,p6,p6,p7, &
         p8,p9,p9,p10,p10,p11,p11,p12,p12,p13,p13,p14], [4,bndrs]))
    edgeid = [ "R ","L ","B ","C ","L1","L2","L3","L4","L5","L6","b1","b2","b3","b4","b5","b6" ]
    potsid = [ "L1","L2","L3","L4","L5","L6","b1","b2","b3","b4","b5","b6" ]
    rileft  = [ n1, n1, n1, n1, n1, n1, n2, n2, n2, n2, n2, n2 ]
    riright = [ n0, n0, n0, n0, n0, n0, n1, n1, n1, n1, n1, n1 ]
  end subroutine initialize2p6l

  subroutine p2l6(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14, &
         l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,k)
    real(wp), dimension(2), intent(inout)    :: p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14, &
         l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12
    real(wp)    :: base1 = 1.0_wp*height/3.0_wp
    real(wp)    :: base2 = 1.0_wp*height/10._wp
    real(wp)    :: tip   = 1.0_wp*height/1.7_wp
    real(wp)    :: tip2  = 0.3_wp*height
    integer(wp), optional    :: k

    ! Top layer
    p1 = [0._wp*width/1.0_wp,base1]
    p2 = [1._wp*width/6.0_wp-skew*0.0_wp,tip-0.0_wp]
    p3 = [2._wp*width/6.0_wp,base1]
    p4 = [3._wp*width/6.0_wp-skew*0.0_wp,tip-0.0_wp]
    p5 = [4._wp*width/6.0_wp,base1]
    p6 = [5._wp*width/6.0_wp-skew*0.0_wp,tip+0.0_wp]
    p7 = [6._wp*width/6.0_wp,base1]
    
    !Bottom layer
    p8  = [0._wp*width/1.0_wp,base2]
    p9  = [1._wp*width/6.0_wp-skew*0.0_wp+0.0_wp,tip2]
    p10 = [2._wp*width/6.0_wp,base2+0.0_wp]
    p11 = [3._wp*width/6.0_wp-skew*0.0_wp+0.0_wp,tip2]
    p12 = [4._wp*width/6.0_wp,base2+0.0_wp]
    p13 = [5._wp*width/6.0_wp-skew*0.0_wp+0.0_wp,tip2]
    p14 = [6._wp*width/6.0_wp,base2]

    ! Lines
    l1   = [p2(1)-p1(1),p2(2)-p1(2)]
    l2   = [p3(1)-p2(1),p3(2)-p2(2)]
    l3   = [p4(1)-p3(1),p4(2)-p3(2)]
    l4   = [p5(1)-p4(1),p5(2)-p4(2)]
    l5   = [p6(1)-p5(1),p6(2)-p5(2)]
    l6   = [p7(1)-p6(1),p7(2)-p6(2)]

    l7   = [p9(1)-p8(1),p9(2)-p8(2)]
    l8   = [p10(1)-p9(1),p10(2)-p9(2)]
    l9   = [p11(1)-p10(1),p11(2)-p10(2)]
    l10  = [p12(1)-p11(1),p12(2)-p11(2)]
    l11  = [p13(1)-p12(1),p13(2)-p12(2)]
    l12  = [p14(1)-p13(1),p14(2)-p13(2)]

    ! Print potential
    open(unit = 1, file = "tpotential1.txt", form='formatted')
    write(1,FMT='(2A15)') "x", "y"
    write(1,*) p1(1), p1(2)
    write(1,*) p2(1), p2(2)
    write(1,*) p3(1), p3(2)
    write(1,*) p4(1), p4(2)
    write(1,*) p5(1), p5(2)
    write(1,*) p6(1), p6(2)
    write(1,*) p7(1), p7(2)
    close(1)

    open(unit = 2, file = "tpotential2.txt", form='formatted')
    write(2,FMT='(2A15)') "x", "y"
    write(2,*) p8(1),  p8(2)
    write(2,*) p9(1),  p9(2)
    write(2,*) p10(1), p10(2)
    write(2,*) p11(1), p11(2)
    write(2,*) p12(1), p12(2)
    write(2,*) p13(1), p13(2)
    write(2,*) p14(1), p14(2)
    close(2)
  end subroutine p2l6

  subroutine deallocatePotential(potential,point,edgecoords,edgeid,potsid,rileft,riright)
    real(wp), allocatable, intent(inout)           :: potential(:,:), point(:,:),edgecoords(:,:)
    real(wp), allocatable, intent(inout)           :: rileft(:), riright(:)
    character(len=4), allocatable, intent(inout)   :: edgeid(:), potsid(:)

    ! Deallocate
    deallocate(potential,point,edgecoords,edgeid,potsid,rileft,riright)
  end subroutine deallocatePotential
end module shapes
