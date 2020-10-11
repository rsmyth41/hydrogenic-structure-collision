    program atom
    use subroutinesMod
    implicit none

    integer :: N1, N2, N, Ncol                                   ! Number of grid steps
    real*8, parameter :: rmin=0.0                   		     ! Minimum radius
    real*8, allocatable :: r(:), v(:)                            ! Vector of positions and potential
    real*8, allocatable :: a1(:), a2(:)               		     ! Vector of Numerov variable
    real*8, allocatable :: u1(:), u2(:), ut(:)      	         ! Vectors of wavefunction of first configuration
    real*8, allocatable :: w1(:), w2(:), wt(:)                   ! Vectors of wavefunction of second configuration
    real*8 :: pn1, pn2, l1, l2, Z                                ! Quantum numbers of 1st, 2nd configs and atomic number
    real*8 :: rmax1, rmax2, rmax                                 ! Max radius of 1st and 2nd configs and the max radius used throughout
    real*8 :: Ein1, Ein2, E1, E2, E, h                           ! Energy (in Hartree atomic units: m=hbar=q=1, c=137) and grid step spacing
    integer :: i, inflex1, inflex2                             ! Integer counter and inflection point array position
    integer :: nodes1, nodes2, nodecount, flag                 ! Number of nodes a plot should have and the number it actually has
    real*8 :: grad11, grad22, grad12, grad21                     ! Gradients from inward and outward integrations about inflex point
    real*8 :: sum1, sum2                                         ! Summations used for modifying energy
    real*8 :: linesum, linestrength, Ediff, wavelen
    real*8 :: A21, f12, g1, g2, temp
    integer*8 :: inflexflag, traflag
    real*8, allocatable :: radial1(:), radial2(:)               !radial range for collision calculation


!!!! ADD E2 EXPRESSIONS FOR A VALUES FROM NIST PAPER !!!!!
!!!! E2 LINE STRENGTH EXPRESSION IS DIFFERENT - CONTAINS POWER OF 2 TERM !!!!
!!!! USE THE VALUES OF 3j SYMBOLS FROM symbol_3j CODE TO DETERMINE THE REDUCED MATRIX ELEMENTS !!!!
!!!! POSSIBLY TRY LOOPING OVER ALL n and l up to approx 10h !!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!! DATA INPUT for nl configs !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    print *, 'Enter the values of (n1 l1) of first configuration'
100 continue
    read *, pn1, l1
    if(l1>pn1 .or. l1==pn1 .or. l1<0.0) then
        print *, 'Error, check the configuration input'
        print *, 'Enter the values of (n1 l1) of first configuration'
        goto 100
    end if
    Ein1=-0.001

    print*, 'Enter the values of (n2 l2) of second configuration'
101 continue
    read *, pn2, l2
    if(l2>pn2 .or. l2==pn2 .or. l2<0.0) then
        print *, 'Error, check the configuration input'
        print *, 'Enter the values of (n2 l2) of first configuration'
        goto 101
    end if
    Ein2=-0.001

    print *, 'Enter the value of Z'
    read *, Z

    !! Calculates MAXIMUM r for both nl configs !!
    rmax1=((3.0*pn1*pn1-l1*(l1+1.0)))*4.0/Z
    rmax2=((3.0*pn2*pn2-l2*(l2+1.0)))*4.0/Z

    !! Calculates the number of nodes each nl config should have !!
    nodes1=pn1-1-l1
    nodes2=pn2-1-l2

    h=0.005  !! Stepsize for calculations

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FIRST CONFIGURATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    N1=real(rmax1/h)
    print *, N1
    E=Ein1

    allocate(r(0:N1))
    allocate(v(0:N1))
    allocate(a1(0:N1))
    allocate(u1(0:N1))
    allocate(u2(0:N1))
    allocate(ut(0:N1))

    do i=1, N1
        r(i)=real(i,8)*h
        v(i)=(-1.0*Z)/r(i)
    end do
    print *, ' '

    v(0)=-100000000.0            ! Large close the origin
    r(0)=0.0                     ! Initial position
    u1(0)=0.0  			         ! Zero at the origin
    u1(1)=0.000001               ! Very small close to the origin
    u2(N1)=0.0                    ! Zero at infinity
    u2(N1-1)=0.0000001            ! Small far from the origin
    a1(0)=100000000.0

2   continue ! Returns from the node section of the code

    !! Calculates ARRAY ELEMENTS of a !!
    do i=1, N1
        a1(i)=2.0*(E-v(i)) - (l1*(l1+1.0))/(r(i)*r(i))
    end do

    !! Calculates INFLECTION POINT from LHS !!
    call inflection(a1, N1, inflex1)
    if(inflex1==-1 .or. inflex1==N1) then
        E=E*10.0
        goto 2
    end if

    !! Uses NUMEROV METHOD !!
    call numerov(h, inflex1, N1, a1, u1, u2)

    !! RESCALING of u2(i) to connect with u1(i) !!
    call scalevector(N1, inflex1, u1, u2)

    !! COUNTS NODES in our plot and SHIFTS INITITAL ENERGY !!
    nodecount=0
    do i=0, inflex1
        if(u1(i)*u1(i+1)<0) then
            nodecount=nodecount+1
        end if
    end do
    do i=N1-1, inflex1, -1
        if(u2(i)*u2(i+1)<0) then
            nodecount=nodecount+1
        end if
    end do
    if(nodecount==nodes1) then
        goto 3
    end if
    if(nodecount/=nodes1 .and. nodecount>nodes1) then
        E=E*1.1
    else if(nodecount/=nodes1 .and. nodecount<nodes1) then
        E=E*0.9
    end if
    goto 2
3   continue

    !! CALCULATING GRADIENTS and MODIFYING ENERGY !!
    call grad(inflex1, N1, u1, u2, r, flag, grad11, grad12)
    if(flag==0) then
        sum1=0.0
        sum2=0.0
        do i=0, inflex1
            sum1=sum1 + 0.5*h*((u1(i)**2) + (u1(i+1)**2))
        end do
        do i=inflex1, N1-1
            sum2=sum2 + 0.5*h*((u2(i)**2) + (u2(i+1)**2))
        end do
        sum1=sum1+sum2
        E=E+((grad11-grad12)*u1(inflex1))/(2*sum1)                !! Approx from Atomic structure theory text !!
        !print *, E
        goto 2
    end if

    !! Constructing and normalising new TOTAL VECTOR ut(i) !!
    call normtotal(inflex1, N1, h, u1, u2, ut)

    print *, 'Inflection point is at', r(inflex1)
    print 1000, 'Energy of n=', pn1, 'l=', l1, 'config is', E   !! Energy output in Ryd rel to lowest !!
    E1=E

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SECOND CONFIGURATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    deallocate(r,v)

    N2=real(rmax2/h)
    print *, N2
    E=Ein2

    allocate(r(0:N2))
    allocate(v(0:N2))
    allocate(a2(0:N2))
    allocate(w1(0:N2))
    allocate(w2(0:N2))
    allocate(wt(0:N2))

    do i=1, N2
        r(i)=real(i,8)*h
        v(i)=(-1.0*Z)/r(i)
    end do
    print *, ' '

    v(0)=-100000000.0            ! Large close the origin
    r(0)=0.0                     ! Initial position
    w1(0)=0.0
    w1(1)=0.000001
    w2(N2)=0.0
    w2(N2-1)=0.0000001
    a2(0)=100000000.0

4   continue

    do i=1, N2
        a2(i)=2.0*(E-v(i)) - (l2*(l2+1.0))/(r(i)*r(i))
    end do

    call inflection(a2, N2, inflex2)
    if(inflex2==-1 .or. inflex2==N) then
        E=E*10.0
    goto 4
    end if
    call numerov(h, inflex2, N2, a2, w1, w2)

    nodecount=0
    do i=0, inflex2
        if(w1(i)*w1(i+1)<0) then
            nodecount=nodecount+1
        end if
    end do
    do i=N2-1, inflex2, -1
        if(w2(i)*w2(i+1)<0) then
            nodecount=nodecount+1
        end if
    end do
    if(nodecount==nodes2) then
        goto 5
    end if
    if(nodecount/=nodes2 .and. nodecount>nodes2) then
        E=E*1.1
    else if(nodecount/=nodes2 .and. nodecount<nodes2) then
        E=E*0.9
    end if
    goto 4
5   continue

    call scalevector(N2, inflex2, w1, w2)
    call grad(inflex2, N2, w1, w2, r, flag, grad21, grad22)
    if(flag==0) then
        sum1=0.0
        sum2=0.0
        do i=0, inflex2
            sum1=sum1 + 0.5*h*((w1(i)**2) + (w1(i+1)**2))
        end do
        do i=inflex2, N2-1
            sum2=sum2 + 0.5*h*((w2(i)**2) + (w2(i+1)**2))
        end do
        sum1=sum1+sum2
        E=E+((grad21-grad22)*w1(inflex2))/(2*sum1)
        goto 4
    end if
  
    call normtotal(inflex2, N2, h, w1, w2, wt)

    print *, 'Inflection point is at', r(inflex2)
    print 1000, 'Energy of n=', pn2, 'l=', l2, 'config is', E  !! Energy output in Ryd rel to lowest !!
    E2=E

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!! TRANSITIONS & A VALUES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(N1>N2) then
        N=N2
    else if(N1<N2) then         !! THIS IF STATEMENT ENSURES THAT WE ONLY INTEGRATE OVER NON-ZERO !!
        N=N1                    !! OF THE "SMALLEST" RADIAL FUNCTION i.e. TO ENSURE WE DONT       !!
    end if                      !! EXCEED THE LIMIT OF THE ARRAY                                  !!

    deallocate(r)
    allocate(r(0:N))
    r(0)=0.0
    do i=1, N
        r(i)=real(i,8)*h
    end do

    traflag=0     !! A FLAG FOR DIFFERENT TYPES OF TRANSITION e.g. E1,E2 etc !!

    if(abs(l1-l2)>1.0 .or. l1-l2==0.0) then
        print *, 'Delta l must equal +/- 1 for E1 transition. No transition calculations carried out.'
        traflag=1
        goto 500
    end if

    linesum=0.0
    do i=0, N-1
        linesum= linesum + 0.5*h*((ut(i)*wt(i)*r(i)) + (ut(i+1)*wt(i+1)*r(i+1)))       !! E1 transition matrix element !!
    end do
    if(l1==l2+1.0) then
        linestrength= 2.0*(l2+1)*linesum*linesum         !! Prefactor of 2 !!
    else if(l1==l2-1.0) then                             !! Expressions here from Atomic Structure Theory text !!
        linestrength= 2.0*l2*linesum*linesum
    end if

    Ediff=abs(E1-E2)                                                           !! Energy difference between levels in au !!
    print *, 'Energy difference (atomic units) =', Ediff
    Ediff=Ediff*4.3597443419*10.0**(-18)                                       !! Energy level difference in J using conversion factor !!
    print *, 'Energy difference (Joules) =', Ediff
    wavelen=(6.62607004081*10.0**(-34)*2.99792458*10.0**8)/(Ediff)             !! Wavelength EM radiation in m. Using wavelen=hc/Edff !!
    wavelen=wavelen*(10.0**(10))                                               !! Converting to angstroms !!
    print *, 'Radiation wavelength (Angstroms) =', wavelen

    if(pn1>pn2) then                    !! Assigning stat weights based on the order they were entered above !!
        g2=(2.0)*(2.0*l1+1.0)           !! These commands ensure that data can be entered in any order !!
        g1=(2.0)*(2.0*l2+1.0)           !! Upper level denoted by g2, lower denoted by g1 !!
    else if(pn2>pn1) then
        g2=(2.0)*(2.0*l2+1.0)
        g1=(2.0)*(2.0*l1+1.0)
    else if(pn1==pn2 .and. l2>l1) then
        g2=(2.0)*(2.0*l2+1.0)
        g1=(2.0)*(2.0*l1+1.0)
    else if(pn1==pn2 .and.l1>l2) then
        g2=(2.0)*(2.0*l1+1.0)
        g1=(2.0)*(2.0*l2+1.0)
    end if

    if(E1>E2) then
        temp=E1                 !! These commands ensure that data can be entered in any order !!
        E1=E2                   !! If energy of lower is greater than energy of upper then reorder energies !!
        E2=temp                 !! E2 now energy of upper level and E1 now energy of lower level !!
    end if

    A21=(2.02613*(10.0**18)*linestrength)/((wavelen**3)*(g2))                    !! Transition probability/rate (upper -> lower) in s-1 !!
    f12=(2.0/3.0)*(E2-E1)*linestrength/g1                                        !! Oscillator strength (upper -> lower) dimensionless !!

    print *, 'Statistical weight of lower level: g1 =', g1
    print *, 'Statistical weight of upper level: g2 =', g2
    print *, 'Oscillator strength: f12 =', f12
    print *, 'Line strength (atomic units) =', linestrength                      !! Linestrength of the transition !!
    print *, 'A_ki values (s^-1): A21 =', A21
    print *, ' '

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FORMATTING & OUTPUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
500 continue

    deallocate(r)
    if(N1<N2) then
        Ncol=N2
        allocate(r(0:Ncol), radial1(0:Ncol), radial2(0:Ncol))
        r(0)=0.0
        do i=0, Ncol
            r(i)=real(i,8)*h
        end do
        do i=0, N1
            radial1(i)=ut(i)
            radial2(i)=wt(i)
        end do
        do i=N1+1, Ncol
            radial1(i)=0.0
            radial2(i)=wt(i)
        end do
    else if(N1>N2) then
        Ncol=N1
        allocate(r(0:Ncol), radial1(0:Ncol), radial2(0:Ncol))
        r(0)=0.0
        do i=0, Ncol
            r(i)=real(i,8)*h
        end do
        do i=0, N2
            radial1(i)=ut(i)
            radial2(i)=wt(i)
        end do
        do i=N2+1, Ncol
            radial1(i)=ut(i)
            radial2(i)=0.0
        end do
    end if

    open(unit=51, file='RADIAL1.DAT')
    open(unit=52, file='RADIAL2.DAT')
    do i=0, Ncol
        write(51,*) r(i), radial1(i)
        write(52,*) r(i), radial2(i)
    end do
    print *, Ncol
    close(51)
    close(52)

    open(unit=20, file='STRUCTURE.OUT')
        write(20,*) ' '
        write(20,*) '-----CONFIGURATIONS & ENERGIES-----'
        write(20,1005) 'Atomic number of ion Z = ', int(Z)
        write(20,1001) 'n', 'l', 'Energy'
        write(20,1002) int(pn1), int(l1), E1
        write(20,1002) int(pn2), int(l2), E2
        write(20,*) ' '
        write(20,*) '-----NUMERICAL GRID-----'
        write(20,1003) 'Number of grid points:', Ncol !The greatest number
        write(20,*) 'Grid step size:', h
        write(20,*) ' '
        write(20,*) '-----TRANSITION CALCULATIONS-----'
        if(traflag==0) then
            write(20,*) 'Radiation wavelength:', wavelen
            write(20,*) 'Line strength:', linestrength
            write(20,*) 'Transition rate:', A21
            write(20,*) 'Oscillator strength:', f12
        else if(traflag==1) then
            write(20,*) 'No transition calculations carried out.'
            goto 501
        end if
501     continue
        write(20,1004) 'Functions u(r) written to "RADIAL1.DAT" & "RADIAL2.DAT".'
    close(20)

    open(unit=30, file='COLLDATA.INP') ! File of input data to be used with collision code
        write(30,*) pn1, l1, E1
        write(30,*) pn2, l2, E2
        write(30,*) f12
        write(30,*) Ncol
        write(30,*) h
    close(30)

1000 format(1x, a, 1x, f3.1, 1x, a, 1x, f3.1, 1x, a, 1x, f15.7, /)
1001 format(1x, a, 3x, a, 4x, a)
1002 format(1x, i2, 2x, i2, 2x, f15.7)
1003 format(1x, a, i6)
1004 format(/, 1x, a, /)
1005 format(1x, a, i2)

    end program atom
