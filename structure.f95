    program hydrogenicAtomIon
    use subroutinesMod
    use variablesMod
    use configMod
    implicit none

    ! Data input for nl configs
    call dataInput(principalNumber1, principalNumber2, angularNumber1, angularNumber2, atomicNumber)

    ! Calculate energy and wavefunctions of nl configurations
    call calculateConfig(principalNumber1, angularNumber1, atomicNumber, E1, ut, N1)
    call calculateConfig(principalNumber2, angularNumber2, atomicNumber, E2, wt, N2)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!! TRANSITIONS & A VALUES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (N1 > N2) then
        N = N2
    else if (N1 < N2) then         !! THIS IF STATEMENT ENSURES THAT WE ONLY INTEGRATE OVER NON-ZERO !!
        N = N1                    !! OF THE "SMALLEST" RADIAL FUNCTION i.e. TO ENSURE WE DONT       !!
    end if                      !! EXCEED THE LIMIT OF THE ARRAY                                  !!

!    deallocate(radialGrid)
    allocate(radialGrid(0: N))
    radialGrid(0) = 0.0
    do i = 1, N
        radialGrid(i) = real(i, 8) * deltaR
    end do

    if(abs(angularNumber1 - angularNumber2) > 1.0 .or. angularNumber1 - angularNumber2 == 0.0) then
        print *, 'Delta l must equal +/- 1 for E1 transition. No transition calculations carried out.'
        traflag = 1
        goto 500
    end if

    linesum=0.0
    do i=0, N-1
        linesum= linesum + 0.5*deltaR*((ut(i)*wt(i)*radialGrid(i)) + (ut(i+1)*wt(i+1)*radialGrid(i+1)))       !! E1 transition matrix element !!
    end do
    if(angularNumber1==angularNumber2+1.0) then
        linestrength= 2.0*(angularNumber2+1)*linesum*linesum         !! Prefactor of 2 !!
    else if(angularNumber1==angularNumber2-1.0) then                             !! Expressions here from Atomic Structure Theory text !!
        linestrength= 2.0*angularNumber2*linesum*linesum
    end if

    Ediff=abs(E1-E2)                                                           !! Energy difference between levels in au !!
    print *, 'Energy difference (atomic units) =', Ediff
    Ediff=Ediff*4.3597443419*10.0**(-18)                                       !! Energy level difference in J using conversion factor !!
    print *, 'Energy difference (Joules) =', Ediff
    wavelen=(6.62607004081*10.0**(-34)*2.99792458*10.0**8)/(Ediff)             !! Wavelength EM radiation in m. Using wavelen=hc/Edff !!
    wavelen=wavelen*(10.0**(10))                                               !! Converting to angstroms !!
    print *, 'Radiation wavelength (Angstroms) =', wavelen

    if(principalNumber1>principalNumber2) then                    !! Assigning stat weights based on the order they were entered above !!
        g2=(2.0)*(2.0*angularNumber1+1.0)           !! These commands ensure that data can be entered in any order !!
        g1=(2.0)*(2.0*angularNumber2+1.0)           !! Upper level denoted by g2, lower denoted by g1 !!
    else if(principalNumber1<principalNumber2) then
        g2=(2.0)*(2.0*angularNumber2+1.0)
        g1=(2.0)*(2.0*angularNumber1+1.0)
    else if(principalNumber1==principalNumber2 .and. angularNumber2>angularNumber1) then
        g2=(2.0)*(2.0*angularNumber2+1.0)
        g1=(2.0)*(2.0*angularNumber1+1.0)
    else if(principalNumber1==principalNumber2 .and.angularNumber2<angularNumber1) then
        g2=(2.0)*(2.0*angularNumber1+1.0)
        g1=(2.0)*(2.0*angularNumber2+1.0)
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FORMATTING & OUTPUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

500 continue

    deallocate(radialGrid)
    if(N1<N2) then
        Ncol=N2
        allocate(radialGrid(0:Ncol), radial1(0:Ncol), radial2(0:Ncol))
        radialGrid(0)=0.0
        do i=0, Ncol
            radialGrid(i)=real(i,8)*deltaR
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
        allocate(radialGrid(0:Ncol), radial1(0:Ncol), radial2(0:Ncol))
        radialGrid(0)=0.0
        do i=0, Ncol
            radialGrid(i)=real(i,8)*deltaR
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
        write(51,*) radialGrid(i), radial1(i)
        write(52,*) radialGrid(i), radial2(i)
    end do
    close(51)
    close(52)

    open(unit=20, file='STRUCTURE.OUT')
        write(20,*) ' '
        write(20,*) '-----CONFIGURATIONS & ENERGIES-----'
        write(20,1005) 'Atomic number of ion Z = ', int(atomicNumber)
        write(20,1001) 'n', 'l', 'Energy'
        write(20,1002) int(principalNumber1), int(angularNumber1), E1
        write(20,1002) int(principalNumber2), int(angularNumber2), E2
        write(20,*) ' '
        write(20,*) '-----NUMERICAL GRID-----'
        write(20,1003) 'Number of grid points:', Ncol !The greatest number
        write(20,*) 'Grid step size:', deltaR
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
        write(30,*) principalNumber1, angularNumber1, E1
        write(30,*) principalNumber2, angularNumber2, E2
        write(30,*) f12
        write(30,*) Ncol
        write(30,*) deltaR
    close(30)

1001 format(1x, a, 3x, a, 4x, a)
1002 format(1x, i2, 2x, i2, 2x, f15.7)
1003 format(1x, a, i6)
1004 format(/, 1x, a, /)
1005 format(1x, a, i2)

    end program hydrogenicAtomIon
