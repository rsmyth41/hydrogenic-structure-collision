    program collision
    implicit none

    real*8, allocatable :: r(:), u(:), w(:), q(:), Eplot(:), sigma(:), totalsum(:)
    integer :: i, j, k, RANGE, RANGEB, N, filestat1, filestat2, collinput
    real*8 :: pn1, l1, pn2, l2, dr, E1, E2, f12, EB1, EB2
    real*8 :: Uex, Eelec, A, c, dE, EMAX, EMIN
    real*8 :: EINIT, EFIN, q1, q2, dq, qsum, EBORN1, EBORN2
    real*8 :: bessel, rtemp1, rtemp2, qtemp1, qtemp2
    real*8 :: EMAXB, dEB
    real*8 :: omegasum1, omegasum2, collstrength, colsum

!!!! USE FIRST BORN APPROXIMATION FROM XeP PAPER !!!!
!!!! BORN APPROXIMATION FROM YOUNGER PAPER !!!!
!!!! USE BETHE-BORN APPROXIMATION FROM SCHRAM PAPER AND COMPARE !!!!
!!!! ADD CODE TO CALCULATE TIME TAKEN FOR SIGMA FILE CREATION !!!!

    open(unit=100, file='COLLDATA.INP', status='old', iostat=collinput)
    if(collinput>0) then
        print *, 'Input file "COLLDATA.INP" could not be opened.'
        stop
    else if(collinput==0) then
        print *, 'File "COLLDATA.INP" successfully opened'
    end if
    read(100,*) pn1, l1, E1
    read(100,*) pn2, l2, E2
    read(100,*) f12
    read(100,*) N
    read(100,*) dr
    close(100)

    allocate(r(0:N), u(0:N), w(0:N), q(0:N))
    allocate(totalsum(1:N-1))

    open(unit=10, file='RADIAL1.DAT', status='old', iostat=filestat1)
    open(unit=20, file='RADIAL2.DAT', status='old', iostat=filestat2)

    if(filestat1>0) then
        print *, 'File "RADIAL1.DAT" could not be opened.'
        stop
    else if(filestat1==0) then
        print *, 'File "RADIAL1.DAT" successfully opened.'
    end if
    if(filestat2>0) then
        print *, 'File "RADIAL2.DAT" could not be opened.'
        stop
    else if (filestat2==0) then
        print *, 'File "RADIAL2.DAT" successfully opened.'
    end if
    print *, ' '

    do i=0, N
        read(10,*) r(i), u(i)
        read(20,*) r(i), w(i)
    end do
    close(10)
    close(20)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!! BETHE-BORN CROSS SECTION !!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    EB1=E1*2.0*13.6057 !convert au -> Ryd -> eV
    EB2=E2*2.0*13.6057
    Uex=abs(EB1-EB2) !Excitation energy
    !print *, Uex !This is 100% fine
    A=4.0*(13.6057**2.0)*f12/Uex  ! OSCILLATOR STRENGTH MAY NEED CHECKED FOR SPECIFIC RUNS !
    c=(1.0/Uex)*((0.984*(pn1+0.08-(0.3/pn1))**(-1.0))*((1.0-(l1/pn1)+0.18*(l1*(l1-1.0))/(pn1)**2.0))**(-1.0))
    !print *, c*Uex

    RANGE=500
    EMAX=100.0
    EMIN=Uex
    dE=abs(EMAX-EMIN)/real(RANGE)

    allocate(Eplot(1:RANGE), sigma(1:RANGE))

    Eelec=EMIN
    open(unit=21, file='SIGMA_BETHE.DAT')
    do i=1, RANGE
        sigma(i)=A*((Eelec-Uex)/(Eelec**2.0))*log(1.0+c*(Eelec-Uex))
        Eplot(i)=Eelec
        write(21,*) Eplot(i), sigma(i)
        Eelec=Eelec+dE
    end do
    print *, 'Bethe-Born cross section written to SIGMA-BETHE.DAT file.'
    close(21)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! BORN APPROXIMATION !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!u(i)*w(i)*bessel(rtemp,qtemp))/(qtemp**(1.5)
!area=area + 0.5*h*((ut(i)**2) + (ut(i+1)**2))

!    EBORN1=E1*2.0*13.6057
!    EBORN2=E2*2.0*13.6057
    EBORN1=E1  !Working in Rydbergs
    EBORN2=E2
    EINIT=abs(EBORN1-EBORN2) !the excitation energy
    EMAXB=100.0/(2.0*13.6)
    RANGEB=50 !Number of energy points (SET AT 100)
    dEB=abs(EMAXB-EINIT)/real(RANGEB)

    open(unit=22, file='SIGMA_BORN.DAT')
    do k=1, RANGEB
        EFIN=EINIT-abs(EBORN1-EBORN2)
        q1=sqrt(2.0*EINIT)
        q2=sqrt(2.0*EFIN)
        dq=(abs(q2+q1)-abs(q2-q1))/real(N)
        qsum=abs(q2-q1)
        do i=0, N
            q(i)=qsum
            qsum=qsum+dq
        end do

        do i=1, N-1
            omegasum1=0.0
            do j=1, N-1
                qtemp1=q(i)
                rtemp1=r(j)
                rtemp2=r(j+1)
                omegasum1=omegasum1+0.5*dr*(u(j)*w(j)*(bessel(rtemp1,qtemp1)/(qtemp1**(1.5))) &
                    +u(j+1)*w(j+1)*bessel(rtemp2,qtemp1)/(qtemp1**(1.5)))
            end do
        totalsum(i)=abs(omegasum1)**(2.0)
        end do

        colsum=0.0
        do i=1, N-2
            colsum=colsum+0.5*dq*(totalsum(i)+totalsum(i+1))
        end do

        collstrength=24.0*colsum                ! Collision Strengths
        collstrength=collstrength/(2.0*EINIT)   ! working out cross-section
        write(22,*) EINIT*2.0*13.6057, collstrength     !Convering energy Ryd -> eV
        EINIT=EINIT+dEB
    end do
    close(22)

    print *, 'Born cross section written to SIGMA-BORN.DAT file.'
    print *, ' '
!    print *, EINIT, collstrength
    end program collision


!!!!!!!!!!!!!! FUNCTION !!!!!!!!!!!!!!!!!!!!
    real*8 function bessel(rtemp, qtemp)
    real*8 :: rtemp, qtemp, var
    var=rtemp*qtemp
    bessel=(sin(var)/(var*var))-cos(var)/var
    end function bessel


!!!!!!!!!!!!!!!!!!!!!!!!
!!!! END OF PROGRAM !!!!
!!!!!!!!!!!!!!!!!!!!!!!!
