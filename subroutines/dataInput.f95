subroutine dataInput(principalNumber1, principalNumber2, angularNumber1, angularNumber2, atomicNumber)
    use variablesMod, only : real_prec
    implicit none

    real(kind = real_prec), intent(out) :: principalNumber1, principalNumber2
    real(kind = real_prec), intent(out) :: angularNumber1, angularNumber2
    real(kind = real_prec), intent(out) :: atomicNumber

    print *, 'Enter the values of (n1 l1) of first configuration'
    100 continue
    read *, principalNumber1, angularNumber1
    if (angularNumber1 > principalNumber1 .or. angularNumber1 == principalNumber1 .or. angularNumber1 < 0.0) then
        print *, 'Error, check the configuration input'
        print *, 'Enter the values of (n1 l1) of first configuration'
        goto 100
    end if

    print *, 'Enter the values of (n2 l2) of second configuration'
    101 continue
    read *, principalNumber2, angularNumber2
    if (angularNumber2 > principalNumber2 .or. angularNumber2 == principalNumber2 .or. angularNumber2 < 0.0) then
        print *, 'Error, check the configuration input'
        print *, 'Enter the values of (n2 l2) of first configuration'
        goto 101
    end if

    print *, 'Enter the value of Z'
    read *, atomicNumber

end subroutine dataInput