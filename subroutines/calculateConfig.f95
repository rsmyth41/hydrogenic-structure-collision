module configMod
    contains
    subroutine calculateConfig(principalNumber1, angularNumber1, atomicNumber, levelEnergy, totalVector, totalNumPoints)
        use subroutinesMod
        use variablesMod, only : real_prec
        implicit none

        real(kind = real_prec), intent(in) :: principalNumber1, angularNumber1, atomicNumber
        real(kind = real_prec), intent(out) :: levelEnergy
        real(kind = real_prec), intent(out), allocatable :: totalVector(:)
        integer, intent(out) :: totalNumPoints

        !! Local Variables
        ! Stepsize for calculations
        real(kind = real_prec), parameter :: deltaR = 0.005
        ! Vector of Numerov variable
        real(kind = real_prec), allocatable :: numerovVector(:)
        ! Vectors of wavefunction of first configuration
        real(kind = real_prec), allocatable :: firstVector(:), secondVector(:)
        ! Integer counter and inflection point array position
        integer :: i, inflectionPoint
        ! Number of nodes a plot should have and the number it actually has
        integer :: nodes, nodecount, flag, nodecountflag
        ! Vector of positions and potential
        real(kind = real_prec), allocatable :: radialGrid(:), potentialVector(:)
        ! Max radius of 1st and 2nd configs and the max radius used throughout
        real(kind = real_prec) :: rmax


        levelEnergy = -0.001
        rmax = ((3.0 * (principalNumber1 ** 2) - angularNumber1 * (angularNumber1 + 1.0))) * 4.0 / atomicNumber
        nodes = real(principalNumber1 - 1 - angularNumber1)

        totalNumPoints = real(rmax / deltaR)

        allocate(radialGrid(0: totalNumPoints))
        allocate(potentialVector(0: totalNumPoints))
        allocate(numerovVector(0: totalNumPoints))
        allocate(firstVector(0: totalNumPoints))
        allocate(secondVector(0: totalNumPoints))
        allocate(totalVector(0: totalNumPoints))

        do i = 1, totalNumPoints
            radialGrid(i) = real(i, 8) * deltaR
            potentialVector(i) = (-1.0 * atomicNumber) / radialGrid(i)
        end do
        print *, ' '

        ! Large repulsive potential close the origin
        potentialVector(0) = -100000000.0

        ! Initial grid position
        radialGrid(0) = 0.0

        ! Zero at the origin
        firstVector(0) = 0.0

        ! Very small close to the origin
        firstVector(1) = 0.000001

        ! Zero at infinity
        secondVector(totalNumPoints) = 0.0

        ! Small far from the origin
        secondVector(totalNumPoints - 1) = 0.0000001

        ! Very large close to the origin
        numerovVector(0) = 100000000.0

        ! Here we return from the node section of the code
2   continue

        !! Calculates array elements of the numerov vector !!
        do i = 1, totalNumPoints
            numerovVector(i) = 2.0 * (levelEnergy - potentialVector(i)) - (angularNumber1 * (angularNumber1 + 1.0)) &
                    & / (radialGrid(i) * radialGrid(i))
        end do

        !! Calculates inflection point from LHS !!
        call determineInflection(numerovVector, totalNumPoints, inflectionPoint)
        if (inflectionPoint == -1 .or. inflectionPoint == totalNumPoints) then
            levelEnergy = levelEnergy * 10.0
            goto 2
        end if

        !! Uses numerov method !!
        call numerovMethod(deltaR, inflectionPoint, totalNumPoints, numerovVector, firstVector, secondVector)

        !! Rescaling of secondVector(i) to connect with firstVector(i) !!
        call scaleVector(totalNumPoints, inflectionPoint, firstVector, secondVector)

        !! Counts nodes in our plot and shifts initial energy !!
        nodecount = 0
        call nodeChecker(nodecount, nodecountflag, inflectionPoint, firstVector, secondVector, nodes, levelEnergy, totalNumPoints)

        if (nodecountflag == 1) goto 2

        !! Calculates gradients and modifies energy !!
        call energyModification(inflectionPoint, totalNumPoints, firstVector, secondVector, radialGrid, flag, deltaR, levelEnergy)
        if (flag == 0) goto 2

        !! Constructing and normalising new total vector totalVector(i) !!
        call normaliseTotal(inflectionPoint, totalNumPoints, deltaR, firstVector, secondVector, totalVector)

        print *, 'Inflection point is at', radialGrid(inflectionPoint)
        !! Energy output in Ryd rel to lowest !!
        print 1000, 'Energy of n =', principalNumber1, 'l =', angularNumber1, 'config is', levelEnergy

        !! DEALLOCATE EVERYTHING HERE
        deallocate(numerovVector, firstVector, secondVector, radialGrid, potentialVector)

1000 format(1x, a, 1x, f3.1, 1x, a, 1x, f3.1, 1x, a, 1x, f15.7, /)

    end subroutine calculateConfig
end module configMod