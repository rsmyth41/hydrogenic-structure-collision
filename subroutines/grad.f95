subroutine grad(inflectionPoint, totalNumPoints, firstVector, secondVector, radialGrid, convergenceFlag, deltaR, levelEnergy)
    use variablesMod, only : real_prec
    implicit none

    integer, intent(in) :: inflectionPoint, totalNumPoints
    integer, intent(out) :: convergenceFlag
    real(kind = real_prec), intent(inout) :: levelEnergy
    real(kind = real_prec), intent(in) :: firstVector(0: totalNumPoints), secondVector(0: totalNumPoints)
    real(kind = real_prec), intent(in) :: radialGrid(0: totalNumPoints), deltaR

    ! Local variables
    integer :: i
    real(kind = real_prec) :: outwardSum, inwardSum, totalSum
    real(kind = real_prec) :: gradient1, gradient2, gradientDifference
    real(kind = real_prec) :: convergenceCondition = 1.0E-9

    gradient1 = (firstVector(inflectionPoint - 1) - firstVector(inflectionPoint)) &
            / (radialGrid(inflectionPoint - 1)- radialGrid(inflectionPoint))
    gradient2 = (secondVector(inflectionPoint + 1) - secondVector(inflectionPoint)) &
            / (radialGrid(inflectionPoint + 1) - radialGrid(inflectionPoint))
    gradientDifference = abs(gradient1 - gradient2)

    convergenceFlag = 0
    if (gradientDifference < convergenceCondition) then
        convergenceFlag = 1
    elseif (gradientDifference > convergenceCondition) then
        convergenceFlag = 0
    end if

    if (convergenceFlag == 0) then
        outwardSum = 0
        inwardSum = 0
        totalSum = 0
        do i = 0, inflectionPoint
            outwardSum = outwardSum + 0.5 * deltaR * ((firstVector(i) ** 2) + (firstVector(i + 1) ** 2))
        end do
        do i = totalNumPoints - 1, inflectionPoint, -1
            inwardSum = inwardSum + 0.5 * deltaR * ((secondVector(i) ** 2) + (secondVector(i + 1) ** 2))
        end do
        totalSum = outwardSum + inwardSum

        !! Approx from Atomic structure theory text !!
        levelEnergy = levelEnergy + ((gradient1 - gradient2) * firstVector(inflectionPoint)) / (2.0 * totalSum)
    end if
end subroutine grad