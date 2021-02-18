subroutine scaleVector(totalNumPoints, inflectionPoint, firstVector, secondVector)
    use variablesMod, only : real_prec
    implicit none

    integer, intent(in) :: inflectionPoint, totalNumPoints
    real(kind = real_prec), intent(in) :: firstVector(0: totalNumPoints)
    real(kind = real_prec), intent(inout) :: secondVector(0: totalNumPoints)

    ! Local variables
    integer :: i
    real(kind = real_prec) :: constant, numerator, denominator

    numerator = firstVector(inflectionPoint)
    denominator = secondVector(inflectionPoint)

    if (denominator == 0) denominator = 1.0E-37

    constant = numerator / denominator
    secondVector(inflectionPoint - 2: totalNumPoints) = secondVector(inflectionPoint - 2: totalNumPoints) * constant
end subroutine scaleVector