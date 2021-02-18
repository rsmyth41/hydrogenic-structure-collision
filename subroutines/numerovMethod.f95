subroutine numerovMethod(deltaR, inflectionPoint, totalNumPoints, numerovVector, firstVector, secondVector)
    use variablesMod, only : real_prec, int_prec
    implicit none

    integer, intent(in) :: inflectionPoint, totalNumPoints
    real(kind = real_prec), intent(in) :: deltaR, numerovVector(0: totalNumPoints)
    real(kind = real_prec), intent(inout) :: firstVector(0: totalNumPoints), secondVector(0: totalNumPoints)

    ! Local variables
    integer :: i
    real(kind = real_prec) :: constant1, constant2, numerator, denominator

    constant1 = 5.0 * (deltaR * deltaR / 6.0)
    constant2 = deltaR * deltaR / 12.0

    ! Outwards integration
    do i = 1, inflectionPoint + 1
        numerator = ((2.0 - constant1 * numerovVector(i)) * firstVector(i) - &
                (1.0 + constant2 * numerovVector(i - 1)) * firstVector(i - 1))
        denominator = (1.0 + constant2 * numerovVector(i + 1))
        firstVector(i + 1) = numerator / denominator
    end do

    ! Inwards integration
    do i = totalNumPoints, inflectionPoint - 1, -1
        numerator = ((2.0 - constant1 * numerovVector(i - 1)) * secondVector(i - 1) - &
                (1.0 + constant2 * numerovVector(i)) * secondVector(i))
        denominator = (1.0 + constant2 * numerovVector(i - 2))
        secondVector(i - 2) = numerator / denominator
    end do
end subroutine numerovMethod