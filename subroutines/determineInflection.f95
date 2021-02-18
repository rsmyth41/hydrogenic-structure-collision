subroutine determineInflection(numerovVector, totalNumPoints, inflectionPoint)
    use variablesMod, only : real_prec
    implicit none

    integer, intent(in) :: totalNumPoints
    integer, intent(out) :: inflectionPoint
    real(kind = real_prec), intent(in) :: numerovVector(0: totalNumPoints)

    ! Local variables
    integer :: i
    real(kind = real_prec) :: checkSignChange

    do i = totalNumPoints - 1, 0, -1
        checkSignChange = numerovVector(i) * numerovVector(i + 1)
        if (checkSignChange == 0.0 .or. checkSignChange < 0.0) then
            inflectionPoint = i
            return
        end if
    end do

    !! setting condition in case point of inflection cannot be determined due to low E !!
    inflectionPoint = -1
end subroutine determineInflection