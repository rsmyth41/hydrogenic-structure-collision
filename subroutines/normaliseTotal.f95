subroutine normaliseTotal(inflectionPoint, totalNumPoints, deltaR, firstVector, secondVector, totalVector)
    use variablesMod, only : real_prec
    implicit none

    integer, intent(in) :: inflectionPoint, totalNumPoints
    real(kind = real_prec), intent(in) :: deltaR, firstVector(0: totalNumPoints), secondVector(0: totalNumPoints)
    real(kind = real_prec), intent(inout) :: totalVector(0: totalNumPoints)

    ! Local variables
    integer :: i
    real(kind = real_prec) :: normalisationConstant, areaUnderCurve

    totalVector(0: inflectionPoint) = firstVector(0: inflectionPoint)
    totalVector(inflectionPoint + 1: totalNumPoints) = secondVector(inflectionPoint + 1: totalNumPoints)

    areaUnderCurve = 0.0
    do i = 0, totalNumPoints - 1
        areaUnderCurve = areaUnderCurve + 0.5 * deltaR * ((totalVector(i) ** 2) + (totalVector(i + 1) ** 2))
    end do

    normalisationConstant = 1.0 / sqrt(areaUnderCurve)
    totalVector = totalVector * normalisationConstant
end subroutine normaliseTotal