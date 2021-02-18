subroutine nodeChecker(nodeCount, nodeCountFlag, inflectionPoint, firstVector, secondVector, numNodes, levelEnergy, totalNumPoints)
    use variablesMod, only : real_prec
    implicit none

    integer, intent(in) :: inflectionPoint, totalNumPoints, numNodes
    integer, intent(inout) :: nodeCount
    integer, intent(out) :: nodeCountFlag
    real(kind = real_prec), intent(inout) :: firstVector(0: totalNumPoints), secondVector(0: totalNumPoints), levelEnergy

    ! Local variables
    integer :: i

    nodeCountFlag = 0
    do i = 0, inflectionPoint
        if (firstVector(i) * firstVector(i + 1) < 0) then
            nodeCount = nodeCount + 1
        end if
    end do

    do i = totalNumPoints - 1, inflectionPoint, -1
        if (secondVector(i) * secondVector(i + 1) < 0) then
            nodeCount = nodeCount + 1
        end if
    end do

    if (nodeCount == numNodes) then
        return
    end if

    if(nodeCount /= numNodes .and. nodeCount > numNodes) then
        levelEnergy = levelEnergy * 1.1
    elseif (nodeCount /= numNodes .and. nodeCount < numNodes) then
        levelEnergy = levelEnergy * 0.9
    end if

    nodeCountFlag = 1
end subroutine nodeChecker
