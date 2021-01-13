subroutine nodeChecker(nodecount, nodecountflag, inflex, u1, u2, nodes, E, N)
    integer, intent(in) :: inflex, N, nodes
    integer, intent(inout) :: nodecount
    integer, intent(out) :: nodecountflag
    real(kind=8), intent(inout) :: u1(0: N), u2(0: N), E

    nodecountflag = 0

    do i = 0, inflex
        if (u1(i) * u1(i + 1) < 0) then
            nodecount = nodecount + 1
        end if
    end do

    do i = N - 1, inflex, -1
        if (u2(i) * u2(i + 1) < 0) then
            nodecount = nodecount + 1
        end if
    end do

    if (nodecount == nodes) then
        goto 1
    end if

    if(nodecount /= nodes .and. nodecount > nodes) then
        E = E * 1.1
    elseif (nodecount /= nodes .and. nodecount < nodes) then
        E = E * 0.9
    end if

    nodecountflag = 1
1   continue

end subroutine nodeChecker
