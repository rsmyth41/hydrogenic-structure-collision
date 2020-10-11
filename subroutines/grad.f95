subroutine grad(inflex, N, u1, u2, r, flag, grad1, grad2)
    integer, intent(in) :: inflex, N
    real*8, intent(in) :: u1(0: N), u2(0: N), r(0: N)
    real*8, intent(out) :: grad1, grad2
    integer, intent(out) :: flag

    grad1 = (u1(inflex - 1) - u1(inflex))/(r(inflex - 1) - r(inflex))
    grad2 = (u2(inflex + 1) - u2(inflex))/(r(inflex + 1) - r(inflex))
    
    flag = 0

    !! setting a small condition of 10**-10 !!
    if (abs(grad1 - grad2) < 0.000000001) then
        flag = 1
    elseif (abs(grad1 - grad2) > 0.000000001) then
        flag = 0
    end if
end subroutine grad