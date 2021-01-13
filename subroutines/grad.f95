subroutine grad(inflex, N, u1, u2, r, flag, grad1, grad2, h, E)
    integer, intent(in) :: inflex, N
    real(kind=8), intent(inout) :: E
    real(kind=8), intent(in) :: u1(0: N), u2(0: N), r(0: N), h
    real(kind=8) ::  sum1, sum2, totalSum
    real(kind=8), intent(out) :: grad1, grad2
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

    sum1 = 0
    sum2 = 0
    totalSum = 0
    if (flag == 0) then
        do i = 0, inflex
            sum1 = sum1 + 0.5 * h * ((u1(i) ** 2) + (u1(i + 1) ** 2))
        end do
        do i = inflex, N - 1
            sum2 = sum2 + 0.5 * h * ((u2(i) ** 2) + (u2(i + 1) ** 2))
        end do
        totalSum = sum1 + sum2
        !! Approx from Atomic structure theory text !!
        E = E + ((grad1 - grad2) * u1(inflex)) / (2 * totalSum)
    end if

end subroutine grad