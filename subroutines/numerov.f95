subroutine numerov(h, inflex, N, a, u1, u2)
    integer, intent(in) :: inflex, N
    real*8, intent(in) :: h, a(0: N)
    real*8, intent(inout) :: u1(0: N), u2(0: N)
    real*8 :: c1, c2, numerator, denominator
    integer*8 :: i

    c1 = 5.0 * (h * h / 6.0)
    c2 = h * h/ 12.0

    ! Forward integration
    do i = 1, inflex + 1   
        numerator = ((2.0 - c1 * a(i)) * u1(i) - (1.0 + c2 *a (i - 1)) * u1(i - 1))
        denominator = (1.0 + c2 * a(i + 1))
        u1(i + 1) = numerator / denominator
    end do

    ! Backwards integration
    do i = N, inflex - 1, -1 
        numerator = ((2.0 - c1 * a(i - 1)) * u2(i - 1)-(1.0+  c2 * a(i)) * u2(i))
        denominator = (1.0 + c2 * a(i - 2))
        u2(i - 2) = numerator / denominator
    end do
end subroutine numerov