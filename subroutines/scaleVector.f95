subroutine scalevector(N, inflex, u1, u2)
    integer, intent(in) :: inflex, N
    real(kind=8), intent(in) :: u1(0: N)
    real(kind=8), intent(inout) :: u2(0: N)
    real(kind=8) :: scale, x, y

    x = u1(inflex)
    y = u2(inflex)

    if (y == 0) then
        y = 1.0E-37
    end if

    scale = x / y

    do i = N, inflex - 2, -1
        u2(i) = u2(i) * scale
    end do
end subroutine scalevector