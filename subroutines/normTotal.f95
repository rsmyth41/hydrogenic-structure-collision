subroutine normtotal(inflex, N, h, u1, u2, ut)
    use variablesMod, only : real_prec

    integer, intent(in) :: inflex, N
    integer :: i
    real(kind = real_prec), intent(in) :: h, u1(0: N), u2(0: N)
    real(kind = real_prec), intent(inout) :: ut(0: N)
    real(kind = real_prec) :: norm, area
    
    do i = 0, inflex
        ut(i) = u1(i)
    end do

    do i = inflex + 1, N
        ut(i) = u2(i)
    end do

    area = 0.0

    do i = 0, N - 1
        area = area + 0.5 * h * ((ut(i) ** 2) + (ut(i+1) ** 2))
    end do

    norm = 1.0 / sqrt(area)

    print *, 'Area = ', area
    print *, 'Normalisation constant = ', norm
    
    ut = ut * norm
end subroutine normtotal