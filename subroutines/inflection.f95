subroutine inflection(a, N, inflex)
    real(kind=8), intent(in) :: a(0: N)
    integer, intent(in) :: N
    integer, intent(out) :: inflex
    integer :: i

    do i = N - 1, 0, -1
        if (a(i) * a(i + 1) == 0.0) then
            inflex = i
            goto 1
        end if

        if (a(i) *a (i + 1) < 0.0) then
            inflex = i
            goto 1
        end if
    end do

    !! setting condition in case point of inflection cannot be determined due to low E !!
    inflex = -1 

1   continue
end subroutine inflection