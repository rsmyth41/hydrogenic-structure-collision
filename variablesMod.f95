module variablesMod
    ! Integer variable precision :: Int*8
    integer, parameter :: int_prec = selected_int_kind(8)

    ! Real variable precision :: Real*8
    integer, parameter :: real_prec = selected_real_kind(8)

    ! Number of grid steps
    integer :: N1, N2, N, Ncol

    ! Integer counter and inflection point array position
    integer :: i

    ! Stepsize for calculations
    real(kind = real_prec), parameter :: deltaR = 0.005

    ! Vector of positions and potential
    real(kind = real_prec), allocatable :: radialGrid(:)

    ! Vectors of wavefunction of first configuration
    real(kind = real_prec), allocatable :: ut(:)

    ! Vectors of wavefunction of second configuration
    real(kind = real_prec), allocatable :: wt(:)

    ! Quantum numbers of 1st, 2nd configs and atomic number
    real(kind = real_prec) :: principalNumber1, principalNumber2, angularNumber1, angularNumber2, atomicNumber

    ! Energy (in Hartree atomic units: m=hbar=q=1, c=137) and grid step spacing
    real(kind = real_prec) :: E1, E2, E

    !radial range for collision calculation
    real(kind = real_prec), allocatable :: radial1(:), radial2(:)

    !variables for transition calculations
    real(kind = real_prec) :: linesum, linestrength, Ediff, wavelen
    real(kind = real_prec) :: A21, f12, g1, g2, temp
    integer(kind = int_prec) :: inflexflag, traflag
end module variablesMod