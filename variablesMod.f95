module variablesMod
    integer, parameter :: int_prec = selected_int_kind(8)
    integer, parameter :: real_prec = selected_real_kind(8)

    ! Number of grid steps
    integer :: N1, N2, N, Ncol

    ! Integer counter and inflection point array position
    integer :: i, inflex1, inflex2

    ! Number of nodes a plot should have and the number it actually has
    integer :: nodes1, nodes2, nodecount, flag, nodecountflag

    ! Minimum radius
    real(kind = real_prec), parameter :: rmin=0.0

    !Stepsize for calculations
    real(kind = real_prec) :: deltaR = 0.005

    ! Vector of positions and potential
    real(kind = real_prec), allocatable :: r(:), v(:)

    ! Vector of Numerov variable
    real(kind = real_prec), allocatable :: a1(:), a2(:)

    ! Vectors of wavefunction of first configuration
    real(kind = real_prec), allocatable :: u1(:), u2(:), ut(:)

    ! Vectors of wavefunction of second configuration
    real(kind = real_prec), allocatable :: w1(:), w2(:), wt(:)

    ! Quantum numbers of 1st, 2nd configs and atomic number
    real(kind = real_prec) :: principalNumber1, principalNumber2, angularNumber1, angularNumber2, atomicNumber

    ! Max radius of 1st and 2nd configs and the max radius used throughout
    real(kind = real_prec) :: rmax1, rmax2, rmax

    ! Energy (in Hartree atomic units: m=hbar=q=1, c=137) and grid step spacing
    real(kind = real_prec) :: Ein1, Ein2, E1, E2, E

    ! Gradients from inward and outward integrations about inflex point
    real(kind = real_prec) :: grad11, grad22, grad12, grad21

    !radial range for collision calculation
    real(kind = real_prec), allocatable :: radial1(:), radial2(:)

    !variables for transition calculations
    real(kind = real_prec) :: linesum, linestrength, Ediff, wavelen
    real(kind = real_prec) :: A21, f12, g1, g2, temp
    integer(kind = int_prec) :: inflexflag, traflag
end module variablesMod