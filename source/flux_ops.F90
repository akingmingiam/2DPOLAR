!===========================================================
! File   : flux_ops.F90
! Module : flux_ops
! Purpose:
!   Assemble numerical fluxes on all faces of the polar mesh
!   using the HLL Riemann solver.
!
!   - compute_radial_fluxes  : faces normal to e_r
!   - compute_angular_fluxes : faces normal to e_theta
!   - compute_all_fluxes     : convenience wrapper
!===========================================================
module flux_ops
    use mesh_types,         only : pi, PolarMesh, EQ, RECON
    use flow_fields,        only : PrimitiveVariables, ConservedVariables
    use fluid_properties,   only : StiffenedGas
    use flux_types,         only : FluxFields
    use reconstruction_ops, only : muscl_reconstruct
    use riemann_solver,     only : flux_solver

    implicit none
    private
    public :: compute_all_fluxes

contains

    !-----------------------------------------------------------
    ! Convenience wrapper: compute both radial and angular fluxes
    !-----------------------------------------------------------
    subroutine compute_all_fluxes(mesh, prim, gas, flux)
        type(PolarMesh),          intent(in)    :: mesh
        type(PrimitiveVariables), intent(in)    :: prim
        type(StiffenedGas),       intent(in)    :: gas
        type(FluxFields),         intent(inout) :: flux

        select case(mesh%params%equation_type)

        case (EQ%CARTESIAN_TRANSFORM)
            call compute_radial_fluxes_metric (mesh, prim, gas, flux)
            call compute_angular_fluxes_rot   (mesh, prim, gas, flux)

        case (EQ%POLAR_EULER)
            call compute_radial_fluxes_polar (mesh, prim, gas, flux)
            call compute_angular_fluxes_polar(mesh, prim, gas, flux)

        case default
            print *, "ERROR: Unknown equation type in compute_all_fluxes"
            stop

        end select

    end subroutine compute_all_fluxes
    
    !-----------------------------------------------------------
    !  Compute numerical fluxes for the metric-form Euler equations.
    !  Indexing:
    !    i = i_lo_phys ... i_hi_phys+1    (faces)
    !    j = j_lo_phys ... j_hi_phys      (cells in theta)
    !-----------------------------------------------------------

    subroutine compute_radial_fluxes(mesh, prim, gas, flux)
        type(PolarMesh),          intent(in)    :: mesh
        type(PrimitiveVariables), intent(in)    :: prim
        type(StiffenedGas),       intent(in)    :: gas
        type(FluxFields),         intent(inout) :: flux

        integer :: i, j
        integer :: iL, iR, jL, jR
        double precision :: nx, ny
        double precision :: rhoL,uL,vL,pL
        double precision :: rhoR,uR,vR,pR
        double precision :: qst(4)
        double precision :: fmass, fmom_u, fmom_v, fenergy

        iL = mesh%params%i_lo_phys
        iR = mesh%params%i_hi_phys
        jL = mesh%params%j_lo_phys
        jR = mesh%params%j_hi_phys

        do j = jL, jR
            do i = iL, iR+1
                ! left cell:  (i-1, j)
                ! right cell: (i,   j)
                nx = cos(mesh%theta_center(j))
                ny = sin(mesh%theta_center(j))

                select case (mesh%params%reconstruction)

                case (RECON%FIRST_ORDER)
                    ! ------------------------------------
                    ! First-order: piecewise constant
                    ! left cell  = (i-1, j)
                    ! right cell = (i,   j)
                    ! ------------------------------------
                    rhoL = prim%density(   i-1, j)
                    uL   = prim%velocity_u(i-1, j)
                    vL   = prim%velocity_v(i-1, j)
                    pL   = prim%pressure(  i-1, j)

                    rhoR = prim%density(   i, j)
                    uR   = prim%velocity_u(i, j)
                    vR   = prim%velocity_v(i, j)
                    pR   = prim%pressure(  i, j)

                case (RECON%MUSCL_MINMOD, RECON%MUSCL_VANLEER)
                    if (mesh%params%guard_cells < 2) then
                        write(*,*) "ERROR: MUSCL needs at least 2 ghost cells."
                        stop
                    end if

                    ! ------------------------------------
                    ! Second-order MUSCL (minmod limiter)
                    ! Interface: i-1/2
                    ! Stencil for cell-centered q:
                    !   q_{i-2}, q_{i-1}, q_i, q_{i+1}
                    ! Near boundaries: fall back to first order.
                    ! ------------------------------------
                    if (i <= iL+1) then
                        ! Not enough neighbors → use first-order
                        rhoL = prim%density(   i-1, j)
                        uL   = prim%velocity_u(i-1, j)
                        vL   = prim%velocity_v(i-1, j)
                        pL   = prim%pressure(  i-1, j)

                        rhoR = prim%density(   i, j)
                        uR   = prim%velocity_u(i, j)
                        vR   = prim%velocity_v(i, j)
                        pR   = prim%pressure(  i, j)
                    else
                        qst(1:4) = prim%density(i-2:i+1, j)
                        call muscl_reconstruct(qst, rhoL, rhoR, mesh%params%reconstruction)
                        qst(1:4) = prim%velocity_u(i-2:i+1, j)
                        call muscl_reconstruct(qst, uL, uR, mesh%params%reconstruction)
                        qst(1:4) = prim%velocity_v(i-2:i+1, j)
                        call muscl_reconstruct(qst, vL, vR, mesh%params%reconstruction)
                        qst(1:4) = prim%pressure(i-2:i+1, j)
                        call muscl_reconstruct(qst, pL, pR, mesh%params%reconstruction)
                    end if

                end select

                call flux_solver(mesh%params%riemann_solver, rhoL,uL,vL,pL, rhoR,uR,vR,pR, &
                                 nx,ny, gas, fmass, fmom_u, fmom_v, fenergy)

                ! multiply by face length
                flux%radial%mass(i,j)   = fmass   * mesh%radial_face_length(i,j)
                flux%radial%mom_u(i,j)  = fmom_u  * mesh%radial_face_length(i,j)
                flux%radial%mom_v(i,j)  = fmom_v  * mesh%radial_face_length(i,j)
                flux%radial%energy(i,j) = fenergy * mesh%radial_face_length(i,j)
                
            end do
        end do

    end subroutine compute_radial_fluxes

    subroutine compute_angular_fluxes(mesh, prim, gas, flux)
        type(PolarMesh),          intent(in)    :: mesh
        type(PrimitiveVariables), intent(in)    :: prim
        type(StiffenedGas),       intent(in)    :: gas
        type(FluxFields),         intent(inout) :: flux

        integer :: i, j
        integer :: iL, iR, jL, jR
        double precision :: nx, ny
        double precision :: rhoL,uL,vL,pL
        double precision :: rhoR,uR,vR,pR
        double precision :: qst(4)
        double precision :: fmass, fmom_u, fmom_v, fenergy

        iL = mesh%params%i_lo_phys
        iR = mesh%params%i_hi_phys
        jL = mesh%params%j_lo_phys
        jR = mesh%params%j_hi_phys

        do j = jL, jR+1
            do i = iL, iR
                ! left cell:  (i, j-1)
                ! right cell: (i, j)

                nx = - sin(mesh%theta_face(j))
                ny =   cos(mesh%theta_face(j))

                select case (mesh%params%reconstruction)

                case (RECON%FIRST_ORDER)
                    ! First-order in theta:
                    ! left cell  = (i, j-1)
                    ! right cell = (i, j)
                    rhoL = prim%density(   i, j-1)
                    uL   = prim%velocity_u(i, j-1)
                    vL   = prim%velocity_v(i, j-1)
                    pL   = prim%pressure(  i, j-1)

                    rhoR = prim%density(   i, j)
                    uR   = prim%velocity_u(i, j)
                    vR   = prim%velocity_v(i, j)
                    pR   = prim%pressure(  i, j)

                case (RECON%MUSCL_MINMOD, RECON%MUSCL_VANLEER)
                    if (mesh%params%guard_cells < 2) then
                        write(*,*) "ERROR: MUSCL needs at least 2 ghost cells."
                        stop
                    end if

                    ! Stencil in theta: j-2, j-1, j, j+1

                    qst(1:4) = prim%density(i, j-2:j+1)
                    call muscl_reconstruct(qst, rhoL, rhoR, mesh%params%reconstruction)
                    qst(1:4) = prim%velocity_u(i, j-2:j+1) 
                    call muscl_reconstruct(qst, uL, uR, mesh%params%reconstruction)
                    qst(1:4) = prim%velocity_v(i, j-2:j+1)
                    call muscl_reconstruct(qst, vL, vR, mesh%params%reconstruction)
                    qst(1:4) = prim%pressure(i, j-2:j+1)
                    call muscl_reconstruct(qst, pL, pR, mesh%params%reconstruction)

                end select

                call flux_solver(mesh%params%riemann_solver, rhoL,uL,vL,pL, rhoR,uR,vR,pR, &
                                 nx,ny, gas, fmass, fmom_u, fmom_v, fenergy)
                
                ! multiply by face length
                flux%angular%mass(i,j)   = fmass   * mesh%angular_face_length(i,j)
                flux%angular%mom_u(i,j)  = fmom_u  * mesh%angular_face_length(i,j)
                flux%angular%mom_v(i,j)  = fmom_v  * mesh%angular_face_length(i,j)
                flux%angular%energy(i,j) = fenergy * mesh%angular_face_length(i,j)

            end do
        end do
    end subroutine compute_angular_fluxes


    subroutine compute_angular_fluxes_rot(mesh, prim, gas, flux)
        type(PolarMesh),          intent(in)    :: mesh
        type(PrimitiveVariables), intent(in)    :: prim
        type(StiffenedGas),       intent(in)    :: gas
        type(FluxFields),         intent(inout) :: flux

        integer :: i, j
        integer :: iL, iR, jL, jR
        double precision :: nx, ny
        double precision :: rhoL,uL,vL,pL
        double precision :: rhoR,uR,vR,pR
        double precision :: qst(4)
        double precision :: fmass, fmom_u, fmom_v, fenergy

        integer :: iL_comp, iR_comp, jL_comp, jR_comp
        double precision :: ctc, stc
        double precision :: urL,utL,urR,utR
        double precision, allocatable :: u_r_cell(:,:), u_t_cell(:,:)

        iL = mesh%params%i_lo_phys
        iR = mesh%params%i_hi_phys
        jL = mesh%params%j_lo_phys
        jR = mesh%params%j_hi_phys
        
        iL_comp = mesh%params%i_lo_comp
        iR_comp = mesh%params%i_hi_comp
        jL_comp = mesh%params%j_lo_comp
        jR_comp = mesh%params%j_hi_comp

        allocate(u_r_cell(iL_comp:iR_comp, jL_comp:jR_comp))
        allocate(u_t_cell(iL_comp:iR_comp, jL_comp:jR_comp))

        do j = jL_comp, jR_comp
            ctc = cos(mesh%theta_center(j))
            stc = sin(mesh%theta_center(j))

            do i = iL_comp, iR_comp
                u_r_cell(i,j) =  prim%velocity_u(i,j) * ctc + prim%velocity_v(i,j) * stc
                u_t_cell(i,j) = -prim%velocity_u(i,j) * stc + prim%velocity_v(i,j) * ctc
            end do
        end do



        do j = jL, jR+1
            do i = iL, iR
                ! left cell:  (i, j-1)
                ! right cell: (i, j)

                nx = - sin(mesh%theta_face(j))
                ny =   cos(mesh%theta_face(j))

                select case (mesh%params%reconstruction)

                case (RECON%FIRST_ORDER)
                    ! First-order in theta:
                    ! left cell  = (i, j-1)
                    ! right cell = (i, j)
                    rhoL = prim%density(   i, j-1)
                    urL   = u_r_cell(i, j-1)
                    utL   = u_t_cell(i, j-1)
                    pL   = prim%pressure(  i, j-1)

                    rhoR = prim%density(   i, j)
                    urR   = u_r_cell(i, j)
                    utR   = u_t_cell(i, j)
                    pR   = prim%pressure(  i, j)

                case (RECON%MUSCL_MINMOD, RECON%MUSCL_VANLEER)
                    if (mesh%params%guard_cells < 2) then
                        write(*,*) "ERROR: MUSCL needs at least 2 ghost cells."
                        stop
                    end if

                    ! Stencil in theta: j-2, j-1, j, j+1

                    qst(1:4) = prim%density(i, j-2:j+1)
                    call muscl_reconstruct(qst, rhoL, rhoR, mesh%params%reconstruction)
                    qst(1:4) = u_r_cell(i, j-2:j+1) 
                    call muscl_reconstruct(qst, urL, urR, mesh%params%reconstruction)
                    qst(1:4) = u_t_cell(i, j-2:j+1)
                    call muscl_reconstruct(qst, utL, utR, mesh%params%reconstruction)
                    qst(1:4) = prim%pressure(i, j-2:j+1)
                    call muscl_reconstruct(qst, pL, pR, mesh%params%reconstruction)

                end select

                ! convert back to (u,v)
                uL = urL * cos(mesh%theta_face(j)) - utL * sin(mesh%theta_face(j))
                vL = urL * sin(mesh%theta_face(j)) + utL * cos(mesh%theta_face(j))
                uR = urR * cos(mesh%theta_face(j)) - utR * sin(mesh%theta_face(j))
                vR = urR * sin(mesh%theta_face(j)) + utR * cos(mesh%theta_face(j))

                call flux_solver(mesh%params%riemann_solver, rhoL,uL,vL,pL, rhoR,uR,vR,pR, &
                                 nx,ny, gas, fmass, fmom_u, fmom_v, fenergy)
                
                ! multiply by face length
                flux%angular%mass(i,j)   = fmass   * mesh%angular_face_length(i,j)
                flux%angular%mom_u(i,j)  = fmom_u  * mesh%angular_face_length(i,j)
                flux%angular%mom_v(i,j)  = fmom_v  * mesh%angular_face_length(i,j)
                flux%angular%energy(i,j) = fenergy * mesh%angular_face_length(i,j)

            end do
        end do

        deallocate(u_r_cell)
        deallocate(u_t_cell)
    end subroutine compute_angular_fluxes_rot


    subroutine compute_radial_fluxes_metric(mesh, prim, gas, flux)
        type(PolarMesh),          intent(in)    :: mesh
        type(PrimitiveVariables), intent(in)    :: prim
        type(StiffenedGas),       intent(in)    :: gas
        type(FluxFields),         intent(inout) :: flux

        integer :: i, j
        integer :: iL, iR, jL, jR
        double precision :: magnitude, dy_dtheta, dx_dtheta
        double precision :: nx, ny
        double precision :: rhoL,uL,vL,pL
        double precision :: rhoR,uR,vR,pR
        double precision :: qst(4)
        double precision :: fmass, fmom_u, fmom_v, fenergy

        iL = mesh%params%i_lo_phys
        iR = mesh%params%i_hi_phys
        jL = mesh%params%j_lo_phys
        jR = mesh%params%j_hi_phys

        do j = jL, jR
            do i = iL, iR+1
                ! left cell:  (i-1, j)
                ! right cell: (i,   j)
                dy_dtheta = (mesh%y_ll(i,j+1) - mesh%y_ll(i,j))
                dx_dtheta = (mesh%x_ll(i,j+1) - mesh%x_ll(i,j))
                magnitude = dsqrt( dx_dtheta**2.d0 + dy_dtheta**2.d0 )
                if (magnitude < 1.0d-12) then
                    nx = 0.0d0
                    ny = 0.0d0
                else
                    nx =   dy_dtheta / magnitude
                    ny =  -dx_dtheta / magnitude
                end if

                select case (mesh%params%reconstruction)

                case (RECON%FIRST_ORDER)
                    ! ------------------------------------
                    ! First-order: piecewise constant
                    ! left cell  = (i-1, j)
                    ! right cell = (i,   j)
                    ! ------------------------------------
                    rhoL = prim%density(   i-1, j)
                    uL   = prim%velocity_u(i-1, j)
                    vL   = prim%velocity_v(i-1, j)
                    pL   = prim%pressure(  i-1, j)

                    rhoR = prim%density(   i, j)
                    uR   = prim%velocity_u(i, j)
                    vR   = prim%velocity_v(i, j)
                    pR   = prim%pressure(  i, j)

                case (RECON%MUSCL_MINMOD, RECON%MUSCL_VANLEER)
                    if (mesh%params%guard_cells < 2) then
                        write(*,*) "ERROR: MUSCL needs at least 2 ghost cells."
                        stop
                    end if

                    ! ------------------------------------
                    ! Second-order MUSCL (minmod limiter)
                    ! Interface: i-1/2
                    ! Stencil for cell-centered q:
                    !   q_{i-2}, q_{i-1}, q_i, q_{i+1}
                    ! Near boundaries: fall back to first order.
                    ! ------------------------------------
                    if (i <= iL+1) then
                        ! Not enough neighbors → use first-order
                        rhoL = prim%density(   i-1, j)
                        uL   = prim%velocity_u(i-1, j)
                        vL   = prim%velocity_v(i-1, j)
                        pL   = prim%pressure(  i-1, j)

                        rhoR = prim%density(   i, j)
                        uR   = prim%velocity_u(i, j)
                        vR   = prim%velocity_v(i, j)
                        pR   = prim%pressure(  i, j)
                    else
                        qst(1:4) = prim%density(i-2:i+1, j)
                        call muscl_reconstruct(qst, rhoL, rhoR, mesh%params%reconstruction)
                        qst(1:4) = prim%velocity_u(i-2:i+1, j)
                        call muscl_reconstruct(qst, uL, uR, mesh%params%reconstruction)
                        qst(1:4) = prim%velocity_v(i-2:i+1, j)
                        call muscl_reconstruct(qst, vL, vR, mesh%params%reconstruction)
                        qst(1:4) = prim%pressure(i-2:i+1, j)
                        call muscl_reconstruct(qst, pL, pR, mesh%params%reconstruction)
                    end if

                end select

                call flux_solver(mesh%params%riemann_solver, rhoL,uL,vL,pL, rhoR,uR,vR,pR, &
                                 nx,ny, gas, fmass, fmom_u, fmom_v, fenergy)

                ! multiply by face length
                flux%radial%mass(i,j)   = fmass   * magnitude
                flux%radial%mom_u(i,j)  = fmom_u  * magnitude
                flux%radial%mom_v(i,j)  = fmom_v  * magnitude
                flux%radial%energy(i,j) = fenergy * magnitude
                
            end do
        end do

    end subroutine compute_radial_fluxes_metric

    subroutine compute_angular_fluxes_metric(mesh, prim, gas, flux)
        type(PolarMesh),          intent(in)    :: mesh
        type(PrimitiveVariables), intent(in)    :: prim
        type(StiffenedGas),       intent(in)    :: gas
        type(FluxFields),         intent(inout) :: flux

        integer :: i, j
        integer :: iL, iR, jL, jR
        double precision :: magnitude, dy_dr, dx_dr
        double precision :: nx, ny
        double precision :: rhoL,uL,vL,pL
        double precision :: rhoR,uR,vR,pR
        double precision :: ur_R, ur_L, utheta_R, utheta_L
        double precision :: qst(4)
        double precision :: fmass, fmom_u, fmom_v, fenergy

        iL = mesh%params%i_lo_phys
        iR = mesh%params%i_hi_phys
        jL = mesh%params%j_lo_phys
        jR = mesh%params%j_hi_phys

        do j = jL, jR+1
            do i = iL, iR
                ! left cell:  (i, j-1)
                ! right cell: (i, j)

                dy_dr = (mesh%y_ll(i+1,j) - mesh%y_ll(i,j))
                dx_dr = (mesh%x_ll(i+1,j) - mesh%x_ll(i,j))
                magnitude = dsqrt( dy_dr**2.d0 + dx_dr**2.d0 )
                nx =  -dy_dr / magnitude
                ny =   dx_dr / magnitude

                select case (mesh%params%reconstruction)

                case (RECON%FIRST_ORDER)
                    ! First-order in theta:
                    ! left cell  = (i, j-1)
                    ! right cell = (i, j)
                    rhoL = prim%density(   i, j-1)
                    uL   = prim%velocity_u(i, j-1)
                    vL   = prim%velocity_v(i, j-1)
                    pL   = prim%pressure(  i, j-1)

                    rhoR = prim%density(   i, j)
                    uR   = prim%velocity_u(i, j)
                    vR   = prim%velocity_v(i, j)
                    pR   = prim%pressure(  i, j)

                case (RECON%MUSCL_MINMOD, RECON%MUSCL_VANLEER)
                    if (mesh%params%guard_cells < 2) then
                        write(*,*) "ERROR: MUSCL needs at least 2 ghost cells."
                        stop
                    end if

                    ! Stencil in theta: j-2, j-1, j, j+1

                    qst(1:4) = prim%density(i, j-2:j+1)
                    call muscl_reconstruct(qst, rhoL, rhoR, mesh%params%reconstruction)
                    qst(1:4) = prim%velocity_u(i, j-2:j+1) 
                    call muscl_reconstruct(qst, uL, uR, mesh%params%reconstruction)
                    qst(1:4) = prim%velocity_v(i, j-2:j+1)
                    call muscl_reconstruct(qst, vL, vR, mesh%params%reconstruction)
                    qst(1:4) = prim%pressure(i, j-2:j+1)
                    call muscl_reconstruct(qst, pL, pR, mesh%params%reconstruction)

                end select

                call flux_solver(mesh%params%riemann_solver, rhoL,uL,vL,pL, rhoR,uR,vR,pR, &
                                 nx,ny, gas, fmass, fmom_u, fmom_v, fenergy)
                
                ! multiply by face length
                flux%angular%mass(i,j)   = fmass   * magnitude
                flux%angular%mom_u(i,j)  = fmom_u  * magnitude
                flux%angular%mom_v(i,j)  = fmom_v  * magnitude
                flux%angular%energy(i,j) = fenergy * magnitude

            end do
        end do
    end subroutine compute_angular_fluxes_metric

    
    !===============================================================
    !  Compute radial numerical fluxes for Polar Euler equations (explicit geometric sources)
    !  Returns *physical* fluxes F_r and F_θ.
    !  Does not include geometric source terms; they are added separately.
    !===============================================================
    subroutine compute_radial_fluxes_polar(mesh, prim, gas, flux)
        type(PolarMesh),          intent(in)    :: mesh
        type(PrimitiveVariables), intent(in)    :: prim
        type(StiffenedGas),       intent(in)    :: gas
        type(FluxFields),         intent(inout) :: flux

        integer :: i, j
        integer :: iL, iR, jL, jR
        double precision :: nx, ny
        double precision :: rhoL,uL,vL,pL
        double precision :: rhoR,uR,vR,pR
        double precision :: qst(4)
        double precision :: fmass, fmom_u, fmom_v, fenergy

        iL = mesh%params%i_lo_phys
        iR = mesh%params%i_hi_phys
        jL = mesh%params%j_lo_phys
        jR = mesh%params%j_hi_phys

        do j = jL, jR
            do i = iL, iR+1
                nx = 1.0d0
                ny = 0.0d0

                select case (mesh%params%reconstruction)

                case (RECON%FIRST_ORDER)
                    ! ------------------------------------
                    ! First-order: piecewise constant
                    ! left cell  = (i-1, j)
                    ! right cell = (i,   j)
                    ! ------------------------------------
                    rhoL = prim%density(   i-1, j)
                    uL   = prim%velocity_u(i-1, j)
                    vL   = prim%velocity_v(i-1, j)
                    pL   = prim%pressure(  i-1, j)

                    rhoR = prim%density(   i, j)
                    uR   = prim%velocity_u(i, j)
                    vR   = prim%velocity_v(i, j)
                    pR   = prim%pressure(  i, j)

                case (RECON%MUSCL_MINMOD, RECON%MUSCL_VANLEER)
                    if (mesh%params%guard_cells < 2) then
                        write(*,*) "ERROR: MUSCL needs at least 2 ghost cells."
                        stop
                    end if

                    ! ------------------------------------
                    ! Second-order MUSCL (minmod limiter)
                    ! Interface: i-1/2
                    ! Stencil for cell-centered q:
                    !   q_{i-2}, q_{i-1}, q_i, q_{i+1}
                    ! Near boundaries: fall back to first order.
                    ! ------------------------------------
                    if (i <= iL+1) then
                        ! Not enough neighbors → use first-order
                        rhoL = prim%density(   i-1, j)
                        uL   = prim%velocity_u(i-1, j)
                        vL   = prim%velocity_v(i-1, j)
                        pL   = prim%pressure(  i-1, j)

                        rhoR = prim%density(   i, j)
                        uR   = prim%velocity_u(i, j)
                        vR   = prim%velocity_v(i, j)
                        pR   = prim%pressure(  i, j)
                    else
                        qst(1:4) = prim%density(i-2:i+1, j)
                        call muscl_reconstruct(qst, rhoL, rhoR, mesh%params%reconstruction)
                        qst(1:4) = prim%velocity_u(i-2:i+1, j)
                        call muscl_reconstruct(qst, uL, uR, mesh%params%reconstruction)
                        qst(1:4) = prim%velocity_v(i-2:i+1, j)
                        call muscl_reconstruct(qst, vL, vR, mesh%params%reconstruction)
                        qst(1:4) = prim%pressure(i-2:i+1, j)
                        call muscl_reconstruct(qst, pL, pR, mesh%params%reconstruction)
                    end if

                end select

                call flux_solver(mesh%params%riemann_solver, rhoL,uL,vL,pL, rhoR,uR,vR,pR, &
                                 nx,ny, gas, fmass, fmom_u, fmom_v, fenergy)
                ! multiply by face length
                flux%radial%mass(i,j)   = fmass   * mesh%radial_face_length(i,j)
                flux%radial%mom_u(i,j)  = fmom_u  * mesh%radial_face_length(i,j)
                flux%radial%mom_v(i,j)  = fmom_v  * mesh%radial_face_length(i,j)
                flux%radial%energy(i,j) = fenergy * mesh%radial_face_length(i,j)
            
            end do
        end do
    end subroutine compute_radial_fluxes_polar

    subroutine compute_angular_fluxes_polar(mesh, prim, gas, flux)
        type(PolarMesh),          intent(in)    :: mesh
        type(PrimitiveVariables), intent(in)    :: prim
        type(StiffenedGas),       intent(in)    :: gas
        type(FluxFields),         intent(inout) :: flux

        integer :: i, j
        integer :: iL, iR, jL, jR
        double precision :: nx, ny
        double precision :: rhoL,uL,vL,pL
        double precision :: rhoR,uR,vR,pR
        double precision :: qst(4)
        double precision :: fmass, fmom_u, fmom_v, fenergy

        iL = mesh%params%i_lo_phys
        iR = mesh%params%i_hi_phys
        jL = mesh%params%j_lo_phys
        jR = mesh%params%j_hi_phys

        do j = jL, jR+1
            do i = iL, iR
                nx = 0.0d0
                ny = 1.0d0

                select case (mesh%params%reconstruction)

                case (RECON%FIRST_ORDER)
                    ! First-order in theta:
                    ! left cell  = (i, j-1)
                    ! right cell = (i, j)
                    rhoL = prim%density(   i, j-1)
                    uL   = prim%velocity_u(i, j-1)
                    vL   = prim%velocity_v(i, j-1)
                    pL   = prim%pressure(  i, j-1)

                    rhoR = prim%density(   i, j)
                    uR   = prim%velocity_u(i, j)
                    vR   = prim%velocity_v(i, j)
                    pR   = prim%pressure(  i, j)

                case (RECON%MUSCL_MINMOD, RECON%MUSCL_VANLEER)
                    if (mesh%params%guard_cells < 2) then
                        write(*,*) "ERROR: MUSCL needs at least 2 ghost cells."
                        stop
                    end if

                    ! Stencil in theta: j-2, j-1, j, j+1
                    qst(1:4) = prim%density(i, j-2:j+1)
                    call muscl_reconstruct(qst, rhoL, rhoR, mesh%params%reconstruction)
                    qst(1:4) = prim%velocity_u(i, j-2:j+1)
                    call muscl_reconstruct(qst, uL, uR, mesh%params%reconstruction)
                    qst(1:4) = prim%velocity_v(i, j-2:j+1)
                    call muscl_reconstruct(qst, vL, vR, mesh%params%reconstruction)
                    qst(1:4) = prim%pressure(i, j-2:j+1)
                    call muscl_reconstruct(qst, pL, pR, mesh%params%reconstruction)

                end select

                call flux_solver(mesh%params%riemann_solver, rhoL,uL,vL,pL, rhoR,uR,vR,pR, &
                                 nx,ny, gas, fmass, fmom_u, fmom_v, fenergy)
                
                ! multiply by face length
                flux%angular%mass(i,j)   = fmass   * mesh%angular_face_length(i,j)
                flux%angular%mom_u(i,j)  = fmom_u  * mesh%angular_face_length(i,j)
                flux%angular%mom_v(i,j)  = fmom_v  * mesh%angular_face_length(i,j)
                flux%angular%energy(i,j) = fenergy * mesh%angular_face_length(i,j)
            end do
        end do
    end subroutine compute_angular_fluxes_polar

end module flux_ops
