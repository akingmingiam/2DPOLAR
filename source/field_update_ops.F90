module field_update_ops
    use mesh_types,       only : PolarMesh, EQ
    use flow_fields,      only : PrimitiveVariables, ConservedVariables
    use fluid_properties, only : StiffenedGas
    use flux_types,       only : FluxFields
    implicit none
    private

    public :: update_conserved_field

contains

subroutine update_conserved_field(mesh, cons, prim, flux, dt)
    type(PolarMesh),          intent(in)    :: mesh
    type(ConservedVariables), intent(inout) :: cons
    type(PrimitiveVariables), intent(in)    :: prim
    type(FluxFields),         intent(in)    :: flux
    double precision,         intent(in)    :: dt

    select case(mesh%params%equation_type)

    case (EQ%CARTESIAN_TRANSFORM)
        call update_conserved_metric(mesh, cons, flux, dt)

    case (EQ%POLAR_EULER)
        call update_conserved_polar(mesh, cons, prim, flux, dt)

    case default
        print *, "ERROR: Unknown equation type in compute_all_fluxes"
        stop

    end select
    
end subroutine update_conserved_field

!===================================================================
!    Update conserved variables using the *metric-form* (transformed)
!    conservation laws:
!
!    ∂(J U)/∂t + ∂(F̂)/∂ξ + ∂(Ĝ)/∂η = 0

!    U_new = U_old - dt / V * (F_out - F_in)
!
!    This corresponds to *Cartesian Euler equations* mapped onto a
!    general curvilinear grid.  All geometric effects are absorbed into
!    the metric coefficients (J, ξ_x, ξ_y, η_x, η_y) and appear inside
!    the transformed fluxes.  **No geometric source terms** appear.
!
!    flux arrays already contain:
!    Flux * Face_Length   
!===================================================================
subroutine update_conserved_metric(mesh, cons, flux, dt)
    type(PolarMesh),          intent(in)    :: mesh
    type(ConservedVariables), intent(inout) :: cons
    type(FluxFields),         intent(in)    :: flux
    double precision,         intent(in)    :: dt

    integer :: i, j
    integer :: iL, iR, jL, jR
    double precision :: Vij

    ! physical domain indices
    iL = mesh%params%i_lo_phys
    iR = mesh%params%i_hi_phys
    jL = mesh%params%j_lo_phys
    jR = mesh%params%j_hi_phys

    do j = jL, jR
        do i = iL, iR

            Vij = mesh%cell_area(i,j)

            !------------------------------
            ! mass conservation
            !------------------------------
            cons%rho(i,j) = cons%rho(i,j) - dt / Vij * &
               ( flux%radial%mass(i+1,j)  - flux%radial%mass(i,j)  &
               + flux%angular%mass(i,j+1) - flux%angular%mass(i,j) )

            !------------------------------
            ! momentum-x conservation
            !------------------------------
            cons%rho_u(i,j) = cons%rho_u(i,j) - dt / Vij * &
               ( flux%radial%mom_u(i+1,j)  - flux%radial%mom_u(i,j)  &
               + flux%angular%mom_u(i,j+1) - flux%angular%mom_u(i,j) )

            !------------------------------
            ! momentum-y conservation
            !------------------------------
            cons%rho_v(i,j) = cons%rho_v(i,j) - dt / Vij * &
               ( flux%radial%mom_v(i+1,j)  - flux%radial%mom_v(i,j)  &
               + flux%angular%mom_v(i,j+1) - flux%angular%mom_v(i,j) )

            !------------------------------
            ! energy conservation
            !------------------------------
            cons%rho_E(i,j) = cons%rho_E(i,j) - dt / Vij * &
               ( flux%radial%energy(i+1,j)  - flux%radial%energy(i,j)  &
               + flux%angular%energy(i,j+1) - flux%angular%energy(i,j) )
            
        end do
    end do
    
end subroutine update_conserved_metric


!===================================================================
!    Update conserved variables using the *physical Euler equations*
!    written in *true polar coordinates*:
!
!    ∂(ρ u_r)/∂t + ... = (ρ u_θ^2 + p) / r
!
!    ∂(ρ u_θ)/∂t + ... = - ρ u_r u_θ / r
!    
!    U_new = U_old - dt / V * (F_out - F_in) + dt * S(U,r)

!    The terms (+p/r, +ρ u_θ² / r, −ρ u_r u_θ / r) are *geometric
!    source terms* and must be added explicitly.

!    Note:
!    We do NOT rename variables here:
!    cons%rho_u <=> ρ u_r
!    cons%rho_v <=> ρ u_θ
!    prim%velocity_u <=> u_r
!    prim%velocity_v <=> u_θ
!===================================================================
subroutine update_conserved_polar(mesh, cons, prim, flux, dt)
    type(PolarMesh),          intent(in)    :: mesh
    type(ConservedVariables), intent(inout) :: cons
    type(PrimitiveVariables), intent(in)    :: prim
    type(FluxFields),         intent(in)    :: flux
    double precision,         intent(in)    :: dt

    integer :: i, j
    integer :: iL, iR, jL, jR
    double precision :: Vij, r
    double precision :: rho, ur, uth, p
    double precision :: S_ur, S_uth

    ! physical domain indices
    iL = mesh%params%i_lo_phys
    iR = mesh%params%i_hi_phys
    jL = mesh%params%j_lo_phys
    jR = mesh%params%j_hi_phys

    do j = jL, jR
        do i = iL, iR

            Vij = mesh%cell_area(i,j)
            r   = mesh%r_center(i)
            rho = prim%density   (i,j)
            ur  = prim%velocity_u(i,j)   
            uth = prim%velocity_v(i,j)  
            p   = prim%pressure  (i,j)

            !--------------------------------------------------
            !  Polar-coordinate geometric source terms:
            !  This implementation treats the geometric source term S(U,r) as if it
            !  were a flux through a fictitious face:
            !  Δ(ρ u_r)  ←  - dt/V * [ ...  - S_r * Δr * Δθ ]
            !  S_r  = (ρ u_θ^2 + p)
            !  S_θ  = -ρ u_r u_θ 
            !  That is,  S_r  is multiplied by the local face-area Δr Δθ first,
            !  then the whole term is divided by the cell volume V = r Δr Δθ.
            !--------------------------------------------------
            if (r > 0.0d0) then
                S_ur  = rho * uth*uth + p
                S_uth = - rho * ur * uth
            else
                ! Special treatment at r = 0 (e.g., force S = 0 or use r_min)
                S_ur  = 0.0d0
                S_uth = 0.0d0
            end if

            cons%rho(i,j) = cons%rho(i,j) - dt / Vij * &
               ( flux%radial%mass(i+1,j)  - flux%radial%mass(i,j)  &
               + flux%angular%mass(i,j+1) - flux%angular%mass(i,j) )

            cons%rho_u(i,j) = cons%rho_u(i,j) - dt / Vij * &
               ( flux%radial%mom_u(i+1,j)  - flux%radial%mom_u(i,j)  &
               + flux%angular%mom_u(i,j+1) - flux%angular%mom_u(i,j) &
               - S_ur * mesh%radial_spacing(i) * mesh%angular_spacing(j))

            cons%rho_v(i,j) = cons%rho_v(i,j) - dt / Vij * &
               ( flux%radial%mom_v(i+1,j)  - flux%radial%mom_v(i,j)  &
               + flux%angular%mom_v(i,j+1) - flux%angular%mom_v(i,j) &
               - S_uth * mesh%radial_spacing(i) * mesh%angular_spacing(j))

            cons%rho_E(i,j) = cons%rho_E(i,j) - dt / Vij * &
               ( flux%radial%energy(i+1,j)  - flux%radial%energy(i,j)  &
               + flux%angular%energy(i,j+1) - flux%angular%energy(i,j) )

        end do
    end do
    
end subroutine update_conserved_polar


! The source term and flux are calculated separately
subroutine update_conserved_polar0(mesh, cons, prim, flux, dt)
    type(PolarMesh),          intent(in)    :: mesh
    type(ConservedVariables), intent(inout) :: cons
    type(PrimitiveVariables), intent(in)    :: prim
    type(FluxFields),         intent(in)    :: flux
    double precision,         intent(in)    :: dt

    integer :: i, j
    integer :: iL, iR, jL, jR
    double precision :: Vij, r
    double precision :: rho, ur, uth, p
    double precision :: S_ur, S_uth

    ! physical domain indices
    iL = mesh%params%i_lo_phys
    iR = mesh%params%i_hi_phys
    jL = mesh%params%j_lo_phys
    jR = mesh%params%j_hi_phys

    do j = jL, jR
        do i = iL, iR

            Vij = mesh%cell_area(i,j)
            r   = mesh%r_center(i)
            rho = prim%density   (i,j)
            ur  = prim%velocity_u(i,j)  
            uth = prim%velocity_v(i,j)  
            p   = prim%pressure  (i,j)

            !--------------------------------------------------
            !  Polar-coordinate geometric source terms:
            !  This version follows the physical Euler equations in polar form:
            !
            !   ∂(ρ u_r)/∂t = ...  + (ρ u_θ² + p)/r
            !
            !  and applies the source term exactly as it appears:
            !
            !        Δ(ρ u_r) = dt * S_ur
            !
            !  where S_ur = (ρ u_θ² + p)/r  is already the contribution per unit volume.
            !--------------------------------------------------
            if (r > 0.0d0) then
                S_ur  = (rho * uth*uth + p) / r
                S_uth = (- rho * ur * uth)  / r
            else
                ! Special treatment at r = 0 (e.g., force S = 0 or use r_min)
                S_ur  = 0.0d0
                S_uth = 0.0d0
            end if

            ! U_new = U_old - dt / V * (F_out - F_in) + dt * S(U,r)
            cons%rho(i,j) = cons%rho(i,j) - dt / Vij * &
               ( flux%radial%mass(i+1,j)  - flux%radial%mass(i,j)  &
               + flux%angular%mass(i,j+1) - flux%angular%mass(i,j) )

            cons%rho_u(i,j) = cons%rho_u(i,j) - dt / Vij * &
               ( flux%radial%mom_u(i+1,j)  - flux%radial%mom_u(i,j)  &
               + flux%angular%mom_u(i,j+1) - flux%angular%mom_u(i,j) ) &
               + dt * S_ur

            cons%rho_v(i,j) = cons%rho_v(i,j) - dt / Vij * &
               ( flux%radial%mom_v(i+1,j)  - flux%radial%mom_v(i,j)  &
               + flux%angular%mom_v(i,j+1) - flux%angular%mom_v(i,j) ) &
               + dt * S_uth

            cons%rho_E(i,j) = cons%rho_E(i,j) - dt / Vij * &
               ( flux%radial%energy(i+1,j)  - flux%radial%energy(i,j)  &
               + flux%angular%energy(i,j+1) - flux%angular%energy(i,j) )

        end do
    end do
    
end subroutine update_conserved_polar0



end module field_update_ops
