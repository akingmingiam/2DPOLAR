!===============================================================
! File   : field_update_ops.F90
! Module : field_update_ops
! Purpose:
!   Update conserved variables using polar FVM fluxes (already
!   multiplied by face length) and then reconstruct primitive
!   variables using stiffened-gas EOS.
!===============================================================
module field_update_ops
    use mesh_types,       only : PolarMesh
    use flow_fields,      only : PrimitiveVariables, ConservedVariables
    use fluid_properties, only : StiffenedGas
    use flux_types,       only : FluxFields
    implicit none
    private

    public :: update_conserved_fields

contains

!===================================================================
! 1) Update CONSERVED variables:
!    U_new = U_old - dt / V * (F_out - F_in)
!
!    flux arrays already contain:
!        Flux * Face_Length   (几何加权通量)
!===================================================================
subroutine update_conserved_fields(mesh, cons, flux, dt)
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
    
end subroutine update_conserved_fields

end module field_update_ops
