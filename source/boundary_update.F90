module boundary_update
    !======================================================================
    !  Module: boundary_update
    !  Purpose:
    !      Apply boundary conditions to primitive and conserved fields
    !      on a polar finite-volume mesh with ghost cells.
    !      Ghost cells:
    !          Radial:
    !              i = i_lo_phys-g ... i_lo_phys-1   (inner ghost)
    !              i = i_hi_phys+1 ... i_hi_phys+g   (outer ghost)
    !
    !          Angular:
    !              j = j_lo_phys-g ... j_lo_phys-1   (lower theta ghost)
    !              j = j_hi_phys+1 ... j_hi_phys+g   (upper theta ghost)
    !======================================================================

    use mesh_types,       only : PolarMesh, BC
    use flow_fields,      only : PrimitiveVariables, ConservedVariables
    use fluid_properties, only : StiffenedGas, specific_total_energy_stiffened
    implicit none
    private

    public :: update_boundary_primitive_conserved

contains

!==================================================================
!  Subroutine: update_boundary_primitive_conserved
!  Purpose   : Top-level BC update for both primitive and conserved
!              variables. Extracts all needed info from `mesh%params`
!              and dispatches to directional BC handlers.
!==================================================================
subroutine update_boundary_primitive_conserved(mesh, prim, cons, gas)
    type(PolarMesh),          intent(in)    :: mesh
    type(PrimitiveVariables), intent(inout) :: prim
    type(ConservedVariables), intent(inout) :: cons
    type(StiffenedGas),       intent(in)    :: gas

    integer :: g
    integer :: i_lo_phys, i_hi_phys
    integer :: j_lo_phys, j_hi_phys
    integer :: bc_r_lo, bc_r_hi
    integer :: bc_t_lo, bc_t_hi

    !-------------------------------------------
    ! Extract index ranges and BC types
    !-------------------------------------------
    g          = mesh%params%guard_cells

    i_lo_phys  = mesh%params%i_lo_phys
    i_hi_phys  = mesh%params%i_hi_phys
    j_lo_phys  = mesh%params%j_lo_phys
    j_hi_phys  = mesh%params%j_hi_phys

    bc_r_lo    = mesh%params%bc_r_lo
    bc_r_hi    = mesh%params%bc_r_hi
    bc_t_lo    = mesh%params%bc_theta_lo
    bc_t_hi    = mesh%params%bc_theta_hi

    !-------------------------------------------
    ! Radial and angular BC application
    !-------------------------------------------
    call apply_bc_r(prim, cons, g, i_lo_phys, i_hi_phys, j_lo_phys, j_hi_phys, &
                    bc_r_lo, bc_r_hi, gas)

    call apply_bc_theta(prim, cons, g, i_lo_phys, i_hi_phys, j_lo_phys, j_hi_phys, &
                        bc_t_lo, bc_t_hi)

end subroutine update_boundary_primitive_conserved



!==================================================================
!  Subroutine: apply_bc_r
!  Purpose   : Apply BC on radial direction (inner/outer radius).
!
!  Inputs:
!      g           : number of ghost layers
!      i_lo_phys   : first physical radial index
!      i_hi_phys   : last  physical radial index
!      j_lo_phys   : first physical angular index
!      j_hi_phys   : last  physical angular index
!      bc_lo       : BC type on inner radius
!      bc_hi       : BC type on outer radius
!==================================================================
subroutine apply_bc_r(prim, cons, g, i_lo_phys, i_hi_phys, j_lo_phys, j_hi_phys, &
                      bc_lo, bc_hi, gas)
    type(PrimitiveVariables), intent(inout) :: prim
    type(ConservedVariables), intent(inout) :: cons
    type(StiffenedGas),       intent(in)    :: gas

    integer, intent(in) :: g
    integer, intent(in) :: i_lo_phys, i_hi_phys
    integer, intent(in) :: j_lo_phys, j_hi_phys
    integer, intent(in) :: bc_lo, bc_hi

    integer :: j, k
    integer :: ig, ii  ! ghost index, interior index

    !========================
    !  INNER radial boundary
    !========================
    select case (bc_lo)

    case (BC%EXTRAPOL)
        ! Zero-gradient: ghost cells copy nearest interior cell i_lo_phys
        do j = j_lo_phys, j_hi_phys
            do k = 1, g
                ig = i_lo_phys - k
                ii = i_lo_phys
                prim%density(   ig,j) = prim%density(   ii,j)
                prim%pressure(  ig,j) = prim%pressure(  ii,j)
                prim%velocity_u(ig,j) = prim%velocity_u(ii,j)
                prim%velocity_v(ig,j) = prim%velocity_v(ii,j)

                cons%rho(  ig,j) = cons%rho(  ii,j)
                cons%rho_u(ig,j) = cons%rho_u(ii,j)
                cons%rho_v(ig,j) = cons%rho_v(ii,j)
                cons%rho_E(ig,j) = cons%rho_E(ii,j)
            end do
        end do

    case (BC%PERIODIC)
        ! Periodic in r (less common, but kept for completeness)
        do j = j_lo_phys, j_hi_phys
            do k = 1, g
                ig = i_lo_phys - k
                ii = i_hi_phys - (k-1)
                prim%density(   ig,j) = prim%density(   ii,j)
                prim%pressure(  ig,j) = prim%pressure(  ii,j)
                prim%velocity_u(ig,j) = prim%velocity_u(ii,j)
                prim%velocity_v(ig,j) = prim%velocity_v(ii,j)

                cons%rho(  ig,j) = cons%rho(  ii,j)
                cons%rho_u(ig,j) = cons%rho_u(ii,j)
                cons%rho_v(ig,j) = cons%rho_v(ii,j)
                cons%rho_E(ig,j) = cons%rho_E(ii,j)
            end do
        end do

    case (BC%AXIS)
         do j = j_lo_phys, j_hi_phys
            do k = 1, g
                ig = i_lo_phys - k
                ii = i_lo_phys
                prim%density(   ig,j) = prim%density(   ii,j)
                prim%pressure(  ig,j) = prim%pressure(  ii,j)
                prim%velocity_u(ig,j) = -prim%velocity_u(ii,j)
                prim%velocity_v(ig,j) = prim%velocity_v(ii,j)

                cons%rho(  ig,j) = cons%rho(  ii,j)
                cons%rho_u(ig,j) = -cons%rho_u(ii,j)
                cons%rho_v(ig,j) = cons%rho_v(ii,j)
                cons%rho_E(ig,j) = cons%rho_E(ii,j)
            end do
        end do

    end select


    !========================
    !  OUTER radial boundary
    !========================
    select case (bc_hi)

    case (BC%EXTRAPOL)
        do j = j_lo_phys, j_hi_phys
            do k = 1, g
                ig = i_hi_phys + k
                ii = i_hi_phys
                prim%density(   ig,j) = prim%density(   ii,j)
                prim%pressure(  ig,j) = prim%pressure(  ii,j)
                prim%velocity_u(ig,j) = prim%velocity_u(ii,j)
                prim%velocity_v(ig,j) = prim%velocity_v(ii,j)

                cons%rho(  ig,j) = cons%rho(  ii,j)
                cons%rho_u(ig,j) = cons%rho_u(ii,j)
                cons%rho_v(ig,j) = cons%rho_v(ii,j)
                cons%rho_E(ig,j) = cons%rho_E(ii,j)
            end do
        end do

    case (BC%PERIODIC)
        do j = j_lo_phys, j_hi_phys
            do k = 1, g
                ig = i_hi_phys + k
                ii = i_lo_phys + (k-1)
                prim%density(   ig,j) = prim%density(   ii,j)
                prim%pressure(  ig,j) = prim%pressure(  ii,j)
                prim%velocity_u(ig,j) = prim%velocity_u(ii,j)
                prim%velocity_v(ig,j) = prim%velocity_v(ii,j)

                cons%rho(  ig,j) = cons%rho(  ii,j)
                cons%rho_u(ig,j) = cons%rho_u(ii,j)
                cons%rho_v(ig,j) = cons%rho_v(ii,j)
                cons%rho_E(ig,j) = cons%rho_E(ii,j)
            end do
        end do

    end select

end subroutine apply_bc_r



!==================================================================
!  Subroutine: apply_bc_theta
!  Purpose   : Apply BC on angular direction (θ).
!==================================================================
subroutine apply_bc_theta(prim, cons, g, i_lo_phys, i_hi_phys, j_lo_phys, j_hi_phys, &
                          bc_lo, bc_hi)
    type(PrimitiveVariables), intent(inout) :: prim
    type(ConservedVariables), intent(inout) :: cons
    integer, intent(in) :: g
    integer, intent(in) :: i_lo_phys, i_hi_phys
    integer, intent(in) :: j_lo_phys, j_hi_phys
    integer, intent(in) :: bc_lo, bc_hi

    integer :: i, k
    integer :: jg, jj  ! ghost index, interior index

    !========================
    !  LOWER theta boundary
    !========================
    select case (bc_lo)

    case (BC%EXTRAPOL)
        do k = 1, g
            do i = i_lo_phys, i_hi_phys
                jg = j_lo_phys - k
                jj = j_lo_phys
                prim%density(   i,jg) = prim%density(   i,jj)
                prim%pressure(  i,jg) = prim%pressure(  i,jj)
                prim%velocity_u(i,jg) = prim%velocity_u(i,jj)
                prim%velocity_v(i,jg) = prim%velocity_v(i,jj)

                cons%rho(  i,jg) = cons%rho(  i,jj)
                cons%rho_u(i,jg) = cons%rho_u(i,jj)
                cons%rho_v(i,jg) = cons%rho_v(i,jj)
                cons%rho_E(i,jg) = cons%rho_E(i,jj)
            end do
        end do

    case (BC%PERIODIC)
        do k = 1, g
            do i = i_lo_phys, i_hi_phys
                jg = j_lo_phys - k
                jj = j_hi_phys - (k-1)
                prim%density(   i,jg) = prim%density(   i,jj)
                prim%pressure(  i,jg) = prim%pressure(  i,jj)
                prim%velocity_u(i,jg) = prim%velocity_u(i,jj)
                prim%velocity_v(i,jg) = prim%velocity_v(i,jj)

                cons%rho(  i,jg) = cons%rho(  i,jj)
                cons%rho_u(i,jg) = cons%rho_u(i,jj)
                cons%rho_v(i,jg) = cons%rho_v(i,jj)
                cons%rho_E(i,jg) = cons%rho_E(i,jj)
            end do
        end do

    end select


    !========================
    !  UPPER theta boundary
    !========================
    select case (bc_hi)

    case (BC%EXTRAPOL)
        do k = 1, g
            do i = i_lo_phys, i_hi_phys
                jg = j_hi_phys + k
                jj = j_hi_phys
                prim%density(   i,jg) = prim%density(   i,jj)
                prim%pressure(  i,jg) = prim%pressure(  i,jj)
                prim%velocity_u(i,jg) = prim%velocity_u(i,jj)
                prim%velocity_v(i,jg) = prim%velocity_v(i,jj)

                cons%rho(  i,jg) = cons%rho(  i,jj)
                cons%rho_u(i,jg) = cons%rho_u(i,jj)
                cons%rho_v(i,jg) = cons%rho_v(i,jj)
                cons%rho_E(i,jg) = cons%rho_E(i,jj)
            end do
        end do

    case (BC%PERIODIC)
        do k = 1, g
            do i = i_lo_phys, i_hi_phys
                jg = j_hi_phys + k
                jj = j_lo_phys + (k-1)
                prim%density(   i,jg) = prim%density(   i,jj)
                prim%pressure(  i,jg) = prim%pressure(  i,jj)
                prim%velocity_u(i,jg) = prim%velocity_u(i,jj)
                prim%velocity_v(i,jg) = prim%velocity_v(i,jj)

                cons%rho(  i,jg) = cons%rho(  i,jj)
                cons%rho_u(i,jg) = cons%rho_u(i,jj)
                cons%rho_v(i,jg) = cons%rho_v(i,jj)
                cons%rho_E(i,jg) = cons%rho_E(i,jj)
            end do
        end do

    end select

end subroutine apply_bc_theta


end module boundary_update
