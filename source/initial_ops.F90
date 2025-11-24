module initial_ops
    use mesh_types,       only : PolarMesh
    use flow_fields,      only : PrimitiveVariables, ConservedVariables
    use fluid_properties, only : StiffenedGas, gas_air
    use field_conversion, only : update_conserved_from_primitive
    implicit none
    private

    public :: initialize_flow_fields

contains

!=====================================================================
!  Subroutine: initial_soln
!  Purpose   : Initialize primitive and conserved flow variables
!              using Stiffened-Gas EOS (here: air).
!
!  Default initial condition:
!     - rho = 1.0
!     - p   = 1.0e5 Pa
!     - u,v = 0
!=====================================================================
subroutine initialize_flow_fields(mesh, prim, cons, gas)
    use boundary_update, only : update_boundary_primitive_conserved
    type(PolarMesh),          intent(in)  :: mesh
    type(PrimitiveVariables), intent(inout) :: prim
    type(ConservedVariables), intent(inout) :: cons
    type(StiffenedGas),       intent(out)   :: gas

    integer :: i, j
    integer :: i_lo_phys, i_hi_phys, j_lo_phys, j_hi_phys
    !-------------------------------------------------------------
    ! Use AIR as the working fluid
    !-------------------------------------------------------------
    gas = gas_air()          ! gamma = 1.4, p_inf = 0
    gas%pressure_floor_ratio = 0.5d0
    gas%pressure_floor_eps   = 0.d0

    !-------------------------------------------------------------
    ! Loop over all cells (except guard cells)
    !-------------------------------------------------------------
    i_lo_phys = mesh%params%i_lo_phys
    i_hi_phys = mesh%params%i_hi_phys
    j_lo_phys = mesh%params%j_lo_phys
    j_hi_phys = mesh%params%j_hi_phys

    do j = j_lo_phys, j_hi_phys
        do i = i_lo_phys, i_hi_phys 

            !---------------------------
            ! Primitive variables
            !---------------------------
            prim%density(i,j)    = 1.0d0
            prim%pressure(i,j)   = 1.0d5
            prim%velocity_u(i,j) = 0.0d0
            prim%velocity_v(i,j) = 0.0d0
        end do
    end do

    call update_conserved_from_primitive(prim, cons, gas, i_lo_phys, i_hi_phys, j_lo_phys, j_hi_phys)
    call update_boundary_primitive_conserved(mesh, prim, cons)

    print *, "Initial solution successfully applied."

end subroutine initialize_flow_fields

end module initial_ops
