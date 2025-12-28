module initial_ops
    use mesh_types,       only : PolarMesh, pi
    use flow_fields,      only : PrimitiveVariables, ConservedVariables
    use fluid_properties, only : StiffenedGas, gas_air
    use field_conversion, only : update_conserved_from_primitive
    implicit none
    private

    public :: initialize_flow_fields
    public :: initialize_sedov

contains

!=====================================================================
!  Subroutine: initial_soln
!  Purpose   : Initialize primitive and conserved flow variables
!              using Stiffened-Gas EOS (here: air).
!=====================================================================
subroutine initialize_flow_fields(mesh, prim, cons, gas)
    use boundary_update, only : update_boundary_primitive_conserved
    type(PolarMesh),          intent(in)    :: mesh
    type(PrimitiveVariables), intent(inout) :: prim
    type(ConservedVariables), intent(inout) :: cons
    type(StiffenedGas),       intent(out)   :: gas

    integer :: i, j
    integer :: i_lo_phys, i_hi_phys, j_lo_phys, j_hi_phys

    double precision :: Ms,G,paf,daf,uaf,pbe,dbe,ube,R,theta,x


    pbe = 1.01325d+5
    dbe = 1.225d0
    ube = 0.d0
    Ms  = 2.670d0
    ! Ms  = 1.5d0
    G   = 1.4d0
    paf = pbe*( 1.d0 + 2.d0*G/( G + 1.d0 )*( Ms*Ms - 1.d0 ) )
    daf = dbe*( Ms*Ms/ ( 1.d0 + (G - 1.d0)/(G+1.d0)*( Ms*Ms - 1.d0 ) ) )
    uaf = 2.d0/(G+1.d0)*( Ms - 1.d0/Ms )* dsqrt(G* pbe/dbe)


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
            ! CARTESIAN_TRANSFORM
            !---------------------------

            ! static
            ! prim%density(i,j)    = 1.0d0
            ! prim%pressure(i,j)   = 1.0d5
            ! prim%velocity_u(i,j) = 0.0d0
            ! prim%velocity_v(i,j) = 0.0d0

            ! Circular shock wave
            R = mesh%r_center(i)
            if(R.le.1.0d0)then
                prim%density(i,j) = (dbe + (daf-dbe)*R/1.d0)
                prim%pressure(i,j) = (pbe + (paf-pbe)*R/1.d0)
                prim%velocity_u(i,j) = ((uaf*R/1.d0) * cos(mesh%theta_center(j)))
                prim%velocity_v(i,j) = ((uaf*R/1.d0) * sin(mesh%theta_center(j)))
            else
                prim%density(i,j) = dbe
                prim%pressure(i,j) = pbe
                prim%velocity_u(i,j) = 0.d0
                prim%velocity_v(i,j) = 0.d0
            endif
            
        end do
    end do

    call update_conserved_from_primitive(prim, cons, gas, i_lo_phys, i_hi_phys, j_lo_phys, j_hi_phys)
    call update_boundary_primitive_conserved(mesh, prim, cons, gas)

    print *, "Initial solution successfully applied."

end subroutine initialize_flow_fields


subroutine initialize_sedov(mesh, prim, cons, gas)
    use boundary_update, only : update_boundary_primitive_conserved
    type(PolarMesh),          intent(in)  :: mesh
    type(PrimitiveVariables), intent(inout) :: prim
    type(ConservedVariables), intent(inout) :: cons
    type(StiffenedGas),       intent(out)   :: gas

    integer :: i, j
    integer :: i_lo_phys, i_hi_phys, j_lo_phys, j_hi_phys
    integer :: i_inner
    double precision :: rho0, p0, ur0, uth0
    double precision :: gamma, E0, V_ring, p_inner, E_total

    !-------------------------------------------------------------
    ! Use AIR as the working fluid
    !-------------------------------------------------------------
    gas = gas_air()          ! gamma = 1.4, p_inf = 0
    gas%pressure_floor_ratio = 0.5d0
    gas%pressure_floor_eps   = 0.d0
    gamma = gas%gamma

    !-------------------------------------------------------------
    ! Loop over all cells (except guard cells)
    !-------------------------------------------------------------
    i_lo_phys = mesh%params%i_lo_phys
    i_hi_phys = mesh%params%i_hi_phys
    j_lo_phys = mesh%params%j_lo_phys
    j_hi_phys = mesh%params%j_hi_phys

    !-------------------------------------------------------------
    ! Background state: (rho0, p0, ur0, uth0)
    !   Sedov paper:  rho0 = 1, p0 very small, u = 0
    !-------------------------------------------------------------
    rho0 = 1.0d0
    p0   = 1.0d-6      ! small ambient pressure
    ur0  = 0.0d0
    uth0 = 0.0d0

    !-------------------------------------------------------------
    ! Step 1: fill background state in all physical cells
    !-------------------------------------------------------------
    do j = j_lo_phys, j_hi_phys
        do i = i_lo_phys, i_hi_phys
            prim%density(i,j)    = rho0
            prim%pressure(i,j)   = p0
            prim%velocity_u(i,j) = ur0    ! radial velocity
            prim%velocity_v(i,j) = uth0   ! azimuthal velocity
        end do
    end do

    !-------------------------------------------------------------
    ! Step 2: deposit Sedov point-blast energy in innermost ring
    !
    !   Choose E0 so that shock radius R=1 at t=1.
    !   In the reference: E0 = 0.244816 (for their non-dimensional setup).
    !   You can reuse that value in your dimensionless test.
    !-------------------------------------------------------------
    ! E0 = 0.244816d0
    E0 = 1.0d0

    ! innermost physical radial index
    i_inner = i_lo_phys

    !-------------------------------------------------------------
    ! Step 2.1: compute total volume (area in 2D) of inner ring
    !           V_ring = sum_j V(i_inner, j)
    !-------------------------------------------------------------
    V_ring = 0.0d0
    do j = j_lo_phys, j_hi_phys
        V_ring  = V_ring + mesh%cell_area(i_inner, j)
    end do
    !-------------------------------------------------------------
    ! Step 2.2: set high pressure Por in inner ring
    !           Por = (gamma - 1) * rho0 * E0 / V_ring
    !-------------------------------------------------------------
    p_inner = (gamma - 1.0d0) * rho0 * E0 / V_ring

    p_inner = (gamma - 1.0d0) * rho0 * E0 / pi / (mesh%radial_spacing(1)**2.0d0)

    do j = j_lo_phys, j_hi_phys
        prim%pressure(i_inner, j) = p_inner
    end do

    call update_conserved_from_primitive(prim, cons, gas, i_lo_phys, i_hi_phys, j_lo_phys, j_hi_phys)
    call update_boundary_primitive_conserved(mesh, prim, cons, gas)


    E_total = 0.d0
    do j = j_lo_phys, j_hi_phys
        do i = i_lo_phys, i_hi_phys
            E_total = E_total + cons%rho_E(i,j) * mesh%cell_area(i,j)
        end do
    end do
    print *, 'E_total at t=0 = ', E_total

    print *, "Initial solution successfully applied."


end subroutine initialize_sedov

end module initial_ops




