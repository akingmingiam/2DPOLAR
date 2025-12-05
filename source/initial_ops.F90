module initial_ops
    use mesh_types,       only : PolarMesh, pi
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

    double precision :: Ms,G,paf,daf,uaf,pbe,dbe,ube,R,theta,x


    pbe = 1.01325d+5
    dbe = 1.225d0
    ube = 0.d0
    ! Ms  = 2.670d0
    Ms  = 1.5d0
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
            ! Primitive variables
            !---------------------------
            ! prim%density(i,j)    = 1.0d0
            ! prim%pressure(i,j)   = 1.0d5
            ! prim%velocity_u(i,j) = 0.0d0
            ! prim%velocity_v(i,j) = 0.0d0


            ! R = mesh%r_center(i)
            ! if(R.le.1.0d0)then
            !     prim%density(i,j) = (dbe + (daf-dbe)*R/1.d0)
            !     prim%pressure(i,j) = (pbe + (paf-pbe)*R/1.d0)
            !     prim%velocity_u(i,j) = (uaf*R/1.d0)
            !     prim%velocity_v(i,j) = 0.d0
            ! else
            !     prim%density(i,j) = dbe
            !     prim%pressure(i,j) = pbe
            !     prim%velocity_u(i,j) = 0.d0
            !     prim%velocity_v(i,j) = 0.d0
            ! endif


            ! theta = mesh%theta_center(j)
            ! if (theta < 1.5d0*pi .and. theta > 0.5d0*pi) then
            !     prim%density(i,j)    = 1.0d0
            !     prim%pressure(i,j)   = 1.0d0
            !     prim%velocity_u(i,j) = 0.0d0     ! u_r = 0
            !     prim%velocity_v(i,j) = 0.0d0     ! u_θ = 0
            ! else
            !     prim%density(i,j)    = 0.125d0
            !     prim%pressure(i,j)   = 0.1d0
            !     prim%velocity_u(i,j) = 0.0d0
            !     prim%velocity_v(i,j) = 0.0d0
            ! end if

            theta = mesh%theta_center(j)
            r = mesh%r_center(i)
            x = r * cos(theta)

            prim%density(i,j)    = 1.0d0
            prim%pressure(i,j)   = 1.0d5
            prim%velocity_u(i,j) =   1.0d0 * cos(theta)
            prim%velocity_v(i,j) = - 1.0d0 * sin(theta)

            ! if (x < -1.0d0) then
            !     prim%density(i,j)    = daf
            !     prim%pressure(i,j)   = paf
            !     prim%velocity_u(i,j) =   uaf * cos(theta)  
            !     prim%velocity_v(i,j) = - uaf * sin(theta)    

            ! else
            !     prim%density(i,j)    = dbe
            !     prim%pressure(i,j)   = pbe
            !     prim%velocity_u(i,j) = ube
            !     prim%velocity_v(i,j) = ube
            ! end if

            ! if (x < -1.0d0) then
            !     prim%density(i,j)    = daf
            !     prim%pressure(i,j)   = paf
            !     prim%velocity_u(i,j) = uaf
            !     prim%velocity_v(i,j) = 0.0d0

            ! else
            !     prim%density(i,j)    = dbe
            !     prim%pressure(i,j)   = pbe
            !     prim%velocity_u(i,j) = ube
            !     prim%velocity_v(i,j) = ube
            ! end if
            
        end do
    end do

    call update_conserved_from_primitive(prim, cons, gas, i_lo_phys, i_hi_phys, j_lo_phys, j_hi_phys)
    call update_boundary_primitive_conserved(mesh, prim, cons, gas)

    print *, "Initial solution successfully applied."

end subroutine initialize_flow_fields

end module initial_ops
