module time_step_control
    use mesh_types,       only : PolarMesh
    use flow_fields,      only : PrimitiveVariables
    use fluid_properties, only : StiffenedGas, sound_speed_stiffened
    implicit none
    private

    public :: compute_global_timestep

contains

!===============================================================================
!  Subroutine: compute_global_timestep
!
!  Purpose:
!      Compute global CFL time step for a single-phase flow using a
!      volume-based formula:
!
!          dt_cfl = CFL * min_{i,j} [ sqrt(V_ij) / (|u_ij| + c_ij) ]
!
!  Arguments:
!      mesh    : [in]  polar mesh (geometry + index ranges + cell_area)
!      prim    : [in]  primitive variables (ρ, p, u, v)
!      gas     : [in]  stiffened-gas parameters (γ, p_inf)
!      CFL     : [in]  CFL number (0 < CFL <= 1)
!      dt_cfl  : [out] global time step from CFL condition
!
!===============================================================================

function compute_global_timestep_cartesian(mesh, prim, gas, CFL, dt_max, dt_min) result(dt)
    type(PolarMesh),          intent(in) :: mesh
    type(PrimitiveVariables), intent(in) :: prim
    type(StiffenedGas),       intent(in) :: gas
    double precision,         intent(in) :: CFL
    double precision,         intent(in) :: dt_max
    double precision,         intent(in) :: dt_min

    double precision :: dt
    double precision :: dt_cfl, dt_cell
    double precision :: rho, p, u, v
    double precision :: vmag, c_sound, a_char
    double precision :: cell_area, h_eff
    double precision, parameter :: tiny_speed = 1.0d-14

    integer :: i, j
    integer :: i_lo, i_hi, j_lo, j_hi

    !-------------------------------------------------------------
    ! Physical domain indices
    !-------------------------------------------------------------
    i_lo = mesh%params%i_lo_phys
    i_hi = mesh%params%i_hi_phys
    j_lo = mesh%params%j_lo_phys
    j_hi = mesh%params%j_hi_phys

    !-------------------------------------------------------------
    ! CFL + initial dt bound
    !-------------------------------------------------------------
    dt_cfl = dt_max

    !-------------------------------------------------------------
    ! Loop over physical cells only
    !-------------------------------------------------------------
    do j = j_lo, j_hi
        do i = i_lo, i_hi

            rho = prim%density(i,j)
            p   = prim%pressure(i,j)
            u   = prim%velocity_u(i,j)
            v   = prim%velocity_v(i,j)

            ! speed magnitude
            vmag = sqrt(u*u + v*v)

            ! local sound speed from stiffened-gas EOS
            c_sound = sound_speed_stiffened(rho, p, gas)

            ! characteristic signal speed
            a_char = vmag + c_sound

            if (a_char <= tiny_speed) cycle    ! static cell → no restriction

            ! area-based cell “volume” (for 2D)
            cell_area = mesh%cell_area(i,j)

            ! effective length scale: sqrt(area)   [for 2D]
            if (cell_area > 0.0d0) then
                h_eff = sqrt(cell_area)
            else
                cycle
            end if

            ! local time-step
            dt_cell = CFL * h_eff / a_char
            if (dt_cell < dt_cfl) dt_cfl = dt_cell
        end do
    end do

    !-------------------------------------------------------------
    ! Return dt with safety lower bound
    !-------------------------------------------------------------
    dt = max(dt_min, dt_cfl)
end function compute_global_timestep_cartesian


function compute_global_timestep(mesh, prim, gas, CFL, dt_max, dt_min) result(dt)
    type(PolarMesh),          intent(in) :: mesh
    type(PrimitiveVariables), intent(in) :: prim
    type(StiffenedGas),       intent(in) :: gas
    double precision,         intent(in) :: CFL, dt_max, dt_min
    double precision :: dt

    integer :: i, j, i_lo, i_hi, j_lo, j_hi
    double precision :: rho, p, ur, uth, c_sound
    double precision :: dr, dtheta, r, dl_r, dl_th
    double precision :: a_r, a_th, dt_cell
    double precision, parameter :: tiny = 1.0d-14

    i_lo = mesh%params%i_lo_phys
    i_hi = mesh%params%i_hi_phys
    j_lo = mesh%params%j_lo_phys
    j_hi = mesh%params%j_hi_phys

    dt = dt_max

    do j = j_lo, j_hi
        do i = i_lo, i_hi

            rho = prim%density(i,j)
            p   = prim%pressure(i,j)
            ur  = prim%velocity_u(i,j)   ! u_r
            uth = prim%velocity_v(i,j)   ! u_θ
            c_sound = sound_speed_stiffened(rho, p, gas)

            dr     = mesh%radial_spacing(i)
            dtheta = mesh%angular_spacing(j)
            r      = mesh%r_center(i)

            dl_r  = dr
            dl_th = max(r * dtheta, tiny)

            a_r  = abs(ur)  + c_sound
            a_th = abs(uth) + c_sound

            if (a_r > tiny)  dt_cell = CFL * dl_r  / a_r
            if (a_th > tiny) dt_cell = min(dt_cell, CFL * dl_th / a_th)

            if (dt_cell < dt) dt = dt_cell

        end do
    end do

    dt = max(dt_min, min(dt, dt_max))
end function compute_global_timestep

end module time_step_control
