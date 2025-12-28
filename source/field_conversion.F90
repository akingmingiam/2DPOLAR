module field_conversion
    !======================================================================
    !  Module: field_conversion
    !  Purpose:
    !      Provide conversion routines between primitive (ρ, p, u, v)
    !      and conserved (ρ, ρu, ρv, ρE) variables for compressible flow.
    !======================================================================

    use flow_fields,      only : PrimitiveVariables, ConservedVariables
    use fluid_properties, only : StiffenedGas, specific_total_energy_stiffened, &
                                 apply_pressure_floor
    implicit none
    private

    public :: update_conserved_from_primitive
    public :: update_primitive_from_conserved

contains

    !==================================================================
    !  Subroutine: update_conserved_from_primitive
    !  Purpose   : Compute conserved variables (ρ, ρu, ρv, ρE)
    !              from primitive variables using Stiffened-Gas EOS:
    !
    !              p = (γ - 1)(ρE - ½ρ(u²+v²)) - γp∞
    !
    !              ⇒ ρE = [p + γp∞]/(γ - 1) + ½ρ(u²+v²)
    !
    !  Inputs :
    !      prim   – primitive fields (ρ, p, u, v)
    !      gas    – fluid parameters (γ, p∞)
    !      i_lo_phys, i_hi_phys, j_lo_phys, j_hi_phys – physical indices
    !
    !  Outputs:
    !      cons   – conserved fields updated accordingly
    !==================================================================
    subroutine update_conserved_from_primitive(prim, cons, gas, &
                                               i_lo_phys, i_hi_phys, &
                                               j_lo_phys, j_hi_phys)
        type(PrimitiveVariables), intent(in)    :: prim
        type(ConservedVariables), intent(inout) :: cons
        type(StiffenedGas),       intent(in)    :: gas
        integer, intent(in) :: i_lo_phys, i_hi_phys, j_lo_phys, j_hi_phys

        integer :: i, j
        double precision :: rho, p, u, v
        ! double precision :: E_int, E_tot

        do j = j_lo_phys, j_hi_phys
            do i = i_lo_phys, i_hi_phys
                rho = prim%density(i,j)
                p   = prim%pressure(i,j)
                u   = prim%velocity_u(i,j)
                v   = prim%velocity_v(i,j)


                cons%rho(i,j)   = rho
                cons%rho_u(i,j) = rho * u
                cons%rho_v(i,j) = rho * v
                cons%rho_E(i,j) = rho * specific_total_energy_stiffened(rho, p, u, v, gas)
            end do
        end do
    end subroutine update_conserved_from_primitive

    !===================================================================
    ! Update PRIMITIVE variables from conserved fields
    ! stiffened-gas EOS:
    !     ρE = (p + γ p_inf)/(γ−1) + 0.5 ρ (u² + v²)
    !     ⇒ p = (γ−1)(ρE − 0.5 ρ (u²+v²)) − γ p_inf
    !===================================================================
    subroutine update_primitive_from_conserved(prim, cons, gas, &
                                       i_lo_phys, i_hi_phys, &
                                       j_lo_phys, j_hi_phys)
        type(PrimitiveVariables), intent(inout) :: prim
        type(ConservedVariables), intent(in)    :: cons
        type(StiffenedGas),       intent(in)    :: gas
        integer, intent(in) :: i_lo_phys, i_hi_phys, j_lo_phys, j_hi_phys

        integer :: i, j
        double precision :: rho, rho_u, rho_v, rho_E
        double precision :: u, v, v2, p
        double precision, parameter :: tiny_rho = 1d-14


        do j = j_lo_phys, j_hi_phys
            do i = i_lo_phys, i_hi_phys

                rho   = cons%rho(i,j)
                rho_u = cons%rho_u(i,j)
                rho_v = cons%rho_v(i,j)
                rho_E = cons%rho_E(i,j)

                ! if (rho <= tiny_rho) then
                !     prim%density(i,j)    = tiny_rho
                !     prim%velocity_u(i,j) = 0d0
                !     prim%velocity_v(i,j) = 0d0
                !     prim%pressure(i,j)   = 0d0
                ! else
                !     u  = rho_u / rho
                !     v  = rho_v / rho
                !     v2 = u*u + v*v

                !     p = (gas%gamma - 1d0)*( rho_E - 0.5d0*rho*v2 ) &
                !         - gas%gamma*gas%p_inf
                !     p = apply_pressure_floor(p, gas)

                !     prim%density(i,j)    = rho
                !     prim%velocity_u(i,j) = u
                !     prim%velocity_v(i,j) = v
                !     prim%pressure(i,j)   = p
                ! end if

                u  = rho_u / rho
                v  = rho_v / rho
                v2 = u*u + v*v

                p = (gas%gamma - 1d0)*( rho_E - 0.5d0*rho*v2 ) &
                    - gas%gamma*gas%p_inf
                prim%density(i,j)    = rho
                prim%velocity_u(i,j) = u
                prim%velocity_v(i,j) = v
                prim%pressure(i,j)   = p

            end do
        end do

    end subroutine update_primitive_from_conserved

end module field_conversion