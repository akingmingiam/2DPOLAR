module fluid_properties
    implicit none
    private
    
    public :: StiffenedGas
    public :: gas_air
    public :: gas_water
    public :: gas_custom
    public :: sound_speed_stiffened
    public :: specific_total_energy_stiffened
    public :: apply_pressure_floor

    !======================================================================
    !  Module: fluid_properties
    !  Purpose:
    !      Provide thermodynamic parameters for fluids using the
    !      Stiffened Gas Equation of State (SG-EOS).
    !
    !  EOS form:
    !      p = (γ - 1) ρ e  -  γ p∞
    !
    !  Parameters:
    !      γ      : adiabatic index
    !      p∞     : stiffness pressure   (Pa)
    !      cv     : specific heat at constant volume (J/kg/K)
    !      q      : reference energy shift (J/kg)
    !
    !  Notes:
    !      All parameters stored inside a derived type StiffenedGas.
    !======================================================================

    !======================================================================
    ! Stiffened Gas EOS parameter structure
    !======================================================================
    type :: StiffenedGas
        double precision :: gamma     ! ratio of specific heats
        double precision :: p_inf     ! stiffness constant (Pa)
        ! --- Pressure regularization parameters ---
        double precision :: pressure_floor_ratio  ! ratio for minimum allowed pressure
        double precision :: pressure_floor_eps    ! tiny positive offset to ensure numerical stability

    end type StiffenedGas


    !======================================================================
    ! Public interface
    !======================================================================

contains
    !======================================================================
    !  Function: gas_air, gas_water
    !  Purpose : return SG-EOS parameters
    !======================================================================
    function gas_air() result(gas)
        type(StiffenedGas) :: gas
        gas%gamma = 1.4d0
        gas%p_inf = 0.0d0                ! ideal gas → no stiffness
    end function gas_air

    function gas_water() result(gas)
        type(StiffenedGas) :: gas
        gas%gamma = 4.4d0
        gas%p_inf = 6.0d8                ! Pa
    end function gas_water
    
    function gas_custom(gamma, p_inf) result(gas)
        double precision, intent(in) :: gamma, p_inf
        type(StiffenedGas) :: gas
        gas%gamma = gamma
        gas%p_inf = p_inf
    end function gas_custom


    !===============================================================================
    !  Function: sound_speed_stiffened
    !
    !  Purpose:
    !      Compute sound speed for stiffened-gas EOS:
    !
    !          p = (γ - 1) ρ e - γ p_inf
    !
    !      gives:
    !
    !          c^2 = γ (p + p_inf) / ρ
    !
    !  Arguments:
    !      rho   : density
    !      p     : pressure
    !      gas   : stiffened-gas parameters (γ, p_inf)
    !
    !===============================================================================
    pure function sound_speed_stiffened(rho, p, gas) result(c)
        double precision, intent(in) :: rho, p
        type(StiffenedGas), intent(in) :: gas
        double precision :: c
        double precision :: arg

        if (rho <= 0.0d0) then
            c = 0.0d0
        else
            arg = gas%gamma * (p + gas%p_inf) / rho
            if (arg > 0.0d0) then
                c = sqrt(arg)
            else
                c = 0.0d0
            end if
        end if

    end function sound_speed_stiffened


    !===============================================================================
    !  Function: specific_internal_energy_stiffened
    !
    !  Purpose:
    !      Compute specific internal energy e for Stiffened Gas EOS:
    !
    !          e = (p + γ p_inf) / (ρ (γ - 1))
    !
    !===============================================================================
    pure function specific_internal_energy_stiffened(rho, p, gas) result(e)
        implicit none
        double precision, intent(in) :: rho, p
        type(StiffenedGas), intent(in) :: gas
        double precision :: e

        if (rho <= 0.0d0) then
            e = 0.0d0
        else
            e = (p + gas%gamma * gas%p_inf) / &
                (rho * (gas%gamma - 1.0d0))
        end if
    end function specific_internal_energy_stiffened


    !===============================================================================
    !  Function: specific_total_energy_stiffened
    !
    !  Purpose:
    !      Compute specific total energy E for Stiffened Gas EOS:
    !
    !         E = e + 0.5 (u^2 + v^2)
    !  NOTE:
    !      This is *specific* total energy (per unit mass).
    !      Conserved total energy density ρE must be obtained as ρ * E.

    !===============================================================================
    pure function specific_total_energy_stiffened(rho, p, u, v, gas) result(E)
        implicit none
        double precision, intent(in) :: rho, p, u, v
        type(StiffenedGas), intent(in) :: gas
        double precision :: E
        double precision :: e_int

        e_int = specific_internal_energy_stiffened(rho, p, gas)
        E  = e_int + 0.5d0 * (u*u + v*v)
    end function specific_total_energy_stiffened

    !===============================================================================
    ! Function: apply_pressure_floor
    !
    ! Purpose:
    !   The SG-EOS gives:
    !       c^2 = γ (p + p_inf) / ρ
    !
    !   For liquids (p_inf > 0), low internal energy can produce negative pressure
    !   which may violate hyperbolicity or cause numerical instability(c^2 < 0 when p < p_inf).
    !
    !   This routine enforces:
    !
    !       p >= - pressure_floor_ratio * p_inf + pressure_floor_eps
    !
    !===============================================================================
    function apply_pressure_floor(p, gas) result(p_fixed)
        double precision, intent(in) :: p
        type(StiffenedGas), intent(in) :: gas
        double precision :: p_fixed
        double precision :: p_min
        
        p_min = - gas%pressure_floor_ratio * gas%p_inf
        
        if (p < p_min) then
            p_fixed = p_min + gas%pressure_floor_eps
        else
            p_fixed = p
        end if
    end function apply_pressure_floor



end module fluid_properties
