!===========================================================
! File   : riemann_solver.F90
! Module : riemann_solver
! Purpose:
!   Riemann solvers for compressible Euler equations
!   (Cartesian velocity components), with stiffened-gas EOS.
!
!   Currently implemented:
!     - HLL Riemann solver: hll_flux_solver
!===========================================================
module riemann_solver
    use fluid_properties, only : StiffenedGas, sound_speed_stiffened, &
                                 specific_total_energy_stiffened
    implicit none
    private

    public :: hll_flux_solver

contains


    !---------------------------------------------------------------
    !  Subroutine: hll_flux_solver
    !
    !  Purpose:
    !     Compute HLL numerical flux for 2D Euler equations
    !     with Cartesian velocity components and stiffened-gas EOS.
    !
    !  Inputs:
    !     (rhoL,uL,vL,pL) : left state (primitive)
    !     (rhoR,uR,vR,pR) : right state (primitive)
    !     (nx, ny)        : unit normal vector (from L -> R)
    !     gas             : stiffened-gas EOS parameters
    !
    !  Outputs:
    !     f_rho, f_rho_u, f_rho_v, f_rho_E : numerical flux components
    !
    !---------------------------------------------------------------
    subroutine hll_flux_solver(rhoL,uL,vL,pL, rhoR,uR,vR,pR, &
                               nx,ny, gas, &
                               f_rho, f_rho_u, f_rho_v, f_rho_E)

        double precision, intent(in)   :: rhoL,uL,vL,pL
        double precision, intent(in)   :: rhoR,uR,vR,pR
        double precision, intent(in)   :: nx,ny
        type(StiffenedGas), intent(in) :: gas
        double precision, intent(out)  :: f_rho, f_rho_u, f_rho_v, f_rho_E

        ! --- local variables ---
        double precision :: unL, unR       ! normal velocities
        double precision :: cL, cR         ! sound speeds
        double precision :: SL, SR, denom
        double precision :: rhoEL, rhoER
        double precision :: U_L(4), U_R(4)
        double precision :: FL(4), FR(4)
        double precision :: FHLL(4)
        double precision, parameter :: tiny = 1.0d-14

        !-----------------------------------------------------------
        ! Left state: conservative variables and physical flux
        !-----------------------------------------------------------
        rhoEL = specific_total_energy_stiffened(rhoL, pL, uL, vL, gas) * rhoL

        U_L(1) = rhoL
        U_L(2) = rhoL * uL
        U_L(3) = rhoL * vL
        U_L(4) = rhoEL

        unL   = uL*nx + vL*ny
        FL(1) = rhoL * unL
        FL(2) = rhoL * uL * unL + pL * nx
        FL(3) = rhoL * vL * unL + pL * ny
        FL(4) = (rhoEL + pL) * unL

        !-----------------------------------------------------------
        ! Right state: conservative variables and physical flux
        !-----------------------------------------------------------
        rhoER = specific_total_energy_stiffened(rhoR, pR, uR, vR, gas) * rhoR

        U_R(1) = rhoR
        U_R(2) = rhoR * uR
        U_R(3) = rhoR * vR
        U_R(4) = rhoER

        unR   = uR*nx + vR*ny
        FR(1) = rhoR * unR
        FR(2) = rhoR * uR * unR + pR * nx
        FR(3) = rhoR * vR * unR + pR * ny
        FR(4) = (rhoER + pR) * unR

        !-----------------------------------------------------------
        ! Wave speed estimates for HLL
        !-----------------------------------------------------------
        cL = sound_speed_stiffened(rhoL, pL, gas)
        cR = sound_speed_stiffened(rhoR, pR, gas)

        SL = min( unL - cL, unR - cR )
        SR = max( unL + cL, unR + cR )

        !-----------------------------------------------------------
        ! HLL flux selection
        !-----------------------------------------------------------
        if (SL >= 0.0d0) then
            FHLL(:) = FL(:)

        else if (SR <= 0.0d0) then
            FHLL(:) = FR(:)

        else
            denom = SR - SL
            if (denom > tiny) then
                FHLL(:) = ( SR*FL(:) - SL*FR(:) + SL*SR*(U_R(:) - U_L(:)) ) / denom
            else
                ! Very rare pathological case: fall back to average flux
                FHLL(:) = 0.5d0 * (FL(:) + FR(:))
            end if
        end if

        f_rho   = FHLL(1)
        f_rho_u = FHLL(2)
        f_rho_v = FHLL(3)
        f_rho_E = FHLL(4)

    end subroutine hll_flux_solver

end module riemann_solver
