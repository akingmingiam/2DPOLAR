module riemann_solver
    use fluid_properties, only : StiffenedGas, sound_speed_stiffened, &
                                 specific_total_energy_stiffened
     use mesh_types,      only : RIEMANN
    implicit none
    private

    public :: flux_solver

contains
   

    subroutine flux_solver(riemann_solver, rhoL,uL,vL,pL, rhoR,uR,vR,pR, &
                        nx,ny, gas, &
                        f_rho, f_rho_u, f_rho_v, f_rho_E)

        integer,            intent(in)   :: riemann_solver
        double precision,   intent(in)   :: rhoL,uL,vL,pL
        double precision,   intent(in)   :: rhoR,uR,vR,pR
        double precision,   intent(in)   :: nx,ny
        type(StiffenedGas), intent(in)   :: gas
        double precision,   intent(out)  :: f_rho, f_rho_u, f_rho_v, f_rho_E

        ! clean default (avoid uninitialized)
        f_rho   = 0.0d0
        f_rho_u = 0.0d0
        f_rho_v = 0.0d0
        f_rho_E = 0.0d0

        select case (riemann_solver)

        case (RIEMANN%HLL_SIMPLE)
        call hll_flux_solver_simple( rhoL,uL,vL,pL, rhoR,uR,vR,pR, &
                                     nx,ny, gas, &
                                     f_rho, f_rho_u, f_rho_v, f_rho_E )

        case (RIEMANN%HLL)
            call hll_flux_solver( rhoL,uL,vL,pL, rhoR,uR,vR,pR, &
                                  nx,ny, gas, &
                                  f_rho, f_rho_u, f_rho_v, f_rho_E )

        case (RIEMANN%HLLC)
            call hllc_flux_solver( rhoL,uL,vL,pL, rhoR,uR,vR,pR, &
                                   nx,ny, gas, &
                                   f_rho, f_rho_u, f_rho_v, f_rho_E )

        case default
            write(*,*) "ERROR: Unknown Riemann solver id = ", riemann_solver
            stop
        end select

    end subroutine flux_solver


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
    subroutine hll_flux_solver_simple(rhoL,uL,vL,pL, rhoR,uR,vR,pR, &
                                      nx,ny, gas, &
                                      f_rho, f_rho_u, f_rho_v, f_rho_E)

        double precision,   intent(in)   :: rhoL,uL,vL,pL
        double precision,   intent(in)   :: rhoR,uR,vR,pR
        double precision,   intent(in)   :: nx,ny
        type(StiffenedGas), intent(in)   :: gas
        double precision,   intent(out)  :: f_rho, f_rho_u, f_rho_v, f_rho_E

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

    end subroutine hll_flux_solver_simple

    
    subroutine hll_flux_solver(rhoL,uL,vL,pL, rhoR,uR,vR,pR, &
                               nx,ny, gas, &
                               f_rho, f_rho_u, f_rho_v, f_rho_E)

        double precision,   intent(in)   :: rhoL,uL,vL,pL
        double precision,   intent(in)   :: rhoR,uR,vR,pR
        double precision,   intent(in)   :: nx,ny
        type(StiffenedGas), intent(in)   :: gas
        double precision,   intent(out)  :: f_rho, f_rho_u, f_rho_v, f_rho_E

        ! --- local variables ---
        double precision :: unL, unR       ! normal velocities
        double precision :: cL, cR         ! sound speeds
        double precision :: SL, SR, denom
        double precision :: rhoEL, rhoER
        double precision :: U_L(4), U_R(4)
        double precision :: FL(4), FR(4)
        double precision :: FHLL(4)
        double precision, parameter :: tiny = 1.0d-14

        ! Roe-averaged quantities (Toro p.328)
        double precision :: sqrt_rhoL, sqrt_rhoR, inv_sum_sqrt
        double precision :: u_roe, v_roe, un_roe
        double precision :: H_L, H_R, H_roe
        double precision :: u2_roe, a_roe
        double precision :: gamma

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

        !--- Roe averages ------------------------------------------------
        sqrt_rhoL    = sqrt(rhoL)
        sqrt_rhoR    = sqrt(rhoR)
        inv_sum_sqrt = 1.0d0 / (sqrt_rhoL + sqrt_rhoR)

        ! Roe-averaged velocity components
        u_roe = (sqrt_rhoL*uL + sqrt_rhoR*uR) * inv_sum_sqrt
        v_roe = (sqrt_rhoL*vL + sqrt_rhoR*vR) * inv_sum_sqrt

        ! Roe-averaged normal velocity
        un_roe = u_roe*nx + v_roe*ny

        ! Left / right total enthalpy H = (rhoE + p) / rho
        H_L = (rhoEL + pL) / rhoL
        H_R = (rhoER + pR) / rhoR

        ! Roe-averaged enthalpy
        H_roe = (sqrt_rhoL*H_L + sqrt_rhoR*H_R) * inv_sum_sqrt

        ! |u_roe|^2
        u2_roe = u_roe*u_roe + v_roe*v_roe

        ! Gas gamma (constant-gamma ideal / stiffened gas)
        gamma = gas%gamma  

        ! Roe-averaged sound speed (ideal-gas-like formula)
        a_roe = sqrt( max( tiny, (gamma - 1.0d0) * (H_roe - 0.5d0*u2_roe) ) )

        ! Einfeldt / Roe wave-speed estimates
        SL = min( unL - cL, un_roe - a_roe )
        SR = max( unR + cR, un_roe + a_roe )

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


    subroutine hllc_flux_solver(rhoL,uL,vL,pL, rhoR,uR,vR,pR, &
                                nx,ny, gas,                   &
                                f_rho, f_rho_u, f_rho_v, f_rho_E)

        implicit none

        !---------------- 输入/输出 ----------------
        double precision,   intent(in)   :: rhoL,uL,vL,pL
        double precision,   intent(in)   :: rhoR,uR,vR,pR
        double precision,   intent(in)   :: nx,ny
        type(StiffenedGas), intent(in)   :: gas
        double precision,   intent(out)  :: f_rho, f_rho_u, f_rho_v, f_rho_E

        !---------------- 局部变量 ----------------
        double precision, parameter :: tiny = 1.0d-14

        ! primitive -> conservative
        double precision :: EL, ER              ! specific total energy
        double precision :: rhoEL, rhoER        ! total energy density
        double precision :: U_L(4), U_R(4)
        double precision :: FL(4), FR(4), F(4)

        ! 声速及法向速度
        double precision :: unL, unR
        double precision :: cL,  cR

        ! 波速
        double precision :: SL, SR, SM          ! left, right, contact wave

        ! Roe 平均
        double precision :: sqrt_rhoL, sqrt_rhoR, inv_sum_sqrt
        double precision :: u_roe, v_roe, un_roe
        double precision :: H_L, H_R, H_roe
        double precision :: u2_roe, a_roe
        double precision :: gamma

        ! 星区域（*）状态
        double precision :: rhoL_star, rhoR_star
        double precision :: uL_star, vL_star
        double precision :: uR_star, vR_star
        double precision :: EL_star, ER_star       ! specific total energy in star
        double precision :: rhoEL_star, rhoER_star ! energy density in star

        ! 方便写公式的中间量
        double precision :: SLm_unL, SRm_unR
        double precision :: denom

        !============================================================
        ! 1. 左右状态：保守量 & 物理通量
        !============================================================
        EL   = specific_total_energy_stiffened(rhoL, pL, uL, vL, gas)
        ER   = specific_total_energy_stiffened(rhoR, pR, uR, vR, gas)
        rhoEL = rhoL * EL
        rhoER = rhoR * ER

        U_L(1) = rhoL
        U_L(2) = rhoL * uL
        U_L(3) = rhoL * vL
        U_L(4) = rhoEL

        U_R(1) = rhoR
        U_R(2) = rhoR * uR
        U_R(3) = rhoR * vR
        U_R(4) = rhoER

        unL = uL*nx + vL*ny
        unR = uR*nx + vR*ny

        FL(1) = rhoL * unL
        FL(2) = rhoL * uL * unL + pL * nx
        FL(3) = rhoL * vL * unL + pL * ny
        FL(4) = (rhoEL + pL) * unL

        FR(1) = rhoR * unR
        FR(2) = rhoR * uR * unR + pR * nx
        FR(3) = rhoR * vR * unR + pR * ny
        FR(4) = (rhoER + pR) * unR

        !============================================================
        ! 2. 声速 & Roe 波速估计
        !============================================================
        cL = sound_speed_stiffened(rhoL, pL, gas)
        cR = sound_speed_stiffened(rhoR, pR, gas)

        ! Roe averages
        sqrt_rhoL    = sqrt(rhoL)
        sqrt_rhoR    = sqrt(rhoR)
        inv_sum_sqrt = 1.0d0 / (sqrt_rhoL + sqrt_rhoR)

        u_roe = (sqrt_rhoL*uL + sqrt_rhoR*uR) * inv_sum_sqrt
        v_roe = (sqrt_rhoL*vL + sqrt_rhoR*vR) * inv_sum_sqrt
        un_roe = u_roe*nx + v_roe*ny

        H_L = (rhoEL + pL) / rhoL
        H_R = (rhoER + pR) / rhoR
        H_roe = (sqrt_rhoL*H_L + sqrt_rhoR*H_R) * inv_sum_sqrt

        u2_roe = u_roe*u_roe + v_roe*v_roe
        gamma  = gas%gamma

        a_roe = sqrt( max(tiny, (gamma - 1.0d0) * (H_roe - 0.5d0*u2_roe) ) )

        ! 左右波速 (Einfeldt / Roe)
        SL = min( unL - cL, un_roe - a_roe )
        SR = max( unR + cR, un_roe + a_roe )

        ! 如果极端退化，直接退成简单 HLL
        if (SR <= SL + tiny) then
            ! HLL 通量
            denom = SR - SL
            if (denom > tiny) then
                F(:) = ( SR*FL(:) - SL*FR(:) + SL*SR*(U_R(:) - U_L(:)) ) / denom
            else
                F(:) = 0.5d0 * (FL(:) + FR(:))
            end if
            f_rho   = F(1)
            f_rho_u = F(2)
            f_rho_v = F(3)
            f_rho_E = F(4)
            return
        end if

        !============================================================
        ! 3. 接触波速度 SM (sstar)
        !    SM = sstar = contact wave speed
        !============================================================
        SLm_unL = SL - unL
        SRm_unR = SR - unR

        denom = rhoL*SLm_unL - rhoR*SRm_unR
        if (abs(denom) < tiny) denom = sign(tiny, denom)

        SM = (pR - pL + rhoL*unL*SLm_unL - rhoR*unR*SRm_unR) / denom

        !============================================================
        ! 4. 构造 * 区域的状态（左星、右星）
        !    先算密度，再算速度，再算能量
        !============================================================
        ! 左星密度
        rhoL_star = rhoL * (SLm_unL) / (SL - SM + tiny)
        ! 右星密度
        rhoR_star = rhoR * (SRm_unR) / (SR - SM + tiny)

        ! 左星速度：保持切向速度不变，只替换法向为 SM
        uL_star = uL + (SM - unL)*nx
        vL_star = vL + (SM - unL)*ny

        ! 右星速度
        uR_star = uR + (SM - unR)*nx
        vR_star = vR + (SM - unR)*ny

        ! 左星能量（specific total energy）
        EL_star = EL + (SM - unL) * ( SM + pL/(rhoL*(SL - unL + tiny)) )
        rhoEL_star = rhoL_star * EL_star

        ! 右星能量
        ER_star = ER + (SM - unR) * ( SM + pR/(rhoR*(SR - unR + tiny)) )
        rhoER_star = rhoR_star * ER_star

        ! 左右星保守量
        ! U*_L
        !  [ rho*, rho* u*, rho* v*, rho*E* ]
        ! U*_R 同理
        !============================================================
        ! 5. HLLC 通量选择
        !============================================================
        if (SL >= 0.0d0) then
            ! 全部由左状态决定
            F(:) = FL(:)

        else if (SR <= 0.0d0) then
            ! 全部由右状态决定
            F(:) = FR(:)

        else if (SM >= 0.0d0) then
            ! 0 在左星扇区内：F = F_L + SL (U_L* - U_L)
            F(1) = FL(1) + SL*(rhoL_star            - U_L(1))
            F(2) = FL(2) + SL*(rhoL_star*uL_star    - U_L(2))
            F(3) = FL(3) + SL*(rhoL_star*vL_star    - U_L(3))
            F(4) = FL(4) + SL*(rhoEL_star           - U_L(4))

        else
            ! SM < 0 < SR：0 在右星扇区内：F = F_R + SR (U_R* - U_R)
            F(1) = FR(1) + SR*(rhoR_star            - U_R(1))
            F(2) = FR(2) + SR*(rhoR_star*uR_star    - U_R(2))
            F(3) = FR(3) + SR*(rhoR_star*vR_star    - U_R(3))
            F(4) = FR(4) + SR*(rhoER_star           - U_R(4))
        end if

        !============================================================
        ! 6. 输出
        !============================================================
        f_rho   = F(1)
        f_rho_u = F(2)
        f_rho_v = F(3)
        f_rho_E = F(4)

    end subroutine hllc_flux_solver



end module riemann_solver
