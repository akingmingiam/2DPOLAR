! Module: output_writer
!
! Purpose:
!     Provide Tecplot ASCII output for primitive flow variables in 2D
!     polar → Cartesian CFD simulations.
!
!===============================================================================

module output_writer
    use mesh_types,    only : PolarMesh, BC, EQ, RECON, RIEMANN, TimeControlParameters, PolarGridParameters
    use flow_fields,   only : PrimitiveVariables
    implicit none
    private
    ! Public interface
    public :: write_primitive_tecplot
    public :: print_time_info
    public :: write_run_info

contains


!===============================================================================
!  Subroutine: write_primitive_tecplot
!
!  Purpose:
!       Export primitive variables (ρ, p, u, v) along with Cartesian coordinates
!       (x, y) to Tecplot ASCII file.
!
!  Output Format:
!       TITLE="Primitive Fields"
!       VARIABLES="X","Y","Density","Pressure","U","V"
!       ZONE I = Nr,  J = Ntheta,  F = POINT
!
!  Note:
!       Only the physical domain [i_lo_phys : i_hi_phys] × [j_lo_phys : j_hi_phys]
!       is exported. Ghost cells are excluded.
!
!===============================================================================
subroutine write_primitive_tecplot(t, step, output_step, dir, mesh, prim)
    implicit none

    !----------------------------
    ! Arguments
    !----------------------------
    double precision,         intent(in) :: t
    integer,                  intent(in) :: step
    integer,                  intent(in) :: output_step
    character(*),             intent(in) :: dir
    type(PolarMesh),          intent(in) :: mesh
    type(PrimitiveVariables), intent(in) :: prim

    !----------------------------
    ! Local variables
    !----------------------------
    integer :: i, j, unit, ios
    integer :: i_lo, i_hi, j_lo, j_hi
    character(len=256) :: errmsg
    character(len=256) :: filename

    !--------------------------------------------
    ! Construct filename: dir + "xxxxxx.plt"
    !--------------------------------------------
    write(filename, '(A,I6.6,".plt")') trim(dir), output_step

    !----------------------------------------
    ! Physical index range
    !----------------------------------------
    i_lo = mesh%params%i_lo_phys
    i_hi = mesh%params%i_hi_phys
    j_lo = mesh%params%j_lo_phys
    j_hi = mesh%params%j_hi_phys

    !----------------------------------------
    ! If theta-high BC is periodic, add 1 layer
    !----------------------------------------
    if (mesh%params%bc_theta_hi == BC%PERIODIC) then
        j_hi = j_hi + 1
    end if

    !----------------------------------------
    ! Open file with error handling
    !----------------------------------------
    unit = 50
    open(unit, file=filename, status="replace", action="write", &
         iostat=ios, iomsg=errmsg)

    if (ios /= 0) then
        print *, "ERROR: Failed to open output file:"
        print *, "  File   :", trim(filename)
        print *, "  Reason :", trim(errmsg)
        return
    end if

    !----------------------------------------
    ! Write Tecplot header
    !----------------------------------------
    write(unit,*) 'TITLE="Primitive Fields"'
    write(unit,*) 'VARIABLES="X","Y","Density","Pressure","Ur","Utheta","U","V"'
    write(unit,'("ZONE I=",I6,", J=",I6,", F=POINT")') &
         (i_hi - i_lo + 1), (j_hi - j_lo + 1)

    !----------------------------------------
    ! Export data
    !----------------------------------------
    do j = j_lo, j_hi
        do i = i_lo, i_hi

            select case(mesh%params%equation_type)

            case (EQ%CARTESIAN_TRANSFORM)
                write(unit,'(8ES26.18)')  &
                    mesh%x_center(i,j), mesh%y_center(i,j), &
                    prim%density(i,j), prim%pressure(i,j), &
                      prim%velocity_u(i,j)*cos(mesh%theta_center(j)) + prim%velocity_v(i,j)*sin(mesh%theta_center(j)), &
                    - prim%velocity_u(i,j)*sin(mesh%theta_center(j)) + prim%velocity_v(i,j)*cos(mesh%theta_center(j)), &
                    prim%velocity_u(i,j), prim%velocity_v(i,j)

            case (EQ%POLAR_EULER)
                write(unit,'(8ES26.18)')  &
                    mesh%x_center(i,j), mesh%y_center(i,j), &
                    prim%density(i,j), prim%pressure(i,j), &
                    prim%velocity_u(i,j), prim%velocity_v(i,j), &
                    prim%velocity_u(i,j)*cos(mesh%theta_center(j)) - prim%velocity_v(i,j)*sin(mesh%theta_center(j)), &
                    prim%velocity_u(i,j)*sin(mesh%theta_center(j)) + prim%velocity_v(i,j)*cos(mesh%theta_center(j))

            case default
                print *, "ERROR: Unknown equation type in compute_all_fluxes"
                stop

            end select

        end do
    end do

    close(unit)

    !----------------------------------------
    ! Print log to terminal (professional)
    !----------------------------------------
    print *, "-----------------------------------------------"
    print '(A,E12.5,A,I8,A,I8,A,A)', &
        "Data output completed for t=", t, " step=", step, &
        " output_step=", output_step, " file=", trim(filename)
    print *, "-----------------------------------------------"
end subroutine write_primitive_tecplot



subroutine print_time_info(step, t, dt)
    implicit none
    integer,          intent(in) :: step
    double precision, intent(in) :: t, dt

    print '(A,I8,A,E12.5,A,E12.5)',  &
          "Step=", step, " | t=", t, " | dt=", dt

end subroutine print_time_info



subroutine write_run_info(output_dir, grid, time)
    implicit none

    !-----------------接口-----------------
    character(len=*),           intent(in) :: output_dir
    type(PolarGridParameters),  intent(in) :: grid
    type(TimeControlParameters),intent(in) :: time

    !-----------------局部变量-----------------
    character(len=256) :: filename
    integer            :: iu, ios

    ! 把整数枚举翻成可读字符串
    character(len=64) :: eq_name, recon_name, riemann_solver
    character(len=64) :: bc_r_lo_name, bc_r_hi_name
    character(len=64) :: bc_th_lo_name, bc_th_hi_name

    ! 注意：这里假设 EQ, RECON, BC 这些 “枚举常量” 已经在某个 module 里定义，
    ! 并且当前文件 use 了那个 module。
    ! 例如：
    !   use equation_system_types,   only: EQ
    !   use reconstruction_types,    only: RECON
    !   use boundary_condition_types,only: BC
    !
    ! 如果你是在同一个 module 里，就不需要再 use。

    !-----------------1. 方程类型名-----------------
    select case (grid%equation_type)
    case (EQ%POLAR_EULER)
        eq_name = 'POLAR_EULER'
    case (EQ%CARTESIAN_TRANSFORM)
        eq_name = 'CARTESIAN_TRANSFORM'
    case default
        eq_name = 'UNKNOWN_EQUATION'
    end select

    !-----------------2. 重构方式名-----------------
    select case (grid%reconstruction)
    case (RECON%FIRST_ORDER)
        recon_name = 'FIRST_ORDER (piecewise constant)'
    case (RECON%MUSCL_MINMOD)
        recon_name = 'MUSCL (minmod limiter)'
    case (RECON%MUSCL_VANLEER)
        recon_name = 'MUSCL (van Leer limiter)'
    case default
        recon_name = 'UNKNOWN_RECONSTRUCTION'
    end select

    !-----------------3. 黎曼求解器名-----------------
    select case (grid%riemann_solver)
    case (RIEMANN%HLL_SIMPLE)
        riemann_solver = 'HLL_SIMPLE'
    case (RIEMANN%HLL)
    riemann_solver = 'HLL'
    case (RIEMANN%HLLC)
        riemann_solver = 'HLLC'
    case default
        riemann_solver = 'UNKNOWN_RIEMANN_SOLVER'
    end select

    !-----------------4. 边界条件名-----------------
    select case (grid%bc_r_lo)
    case (BC%AXIS)
        bc_r_lo_name = 'AXIS'
    case (BC%EXTRAPOL)
        bc_r_lo_name = 'EXTRAPOLATION'
    case default
        bc_r_lo_name = 'UNKNOWN_BC'
    end select

    select case (grid%bc_r_hi)
    case (BC%AXIS)
        bc_r_hi_name = 'AXIS'
    case (BC%EXTRAPOL)
        bc_r_hi_name = 'EXTRAPOLATION'
    case default
        bc_r_hi_name = 'UNKNOWN_BC'
    end select

    select case (grid%bc_theta_lo)
    case (BC%PERIODIC)
        bc_th_lo_name = 'PERIODIC'
    case default
        bc_th_lo_name = 'UNKNOWN_BC'
    end select

    select case (grid%bc_theta_hi)
    case (BC%PERIODIC)
        bc_th_hi_name = 'PERIODIC'
    case default
        bc_th_hi_name = 'UNKNOWN_BC'
    end select

    !-----------------5. 打开 result/run_info.txt-----------------
    filename = trim(output_dir)//'/run_info.txt'

    open(newunit=iu, file=filename, status='replace', action='write', &
         iostat=ios)
    if (ios /= 0) then
        write(*,*) 'WARNING: cannot open ', trim(filename), ' iostat=', ios
        return
    end if

    !-----------------6. 写入说明-----------------
    write(iu,'(A)') '=============================================='
    write(iu,'(A)') '  Simulation run information'
    write(iu,'(A)') '=============================================='
    write(iu,'(A)') ''

    ! --- 网格参数 ---
    write(iu,'(A)') '--- Polar grid parameters ---'
    write(iu,'(A,I8)')  'n_r         = ', grid%n_r
    write(iu,'(A,I8)')  'n_theta     = ', grid%n_theta
    write(iu,'(A,I8)')  'guard_cells = ', grid%guard_cells
    write(iu,'(A,1X,ES16.8)') 'r_min       = ', grid%r_min
    write(iu,'(A,1X,ES16.8)') 'r_max       = ', grid%r_max
    write(iu,'(A,1X,ES16.8)') 'theta_min   = ', grid%theta_min
    write(iu,'(A,1X,ES16.8)') 'theta_max   = ', grid%theta_max
    write(iu,'(A)') ''

    ! --- 索引范围（方便调试） ---
    write(iu,'(A)') '--- Index ranges ---'
    write(iu,'(A,2(I8,1X))') 'i_lo_phys, i_hi_phys   = ', grid%i_lo_phys, grid%i_hi_phys
    write(iu,'(A,2(I8,1X))') 'j_lo_phys, j_hi_phys   = ', grid%j_lo_phys, grid%j_hi_phys
    write(iu,'(A,2(I8,1X))') 'i_lo_comp, i_hi_comp   = ', grid%i_lo_comp, grid%i_hi_comp
    write(iu,'(A,2(I8,1X))') 'j_lo_comp, j_hi_comp   = ', grid%j_lo_comp, grid%j_hi_comp
    write(iu,'(A)') ''

    ! --- 方程 & 重构 & 边界 ---
    write(iu,'(A)') '--- Equation & reconstruction ---'
    write(iu,'(A,A)') 'equation_type   = ', trim(eq_name)
    write(iu,'(A,A)') 'reconstruction  = ', trim(recon_name)
    write(iu,'(A,A)') 'riemann_siolver = ', trim(riemann_solver)
    write(iu,'(A)') ''

    write(iu,'(A)') '--- Boundary conditions ---'
    write(iu,'(A,A)') 'bc_r_lo      = ', trim(bc_r_lo_name)
    write(iu,'(A,A)') 'bc_r_hi      = ', trim(bc_r_hi_name)
    write(iu,'(A,A)') 'bc_theta_lo  = ', trim(bc_th_lo_name)
    write(iu,'(A,A)') 'bc_theta_hi  = ', trim(bc_th_hi_name)
    write(iu,'(A)') ''

    ! --- 时间控制 ---
    write(iu,'(A)') '--- Time control parameters ---'
    write(iu,'(A,1X,ES16.8)') 'CFL          = ', time%CFL
    write(iu,'(A,1X,ES16.8)') 't_start      = ', time%t
    write(iu,'(A,1X,ES16.8)') 't_end        = ', time%t_end
    write(iu,'(A,1X,ES16.8)') 'dt_init      = ', time%dt
    write(iu,'(A,1X,ES16.8)') 'dt_max       = ', time%dt_max
    write(iu,'(A,1X,ES16.8)') 'dt_min       = ', time%dt_min
    write(iu,'(A,I8)')        'step_start   = ', time%step
    write(iu,'(A,I8)')        'max_steps    = ', time%max_steps
    write(iu,'(A,I8)')        'output_step(start) = ', time%output_step
    write(iu,'(A,I8)')        'max_output_steps   = ', time%max_output_steps
    write(iu,'(A)') ''
    write(iu,'(A)') 'Note: output times are stored in time%time_outputs(:)'
    write(iu,'(A)') ''

    close(iu)

end subroutine write_run_info





end module output_writer
