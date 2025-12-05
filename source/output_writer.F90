! Module: output_writer
!
! Purpose:
!     Provide Tecplot ASCII output for primitive flow variables in 2D
!     polar → Cartesian CFD simulations.
!
!===============================================================================

module output_writer
    use mesh_types,    only : PolarMesh, BC, EQ
    use flow_fields,   only : PrimitiveVariables
    implicit none
    private
    ! Public interface
    public :: write_primitive_tecplot
    public :: print_time_info

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




end module output_writer
