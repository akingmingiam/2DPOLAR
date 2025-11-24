! subroutine write_primitive_tecplot(filename, mesh, prim)
!     !----------------------------
!     ! Arguments
!     !----------------------------
!     character(*),             intent(in) :: filename
!     type(PolarMesh),          intent(in) :: mesh
!     type(PrimitiveVariables), intent(in) :: prim

!     !----------------------------
!     ! Local variables
!     !----------------------------
!     integer :: i, j, unit
!     integer :: i_lo, i_hi
!     integer :: j_lo, j_hi
!     integer :: ios                      ! I/O status code
!     character(len=256) :: errmsg        ! I/O error message

!     !----------------------------------------
!     ! Physical index range
!     !----------------------------------------
!     i_lo = mesh%params%i_lo_phys
!     i_hi = mesh%params%i_hi_phys
!     j_lo = mesh%params%j_lo_phys
!     j_hi = mesh%params%j_hi_phys

!     !--------------------------------------------
!     ! If theta-high BC is periodic, add 1 layer
!     !--------------------------------------------
!     if (mesh%params%bc_theta_hi == BC_PERIODIC) then
!         j_hi = j_hi + 1
!     end if

!     !----------------------------------------
!     ! Try opening the file with error handling
!     !----------------------------------------
!     unit = 50
!     open(unit, file=filename, status="replace", action="write", &
!          iostat=ios, iomsg=errmsg)

!     if (ios /= 0) then
!         print *, "ERROR: Failed to open output file:"
!         print *, "  File: ", trim(filename)
!         print *, "  Reason:", trim(errmsg)
!         print *, "  ==> Check if directory exists or if permission is denied."
!         return
!     end if

!     !----------------------------------------
!     ! Write Tecplot header
!     !----------------------------------------
!     write(unit,*) 'TITLE="Primitive Fields"'
!     write(unit,*) 'VARIABLES="X","Y","Density","Pressure","U","V"'
!     write(unit,'("ZONE I=",I5,", J=",I5,", F=POINT")') &
!          (i_hi - i_lo + 1), (j_hi - j_lo + 1)

!     !----------------------------------------
!     ! Export data
!     !----------------------------------------
!     do j = j_lo, j_hi
!         do i = i_lo, i_hi
!             write(unit,'(6E20.10)')  &
!                 mesh%x_center(i,j), mesh%y_center(i,j), &
!                 prim%density(i,j), prim%pressure(i,j), &
!                 prim%velocity_u(i,j), prim%velocity_v(i,j)
!         end do
!     end do

!     close(unit)
! end subroutine write_primitive_tecplot





! subroutine compute_global_timestep(mesh, prim, gas)
!     type(PolarMesh),          intent(in)  :: mesh
!     type(PrimitiveVariables), intent(in)  :: prim
!     type(StiffenedGas),       intent(in)  :: gas

!     ! --- local variables ---
!     integer :: i, j
!     integer :: i_lo, i_hi, j_lo, j_hi
!     double precision :: rho, p, u, v
!     double precision :: vmag, c_sound, a_char
!     double precision :: cell_area, h_eff
!     double precision :: CFL, dt_cfl, dt_cell
!     double precision, parameter :: tiny_speed = 1.0d-14

!     ! physical index range (without ghost cells)
!     i_lo = mesh%params%i_lo_phys
!     i_hi = mesh%params%i_hi_phys
!     j_lo = mesh%params%j_lo_phys
!     j_hi = mesh%params%j_hi_phys

!     ! initialize CFL and dt_cfl
!     CFL    = mesh%time%CFL
!     dt_cfl = mesh%time%dt_max


!     ! loop over physical cells
!     do j = j_lo, j_hi
!         do i = i_lo, i_hi

!             rho = prim%density(i,j)
!             p   = prim%pressure(i,j)
!             u   = prim%velocity_u(i,j)
!             v   = prim%velocity_v(i,j)

!             ! velocity magnitude (Cartesian)
!             vmag = sqrt(u*u + v*v)

!             ! local sound speed
!             c_sound = sound_speed_stiffened(rho, p, gas)

!             ! characteristic speed scale
!             a_char = vmag + c_sound

!             if (a_char <= tiny_speed) cycle   ! static cell, skip contribution

!             ! cell "volume": in 2D we use cell_area as volume per unit depth
!             cell_area = mesh%cell_area(i,j)

!             ! effective length scale: cbrt(V)
!             if (cell_area > 0.0d0) then
!                 h_eff = sqrt(cell_vol)
!             else
!                 cycle
!             end if

!             ! local time step
!             dt_cell = CFL * h_eff / a_char
!             if (dt_cell < dt_cfl) dt_cfl = dt_cell
!         end do
!     end do
!     mesh%time%dt = max(mesh%time%dt_min, dt_cfl)

! end subroutine compute_global_timestep
