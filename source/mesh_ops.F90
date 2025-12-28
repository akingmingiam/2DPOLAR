module mesh_ops
    use mesh_types,     only : pi, PolarMesh, TimeControlParameters, BC, EQ, RECON, RIEMANN, allocate_mesh_arrays
    use flow_fields,    only : PrimitiveVariables, ConservedVariables, allocate_primitive_fields, allocate_conserved_fields
    implicit none
contains

!=======================================================================
!  Subroutine: mesh_setup
!  Purpose   : Construct a uniform polar mesh with guard cells.
!              This includes computing cell centers, cell sizes
!              (Δr, Δθ), and cell areas in polar coordinates.
!=======================================================================
subroutine mesh_setup(mesh, time, prim, cons)
    type(PolarMesh), intent(out)             :: mesh
    type(TimeControlParameters), intent(out) :: time
    type(PrimitiveVariables),intent(out)     :: prim
    type(ConservedVariables),intent(out)     :: cons

    integer :: n_r, n_theta, g
    integer :: i, j
    double precision :: dr, dtheta, dt_out

    !---------------------------------------------------------------
    ! Assign mesh parameters (fixed, not modifiable externally)
    !---------------------------------------------------------------
    mesh%params%n_r         = 100
    mesh%params%n_theta     = 200
    mesh%params%guard_cells = 2
    mesh%params%r_min       = 0.0d0
    mesh%params%r_max       = 2.0d0
    mesh%params%theta_min   = 0.0d0
    mesh%params%theta_max   = 2.0d0 * pi

    !---------------------------------------------------------------
    ! Assign boundary condition types & equation type & reconstrucion
    !---------------------------------------------------------------
    mesh%params%equation_type  = EQ%CARTESIAN_TRANSFORM
    mesh%params%reconstruction = RECON%MUSCL_MINMOD
    mesh%params%riemann_solver = RIEMANN%HLLC
    mesh%params%bc_r_lo        = BC%EXTRAPOL   
    mesh%params%bc_r_hi        = BC%EXTRAPOL       
    mesh%params%bc_theta_lo    = BC%PERIODIC
    mesh%params%bc_theta_hi    = BC%PERIODIC

    !---------------------------------------------------------------
    ! Time control
    !---------------------------------------------------------------
    time%CFL               = 0.4d0
    time%t                 = 0.0d0
    time%t_end             = 4.0d-3
    time%dt                = 1.0d0  
    time%dt_max            = 1.0d0     
    time%dt_min            = 1.d-14 
    time%step              = 0
    time%max_steps         = 10000000
    time%output_step       = 0
    time%max_output_steps  = 10
    if (allocated(time%time_outputs)) deallocate(time%time_outputs)
    allocate(time%time_outputs(0:time%max_output_steps))
    dt_out = (time%t_end - time%t) / dble(time%max_output_steps)
    do i = 0, time%max_output_steps
        time%time_outputs(i) = time%t + i * dt_out
    end do

    !---------------------------------------------------------------
    ! Derive index ranges (for physical / computational domains)
    !---------------------------------------------------------------
    mesh%params%i_lo_phys = 1
    mesh%params%i_hi_phys = mesh%params%n_r
    mesh%params%j_lo_phys = 1
    mesh%params%j_hi_phys = mesh%params%n_theta

    mesh%params%i_lo_comp = 1 - mesh%params%guard_cells
    mesh%params%i_hi_comp = mesh%params%n_r     + mesh%params%guard_cells
    mesh%params%j_lo_comp = 1 - mesh%params%guard_cells
    mesh%params%j_hi_comp = mesh%params%n_theta + mesh%params%guard_cells

    !---------------------------------------------------------------
    ! Local shorthand for readability
    !---------------------------------------------------------------
    n_r     = mesh%params%n_r
    n_theta = mesh%params%n_theta
    g       = mesh%params%guard_cells

    !---------------------------------------------------------------
    ! Allocate mesh arrays on extended domain [1−g : n+g]
    !---------------------------------------------------------------
    call allocate_mesh_arrays(mesh)

    !---------------------------------------------------------------
    ! Uniform spacing
    !---------------------------------------------------------------
    dr     = (mesh%params%r_max - mesh%params%r_min)         / dble(n_r)
    dtheta = (mesh%params%theta_max - mesh%params%theta_min) / dble(n_theta)

    !---------------------------------------------------------------
    ! Calculate geometric parameters
    !---------------------------------------------------------------
    do i = 1-g, n_r+g
        mesh%radial_spacing(i) = dr
        mesh%r_center(i) = mesh%params%r_min + dr * dble(i - 0.5d0)
        mesh%r_face(i) = mesh%r_center(i) - mesh%radial_spacing(i) / 2.d0
        mesh%r_centroid(i) = mesh%r_center(i) + mesh%radial_spacing(i) / (12.d0 * mesh%r_center(i))
    end do

    do j = 1-g, n_theta+g
        mesh%angular_spacing(j) = dtheta
        mesh%theta_center(j) = mesh%params%theta_min + dtheta * dble(j - 0.5d0)
        mesh%theta_face(j) = mesh%theta_center(j) - mesh%angular_spacing(j) / 2.d0
        mesh%theta_centroid(j) = mesh%theta_center(j)
    end do

    do j = 1-g, n_theta+g
        do i = 1-g, n_r+g
            mesh%cell_area(i,j) = mesh%r_center(i) * dr * dtheta
            mesh%x_center(i,j) = mesh%r_center(i) * cos(mesh%theta_center(j))
            mesh%y_center(i,j) = mesh%r_center(i) * sin(mesh%theta_center(j))

            mesh%x_ll(i,j) = mesh%r_face(i) * cos(mesh%theta_face(j))
            mesh%y_ll(i,j) = mesh%r_face(i) * sin(mesh%theta_face(j))

            mesh%radial_face_length(i,j) = mesh%r_face(i) * mesh%angular_spacing(j)
            mesh%angular_face_length(i,j) = mesh%radial_spacing(i)
        end do
    end do
               

    !---------------------------------------------------------------
    ! Allocate PRIMITIVE variables (ρ, p, u, v), CONSERVED variables (ρ, ρu, ρv, ρE)
    !---------------------------------------------------------------
    call allocate_primitive_fields( prim, 1-g, n_r+g, 1-g, n_theta+g )
    call allocate_conserved_fields( cons, 1-g, n_r+g, 1-g, n_theta+g )

    ! Print summary
    call print_mesh_summary(mesh)
end subroutine mesh_setup


!=======================================================================
!  Subroutine: print_mesh_summary
!  Purpose   : Print formatted information about the mesh.
!=======================================================================
subroutine print_mesh_summary(mesh)
    type(PolarMesh), intent(in) :: mesh
    double precision :: dr, dtheta

    dr     = (mesh%params%r_max - mesh%params%r_min) / dble(mesh%params%n_r)
    dtheta = (mesh%params%theta_max - mesh%params%theta_min) / dble(mesh%params%n_theta)

    write(*,*)
    write(*,'(a)') '--- Polar Mesh Configuration ---'

    write(*,'(a,1x,i4,1x,a,1x,i4)') '  Cells           : ', &
                                    mesh%params%n_r, 'x', mesh%params%n_theta

    write(*,'(a,1x,f12.6,1x,a,1x,f12.6)') '  r-range         :', &
                                          mesh%params%r_min, '→', mesh%params%r_max

    write(*,'(a,1x,f12.6,1x,a,1x,f12.6)') '  theta-range     :', &
                                          mesh%params%theta_min, '→', mesh%params%theta_max

    write(*,'(a,1x,i5)') '  Guard cells     :', mesh%params%guard_cells

    write(*,'(a,1x,f12.6)') '  Δr (uniform)    :', dr
    write(*,'(a,1x,f12.6)') '  Δθ (uniform)    :', dtheta

    write(*,'(a)') '--------------------------------'
    write(*,*)
end subroutine print_mesh_summary


end module mesh_ops
