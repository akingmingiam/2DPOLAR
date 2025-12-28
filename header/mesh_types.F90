module mesh_types
    implicit none
    private

    public :: pi
    public :: PolarGridParameters
    public :: TimeControlParameters
    public :: PolarMesh
    public :: allocate_mesh_arrays
    public :: BC, EQ, RECON, RIEMANN

    !=======================================================================
    !  Mathematical constants
    !=======================================================================
    ! real(kind=8), parameter :: pi = 3.1415926535897932384626433832795d0
    real(kind=8), parameter :: pi = acos(-1.0d0)

    !=======================================================================
    !  Boundary condition, Equation type, Reconstruction
    !=======================================================================
    type :: BoundaryConditionTypes
        integer :: SYMMETRY  = 1   ! reflective wall
        integer :: EXTRAPOL  = 2   ! zero-gradient
        integer :: PERIODIC  = 3   ! periodic
        integer :: AXIS      = 4   ! r = 0 axis
    end type BoundaryConditionTypes

    type :: EquationSystemTypes
        integer :: CARTESIAN_TRANSFORM = 1   ! Cartesian Euler + mapping (no geometric sources)
        integer :: POLAR_EULER         = 2   ! Polar Euler equations (explicit geometric sources)
    end type EquationSystemTypes

    type :: ReconstructionTypes
        integer :: FIRST_ORDER    = 1   ! piecewise constant: q_L = q_{i-1}, q_R = q_i
        integer :: MUSCL_MINMOD   = 2   ! MUSCL reconstruction with minmod limiter
        integer :: MUSCL_VANLEER  = 3   ! MUSCL reconstruction with van Leer limiter
    end type ReconstructionTypes

    type :: RiemannSolverTypes
        integer :: HLL_SIMPLE  = 1   ! HLL solver 
        integer :: HLL         = 2   ! HLL solver (Roe)
        integer :: HLLC        = 3   ! HLLC solver (Roe)
    end type RiemannSolverTypes


    type(BoundaryConditionTypes), parameter :: BC = BoundaryConditionTypes()
    type(EquationSystemTypes),    parameter :: EQ = EquationSystemTypes()
    type(ReconstructionTypes),    parameter :: RECON = ReconstructionTypes()
    type(RiemannSolverTypes),     parameter :: RIEMANN = RiemannSolverTypes()
     
    !=======================================================================
    !  Type: PolarGridParameters
    !  Purpose:
    !      Store configuration values defining the computational mesh:
    !      resolution (Nr, Nθ), ghost layers, radial/azimuthal ranges.
    !=======================================================================
    type :: PolarGridParameters

        ! Equation type
        integer :: equation_type

        ! Reconstruction
        integer :: reconstruction

        ! Riemann Solver
        integer :: riemann_solver

        ! Grid resolution
        integer :: n_r          ! Number of radial cells (physical domain)
        integer :: n_theta      ! Number of angular cells (physical domain)
        integer :: guard_cells  ! Number of ghost cells on each side

        ! Physical domain extents
        double precision :: r_min      ! Inner radius  (m)
        double precision :: r_max      ! Outer radius  (m)
        double precision :: theta_min  ! Minimum angle (rad)
        double precision :: theta_max  ! Maximum angle (rad)

        ! Boundary condition types
        integer :: bc_r_lo      ! Type of boundary at r = r_min
        integer :: bc_r_hi      ! Type of boundary at r = r_max
        integer :: bc_theta_lo  ! Type of boundary at theta = theta_min
        integer :: bc_theta_hi  ! Type of boundary at theta = theta_max

        ! Index ranges (for convenience)
        ! Physical domain indices
        integer :: i_lo_phys, i_hi_phys
        integer :: j_lo_phys, j_hi_phys

        ! Computational domain indices (including ghost cells)
        integer :: i_lo_comp, i_hi_comp
        integer :: j_lo_comp, j_hi_comp
    end type PolarGridParameters

    !=======================================================================
    !  Type: TimeControlParameters
    !=======================================================================
    type :: TimeControlParameters
        double precision :: CFL                            ! Courant-Friedrichs-Lewy number
        double precision :: t                              ! simulation start time
        double precision :: t_end                          ! simulation end time
        double precision :: dt                             ! time
        double precision :: dt_max                         ! maximum timestep
        double precision :: dt_min                         ! minimum timestep

        integer :: step                                    ! starting step index
        integer :: max_steps                               ! hard limit on step count
        integer :: output_step                             ! output step
        integer :: max_output_steps                        ! max number of output frames
        double precision, allocatable :: time_outputs(:)   ! output times (user-defined)
    end type TimeControlParameters




    !=======================================================================
    !  Type: PolarMesh
    !  Indexing convention:
    !      Arrays are allocated on [1−g : N+g] to include ghost cells.
    !
    !      Physical domain:
    !          i = 1 ... n_r
    !          j = 1 ... n_theta
    !
    !      Ghost layers:
    !          i = 1-g ... 0        and     n_r+1 ... n_r+g
    !          j = 1-g ... 0        and     n_theta+1 ... n_theta+g
    !
    !=======================================================================
    type :: PolarMesh

        ! Grid configuration (resolution and geometric limits)
        type(PolarGridParameters)     :: params
        double precision, allocatable :: r_center(:)              ! radial centers
        double precision, allocatable :: theta_center(:)          ! angular centers
        double precision, allocatable :: r_centroid(:)            ! cell centroid radial coordinate
        double precision, allocatable :: theta_centroid(:)        ! cell centroid angular coordinate
        double precision, allocatable :: r_face(:)                ! radial face locations for i - 1/2
        double precision, allocatable :: theta_face(:)            ! angular face locations for j - 1/2
        double precision, allocatable :: radial_spacing(:)        ! Δr
        double precision, allocatable :: angular_spacing(:)       ! Δθ
        double precision, allocatable :: cell_area(:,:)           ! cell volumes/areas
        double precision, allocatable :: x_center(:,:)            ! Cartesian coordinates of cell centers
        double precision, allocatable :: y_center(:,:)            ! Cartesian coordinates of cell centers
        double precision, allocatable :: radial_face_length(:,:)  ! length of radial faces  rΔθ for (i - 1/2, j)
        double precision, allocatable :: angular_face_length(:,:) ! length of angular faces Δr  for (i, j - 1/2)
        double precision, allocatable :: x_ll(:,:)                ! left-lower half index of x (i - 1/2, j - 1/2)
        double precision, allocatable :: y_ll(:,:)                ! left-lower half index of y (i - 1/2, j - 1/2)


    end type PolarMesh

contains
    !=======================================================================
    !  Subroutine: allocate_mesh_arrays + deallocate_mesh_arrays
    !=======================================================================
    subroutine allocate_mesh_arrays(mesh)
        type(PolarMesh), intent(inout) :: mesh
        integer :: n_r, n_theta, g

        n_r     = mesh%params%n_r
        n_theta = mesh%params%n_theta
        g       = mesh%params%guard_cells

        call deallocate_mesh_arrays(mesh)

        allocate(mesh%r_center(1-g:n_r+g))
        allocate(mesh%theta_center(1-g:n_theta+g))
        allocate(mesh%r_centroid(1-g:n_r+g))
        allocate(mesh%theta_centroid(1-g:n_theta+g))
        allocate(mesh%r_face(1-g:n_r+g))
        allocate(mesh%theta_face(1-g:n_theta+g))
        allocate(mesh%radial_spacing(1-g:n_r+g))
        allocate(mesh%angular_spacing(1-g:n_theta+g))
        allocate(mesh%cell_area(1-g:n_r+g, 1-g:n_theta+g))
        allocate(mesh%x_center(1-g:n_r+g, 1-g:n_theta+g))
        allocate(mesh%y_center(1-g:n_r+g, 1-g:n_theta+g))
        allocate(mesh%radial_face_length(1-g:n_r+g, 1-g:n_theta+g))
        allocate(mesh%angular_face_length(1-g:n_r+g, 1-g:n_theta+g))
        allocate(mesh%x_ll(1-g:n_r+g, 1-g:n_theta+g))
        allocate(mesh%y_ll(1-g:n_r+g, 1-g:n_theta+g))

    end subroutine allocate_mesh_arrays

    subroutine deallocate_mesh_arrays(mesh)
        type(PolarMesh), intent(inout) :: mesh

        if (allocated(mesh%r_center)) deallocate(mesh%r_center)
        if (allocated(mesh%theta_center)) deallocate(mesh%theta_center)
        if (allocated(mesh%x_center)) deallocate(mesh%r_centroid)
        if (allocated(mesh%y_center)) deallocate(mesh%theta_centroid)
        if (allocated(mesh%r_face)) deallocate(mesh%r_face)
        if (allocated(mesh%theta_face)) deallocate(mesh%theta_face)
        if (allocated(mesh%radial_spacing)) deallocate(mesh%radial_spacing)
        if (allocated(mesh%angular_spacing)) deallocate(mesh%angular_spacing)
        if (allocated(mesh%cell_area)) deallocate(mesh%cell_area)
        if (allocated(mesh%x_center)) deallocate(mesh%x_center)
        if (allocated(mesh%y_center)) deallocate(mesh%y_center)
        if (allocated(mesh%radial_face_length)) deallocate(mesh%radial_face_length)
        if (allocated(mesh%angular_face_length)) deallocate(mesh%angular_face_length)
        if (allocated(mesh%x_ll)) deallocate(mesh%x_ll)
        if (allocated(mesh%y_ll)) deallocate(mesh%y_ll)
    end subroutine deallocate_mesh_arrays

end module mesh_types
