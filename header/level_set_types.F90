module level_set_types
    use mesh_types, only : PolarMesh
    implicit none
    private

    public :: LevelSetFields
    public :: allocate_level_set_fields
    public :: deallocate_level_set_fields

    !======================================================================
    !  Module: level_set_types
    !  Purpose:
    !      Store cell-centered data required by a level-set solver,
    !      including the signed-distance function and derived geometric
    !      quantities used for interface capturing and reinitialization.
    !======================================================================

    !------------------------------------------------------------------
    !  Type: LevelSetFields
    !  Cell-centered level-set quantities defined on the computational
    !  domain, including ghost cells.
    !------------------------------------------------------------------
    type :: LevelSetFields
        real(kind=8), allocatable :: phi           (:,:)  ! signed distance function φ
        real(kind=8), allocatable :: gradient_r    (:,:)  ! ∂φ/∂r (radial derivative)
        real(kind=8), allocatable :: gradient_theta(:,:)  ! ∂φ/∂θ (angular derivative)
        real(kind=8), allocatable :: curvature     (:,:)  ! κ = ∇·(∇φ / |∇φ|)
    end type LevelSetFields

contains

    !------------------------------------------------------------------
    !  Allocate level-set arrays on computational domain indices.
    !------------------------------------------------------------------
    subroutine allocate_level_set_fields(ls, mesh)
        type(LevelSetFields), intent(inout) :: ls
        type(PolarMesh),      intent(in)    :: mesh
        integer :: iL, iR, jL, jR

        iL = mesh%params%i_lo_comp
        iR = mesh%params%i_hi_comp
        jL = mesh%params%j_lo_comp
        jR = mesh%params%j_hi_comp

        call deallocate_level_set_fields(ls)

        allocate(ls%phi           (iL:iR, jL:jR))
        allocate(ls%gradient_r    (iL:iR, jL:jR))
        allocate(ls%gradient_theta(iL:iR, jL:jR))
        allocate(ls%curvature     (iL:iR, jL:jR))
    end subroutine allocate_level_set_fields

    !------------------------------------------------------------------
    !  Deallocate level-set arrays
    !------------------------------------------------------------------
    subroutine deallocate_level_set_fields(ls)
        type(LevelSetFields), intent(inout) :: ls

        if (allocated(ls%phi           )) deallocate(ls%phi)
        if (allocated(ls%gradient_r    )) deallocate(ls%gradient_r)
        if (allocated(ls%gradient_theta)) deallocate(ls%gradient_theta)
        if (allocated(ls%curvature     )) deallocate(ls%curvature)
    end subroutine deallocate_level_set_fields

end module level_set_types
