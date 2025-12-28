!===========================================
! File   : flux_types.F90
! Module : flux_types
! Purpose:
!   - Define storage for numerical fluxes
!     on radial and angular faces.
!===========================================
module flux_types
    use mesh_types, only : PolarMesh
    implicit none
    private

    public :: FaceFlux2D
    public :: FluxFields
    public :: allocate_flux_fields
    public :: deallocate_flux_fields

    !---------------------------------------
    ! Flux on a 2D face family
    !---------------------------------------
    type :: FaceFlux2D
        real(kind=8), allocatable :: mass   (:,:)   ! ρ   flux
        real(kind=8), allocatable :: mom_u  (:,:)   ! ρu  flux
        real(kind=8), allocatable :: mom_v  (:,:)   ! ρv  flux
        real(kind=8), allocatable :: energy (:,:)   ! ρE  flux
    end type FaceFlux2D

    !---------------------------------------
    ! Fluxes on radial and angular faces
    !---------------------------------------
    type :: FluxFields
        type(FaceFlux2D) :: radial   ! faces normal to e_r
        type(FaceFlux2D) :: angular  ! faces normal to e_θ
    end type FluxFields

contains

    !===============================================================
    ! Allocate flux arrays according to mesh physical index ranges.
    !
    !  Radial faces  (F_r): i = i_lo_phys ... i_hi_phys+1, j = j_lo_phys ... j_hi_phys
    !  Angular faces (F_θ): i = i_lo_phys ... i_hi_phys,   j = j_lo_phys ... j_hi_phys+1
    !===============================================================
    subroutine allocate_flux_fields(flux, mesh)
        type(FluxFields), intent(inout) :: flux
        type(PolarMesh),  intent(in)    :: mesh

        integer :: iL, iR, jL, jR

        iL = mesh%params%i_lo_phys
        iR = mesh%params%i_hi_phys
        jL = mesh%params%j_lo_phys
        jR = mesh%params%j_hi_phys

        call deallocate_flux_fields(flux)
        !---- radial faces: size (Nr+1, Nθ)
        allocate(flux%radial%mass  (iL:iR+1, jL:jR))
        allocate(flux%radial%mom_u (iL:iR+1, jL:jR))
        allocate(flux%radial%mom_v (iL:iR+1, jL:jR))
        allocate(flux%radial%energy(iL:iR+1, jL:jR))

        !---- angular faces: size (Nr, Nθ+1)
        allocate(flux%angular%mass  (iL:iR, jL:jR+1))
        allocate(flux%angular%mom_u (iL:iR, jL:jR+1))
        allocate(flux%angular%mom_v (iL:iR, jL:jR+1))
        allocate(flux%angular%energy(iL:iR, jL:jR+1))
    end subroutine allocate_flux_fields


    subroutine deallocate_flux_fields(flux)
        type(FluxFields), intent(inout) :: flux

        if (allocated(flux%radial%mass  )) deallocate(flux%radial%mass)
        if (allocated(flux%radial%mom_u )) deallocate(flux%radial%mom_u)
        if (allocated(flux%radial%mom_v )) deallocate(flux%radial%mom_v)
        if (allocated(flux%radial%energy)) deallocate(flux%radial%energy)

        if (allocated(flux%angular%mass  )) deallocate(flux%angular%mass)
        if (allocated(flux%angular%mom_u )) deallocate(flux%angular%mom_u)
        if (allocated(flux%angular%mom_v )) deallocate(flux%angular%mom_v)
        if (allocated(flux%angular%energy)) deallocate(flux%angular%energy)
    end subroutine deallocate_flux_fields

end module flux_types
