module flow_fields
    implicit none
    private
    !======================================================================
    !  Module: flow_fields
    !  Purpose:
    !      Store primitive and conserved flow field variables for a
    !      2D CFD solver (cell-centered fields).
    !
    !      - PrimitiveVariables : (ρ, p, u, v)
    !      - ConservedVariables : (ρ, ρu, ρv, ρE)
    !======================================================================


    !======================================================================
    !  Type: PrimitiveVariables
    !  Cell-centered primitive variables:
    !======================================================================
    type, public :: PrimitiveVariables
        real(kind=8), allocatable :: density    (:,:)  ! ρ
        real(kind=8), allocatable :: pressure   (:,:)  ! p
        real(kind=8), allocatable :: velocity_u (:,:)  ! u-component Cartesian component
        real(kind=8), allocatable :: velocity_v (:,:)  ! v-component Cartesian component
    end type PrimitiveVariables


    !======================================================================
    !  Type: ConservedVariables
    !  Cell-centered conserved variables for finite-volume schemes
    !======================================================================
    type, public :: ConservedVariables
        real(kind=8), allocatable :: rho   (:,:)  ! mass density
        real(kind=8), allocatable :: rho_u (:,:)  ! momentum in u-direction
        real(kind=8), allocatable :: rho_v (:,:)  ! momentum in v-direction
        real(kind=8), allocatable :: rho_E (:,:)  ! total energy density
    end type ConservedVariables


    !======================================================================
    !  Public helper interfaces
    !======================================================================
    public :: allocate_primitive_fields
    public :: allocate_conserved_fields
    public :: deallocate_primitive_fields
    public :: deallocate_conserved_fields

contains

    !------------------------------------------------------------------
    !  Allocate variables on [i_lo:i_hi, j_lo:j_hi]
    !------------------------------------------------------------------
    subroutine allocate_primitive_fields(prim, i_lo, i_hi, j_lo, j_hi)
        type(PrimitiveVariables), intent(inout) :: prim
        integer, intent(in) :: i_lo, i_hi, j_lo, j_hi
       
        call deallocate_primitive_fields(prim)
        allocate(prim%density   (i_lo:i_hi, j_lo:j_hi))
        allocate(prim%pressure  (i_lo:i_hi, j_lo:j_hi))
        allocate(prim%velocity_u(i_lo:i_hi, j_lo:j_hi))
        allocate(prim%velocity_v(i_lo:i_hi, j_lo:j_hi))
    end subroutine allocate_primitive_fields

    subroutine allocate_conserved_fields(cons, i_lo, i_hi, j_lo, j_hi)
        type(ConservedVariables), intent(inout) :: cons
        integer, intent(in) :: i_lo, i_hi, j_lo, j_hi
        
        call deallocate_conserved_fields(cons)
        allocate(cons%rho  (i_lo:i_hi, j_lo:j_hi))
        allocate(cons%rho_u(i_lo:i_hi, j_lo:j_hi))
        allocate(cons%rho_v(i_lo:i_hi, j_lo:j_hi))
        allocate(cons%rho_E(i_lo:i_hi, j_lo:j_hi))
    end subroutine allocate_conserved_fields


    !------------------------------------------------------------------
    !  Deallocate variables
    !------------------------------------------------------------------
    subroutine deallocate_primitive_fields(prim)
        type(PrimitiveVariables), intent(inout) :: prim

        if (allocated(prim%density   )) deallocate(prim%density)
        if (allocated(prim%pressure  )) deallocate(prim%pressure)
        if (allocated(prim%velocity_u)) deallocate(prim%velocity_u)
        if (allocated(prim%velocity_v)) deallocate(prim%velocity_v)
    end subroutine deallocate_primitive_fields

    subroutine deallocate_conserved_fields(cons)
        type(ConservedVariables), intent(inout) :: cons

        if (allocated(cons%rho  )) deallocate(cons%rho)
        if (allocated(cons%rho_u)) deallocate(cons%rho_u)
        if (allocated(cons%rho_v)) deallocate(cons%rho_v)
        if (allocated(cons%rho_E)) deallocate(cons%rho_E)
    end subroutine deallocate_conserved_fields


end module flow_fields
