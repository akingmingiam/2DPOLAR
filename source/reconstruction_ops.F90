module reconstruction_ops
    use mesh_types, only : RECON

    implicit none
    private

    public :: muscl_reconstruct

contains

    !----------------------------------------------------------
    !  MUSCL reconstruction with minmod limiter (uniform grid)
    !
    !  Input:
    !    q_stencil(1:4) =
    !        q_stencil(1) = q_{i-2}
    !        q_stencil(2) = q_{i-1}
    !        q_stencil(3) = q_{i}
    !        q_stencil(4) = q_{i+1}
    !
    !  Output:
    !    qL : left state at interface i-1/2, from cell (i-1)
    !    qR : right state at interface i-1/2, from cell (i)
    !----------------------------------------------------------
    subroutine muscl_reconstruct(q_stencil, qL, qR, limiter_type)
        double precision, intent(in)  :: q_stencil(4)
          integer,        intent(in)  :: limiter_type
        double precision, intent(out) :: qL, qR

        double precision             :: dq_m, dq_p  ! slopes
        double precision             :: r, phi      ! limited
        double precision, parameter  :: eps = 1.0d-16

        ! Slopes around cell i-1:
        dq_m = q_stencil(2) - q_stencil(1)
        dq_p = q_stencil(3) - q_stencil(2)
        r    = dq_m / (dq_p + eps)
        phi = compute_limiter(r, limiter_type)
        qL = q_stencil(2) + 0.5d0 * phi * dq_p

        ! Slopes around cell i:
        !   dq_m = q_i       - q_{i-1}
        !   dq_p = q_{i+1}   - q_i
        dq_m = q_stencil(3) - q_stencil(2)
        dq_p = q_stencil(4) - q_stencil(3)
        r    = dq_m / (dq_p + eps)
        phi = compute_limiter(r, limiter_type)
        qR = q_stencil(3) - 0.5d0 * phi * dq_p

    end subroutine muscl_reconstruct


    function compute_limiter(r, limiter_type) result(phi)
        double precision, intent(in) :: r
        integer, intent(in)          :: limiter_type
        double precision             :: phi

        select case (limiter_type)
        case (RECON%MUSCL_MINMOD)
            phi = minmod(1.0d0, r)

        case (RECON%MUSCL_VANLEER)
            phi = vanleer_limiter(r)

        case default
            stop "Unknown reconstruction type!"
        end select

    end function compute_limiter


    !----------------------------------------------------------
    !  Minmod limiter: mm = minmod(a, b)
    !
    !  minmod(a,b) =
    !     0                      if a*b <= 0
    !     sign(a)*min(|a|,|b|)   if a*b  > 0
    !----------------------------------------------------------
    pure function minmod(a, b) result(mm)
        double precision, intent(in) :: a, b
        double precision             :: mm

        if (a*b <= 0.0d0) then
            mm = 0.0d0
        else
            if (abs(a) < abs(b)) then
                mm = a
            else
                mm = b
            end if
        end if
    end function minmod

     !----------------------------------------------------------
    !  Van Leer limiter:
    !
    !      phi(r) = (r + |r|) / (1 + |r|),   r > 0
    !             = 0,                       r <= 0
    !----------------------------------------------------------
    pure function vanleer_limiter(r) result(phi)
        double precision, intent(in) :: r
        double precision             :: phi
        double precision, parameter  :: eps = 1.0d-16

        if (r <= 0.0d0) then
            phi = 0.0d0
        else
            phi = (r + abs(r)) / (1.0d0 + abs(r) + eps)
        end if
    end function vanleer_limiter

end module reconstruction_ops
