module ScalableMatrix
  use Globals
  use ReadableMatrix
  implicit none
  private 
  public :: scalable_mat, scale_by_diagonal ,scale_system, scale_back
  type, abstract, extends(readable_mat) :: scalable_mat
     ! empty type
   contains
     !>-----------------------------------------------------------
     !> Procedure to scale a matrix M with a diagonal matrix D obtaing
     !>                M(out) = D M D
     !> Procedure restriced to square matrix
     !>-----------------------------------------------------------
     procedure(scale_by_diagonal), deferred :: diagonal_scale
  end type scalable_mat
  abstract interface
     subroutine scale_by_diagonal(this,lun_err,diagonal)
       use Globals
       import scalable_mat
       class(scalable_mat),  intent(inout) :: this
       integer,           intent(in   ) :: lun_err
       real(kind=double), intent(in   ) :: diagonal(this%ncol)
     end subroutine scale_by_diagonal
  end interface
contains
  !>-----------------------------------------------------------------
    !> Procedure to scale linear system 
    !>             $A x = b $
    !> trasforming it into
    !>             $M y = c $
    !> with 
    !>    $ M = D(A) A D(A)              $
    !>    $ y = D(A)^{-1}  x             $
    !>    $ c = D(A)       b             $
    !>----------------------------------------------------------------
    subroutine scale_system(lun_err,matrix, rhs, sol, diagonal)
      use Globals
      use Matrix
      integer,             intent(in   ) :: lun_err
      class(scalable_mat),    intent(inout) :: matrix
      real(kind=double),   intent(inout) :: rhs(matrix%nrow)
      real(kind=double),   intent(inout) :: sol(matrix%ncol)
      real(kind=double),   intent(in   ) :: diagonal(matrix%ncol)
            
   
      ! Diag(M)^{-1/2} M Diag(M)^{-1/2}
      call matrix%diagonal_scale(lun_err,diagonal)
      
      ! sol = Diag(M)^{1/2} sol 
      sol = sol / diagonal
      ! rhs = Diag(M)^{-1/2} rhs
      rhs = rhs * diagonal

    end subroutine scale_system

    !>-----------------------------------------------------------------
    !> Procedure to compute the solution x linear system 
    !>             $A x = b $
    !> after a diagonal scaling
    !>----------------------------------------------------------------
    subroutine scale_back(ncol, sol, inv_sqrt_diagonal)
      use Globals
      use Matrix
      integer,           intent(in   ) :: ncol
      real(kind=double), intent(inout) :: sol(ncol)
      real(kind=double), intent(inout) :: inv_sqrt_diagonal(ncol)

      sol = sol * inv_sqrt_diagonal

    end subroutine scale_back

end module ScalableMatrix
