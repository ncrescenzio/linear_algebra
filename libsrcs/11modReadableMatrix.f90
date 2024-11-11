module ReadableMatrix
  use Globals
  use Matrix
  implicit none
  private 
  public :: readable_mat, extract_diagonal 
  type, abstract, extends(abs_matrix) :: readable_mat
     ! empty type
   contains
     !>-----------------------------------------------------------
     !> Procedure to extract the diagonal from a matrix
     procedure(extract_diagonal), deferred :: get_diagonal
  end type readable_mat
  abstract interface
     subroutine extract_diagonal(this,diagonal)
       use Globals
       import readable_mat
       class(readable_mat), intent(in)  :: this
       real(kind=double), intent(inout) :: diagonal(min(this%nrow,this%ncol))
     end subroutine extract_diagonal
  end interface
end module ReadableMatrix
