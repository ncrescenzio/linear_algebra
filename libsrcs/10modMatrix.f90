module Matrix
  use Globals
  use LinearOperator
  implicit none
  private
  !> Abstract type definition for matrix-type, which represents
  !> any linear operator for REAL^{NCOL} to REAL^{NROW}
  !> Additional label are included to described the matrix structure
  !> that will limit or active matrix operations. For example only
  !> symmetric matrices can be passed to PCG procedure, or can have 
  !> a Cholesky decomposition. While triangular matric can be easily solved.
  !> The Abstract inferface requires that any type extdending the 
  !> class abs_matrix must have a procedure defining the operation
  !> matrix-times-vector ( being extension of the abs_linop type)
  !> and matrix(transpose)-times-vector 
  public :: abs_matrix,multiplication_transpose, TransposeMatrix
  type, extends(abs_linop), abstract :: abs_matrix
   contains
     !> Publlic procedure to set the properties
     !> adjoint operator AT (nrow ncol etc) 
     !> for origianl operator A
     procedure, public, pass:: set_adj_properties
     !> Procedure for computation of
     !> (matrix-transpose) times (vector)
     !> Each type extending the abs_matrix type MUST
     !> contain this procedure with this name
     !> (it will not compile otherwise)
     procedure(multiplication_transpose), deferred :: matrix_transpose_times_vector
     !> Short-hand version of matrix_transpose_times_vector 
     !> without info and lun_err
     procedure, public, pass:: MTxv => mtxv_abs_matrix
  end type abs_matrix
  
  
  abstract interface
     !>-------------------------------------------------------------
     !> Abstract procedure defining the interface for a general
     !> matrix(transpose)-vector multiplication
     !>         vec_out = (M)^T times (vec_in)
     !> (public procedure for class abs_matrix)
     !> 
     !> usage:
     !>     call 'var'%MTxv(vec_in,vec_out,[info])
     !>
     !> where 
     !> \param[in   ] vec_in  -> real, dimension('var'%nrow)
     !>                          vector to be multiplied
     !> \param[inout] vec_out -> real, dimension('var'%ncol)
     !>                          vector (M) times (vec_in) 
     !> \param[in   ] info    -> integer. Info number
     !>                          in case of error   
     !<-------------------------------------------------------------
     subroutine multiplication_transpose(this,vec_in,vec_out,info,lun_err)
       use Globals
       import abs_matrix
       implicit none
       class(abs_matrix), intent(inout) :: this
       real(kind=double), intent(in   ) :: vec_in(this%nrow)
       real(kind=double), intent(inout) :: vec_out(this%ncol)
       integer,           intent(inout) :: info
       integer,           intent(in   ) :: lun_err

     end subroutine multiplication_transpose
  end interface

  !>------------------------------------------------------------------
  !> Derived type used to create list of pointers to abstract matrices
  !> (at the best of my knowledge Fortan does not allow to
  !>  create arrays of pointers)
  !> E.g. in a Block-matrix
  !> M = (A D) with A and D of different type, we can create
  !> an array type that points toward A and D.
  !> Declarations: 
  !>  type(spmat)     :: A
  !>  type(diagmat  ) :: D
  !>  type(array_mat) :: list(2)
  !> Assignment
  !>  list(1)%mat => A
  !>  list(2)%mat => D
  !> See BlockMatrix module for a concrete example
  type, public :: array_mat
     !> Dimension (nmats)
     !> Array that contains the non-zero blocks
     class(abs_matrix), pointer :: mat
  end type array_mat  

  !>----------------------------------------------------------------------
  !> Structure variable containg member storing sparse real matrices
  !> in csr ( compress sparse row ) and ssr (symmetric sparse row format )
  !>----------------------------------------------------------------------
  type, extends(abs_linop), public :: TransposeMatrix
     class(abs_matrix), pointer :: original
   contains
     !> Procedure to compute 
     !>         y = M * x 
     !> (public procedure for type spmat)
     procedure, public, pass :: matrix_times_vector => TransposeMatrixxv
     !> Procedure to compute 
     !>         y = M * x 
     !> (public procedure for type spmat)
     procedure, public, pass :: init => init_TransposeMatrix
  end type TransposeMatrix

contains

  !>-------------------------------------------------------------
  !> Procedure for matrix transpose -vector multiplication
  !>         vec_out = (M^T) times (vec_in)
  !> (public procedure for class abs_linop)
  !> It allows a simpler interface with matrix_times_vector
  !> procedure
  !>
  !> usage:
  !>     call 'var'%Mxv(vec_in,vec_out,[lun_err,info])
  !>
  !> where 
  !> \param[in   ] vec_in            -> real, dimension('var'%ncol)
  !>                                   vector to be multiplied
  !> \param[inout] vec_out           -> real, dimension('var'%nrow)
  !>                                    vector (M) times (vec_in) 
  !> \param[in   ] (optional) info    -> integer. Info number
  !>                                    in case of error 
  !> \param[in   ] (optional) lun_err -> integer. Info number
  !>                                    in case of error 
  !<-------------------------------------------------------------
  subroutine mtxv_abs_matrix(this,vec_in,vec_out,info,lun_err)
    implicit none
    class(abs_matrix),  intent(inout) :: this
    real(kind=double), intent(in   ) :: vec_in(this%nrow)
    real(kind=double), intent(inout) :: vec_out(this%ncol)
    integer, optional, intent(inout) :: info
    integer, optional, intent(in   ) :: lun_err
    ! local
    logical :: rc
    integer :: info_loc,lun_loc=0

    !
    ! select logical unit for errors
    !
    if ( present(lun_err) ) lun_loc = lun_err

    call this%matrix_transpose_times_vector(vec_in,vec_out,info_loc,lun_loc)

    !
    ! if any error occured:
    ! - print warning if info was passed
    ! - stop otherwise
    !
    if (info_loc .ne.0 ) then
       if ( present(info) ) then
          info=info_loc
          rc = IOerr(lun_loc, wrn_val, 'mtxv_abs_matrix', &
               ' Error in matrix transpose times vector multiplication')   
       else
          rc = IOerr(lun_loc, err_val, 'mtxv_abs_matrix' , &
               ' Error in matrix transpose times vector multiplication')   
       end if
    end if
  end subroutine mtxv_abs_matrix
  

  subroutine init_TransposeMatrix(this,matrix,info)
    implicit none
    class(TransposeMatrix), intent(inout) :: this
    class(abs_matrix),  target,   intent(in   ) :: matrix
    integer, optional,      intent(inout) :: info
    
    ! set general properties transpose operator
    call matrix%set_adj_properties(this)

    ! assing pointer to original matrix
    this%original => matrix

  end subroutine init_TransposeMatrix

  !>-------------------------------------------------------------
  !> Procedure to compute Matrix vector product
  !>         vec_out = (M) times (vec_in)
  !> (public procedure for type spmat)
  !> 
  !> usage:
  !>     call 'var'%Mxv(vec_in,vec_out,[info])
  !>
  !> where 
  !> \param[in   ] vec_in          -> real. dimension('var'%ncol)
  !>                                  vector to be multiplied
  !> \param[inout] vec_out         -> real. dimension('var'%nrow)
  !>                                  vector (M) times (vec_in) 
  !> \param[in   ] info -> integer. Info number
  !>                                  in case of error   
  !<-------------------------------------------------------------
  subroutine TransposeMatrixxv(this,vec_in,vec_out, info,lun_err)
    use Globals
    implicit none
    class(TransposeMatrix),    intent(inout) :: this
    real(kind=double), intent(in   ) :: vec_in(this%ncol)
    real(kind=double), intent(inout) :: vec_out(this%nrow)
    integer,           intent(inout) :: info
    integer,           intent(in   ) :: lun_err

    !
    ! use the matrix transpose x vector operator of the original matrix
    !
    call this%original%matrix_transpose_times_vector(vec_in,vec_out,info,lun_err)

  end subroutine TransposeMatrixxv

!!$  function TRANS(matrix) result(this) 
!!$    implicit none
!!$    class(abs_matrix), target,  intent(in   ) :: matrix
!!$    class(TransposeMatrix) :: this
!!$    ! local
!!$    logical :: rc
!!$    integer :: info_loc
!!$
!!$    call this%init(matrix,info_loc)
!!$    if ( info_loc .ne. 0)  rc = IOerr(0, wrn_val, 'TRANS', &
!!$         '  type t_spmat ',info_loc)
!!$
!!$
!!$  end function TRANS

  !>------------------------------------------------------
  !> Subruotine to define all proprieties of
  !> adjoint linear operator. Useful in the 
  !> initialiazionion of concrete type extending the
  !> adj_linop class.
  !> (public procedure for class adj_linop)
  !> 
  !> usage:
  !>     call 'var'%set_adj_properties(matrix)
  !<-------------------------------------------------------------
  subroutine set_adj_properties(matrix,this)
    implicit none
    class(abs_matrix),      intent(in   ) :: matrix
    class(TransposeMatrix), intent(inout) :: this
    

    this%nrow         = matrix%ncol
    this%ncol         = matrix%nrow
    this%is_symmetric = matrix%is_symmetric
    this%unitary_diag = matrix%unitary_diag
    if ( SCAN(matrix%name,'?').eq. 0 ) then
       this%name=etb('('//etb(this%name)//')^T')
    end if

    ! set properties of transpose matrix
    select case (matrix%triangular)
    case ('N')
       this%triangular = 'N'
    case ('U')
       this%triangular = 'L'
    case ('L')
       this%triangular = 'U'
    end select


  end subroutine set_adj_properties

  

  
end module Matrix




