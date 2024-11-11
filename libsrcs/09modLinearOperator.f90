!> This module contains the structure varibles for 
!> abstract definition of linear operator (matrix)
!> and some simple extension, like linear combination
!> M=\sum_i a_i M_i
!> or composition
!> M = M1 * M2 * ...* Mn
!> or an combination of the two, like
!> M = A -2 * B * C
!> Abstract type contains minimal informations to describe
!> a linear operator :: domain and codomain dimension, 
!> symmetry or other properties.
!> The most important part of the module is the definition of the
!> abstact procedure for linear operator application that computes
!>
!> vec_out = linear_operator(vec_in)
!>
module LinearOperator
  use Globals
  use Scratch
  implicit none
  !> Few globlas parameters
  !> Set to .True. to trigger debugging at all level
  logical, parameter :: full_debug=.False.
  !> Maximum lenght of matrix name
  integer, parameter :: max_length_name=70
  private
  public :: abs_linop,multiplication,info_abs_linop,array_linop
  !> Abstract type definition for linop-type, which
  !> represents any linear operator from REAL^{NCOL} to
  !> REAL^{NROW}.
  !> Additional label are included to described
  !> the linerar operator structure that will limit or
  !> active matrix operations (only symmetric matrices can
  !> be passed to PCG procedure or can possibly have a
  !> Cholesky decomposition). While triangular matrices can be
  !> easily solved.
  !> The Abstract inferface requires that any
  !> type extdending the class abs_linop must have a
  !> procedure defining the operation matrix-times-vector
  !>
  !> ISSUES 1: "replace intent(in) or intent(inout)?"
  !>         Now is intent(inout)
  !> Note that we set intent(inout) for abs_linop in the
  !> muliptlication procedure while, in theory, application
  !> of linear operator should not affects its components.
  !> This is done due to the fact that scratch arrays
  !> included in the concrete type extendining abstract type
  !> linop may be used. For example matrix-vector operation
  !> $M \times v$ with $M =A^T A$ a scratch array is
  !> required. This can be removed passing an optional
  !> scratch arrays properly sized.
  !>
  !> PRO   : intent(in) used, avoid modification 
  !> CONTRA: not self-contained, each time you need to pass
  !> or create scratch array with auxilary (type scrt in
  !> module Scratch may help) 
  !>
  !> ISSUES 2: "
  !> "Use:
  !>
  !>  vec_out = alpha*vec_out + beta *M*vec_in
  !>
  !>  instead of
  !>
  !>  vec_out = M * Vec_in"
  !>
  !> The current impletation return  vec_out = M *vec_in.
  !> However, in practically all concrete types we built,
  !> we first initialize vec_out to 
  !> zero and that make all computations. If insteand we add 
  !> an new arguments, possibly optional, that set vec_out to zero
  !> (letting the user to define it) we could compute
  !> vec_out = vec_out +M *vec_in with no additional effort.
  !> It will be a sort of daxpy procedure.
  !> PRO   : save computation time, 
  !>         save storage requirement, for example for linear combination
  !> CONTRA: more arguments, less clean procedure
  !>
  !> POSSIBLE SOLUTION :  Use
  !> vec_out = alpha*vec_out + beta *M*vec_in
  !> in a subroutine (say aypbMx), the only use we really have
  !> to write. No optional arguments would be used
  !> (it is always painfull to handle the presence or not
  !> of optional arguments info and lun_err).
  !> We introduce a new subroutine
  !> computing only the matrix vector productor y=M*x.
  !> This would be done one for all in the abs_linop
  type, abstract :: abs_linop
     !> Number of rows
     integer :: nrow=0
     !> Number of columns
     integer :: ncol=0
     !> Logical flag for symmetry
     logical :: is_symmetric=.false.
     !> Character indicating if the matrix is 
     !> lower ('L')  or upper ('U') triangular or
     !> it is not ('N'), which is the default case
     character(len=1) :: triangular='N'
     !> Logical flag matrix with one on the diagonal
     ! TODO: Used only in SparseMatrix Module remove?
     logical :: unitary_diag=.false.
     !> Matrix Name
     character(len=70) :: name='?'
     !> Real estimating the linearity of the
     !> operator i.e.
     !> linearity = |OP(x)+OP(y)-OP(x+y)|/(max(|x|,|y|,|x+y|)
     !> It is zero for real linear operators
     !> (up to machine precision).
     !> It may be greater that zero for quasi linear
     !> operator like linear solver that acts
     !> like approximate inverse.
     !> See inverse type for examples.
     real(kind=double) :: linearity = zero
     !> Logical to check inializaztion
     logical :: is_initialized=.False.
   contains
     !> Procedure for computation of (matrix) times (vector)
     !> Each type extending the abs_linop type MUST
     !> contain this procedure with this name
     !> (it will not compile otherwise)
     procedure(multiplication), deferred :: matrix_times_vector
     !> Short-hand version of matrix_times_vector 
     !> without info and lun_err
     procedure, public, pass:: Mxv => Mxv_abs_linop
     !> Procedure printing info of the
     !> linear operator. Dimensions, symmetric etc.
     !> It can be overriden 
     procedure, public, pass:: info => info_abs_linop
     !> Procedure to copy properties from one linear operator
     !> into another.
     procedure, public, pass:: is_like
     !> Procedure to reset properties to default.
     !> (zero dimensions, no name, etc.)
     procedure, public, pass:: to_default
  end type abs_linop
     
  abstract interface
     !>-------------------------------------------------------------
     !> Abstract procedure defining the interface for a general
     !> matrix-vector multiplication
     !>         vec_out = (M) times (vec_in)
     !> (public procedure for class abs_linop)
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
     subroutine multiplication(this,vec_in,vec_out,info,lun_err)
       use Globals
       import abs_linop
       implicit none
       class(abs_linop),  intent(inout) :: this
       real(kind=double), intent(in   ) :: vec_in(this%ncol)
       real(kind=double), intent(inout) :: vec_out(this%nrow)
       integer,           intent(inout) :: info
       integer,           intent(in   ) :: lun_err
     end subroutine multiplication
  end interface

  !>---------------------------------------------------------
  !> Derived type used to create list of pointers to abstract 
  !> Linera Operators, that can be allocatable.
  !> (Fortan does not allow to create arrays of pointers)
  !> E.g. in a Block-Matrix
  !> M = (A D) with A and D of different type, we can create
  !> an array type that points toward A and D.
  !> Declarations: 
  !>  type(spmat)     :: A
  !>  type(diagmat  ) :: D
  !>  type(array_mat) :: list(2)
  !> Assignment
  !>  list(1)%linop => A
  !>  list(2)%linop => D
  !> See BlockMatrix module for a concrete example
  !>---------------------------------------------------------
  type, public :: array_linop
     !> Dimension (nmats)
     !> Array that contains the non-zero blocks
     class(abs_linop), pointer :: linop=> null()
  end type array_linop
  
  !>---------------------------------------------------------
  !> Structure variable containing the variables and the procedures
  !> to handle a block matrices e.g.
  !> ( A    0 B  0)
  !> ( 0   -C 0  D)
  !> ( B    0 D  0)
  !> and perform Matrix-vector and Matrix(transpose)-vector operations
  type, extends(abs_linop), public :: block_linop
     !> Number of block matrices
     !> It is less or equal to the number of blocks
     !> In the example nlinops = 4
     integer :: nlinops=0
     !> Dimension (nlinops)
     !> Array that contains the pointer to non-zero blocks
     !> In the example  mmdiff
     !> linops(1)%linop => A
     !> linops(2)%linop => B
     !> linops(3)%linop => C
     !> linops(4)%linop => D
     type(array_linop), allocatable :: linop_list(:)
     !> Number of blocks in each row
     !> In the example nblock_row=3
     integer :: nblock_row=0
     !> Number of blocks in each column
     !> In the example nblock_col=4
     integer :: nblock_col=0
     !> Number of non-zero blocks
     !> In the example nnzblock=6
     integer :: nnzblock=0
     !> Dimension (3,nnzblock)
     !> It contains the list of the matrices used
     !> reading the Block Matrix line by line.
     !> Given i in 1, nnzblock
     !> - imat = block_structure(1,i)   
     !>   imat = index of the i^th-block in linops arrays
     !> - block_structure(2,i) = row    index of the i^th block 
     !> - block_structure(3,i) = column index of the i^th block 
     !> In the example
     !> block_structure(:,1)=(1,1,1)
     !> block_structure(:,2)=(2,1,3)
     !> block_structure(:,3)=(3,2,2)
     !> block_structure(:,4)=(4,2,4)
     !> block_structure(:,5)=(2,3,1)
     !> block_structure(:,6)=(4,3,3)     
     integer, allocatable :: block_structure(:,:)
     !> Dimension(nnzblock)
     !> Real coefficient multiplying each non zero operator
     !> In the example :
     !> alpha(1)=1
     !> alpha(2)=1
     !> alpha(3)=1
     !> alpha(4)=-1
     !> alpha(5)=1
     !> alpha(6)=1
     real(kind=double), allocatable :: alphas(:)
     !> Dimension(nblock_row)
     !> Row dimension 
     !> In the example
     !> nrow_vec(1) = matrix_A%nrow ( must match matrix_B%nrow)
     !> nrow_vec(2) = matrix_C%nrow ( must match matrix_D%nrow)
     !> nrow_vec(3) = matrix_B%ncol ( must match matrix_D%nrow)
     integer, allocatable :: nrow_vec(:)
     !> Dimension(nblock_row)
     !> Row dimension 
     !> In the example
     !> ncol_vec(1) = matrix_A%ncol ( must match matrix_B%nrow)
     !> ncol_vec(2) = matrix_C%ncol 
     !> ncol_vec(3) = matrix_B%ncol ( must match matrix_D%ncol)
     !> ncol_vec(4) = matrix_D%ncol 
     integer, allocatable :: ncol_vec(:)     
     !> Dimension(nrow)
     !> Scratch array for Matrix Vector application
     real(kind=double), allocatable :: scr_nrow(:)
     !> Dimension(ncol)
     !> Scratch array for Matrix Vector application
     real(kind=double), allocatable :: scr_ncol(:)
     !>----------------------------------------------
     !> SADDLE POINT VARIABLES
     !>----------------------------------------------
     !> There is a dedicate subroutine (see 
     !> saddle_point => init_saddle_point )
     !> to initialize this particular
     !> type of block matrices.
     !> ---------------------------------------------
     !> Logical to identify saddle point matrix i.e.,
     !> M=(A  B1T)
     !>   (B2 -C )
     logical :: is_saddle_point=.False.
     !> Logical telling if B1 is equal to B2
     logical :: B1equalB2=.False.
   contains
     !> Static constructor
     !> (procedure public for type block_linop)
     procedure, public, pass :: init => init_block_linop
     !> Static constructor
     !> (procedure public for type block_linop)
     procedure, public, pass :: block_diagonal => init_block_diagonal
     !> Static constructor
     !> (procedure public for type block_linop)
     procedure, public, pass :: saddle_point => init_saddle_point
     !> Static destructor
     !> (procedure public for type block_linop)
     procedure, public, pass :: kill => kill_block_linop
     !> Procedure to compute 
     !>         y = M * x 
     !> (public procedure for type block_linop)
     procedure, public, pass :: matrix_times_vector => Mxv_block
     !> Procedure to build transpose structure 
     !> Given a block matrix return the array
     !> describing the transposed structure
     !> (public procedure for type block_linop)
     procedure, public, pass :: structure_transpose
  end type block_linop

  !>---------------------------------------------------
  !> Structure variable to build implicitely matrices
  !> given as linear combination or compositions  
  !> of linear operator.
  !> For example given two matrices A and B 
  !> it is possible to define an new linear operator
  !> M such that
  !> M * v = A*v + B*v
  !> or 
  !> M*v = A*B*v
  !> This type is the building block for type new_linop
  !>---------------------------------------------------
  public :: pair_linop
  type, extends(abs_linop) :: pair_linop
     !> Number of linear operators
     integer :: nlinops=0
     !> List of linwear operator
     type(array_linop), allocatable :: linop_list(:)
     !> Reals scaling each linear operator
     real(kind=double), allocatable :: alphas(:)
     !> String (2 characters) to tell if current pair
     !> describes a linear combinations of operators (LC)
     !> or a product (LP)
     character(len=2)  :: type='??'! LC,LP
     !> String to specify that some particular pattern was
     !> of operators was used. For example:
     !> Least square matrix : 'ATA'
     !> Weighted least square matrix : 'ATWA'
     !> They will be used when explicit forms of the operator
     !> will be formed
     character(len=70) :: pattern='none'
     !> Leading real number used to scale matrix or change sign
     real(kind=double) :: alpha=one
     !> Integer to store position in array of pairs.
     !> Used by new_linop.
     integer :: index=0
     !> Size of scr space
     integer :: nscr=0
     !> Scratch array
     real(kind=double), allocatable :: scr(:)
   contains
     !> Procedure to compute 
     !>         y = M * x 
     !> (public procedure for type pair_linop)
     procedure, public , pass:: matrix_times_vector => add_Mxv
     !> Procedure to compute 
     !>         y = M * x 
     !> (public procedure for type pair_linop)
     procedure, public , pass:: kill => kill_pair_linop
     !> (public procedure for type pair_linop)
     procedure, public , pass:: init => init_pair_multiple
     !> (public procedure for type pair_linop)
     procedure, public , pass:: info => info_pair_linop
     !> (public procedure for type pair_linop)
     !procedure, public , pass:: kill => kill_pair_linop
     !> (public procedure for type pair_linop)
     procedure, public , pass:: init_pair_times_scalar
  end type pair_linop

  public ::  operator(-) ,operator(+),operator(*)
  public :: assignment(=)

  !>---------------------------------------------------
  !> Structure variable to define implicitely
  !> any linear combination or productor or linear
  !> operator. 
  !> For example given the matrices A,B,C 
  !> it is possible to define an new linear operator
  !> M such that
  !> M * v = A*v + B*v
  !> writing
  !> M = A + B        ( M*v = A*v+B*v         \forall v )
  !> M = A - B        ( M*v = A*v-B*v         \forall v )
  !> M = A * B        ( M*v = A*B*v           \forall v )
  !> M = 2*A*B - 3*C  ( M*v = 2*A*B*v-3*C*v   \forall v )
  !> The operators (+,-,*) creates automatically a
  !> structure made by list of pair_linop that coincides
  !> with the expression given. For example, writing
  !>
  !> M = A*B - C
  !>
  !> two pair_linop operators
  !>
  !> N = A*B
  !> M = N-C
  !>
  !> are created in the array pair_list.
  !> With the current algorithm the
  !> operator M is described by the last pair_linop
  !> varible in pair_list (N-C in the example).
  !> (WARNING tested with gfortran only)
  !> This type work with a delicate handle of
  !> fortran pointers.
  !>-----------------------------------------------------
  type, extends(abs_linop), public :: new_linop
     !> Number of pairs
     integer :: npairs=0
     !> List of pairs 
     type(pair_linop), allocatable :: pair_list(:)
   contains
     procedure, public , pass:: matrix_times_vector => new_linop_Mxv
     procedure, public :: info => info_new_linop
     procedure, public :: kill => kill_new_linop
     procedure, public, nopass ::  copy_new_linop
     procedure, public , pass:: append_pair
  end type new_linop

  INTERFACE ASSIGNMENT (=)
     MODULE PROCEDURE copy_new_linop
  END INTERFACE ASSIGNMENT (=)

  INTERFACE OPERATOR (+)
     module PROCEDURE mmsum
  END INTERFACE OPERATOR (+)

  INTERFACE OPERATOR (-)
     module PROCEDURE mmdiff
  END INTERFACE OPERATOR (-)

  INTERFACE OPERATOR (*)
     PROCEDURE mmprod ! matrix  * matrix
     PROCEDURE rmprod ! real    * matrix
     PROCEDURE mrprod ! matrix  * real 
     PROCEDURE improd ! integer * matrix
     PROCEDURE miprod ! matrix  * integer
  END INTERFACE OPERATOR (*)
  
contains
  !>------------------------------------------------------
  !> Subruotine to copy all proprieties of
  !> adjoint linear operator. 
  !> (public procedure for class adj_linop)
  !> 
  !> usage:
  !>     call 'var'%is_like(matrix)
  !<------------------------------------------------------
  subroutine is_like(this,matrix)
    implicit none
    class(abs_linop),  intent(inout) :: this
    class(abs_linop),  intent(in   ) :: matrix

    this%nrow         = matrix%nrow
    this%ncol         = matrix%ncol
    this%is_symmetric = matrix%is_symmetric
    this%unitary_diag = matrix%unitary_diag
    this%triangular   = matrix%triangular
    this%name         = matrix%name
    this%linearity    = matrix%linearity

  end subroutine is_like

  !>------------------------------------------------------
  !> Subruotine to rest  all proprieties of
  !>  linear operator. Useful in the 
  !> initialiazionion of concrete type extending the
  !> adj_linop class.
  !> (public procedure for class adj_linop)
  !> 
  !> usage:
  !>     call 'var'%is_like(matrix)
  !<-------------------------------------------------------------
  subroutine to_default(this)
    implicit none
    class(abs_linop),  intent(inout) :: this

    this%nrow           = 0
    this%ncol           = 0
    this%is_symmetric   = .False.
    this%unitary_diag   = .False.
    this%triangular     = 'N'
    this%name           = '?'
    this%is_initialized = .False.
    
  end subroutine to_default
  

  !>------------------------------------------------------
  !> Subruotine to print linear operator properties
  !> (public procedure for class adj_linop)
  !> 
  !> usage:
  !>     call 'var'%info(lun)
  !> where :
  !>
  !> \param[in   ] lun -> integer. I/O logical unit 
  !<-------------------------------------------------------------
  subroutine info_abs_linop(this,lun)
    implicit none
    class(abs_linop),  intent(in   ) :: this
    integer,           intent(in   ) :: lun
    !local 
    character(len=256) :: msg

       
    write(msg,*) 'NROW= ', this%nrow, ' NCOL= ', this%ncol
    if ( this%is_symmetric) then
       msg=etb(msg)//'| SYMMETRIC' 
    else
       if ( this%triangular .ne. 'N') then
          if ( this%triangular .eq. 'U')  then
             msg=etb(msg)//' | UPPER TRIANGULAR' 
          else
             msg=etb(msg)//' | LOWER TRIANGULAR' 
          end if
       else
          msg=etb(msg)//' | NON SYMMETRIC' 
       end if
    end if
    
    if ( this%name .ne. '?') msg=etb(etb(this%name)//' :'//etb(msg))

    write(lun,*) etb(msg)

  end subroutine info_abs_linop
  
  !>-------------------------------------------------------------
  !> Procedure for matrix-vector multiplication
  !>         vec_out = (M) times (vec_in)
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
  recursive subroutine mxv_abs_linop(this,vec_in,vec_out,info,lun_err)
    implicit none
    class(abs_linop),  intent(inout) :: this
    real(kind=double), intent(in   ) :: vec_in(this%ncol)
    real(kind=double), intent(inout) :: vec_out(this%nrow)
    integer, optional, intent(inout) :: info
    integer, optional, intent(in   ) :: lun_err
    ! local
    logical :: rc
    integer :: info_loc=0,lun_loc=0
    character(len=256):: str=''
    !
    ! select logical unit for errors
    !
    if ( present(lun_err) ) lun_loc = lun_err

    call this%matrix_times_vector(vec_in,vec_out,info_loc,lun_loc)

    !
    ! if any error occured:
    ! - print warning if info was passed
    ! - stop otherwise
    !
    if (info_loc .ne.0 ) then
       if (this%name.ne.'?') str=this%name
       if ( present(info) ) then
          info=info_loc
          rc = IOerr(lun_loc, wrn_val, 'Mxv_abs_linop', &
               ' Error in matrix '//etb(str)//&
               'vector multiplication info=',info_loc)   
       else
          rc = IOerr(lun_loc, err_val, 'Mxv_abs_linop' , &
               ' Error in matrix  '//etb(str)//&
               ' vector multiplication info=',info_loc)   
       end if
    end if
  end subroutine mxv_abs_linop
    
    

  !>------------------------------------------------------
  !> Function the overides the "=" assignment.
  !> (public procedure for class new_linop)
  !> 
  !> usage:
  !>     'var' = A1
  !<------------------------------------------------------
  subroutine copy_new_linop(copy,original)
    implicit none
    type(new_linop), target, intent(inout) :: copy
    type(new_linop), target, intent(in   ) :: original 

    !
    ! free memort to rewrite
    !
    if (copy%is_initialized) call copy%kill(0)

    !
    ! copy abs_linop properties (nrow, ncol ,symmetry name)
    !
    call copy%is_like(original)


    !
    ! copy new_linop structure
    !
    copy%npairs=original%npairs
    allocate(copy%pair_list(copy%npairs))
    call copy_with_pointers(&
         original%npairs,&
         original%pair_list,&
         copy%pair_list(1:original%npairs))

    !
    ! set initialization flag to true
    !
    copy%is_initialized=.True.
    
  end subroutine copy_new_linop
  

 
  !>------------------------------------------------------
  !> Function that give two linear operator A1 and A2
  !> defines, implicitely, the linear operator
  !> A=A1+A2
  !> (public procedure for class pair_linop)
  !> 
  !> usage:
  !>     'var' = A1 + A2
  !<-------------------------------------------------------------
  recursive function mmprod(matrix_1,matrix_2) result(this)
    use Globals
    implicit none
    class(abs_linop), target, intent(in) :: matrix_1
    class(abs_linop), target, intent(in) :: matrix_2
    type(new_linop), target :: this
    ! local
    logical :: rc
    integer :: res,npairs,i,j,info,lun_err=0
    type(pair_linop), allocatable :: temp(:)

    logical :: debug=.False.
    character(len=70) :: str1=' ',str2=' '

    if (full_debug)  debug=.true.
    !debug=.True.


    if(debug) then
       write(0,*)'**********************************'
       write(0,*) ' BEGIN MMPROD: TODO ',etb(matrix_1%name),' * ',etb(matrix_2%name)
    end if

    if (debug) then
       write(0,*) 'this%is_initialized',this%is_initialized
    end if

    !
    ! check consistency
    !
    if (matrix_1%ncol .ne. matrix_2%nrow)  then
       write(str1,*) 'LHS ncol=',matrix_1%ncol
       write(str2,*) 'RHS ncol=',matrix_2%nrow
       rc = IOerr(0, err_val, 'mmprod', &
            ' dimension must agree. Passed:'//etb(str1)//etb(str2))
       
    end if

    !
    ! init 
    !
    call append_pair(this,matrix_1,matrix_2,'*')


    if(debug) then
       write(0,*) ' END MMPROD'
       call this%info(6)
       write(0,*)'**********************************'
    end if
    
  end function mmprod

  !>------------------------------------------------------
  !> Function that give two linear operator A1 and A2
  !> defines, implicitely, the linear operator
  !> A=A1+A2
  !> (public procedure for class pair_linop)
  !> 
  !> usage:
  !>     'var' = A1 + A2
  !<-------------------------------------------------------------
  recursive function rmprod(factor,matrix) result(this)
    use Globals
    implicit none
    real(kind=double),        intent(in) :: factor
    class(abs_linop), target, intent(in) :: matrix
    type(new_linop), target :: this
    ! local
    logical :: rc
    integer :: res,npairs,i,info,lun_err=0
    integer :: nscr_required
    logical :: debug=.False.


    if (full_debug)  debug=.true.
    !debug=.True.
    if(debug) then
       write(0,*) ' '
       write(0,*)'**********************************'
       write(0,*) ' BEGIN RMPROD'
    end if

    call factor_times_matrix(factor,matrix,this)
    
    if(debug) then
       write(0,*) ' END RMPROD'
       call this%info(6)
       write(0,*)'**********************************'

    end if
       
  end function rmprod

   recursive function mrprod(matrix,factor) result(this)
    use Globals
    implicit none
    class(abs_linop), target, intent(in) :: matrix
    real(kind=double),        intent(in) :: factor
    type(new_linop), target :: this
    ! local
    logical :: rc
    integer :: res,npairs,i,info,lun_err=0
    integer :: nscr_required
    logical :: debug=.False.

    if (full_debug)  debug=.true.
    !debug=.True.

    if(debug) then
       write(0,*) ' '
       write(0,*)'**********************************'
       write(0,*) ' BEGIN MRPROD'
    end if

    call factor_times_matrix(factor,matrix,this)

    if(debug) then
       write(0,*) ' END MRPROD'
       call this%info(6)
       write(0,*)'**********************************'
       
    end if

  end function mrprod

  recursive function improd(factor,matrix) result(this)
    use Globals
    implicit none
    integer          ,        intent(in) :: factor
    class(abs_linop), target, intent(in) :: matrix
    type(new_linop), target :: this
    ! local
    logical :: rc
    integer :: res,npairs,i,info,lun_err=0
    integer :: nscr_required
    logical :: debug=.False.
    real(kind=double) :: rfactor


    if (full_debug)  debug=.true.
    !debug=.True.


    if(debug) then
       write(0,*) ' '
       write(0,*)'**********************************'
       write(0,*) ' BEGIN IMPROD'
    end if

    !
    ! init structure
    !
    rfactor=one*factor
    call factor_times_matrix(rfactor,matrix,this)
    
    if(debug) then
       write(0,*) ' END IMPROD'
       call this%info(6)
       write(0,*)'**********************************'

    end if

  end function improd

   recursive function miprod(matrix,factor) result(this)
    use Globals
    implicit none
    class(abs_linop), target, intent(in) :: matrix
    integer          ,        intent(in) :: factor
    type(new_linop), target :: this
    ! local
    logical :: rc
    integer :: res,npairs,i,info,lun_err=0
    logical :: debug=.False.
    real(kind=double) :: rfactor



    if (full_debug)  debug=.true.
    !debug=.True.

    if(debug) then
       write(0,*) ' '
       write(0,*)'**********************************'
       write(0,*) ' BEGIN MIPROD'
    end if

    !
    ! initialize matrix
    !
    rfactor=one*factor
    call factor_times_matrix(rfactor,matrix,this)
    
    
    if(debug) then
       write(0,*) ' END MIPROD'
        call this%info(6)
       write(0,*)'**********************************'      
    end if

  end function miprod

  
  recursive subroutine factor_times_matrix(factor,matrix,this)
    use Globals
    implicit none
    real(kind=double),        intent(in   ) :: factor
    class(abs_linop), target, intent(in   ) :: matrix
    type(new_linop), target,  intent(inout):: this
    ! local
    logical :: rc
    integer :: res,npairs,i,j,info,lun_err=0
    logical :: debug=.False.
    character(len=1) :: operator

    if (full_debug)  debug=.true.
    !debug=.True.
    if(debug) then
       write(0,*) ' '
       write(0,*)'**********************************'
       write(0,*) ' factor_times_matrix'
    end if

    !
    ! free memory to rewrite
    !
    if (this%is_initialized) call this%kill(0)

    !
    ! copy properties 
    !
    call this%is_like(matrix)
    

    !
    ! If a pairs exists scale the constant alpha
    ! of the last pair. Otherwise initialize
    ! a pair.
    !
    select type ( matrix)
    type is ( new_linop )
       this%npairs=matrix%npairs
       allocate(this%pair_list(matrix%npairs))
       call copy_with_pointers(&
            matrix%npairs,&
            matrix%pair_list,this%pair_list(1:matrix%npairs))
       this%pair_list(this%npairs)%alpha=&
            this%pair_list(this%npairs)%alpha*factor
       write(this%name,'(a,1pe8.2,a,a)') '(',&
            this%pair_list(this%npairs)%alpha,&
            ')*(',etb(this%pair_list(this%npairs)%name)
    class default
       this%npairs=1
       allocate(this%pair_list(this%npairs))
       call this%pair_list(this%npairs)%init_pair_times_scalar(&
            factor,matrix)
       this%pair_list(this%npairs)%index=1
       this%name=etb(this%pair_list(this%npairs)%name)
    end select

    this%is_initialized=.True.
    
    if(debug) then
       do i=1,this%npairs
          write(*,'(a,I3,1pe9.2,I16,I16,I16)') etb(this%pair_list(i)%name),&
               this%pair_list(i)%index,&
               this%pair_list(i)%alpha,&
               loc(this%pair_list(i))
          do j=1,this%pair_list(i)%nlinops
             write(*,'(a,I16)') etb(this%pair_list(i)%linop_list(j)%linop%name),&
                  loc(this%pair_list(i)%linop_list(j)%linop)
          end do
       end do

       write(0,*) ' END FACTOR_TIMES_MATRIX'
       write(0,*)'**********************************'

    end if
    
       
  end subroutine factor_times_matrix


  
 


  
  !>------------------------------------------------------
  !> Function that give two linear operator A1 and A2
  !> defines, implicitely, the linear operator
  !> A=A1+A2
  !> (public procedure for class pair_linop)
  !> 
  !> usage:
  !>     'var' = A1 + A2
  !<------------------------------------------------------
  recursive function mmsum(matrix_1,matrix_2) result(this)
    use Globals
    implicit none
    class(abs_linop), target, intent(in) :: matrix_1
    class(abs_linop), target, intent(in) :: matrix_2
    type(new_linop), target :: this
    ! local
    logical :: rc
    integer :: res,npairs,i,j,info,lun_err=0
    type(pair_linop), allocatable :: temp(:)
    logical :: debug=.False.

    if (full_debug)  debug=.true.
    !debug=.True.

    
    if(debug) then
       write(0,*) ' '
       write(0,*)'**********************************'
       write(0,*) ' BEGIN MMSUM'
    end if

    !
    ! check dimensions
    !
    if (matrix_1%nrow .ne. matrix_2%nrow)  &
         rc = IOerr(0, err_val, 'mmsum', &
         ' nrow not consistent'//&
         etb(matrix_1%name)//' + '//&
         etb(matrix_2%name))
    
    if (matrix_1%ncol .ne. matrix_2%ncol)  &
         rc = IOerr(0, err_val, 'mmsum', &
         ' ncol not consistent'//&
         etb(matrix_1%name)//' + '//&
         etb(matrix_2%name))

    !
    ! init. matrix
    !
    call append_pair(this,matrix_1,matrix_2,'+')
    
    if(debug) then
       write(0,*) ' END MMSUM'
       call this%info(6)
       write(0,*)'**********************************'
    end  if

  end function mmsum

  !>------------------------------------------------------
  !> Function that given two linear operator A1 and A2
  !> defines, implicitely, the linear operator
  !> A=A1-A2
  !> (public procedure for class pair_linop)
  !> 
  !> usage:
  !>     'var' = A1 - A2
  !<-------------------------------------------------------------
  recursive function mmdiff(matrix_1,matrix_2) result(this)
    use Globals
    implicit none
    class(abs_linop), target, intent(in) :: matrix_1
    class(abs_linop), target, intent(in) :: matrix_2
    type(new_linop), target :: this
    ! local
    logical :: rc
    integer :: res,npairs,i,j,info,lun_err=0
    type(pair_linop), allocatable :: temp(:)
    logical :: debug=.False.

    if (full_debug)  debug=.true.
    !debug=.True.

    if(debug) then
       write(0,*) ' '
       write(0,*)'**********************************'
       write(0,*) ' BEGIN MMDIFF'
    end if

    !
    ! Check consistency
    !
    if (matrix_1%nrow .ne. matrix_2%nrow)  &
         rc = IOerr(0, err_val, 'mmdiff', &
         ' dimension must agree ')

    if (matrix_1%ncol .ne. matrix_2%ncol)  &
         rc = IOerr(0, err_val, 'mmdiff', &
         ' dimension must agree ')

    !
    ! init 
    !
    call append_pair(this,matrix_1,matrix_2,'-')

    if(debug) then
       write(0,*) ' END MMDIFF',this%nrow,this%ncol
       call this%info(6)
       write(0,*)'**********************************'
    end if
    
  end function mmdiff


  !>------------------------------------------------------
  !> Subroutine handling the initializiation of
  !> operaotr A= A1 (.operator.) A2
  !> with .operator.=(+,-.*).
  !> (public procedure for class new_linop)
  !> 
  !> usage:
  !>     call this%(matrix_1,matrix_2,operator)
  !<-------------------------------------------------------------
  recursive subroutine append_pair(this,matrix_1,matrix_2,&
       operator)
    use Globals
    implicit none
    class(new_linop), target, intent(inout) :: this
    class(abs_linop), target, intent(in) :: matrix_1
    class(abs_linop), target, intent(in) :: matrix_2
    character(len=1), intent(in ) :: operator
    !local
    logical :: rc
    integer :: res
    integer :: i,j,start
    integer :: lun_err=0
    integer :: info,npairs,nlinops,nlinops1,nlinops2
    logical :: debug=.False.,try2simplify=.FALSE.
    character(len=2) :: type
    type(array_linop), allocatable :: temp(:),temp2(:)
    real(kind=double), allocatable :: alphas(:),alphas2(:)
    real(kind=double) :: alpha
    type(pair_linop), pointer  :: matrix
    
    if (full_debug)  debug=.true.
    !debug=.True.


    
    if (debug) then
       write(0,*)'**********************************'
       write(0,*) ' BEGIN APPEND PAIR'
       write(0,*) 'TODO ',etb(matrix_1%name),' (',operator,') ',etb(matrix_2%name)
    end if

    !
    ! free memory to rewrite
    !
    if (debug) then
       write(0,*) 'this%is_initialized',this%is_initialized
    end if
    if (this%is_initialized) call this%kill(0)
    if ( &
         (operator .eq. '+') .or. &
         (operator .eq. '-') &
         ) then
       type='LC'
    else
       type='LP'
    end if
    if (debug) write(0,*) type
    

    !
    ! different action is one of the the matrices
    ! is a new_linop tpye
    !
    select type (matrix_1)
    type is (new_linop)
       select type (matrix_2)
       type is (new_linop)
          if (debug) then
             write(0,*) ' merge 2 new linop'
          end if
          !
          ! merge the two list
          !
          if (debug) then
             write(0,*) 'matrix1'
             call matrix_1%info(6)
             write(0,*) 'matrix2'
             call matrix_2%info(6)
          end if

          if ( &
               (operator .eq. '+') .or. &
               (operator .eq. '-') &
               ) then
             type='LC'
          else
             type='LP'
          end if

          if ( ( matrix_1%pair_list(matrix_1%npairs)%type .eq. type) .and. &
               ( matrix_2%pair_list(matrix_2%npairs)%type .eq. type) ) then
             !
             ! allocate space
             !
             this%npairs=matrix_1%npairs+matrix_2%npairs-1
             allocate(this%pair_list(this%npairs), stat=res)
             if (res .ne. 0)   &
                  rc = IOerr(0, err_alloc, 'append_pair', &
                  'type new_linop member pairs list ')

!!$             !
!!$             ! copy pairs and reset indeces
!!$             !
!!$             call copy_with_pointers(matrix_1%npairs-1,matrix_1%pair_list,this%pair_list)
!!$             call copy_with_pointers(matrix_2%npairs-1,matrix_2%pair_list,&
!!$                  this%pair_list(matrix_1%npairs:matrix_1%npairs-matrix_2%npairs-2))
!!$             do i=1,this%npairs
!!$                this%pair_list(i)%index=i
!!$             end do
             
             !
             ! copy pairs and reset indeces, pointers are worng but
             ! redirect pointer will fix this issues
             !
             do i=1,matrix_1%npairs-1
                this%pair_list(i)=matrix_1%pair_list(i)   
             end do
             do i=1,matrix_2%npairs-1
                this%pair_list((matrix_1%npairs-1)+i)=&
                     matrix_2%pair_list(i)
             end do
             do i=1,this%npairs
                this%pair_list(i)%index=i
             end do

             
             !
             ! create new operator ( M1%list, M2%list)
             !
             nlinops1=matrix_1%pair_list(matrix_1%npairs)%nlinops
             nlinops2=matrix_2%pair_list(matrix_2%npairs)%nlinops
             nlinops = nlinops1 + nlinops2
             allocate( temp(nlinops) , alphas(nlinops), stat=res)
             if (res .ne. 0)   &
                  rc = IOerr(0, err_alloc, 'append_pair', &
                  ' work arrays temp, alphas ')
             ! pointers
             do i=1,nlinops1
                temp(i) = matrix_1%pair_list(matrix_1%npairs)%&
                  linop_list(i)
             end do
             do i=1,nlinops2
                temp(nlinops1+i) = &
                     matrix_2%pair_list(matrix_2%npairs)%&
                     linop_list(i)
                select type ( matrix=> temp(nlinops1+i)%linop )
                type is (pair_linop)
                   write(0,*) matrix%index
                   matrix%index=matrix%index+nlinops1-1
                   write(0,*) matrix%index
                end select
                

             end do

                    
             ! alphas
             select case (operator)
             case ('+')
                alphas(1:nlinops1)=&
                     matrix_1%pair_list(matrix_1%npairs)%alpha*&
                     matrix_1%pair_list(matrix_1%npairs)%alphas(1:nlinops1)
                alphas(1+nlinops1:nlinops)=&
                     matrix_2%pair_list(matrix_2%npairs)%alpha*&
                     matrix_2%pair_list(matrix_2%npairs)%alphas(1:nlinops2)
                alpha=one
             case ('-')
                alphas(1:nlinops1)=matrix_1%pair_list(matrix_1%npairs)%alpha*&
                     matrix_1%pair_list(matrix_1%npairs)%alphas(1:nlinops1)
                alphas(1+nlinops1:nlinops)=-&
                     matrix_2%pair_list(matrix_2%npairs)%alpha*&
                     matrix_2%pair_list(matrix_2%npairs)%alphas(1:nlinops2)
                alpha=one
             case ('*')
                alphas(:)=one
                alpha=&
                     matrix_1%pair_list(matrix_1%npairs)%alpha*&
                     matrix_2%pair_list(matrix_2%npairs)%alpha
             end select

             if (try2simplify) call simply_list_alphas(type,nlinops, temp, alphas)

             call this%pair_list(this%npairs)%&
                  init(info,lun_err,&
                  nlinops,temp, type,alphas,alpha)
             this%pair_list(this%npairs)%index=this%npairs
             !
             ! redirect p
             !
             call redirect_pointers(this%npairs, this%pair_list)

          else
             !
             ! merge the pairs in matrix_1 and matrix 2 
             !
             npairs=matrix_1%npairs+matrix_2%npairs
             this%npairs=npairs+1
             allocate( this%pair_list(this%npairs), stat=res)
             if (res .ne. 0)   &
                  rc = IOerr(0, err_alloc, 'append_pair', &
                  'type new_linop member pairs list ')
             
             
             call copy_with_pointers(&
                  matrix_1%npairs,&
                  matrix_1%pair_list(1:matrix_1%npairs),&
                  this%pair_list(1:matrix_1%npairs))
             call copy_with_pointers(&
                  matrix_2%npairs,&
                  matrix_2%pair_list(1:matrix_2%npairs),&   ! original
                  this%pair_list(matrix_1%npairs+1:matrix_1%npairs+matrix_2%npairs)) ! copy
             do i=1,matrix_2%npairs
                this%pair_list(matrix_1%npairs+i)%index = &
                     this%pair_list(matrix_1%npairs+i)%index+&
                     matrix_1%npairs
             end do

             
             !
             ! set structure of the last pair to append M1 .oper. M2
             !
             nlinops=2
             allocate(temp(nlinops),alphas(nlinops))
             temp(1)%linop=> this%pair_list(matrix_1%npairs)
             temp(2)%linop=> this%pair_list(npairs)
             select case (operator)
             case ('+')
                alphas(:) = one
                alpha=one
             case ('-')
                alphas(1)=one
                alphas(2)=-one
                alpha=one
             case ('*')
                alphas(:)=one
                alpha=one
             end select
             
             if (try2simplify) call simply_list_alphas(type,nlinops, temp, alphas)
             

             !
             ! include append pair matrix1 +  matrix2 to list
             !
             !
             ! init last pair
             !
             call this%pair_list(this%npairs)%init(info,lun_err,&
                  nlinops,temp, type,alphas,&
                  alpha)           
             this%pair_list(this%npairs)%index = this%npairs

          end if
       class default
          if (debug) then
             write(0,*) ' '
             write(0,*) 'matrix1 new linop matrix 2 general'
             write(0,*) '***matrix1'
             call matrix_1%info(6)
             write(0,*) '***matrix2'
             call matrix_2%info(6)
          end if
          
          if ( &
               (operator .eq. '+') .or. &
               (operator .eq. '-') &
               ) then
             type='LC'
          else
             type='LP'
          end if

          if ( matrix_1%pair_list(matrix_1%npairs)%type .eq. type) then
             if (debug) then
                write(0,*) ' SAME TYPE APPENDING TO LAST PAIR'
             end if

             !
             ! allocate space
             !
             this%npairs=matrix_1%npairs
             allocate(this%pair_list(this%npairs), stat=res)
             if (res .ne. 0)   &
                  rc = IOerr(0, err_alloc, 'append_pair', &
                  'type new_linop member pairs list ')
             do i=1,matrix_1%npairs-1
                this%pair_list(i)=matrix_1%pair_list(i)
             end do
             
             !
             ! create new operator ( M1%list, M2)
             !            
             nlinops1 = matrix_1%pair_list(matrix_1%npairs)%nlinops
             nlinops = nlinops1+1
             allocate( temp(nlinops) , alphas(nlinops), stat=res)
             if (res .ne. 0)   &
                  rc = IOerr(0, err_alloc, 'append_pair', &
                  ' work arrays temp alphas ')
             ! pointers
             do i=1,nlinops1
                temp(i) = matrix_1%pair_list(matrix_1%npairs)%&
                  linop_list(i)
             end do
             temp(nlinops)%linop => matrix_2

             ! alphas
             select case (operator)
             case ('+')
                alphas(1:nlinops1)= matrix_1%pair_list(matrix_1%npairs)%alphas(1:nlinops1)
                alphas(nlinops)   = one
                alpha = one
             case ('-')
                alphas(1:nlinops1)=matrix_1%pair_list(matrix_1%npairs)%alphas(1:nlinops)
                alphas(nlinops)   =-one
                alpha = one
             case ('*')
                alphas(:)=one
                alpha=matrix_1%pair_list(matrix_1%npairs)%alpha
             end select

             if (try2simplify) call simply_list_alphas(type,nlinops, temp, alphas)

             call this%pair_list(matrix_1%npairs)%&
                  init(info,lun_err,&
                  nlinops,temp, type,alphas,&
                  alpha)
             this%pair_list(matrix_1%npairs)%index=matrix_1%npairs
             !
             ! redirect p
             !
             call redirect_pointers(this%npairs, this%pair_list)
             
          else
             !
             ! copy structure of matrix1 into first (matrix_1%npairs)-pairs
             ! of this
             ! 
             this%npairs=matrix_1%npairs+1
             allocate(this%pair_list(this%npairs))
             call copy_with_pointers(&
                  matrix_1%npairs,&
                  matrix_1%pair_list(1:matrix_1%npairs),&
                  this%pair_list(1:matrix_1%npairs))

             !
             ! set structure of the last pair to append M1 .oper. M2
             !
             allocate(temp(2),alphas(2))
             temp(1)%linop=> this%pair_list(matrix_1%npairs)
             temp(2)%linop=> matrix_2
             select case (operator)
             case ('+')
                alphas(:) = one
                alpha=one
             case ('-')
                alphas(1)=one
                alphas(2)=-one
                alpha=one
             case ('*')
                alphas(:)=one
                alpha=one
             end select

             !
             ! init last pair
             !
             call this%pair_list(this%npairs)%init(info,lun_err,&
                  nlinops,temp, type,alphas,&
                  alpha)           
             this%pair_list(this%npairs)%index = this%npairs

          end if
       end select
       
    class default
       select type (matrix_2)
          
       type is (new_linop) ! A + B*C
          if (debug) then
             write(0,*) ' matrix 1=general, matrix 2 = new_linop'
             write(0,*) 'matrix1'
             call matrix_1%info(6)
             write(0,*) 'matrix1 list'

             write(0,*) 'matrix2'
             call matrix_2%info(6)
          end if

          if ( &
               (operator .eq. '+') .or. &
               (operator .eq. '-') &
               ) then
             type='LC'
          else
             type='LP'
          end if

          if ( matrix_2%pair_list(matrix_2%npairs)%type .eq. type) then
             !
             ! allocate space
             !
             this%npairs=matrix_2%npairs
             allocate(this%pair_list(this%npairs))
             do i=1,matrix_2%npairs-1
                this%pair_list(i)=matrix_2%pair_list(i)
             end do
             
             !
             ! create new operator ( M1%list, M2)
             !            
             nlinops2 = matrix_2%pair_list(matrix_2%npairs)%nlinops
             nlinops  = nlinops2+1
             allocate( temp(nlinops) , alphas(nlinops))
             ! pointers
             temp(1)%linop => matrix_1
             do i=1,nlinops2
                temp(1+i) = matrix_2%pair_list(matrix_2%npairs)%&
                     linop_list(i)
             end do
             
             ! alphas
             
             select case (operator)
             case ('+')
                alphas(1)=one
                alphas(2:nlinops)=&
                     matrix_2%pair_list(matrix_2%npairs)%alpha*&
                     matrix_2%pair_list(matrix_2%npairs)%alphas(1:nlinops2)
                alpha=one
             case ('-')
                alphas(1)=one
                alphas(2:nlinops)=&
                     -matrix_2%pair_list(matrix_2%npairs)%alpha*&
                     matrix_2%pair_list(matrix_2%npairs)%alphas(1:nlinops2)
                alpha=one
             case ('*')
                alphas(:)=one
                alpha=matrix_2%pair_list(matrix_2%npairs)%alpha
             end select


             if (try2simplify) call simply_list_alphas(type,nlinops, temp, alphas)


             call this%pair_list(matrix_2%npairs)%&
                  init(info,lun_err,&
                  nlinops,temp, type,alphas,&
                  alpha)
             this%pair_list(matrix_2%npairs)%index=matrix_2%npairs
             !
             ! redirect p
             !
             call redirect_pointers(this%npairs, this%pair_list)
          else
             !
             ! copy structure of matrix 2 into first
             ! pairs of this
             !
             this%npairs=matrix_2%npairs+1
             allocate(this%pair_list(this%npairs))
             do i=1,matrix_2%npairs
                this%pair_list(i) = matrix_2%pair_list(i)
             end do
             call copy_with_pointers(matrix_2%npairs,&
                  matrix_2%pair_list,&
                  this%pair_list)

             !
             ! set pair M1 .oper. M2 in temp
             !
             nlinops=2
             allocate(temp(nlinops),alphas(nlinops), stat=res)
             if (res .ne. 0)   &
                  rc = IOerr(0, err_alloc, 'append_pair', &
                  ' work array temp')
             temp(1)%linop=> matrix_1
             temp(2)%linop=> this%pair_list(matrix_2%npairs)
             select case (operator)
             case ('+')
                alphas(:) = one
                alpha=one
             case ('-')
                alphas(1)=one
                alphas(2)=-one
                alpha=one
             case ('*')
                alphas(:)=one
                alpha=one
             end select


             !
             ! init last pair
             !
             call this%pair_list(this%npairs)%init(info,lun_err,&
                  nlinops,temp, type,alphas,&
                  alpha)           
             this%pair_list(this%npairs)%index = this%npairs
             
          end if
       class default
          !
          ! set temporary arrays
          !
          allocate(temp(2),alphas(2), stat=res)
          if (res .ne. 0)   &
               rc = IOerr(0, err_alloc, 'append_pair', &
               ' work array temp')

          temp(1)%linop=> matrix_1
          temp(2)%linop=> matrix_2

          select case (operator)
          case ('+')
             alphas(:) = one
             alpha=one
          case ('-')
             alphas(1)=one
             alphas(2)=-one
             alpha=one
          case ('*')
             alphas(:)=one
             alpha=one
          end select
          !
          ! init new linop with just one pair
          !
          this%npairs=1
          allocate(this%pair_list(this%npairs))
          call this%pair_list(1)%&
               init(info,lun_err,&
               2,temp, type,alphas,&
               alpha)          
          this%pair_list(this%npairs)%index=1

          deallocate(temp,alphas, stat=res)
          if (res .ne. 0)   &
               rc = IOerr(0, err_dealloc, 'append_pair', &
               ' work array temp')
       end select
    end select

    !
    ! the resulting operator has the
    ! same properties of the last pairs
    ! (it is, matematically speaking, the last pairs)
    !
    call this%is_like(this%pair_list(this%npairs))
    this%is_initialized=.True.

    
    if (debug) write(0,*) 'out ', etb(this%name)

    if (debug) then
       write(0,*) ' END APPEND PAIR'
       write(0,*)'**********************************'
    end if

  end subroutine append_pair

  !>------------------------------------------------------
  !> Procedure aimed to symply the structure of
  !> linear combinatio and product.
  !> NOT WORKING YET
  !<-------------------------------------------------------------
  subroutine simply_list_alphas(type,nlinops, temp, alphas)
    implicit none    
    character(len=2),               intent(in   ) :: type
    integer,                        intent(inout) :: nlinops
    type(array_linop), allocatable, intent(inout) :: temp(:)
    real(kind=double), allocatable, intent(inout) :: alphas(:)
    
    !
    integer :: i,j,start, nsimplified
    logical :: debug,simplified
    type(array_linop), allocatable :: temp2(:)
    real(kind=double), allocatable :: alphas2(:)

    if (full_debug)  debug=.true.
    !debug=.True.
    
    if (debug) write(0,*) ' BEGIN SIMPLIFY ', nlinops
    
    nsimplified=nlinops
    if (debug) write(0,*) nsimplified
    simplified=.False.
    do i=1,nlinops
       select type( matrix => temp(i)%linop )
       type is (pair_linop)
          if ( (matrix%type .eq. 'LC') .and. ( type .eq.'LC') ) then
             nsimplified = nsimplified - 1 + matrix%nlinops
             simplified=.TRUE.
          end if
          if ( ( matrix%type    .eq. 'LP') .and. &
               ( type           .eq. 'LC') .and. &
               ( matrix%nlinops .eq. 1   ) ) then
             simplified=.TRUE.
          end if
       end select
    end do
    if ( simplified) then
       allocate(temp2(nsimplified),alphas2(nsimplified))
       start=1
       do i=1,nlinops
          select type( matrix => temp(i)%linop )
          type is (pair_linop)
             if ( (matrix%type .eq. 'LC') .and. ( type .eq. 'LC') ) then
                do j=1,matrix%nlinops
                   temp2(start+j-1)%linop  => matrix%linop_list(i)%linop
                   alphas(start+j-1) = alphas(i)*matrix%alpha*matrix%alphas(j) 
                end do
                start=start+matrix%nlinops
             end if
             if ( ( matrix%nlinops .eq. 1   ) .and. &
                  ( matrix%type    .eq. 'LP')  .and. &
                  ( type           .eq. 'LC')  ) then
                temp2(start)%linop  => matrix%linop_list(1)%linop
                alphas(start) = alphas(i)*matrix%alpha*matrix%alphas(1) 
                start=start+matrix%nlinops
             end if
          class default
             temp2(start)  = temp(i)
             alphas(start) = alphas(i)
             start=start+1
          end select
       end do
       deallocate(temp,alphas)
       allocate(temp(nsimplified),alphas(nsimplified))
       temp=temp2
       alphas=alphas2
       nlinops=nsimplified
    end if
    if (debug) write(0,*) 
    if (debug) write(0,*) ' END SIMPLIFY ', nlinops
  end subroutine simply_list_alphas

  !>------------------------------------------------------
  !> Procedure handling the proper reassignment
  !> of pointers describing the operator defined in new_linop.
  !> For example the operator
  !>
  !> M=A+B*C-D
  !>
  !> is decomposed in as
  !>
  !> N=(B,*,C,1)
  !> P=(A,+,N,2)
  !> Q=(P,-,D,3) (this coincides with M)
  !>
  !> But if we write M'=M, which copy each component of
  !> M we will get:
  !> 
  !> N'=N  -> N'=(B,*,C,1) right
  !> P'=P  -> P'=(A,+,N,2) N is wrong, should be N'
  !> Q'=Q  -> Q'=(P,+,D,2) P is wrong, should be P'
  !>
  !> This subrotuine will instead creates
  !> N'=(B ,*,C ,1)
  !> P'=(A ,+,N',2)
  !> Q'=(P',+,D, 2) ( this is M' )
  !>
  !> This is fundamental in append_pair procedure.
  !<-------------------------------------------------------------
  subroutine copy_with_pointers(npairs, original, copy)
    implicit none
    integer :: npairs
    type(pair_linop),target :: original(npairs)
    type(pair_linop),target :: copy(npairs)
    !
    integer :: i,j,index 
    logical :: debug=.False.,found
    class(abs_linop), pointer :: matrix


    if (full_debug)  debug=.true.
    debug=.False.


    if (debug) then
       write(0,*)'**********************************'
       write(0,*) ' BEGIN COPY POINTER'
    end if

    !
    ! we copy the 
    !
    copy=original


    !
    ! reassing pointer in copy to itself components
    !
    do i=1, npairs
       ! i=2 mat=(A,+,N,2)
       ! i=3 mat=(P,+,D,2)
       if (debug) then
          write(0,*) 'pair in ',i
          call copy(i)%info(6)
       end if

       do j=1,copy(i)%nlinops
          select type (matrix=>copy(i)%linop_list(j)%linop)
          class is ( pair_linop )
             
             if (debug) write(0,*) 'i,j,index',i,j,matrix%index,matrix%name
             ! i=2 mat=(A,+,N,2)
             index=matrix%index

             ! index=2
             if (debug) write(0,*) loc(matrix),&
                  loc(copy(i)%linop_list(j)%linop),&
                  loc(copy(index))       
             copy(i)%linop_list(j)%linop => copy(index)

             if (debug) write(0,*) loc(matrix),&
                  loc(copy(i)%linop_list(j)%linop),&
                  loc(copy(index))
             !mat=(A,+,N,2)

             ! other classes (class default)
             ! means that current_operator in list points
             ! toward linear perators outside variable
             ! "outside" like A,B,C,D 
             
          end select
          
       end do       
       if (debug) then
          write(0,*) 'pair out ',i
          call copy(i)%info(6)
       end if
    end do

    if (debug) then
       write(0,*) ' END COPY POINTER'
       write(0,*)'**********************************'

    end if

  end subroutine copy_with_pointers

  !>------------------------------------------------------
  !> Procedure handling the proper reassignment
  !> of pointers describing the operator defined in new_linop.
  !> Works like copy_with_pointers, reassigning pointers 
  !> according to index. Only one should remains
  !>---------------------------------------------------------
  subroutine redirect_pointers(npairs, copy)
    implicit none
    integer :: npairs
    type(pair_linop),target :: copy(npairs)
    !
    integer :: i,j,index 
    logical :: debug=.False.,found
    class(abs_linop), pointer :: matrix

    if (full_debug)  debug=.true.
    debug=.False.


    if (debug) then
       write(0,*)'**********************************'
       write(0,*) ' BEGIN REDIRECT POINTER'
    end if

    !
    ! reassing pointer in copy to itself components
    !
    do i=1, npairs
       ! i=2 mat=(A,+,N,2)
       ! i=3 mat=(P,+,D,2)
       if (debug) then
          write(0,*) 'pair in ',i
          call copy(i)%info(6)
       end if

       do j=1,copy(i)%nlinops
          select type (matrix=>copy(i)%linop_list(j)%linop)
          class is ( pair_linop )
             
             if (debug) write(0,*) 'i,j,index',i,j,matrix%index,matrix%name
             ! i=2 mat=(A,+,N,2)
             index=matrix%index

             ! index=2
             if (debug) write(0,*) loc(matrix),&
                  loc(copy(i)%linop_list(j)%linop),&
                  loc(copy(index))       
             copy(i)%linop_list(j)%linop => copy(index)

             if (debug) write(0,*) loc(matrix),&
                  loc(copy(i)%linop_list(j)%linop),&
                  loc(copy(index))
             !mat=(A,+,N,2)
          end select
          
       end do

       ! class default means that rhs points
       ! toward linear perators outside variable
       ! "outside" like A,B,C,D
       if (debug) then
          write(0,*) 'pair out ',i
          call copy(i)%info(6)
       end if
    end do

    if (debug) then
       write(0,*) ' END REDIRECT POINTER'
       write(0,*)'**********************************'

    end if

  end subroutine redirect_pointers

  subroutine init_pair_multiple(this, info,lun_err,&
       nlinops,list_linop, type,alphas,alpha)
    implicit none
    class(pair_linop), target,   intent(inout) :: this
    integer,                     intent(in   ) :: info,lun_err
    integer,                     intent(in   ) :: nlinops
    type(array_linop), target,   intent(in   ) :: list_linop(nlinops)
    character(len=2),            intent(in   ) :: type
    real(kind=double), optional, intent(in   ) :: alphas(nlinops)
    real(kind=double), optional, intent(in   ) :: alpha
    
    !local
    integer :: res
    logical :: rc,fine,is_lower_triangular ,is_upper_triangular 
    integer :: nrow, ncol,nin,nout,i,imat
    character(len=1)  :: operator
    character(len=max_length_name) :: str
    logical :: debug
    real(kind=double) :: factor

    if (full_debug)  debug=.true.
    !debug=.True.

    
    if (debug) then
       write(0,*) 'BEGIN init_pair_multiple ' , type ,nlinops
    end if
    
    if (present(alpha)) then
       this%alpha=alpha
    else
       this%alpha=one
    end if
    
    
    this%type=type
    
    if ( type .eq. 'LC') then
       !
       ! check constistency
       !
       nrow=list_linop(1)%linop%nrow
       ncol=list_linop(1)%linop%ncol
       fine=.True.
       i=2
       !write(0,*) nlinops,nrow,ncol
       do while ( (fine) .and. (i .le. nlinops) )
          !write(0,*) list_linop(i)%linop%nrow,list_linop(i)%linop%ncol
          fine=(nrow .eq. list_linop(i)%linop%nrow)  .and. &
               (ncol .eq. list_linop(i)%linop%ncol )
          i=i+1
       end do
       if (.not. fine) then 
          rc = IOerr(lun_err, wrn_val,&
               'init_pair_multiple', &
               'mismatch dimension in operator', i)
          return
       end if

       !
       ! set structure
       !
       this%nlinops=nlinops
       allocate( &
            this%linop_list(this%nlinops),&
            this%alphas(this%nlinops),&
            stat=res)
       if (res .ne. 0)   &
            rc = IOerr(0, err_alloc, 'init_pair_multiple', &
            ' pair linop member linop_list, alphas')
       !
       ! copy pointers
       !
       this%linop_list=list_linop
       if (present(alphas)) then
          this%alphas=alphas
       else
          this%alphas=one
       end if
       
       !
       ! set properties
       !
       this%nrow = list_linop(1)%linop%nrow
       this%ncol = list_linop(1)%linop%ncol
       this%is_symmetric=.True.
       this%triangular='N'
       is_lower_triangular=.True.
       is_upper_triangular=.True.
       this%name='('
       do i=1,nlinops
          this%is_symmetric=&
               this%is_symmetric .and.  &
               list_linop(i)%linop%is_symmetric
          is_lower_triangular = &
               is_lower_triangular .and. &
               (list_linop(i)%linop%triangular .eq. 'L')
          is_upper_triangular = &
               is_upper_triangular .and. &
               (list_linop(i)%linop%triangular .eq. 'U')
          if ( this%alphas(i) .ge. zero) then
             if (this%alphas(i) .eq. one) then
                str='+'
                if (i.eq.1) str=''
             else
                write(str,'(a,1pe8.2,a)') '+',this%alphas(i),'*'
             end if             
          else
             if (this%alphas(i) .eq. -one) then
                str='-'
             else
                write(str,'(a,1pe8.2,a)') '-',abs(this%alphas(i)),'*'
             end if
          end if
          this%name=etb(etb(this%name)//etb(str)//etb(list_linop(i)%linop%name))
       end do
       this%name=etb(etb(this%name)//')')
       if (is_lower_triangular) this%triangular = 'L'
       if (is_upper_triangular) this%triangular = 'U'

       !
       ! set work space required
       !
       this%nscr = list_linop(1)%linop%nrow
    else if ( type .eq. 'LP') then
       !
       ! check dimensions and set dimension  of workspace
       !
       
       ! this makes sense, check how Mxv works 
       this%nscr=list_linop(nlinops)%linop%ncol
       
       do imat = nlinops,2,-1
          !call list_linop(imat)%linop%info(6)
          !call list_linop(imat-1)%linop%info(6)
          nout = list_linop(imat)%linop%nrow
          
          this%nscr=max(this%nscr,nout)
          !write(0,*) nout,this%nscr
          nin = list_linop(imat-1)%linop%ncol 
          this%nscr=max(this%nscr,nin)
          !write(0,*) nin,this%nscr
          if ( nout .ne. nin ) &
               rc = IOerr(lun_err, err_val, &
               ' init_pair_multiple ', &
               ' Dimension mismatch in matrix =', imat)
       end do
       this%nscr=2*this%nscr

       !
       ! set structure
       !
       this%nlinops=nlinops
       allocate( &
            this%linop_list(this%nlinops),&
            this%alphas(this%nlinops),& ! useless 
            stat=res)
       if (res .ne. 0)   &
            rc = IOerr(0, err_alloc, 'mmprod', &
            etb(this%name)//' pair linop member linop_list alphas')
       

       !
       ! set properties
       !
       this%alphas=one
       this%linop_list =  list_linop
       this%nrow = list_linop(1)%linop%nrow
       this%ncol = list_linop(nlinops)%linop%ncol
       this%is_symmetric=.False.
       this%triangular='N'

       this%name=' '
       this%name='('//etb(list_linop(1)%linop%name)
       do i=2,nlinops
          this%name=etb(etb(this%name)//'*'//etb(list_linop(i)%linop%name))
       end do
       this%name=etb(etb(this%name)//')')
    end if

    !
    ! set name and initialization flag to .true.
    !
    if ( this%alpha .ne. one) then
       write(str,'(a,1pe8.2,a,a)') '(',&
            this%alpha,&
            ')*',etb(this%name)
       this%name=etb(str)
    end if
    this%is_initialized=.True.
    
    !
    ! allocate work space
    !
    allocate(this%scr(this%nscr), stat=res)
    if (res .ne. 0)   &
         rc = IOerr(0, err_alloc, 'mmprod', &
         ' work array scr of type pair_linop')

    if (debug) then
       write(0,*) 'END init_pair_multiple' , etb(this%name),this%alpha
    end if

    
  end subroutine init_pair_multiple
  

  subroutine info_pair_linop(this,lun)
    implicit none
    class(pair_linop),  intent(in   ) :: this
    integer,           intent(in   ) :: lun
    character(len=256) :: str1,str2
    !local
    integer :: i,j
    character(len=3) :: str

   
    write(*,'(a,a,1x,a,I3,a,1pe9.2,1x,a,a,a,I3,a,I16)') &
         '  PAIR',etb(this%name),&
         'index=',this%index,&
         ' alpha=',this%alpha,&
         ' type= ',&
         etb(this%type), &
         ' nlinops= ',this%nlinops,&
         ' memory loc. =',loc(this) 
    
         
    do j=1,this%nlinops
       write(str,'(I3)') j
       write(*,'(a2,a,a,1pe9.2,a,a,1x,I16)') &            
            '  a_',etb(str),'=',this%alphas(j),&
            ' * ',etb(this%linop_list(j)%linop%name),&
            loc(this%linop_list(j)%linop)
    end do

    
  end subroutine info_pair_linop

  subroutine info_new_linop(this,lun)
    implicit none
    class(new_linop),  intent(in   ) :: this
    integer,           intent(in   ) :: lun
    ! local
    character(len=256) :: str1,str2
    integer :: j

    write(lun,*) 'NEW_LINOP :', etb(this%name), 'npairs', this%npairs, ' loc ', loc(this), this%is_initialized
    if (allocated( this%pair_list) ) then
       do j=1,this%npairs
          call this%pair_list(j)%info(6)
       end do
    end if
    
  end subroutine info_new_linop

  !>----------------------------------------------------
  !> Static destructor.
  !> (procedure public for type new_linop)
  !> Free memory and diassociate pointers
  !>
  !> usage:
  !>     call 'var'%kill(lun_err)
  !>
  !> where:
  !> \param[in] lun_err -> integer. Error logical unit
  !<----------------------------------------------------
  subroutine kill_new_linop(this,lun_err)
    implicit none
    class(new_linop),  intent(inout) :: this
    integer,           intent(in   ) :: lun_err
    ! local
    character(len=256) :: str1,str2
    integer :: j
    logical :: rc
    integer :: res

    !
    ! free memory
    !
    if (allocated(this%pair_list) ) then       
       do j=this%npairs,1,-1
          !write(0,*) 'killing', j,this%npairs
          call this%pair_list(j)%kill(lun_err)
       end do
       this%npairs=0
       deallocate(this%pair_list,stat=res)
       if (res .ne. 0)  &
            rc = IOerr(0, err_dealloc,&
            'kill_new_linop', &
            'member pair_list')
    end if

    !
    ! reset proprieties
    !
    call this%to_default()    
  end subroutine kill_new_linop
  
    
  subroutine init_pair_times_scalar(this, factor,matrix)
    implicit none
    class(pair_linop), target, intent(inout) :: this
    real(kind=double),         intent(in   ) :: factor
    class(abs_linop), target,  intent(in   ) :: matrix

    !
    ! exactly same approach
    !
    call this%is_like(matrix)
    write(this%name,'(a,1pe8.2,a,a)') '(',factor,')*',etb(matrix%name)

    this%nlinops=1
    allocate(this%linop_list(1),this%alphas(1))

    this%type='LP'
    this%linop_list(1)%linop  => matrix
    this%alpha  = factor
    this%alphas = one

  end subroutine init_pair_times_scalar

  !>----------------------------------------------------
  !> Static destructor.
  !> (procedure public for type pair_linop)
  !> Free memory and deassociate pointers
  !>
  !> usage:
  !>     call 'var'%kill(lun_err)
  !>
  !> where:
  !> \param[in] lun_err -> integer. Error logical unit
  !<----------------------------------------------------
  subroutine kill_pair_linop(this,lun_err)
    use Globals
    implicit none
    class(pair_linop),intent(inout) :: this
    integer,         intent(in   )  :: lun_err
    ! local
    logical :: rc
    integer :: res,i

    !
    ! free memory
    !
    if ( allocated(this%linop_list) ) then
       !
       ! deassociate pointers
       !
       do i=1,this%nlinops
          this%linop_list(i)%linop => null()
       end do
       deallocate(this%linop_list,this%alphas,stat=res)
       if (res .ne. 0)  &
            rc = IOerr(0, err_dealloc,&
            'kill_pair_linop', &
            'member linop_list')
    end if
    if ( allocated(this%scr)) then
       deallocate(this%scr,stat=res)
       if (res .ne. 0)  &
            rc = IOerr(0, err_dealloc,&
            'kill_pair_linop', &
            'member scr')
    end if
    this%nscr = 0

    
    

    !
    ! set to defualt properties
    !
    call this%to_default()    
  end subroutine kill_pair_linop

  !>-------------------------------------------------------------
  !> Concrete procedure procedure Mxv of parent type abs_linop.
  !> It performes the matrix-vector multiplication for type
  !> pair_linop
  !> 
  !>         vec_out = (M) times (vec_in)
  !> (public procedure for class abs_linop)
  !> 
  !> usage:
  !>     call 'var'%matrix_times_vector(vec_in,vec_out,[info,lun_err])
  !>
  !> where 
  !> \param[in]          vec_in  -> real, dimension('var'%ncol)
  !>                                 vector to be multiplied
  !> \param[inout]       vec_out -> real, dimension('var'%nrow)
  !>                                 vector (M) times (vec_in) 
  !> \param[in,optional] info    -> integer. Info number
  !>                                 in case of error
  !> \param[in,optional] lun_err -> integer. Info number
  !>                                 in case of error  
  !<-------------------------------------------------------------
  recursive subroutine add_Mxv(this,vec_in,vec_out,info,lun_err)
    use Globals
    implicit none
    class(pair_linop), intent(inout) :: this
    real(kind=double), intent(in   ) :: vec_in(this%ncol)
    real(kind=double), intent(inout) :: vec_out(this%nrow)
    integer,           intent(inout) :: info
    integer,           intent(in   ) :: lun_err
    !
    logical :: debug=.False.
    integer :: nin,nout,offset,imat,lun

    info=0

    if (full_debug)  debug=.true.
    !debug=.True.

    if (debug) then
       write(0,*)'**********************************'
       write(0,*) ' add pair mxv'
    end if
    if (debug) call this%info(6)
    select case (this%type)
    case ('LC')
       info = 0
       vec_out=zero
       do imat = 1,this%nlinops      
          ! yi = ai Mi x 
          call this%linop_list(imat)%linop%matrix_times_vector(&
               vec_in,this%scr(1:this%nrow),info,lun)
          if ( info .ne. 0 ) return
          
          ! sum yi to y
          call daxpy(this%nrow,this%alphas(imat),this%scr,1,vec_out,1)
       end do
       call dscal(this%nrow,this%alpha,vec_out,1)

    case ('LP')
       !
       ! copy vec_in into scr
       !
       if (this%nlinops .eq. 1) then
          call this%linop_list(1)%linop%matrix_times_vector( &
               vec_in,&
               vec_out,&
               info=info,&
               lun_err=lun)
          if ( info .ne. 0  ) return
       else
          this%scr(1:this%ncol) = vec_in(1:this%ncol)
          offset=this%nscr/2

          !
          ! multiply from the last to the second
          ! this section is ignore if this%nlinops=1
          !
          ! write(0,*) offset,allocated(this%scr),size(this%scr)
          do imat  = this%nlinops,2,-1     
             nin  = this%linop_list(imat)%linop%ncol
             nout = this%linop_list(imat)%linop%nrow
             call this%linop_list(imat)%linop%matrix_times_vector( &
                  this%scr(1:nin),&
                  this%scr(offset+1:offset+nout),&
                  info=info,&
                  lun_err=lun)
             if ( info .ne. 0 ) return
             this%scr(1:nout) = this%scr(offset+1:offset+nout)
          end do

          !
          ! multiply by first operator
          !
          nin  = this%linop_list(1)%linop%ncol
          nout = this%linop_list(1)%linop%nrow
          call this%linop_list(1)%linop%matrix_times_vector( &
               this%scr(1:nin),&
               vec_out(1:nout),&
               info=info,&
               lun_err=lun)

          if ( info .ne. 0) return
       end if
       !
       ! multiply by leading scalar alpha
       ! dscal shoul avoid useless operations
       !
       call dscal(this%nrow,this%alpha,vec_out,1)

    end select

    !if (debug) write(0,*)  'ncol', this%nrow,'vec out', vec_out
    if (debug) then
       write(0,*) ' END pair mxv'
       write(0,*)'**********************************'
    end if

    
  end subroutine add_Mxv

  !>-------------------------------------------------------------
  !> Concrete procedure procedure Mxv of parent type abs_linop.
  !> It performes the matrix-vector multiplication for type
  !> new_linop
  !> 
  !>         vec_out = (M) times (vec_in)
  !> (public procedure for class abs_linop)
  !> 
  !> usage:
  !>     call 'var'%matrix_times_vector(vec_in,vec_out,[info,lun_err])
  !>
  !> where 
  !> \param[in]          vec_in  -> real, dimension('var'%ncol)
  !>                                 vector to be multiplied
  !> \param[inout]       vec_out -> real, dimension('var'%nrow)
  !>                                 vector (M) times (vec_in) 
  !> \param[in,optional] info    -> integer. Info number
  !>                                 in case of error
  !> \param[in,optional] lun_err -> integer. Info number
  !>                                 in case of error  
  !<-------------------------------------------------------------
  recursive subroutine new_linop_Mxv(this,vec_in,vec_out,info,lun_err)
    use Globals
    implicit none
    class(new_linop),  intent(inout) :: this
    real(kind=double), intent(in   ) :: vec_in(this%ncol)
    real(kind=double), intent(inout) :: vec_out(this%nrow)
    integer,           intent(inout) :: info
    integer,           intent(in   ) :: lun_err
    !
    integer :: i
    logical ::debug=.False.

    if (full_debug)  debug=.true.
    !debug=.True.


    
    if (debug) then
       write(0,*) '*************************************** '
       write(0,*) 'BEGIN new_linop Mxv'
       write(0,*) 'this%npairs',this%npairs
       do i=1,this%npairs
          call  this%pair_list(i)%info(6)
       end do
    end if
    call this%pair_list(this%npairs)%matrix_times_vector(&
         vec_in, vec_out,info,lun_err)

    if (debug) then
       write(0,*) 'END new_linop Mxv'
       write(0,*) '*************************************** '
    end if

  end subroutine new_linop_Mxv

  

  

  !>-------------------------------------------------------------
  !> Static constructor.
  !> (procedure public for type block_linop)
  !> Instantiate variable of type block_linop
  !>
  !> usage:
  !>     call 'var'%init(lun_err, nrow, nnzblk )
  !>
  !> where:
  !> \param[in] lun_err               -> integer. Error logical unit
  !> \param[in] nlinops               -> integer. Number of linear operators
  !> \param[in] linops                -> type(array_linop) Array containg
  !>                                       the pointers to lin.op. used.
  !>                                       Dimension=nlinops
  !> \param[in] nrow_block            -> integer. Number of columns
  !> \param[in] ncol_block            -> integer. Number of columns
  !> \param[in] nnzblk                -> integer. Number of non-zero term
  !>                                     stored in the matrix
  !> \param[in] block_structure       -> integer. dimension(3,nnzblock)
  !>                                     Table describing the block structure
  !>                                     reading it by row
  !>                                     block_structure(1,i) = index of
  !>                                         linear operator in array linops 
  !>                                     block_structure(2,i) = row index of
  !>                                         linear operator 
  !>                                     block_structure(3,i) = column index of
  !>                                         linear operator 
  !> \param[in] is_symmetric          -> Logical. Flag for symmetric or not matrix
  !> \param[in] is_symmetric          -> integer. Dimension of the kernel
  !> \param[in,optinal] alphas        -> real.(dimensions=nnzblk)
  !>                                     Scalars premultiplynong each non
  !>                                     zeros lin. op. Default = 1.0
  !> Example
  !> ( A    0 B  0)
  !> ( 0   -C 0  D)
  !> ( B    0 D  0)
  !> 
  !> nliops=4
  !> nrow_block= 3
  !> ncol_block= 4
  !>
  !> Linops(1)%linop=> A
  !> Linops(2)%linop=> B
  !> Linops(3)%linop=> C
  !> Linops(4)%linop=> D
  !>
  !> nnz_blcok=6
  !>
  !> Block Structure(:,1)=(1,1,1)
  !> Block Structure(:,2)=(2,1,2)
  !> Block Structure(:,3)=(3,2,2)
  !> Block Structure(:,4)=(4,2,4)
  !> Block Structure(:,5)=(2,3,1)
  !> Block Structure(:,6)=(4,3,3)
  !>
  !> alphas=(/one,one,-one,one,one,one/)
  !<-------------------------------------------------------------
  ! TODO automatic detection of symmetry or unsymmetric 
  subroutine init_block_linop(this, lun_err, &
       nlinops, linops,&
       nblock_row, nblock_col, &
       nnzblock, block_structure, &
       is_symmetric,&
       alphas)
    use Globals
    implicit none
    !var
    class(block_linop),          intent(inout) :: this
    integer,                     intent(in   ) :: lun_err
    integer,                     intent(in   ) :: nlinops
    type(array_linop),           intent(in   ) :: linops(nlinops)
    integer,                     intent(in   ) :: nblock_row
    integer,                     intent(in   ) :: nblock_col
    integer,                     intent(in   ) :: nnzblock
    integer,                     intent(in   ) :: block_structure(3,nnzblock)
    logical,  optional,          intent(in   ) :: is_symmetric
    real(kind=double), optional, intent(in   ) :: alphas(nnzblock)
    ! local vars
    logical :: rc
    integer :: res
    integer :: imat, nrow_loc, ncol_loc,irow,icol,iterm,ilinop
    character(len=256) :: msg

    !
    ! free memory in case of re-initialization
    !
    if(this%is_initialized) call this%kill(lun_err)

    !
    ! copy block structure infos
    !
    this%nlinops    = nlinops
    this%nblock_row = nblock_row
    this%nblock_col = nblock_col
    this%nnzblock   = nnzblock
    allocate(&
         this%linop_list(nlinops),&
         this%nrow_vec(nblock_row),&
         this%ncol_vec(nblock_col),&
         this%block_structure(3,nnzblock),&
         this%alphas(nnzblock),&
         stat=res)
    if (res .ne. 0) &
    rc = IOerr(lun_err, err_alloc, 'init_block_linop', &
    ' type block_linop member linops,'//&
    ' nrow_vec, ncol_vec'//&
    ' block_structure scrt_nrow, scr_ncol')

    this%linop_list          = linops
    this%block_structure = block_structure

    !
    ! check dimensions consistency and
    ! define blocks' dimensions
    !
    this%nrow_vec=0
    this%ncol_vec=0
    do iterm = 1, nnzblock
       ilinop = block_structure(1,iterm)
       irow   = block_structure(2,iterm)
       icol   = block_structure(3,iterm)

       !
       ! this code is left for debug porpuse
       ! write(0,*) ilinop,irow, icol
       ! call linops(ilinop)%linop%info(6)

       !
       ! nrow assignment and checks
       !      
       nrow_loc = linops(ilinop)%linop%nrow
       if (  nrow_loc .le. 0 ) then
          write(msg,*) 'matrix ', ilinop, &
               'at block ', irow, icol,&
               'has nrow .le. 0 '       
          rc = IOerr(lun_err, err_val, 'init_block_linop', &
               etb(msg))
       end if
          

       if ( this%nrow_vec(irow) .eq. 0 ) then
          ! set row dimension
          this%nrow_vec(irow) = nrow_loc
       else
          ! check row dimension match
          if (this%nrow_vec(irow) .ne. nrow_loc) then
             rc = IOerr(lun_err, err_val, 'init_block_linop', &
                  '  not consistent dimensions for rows value')
          end if
       end if
       
       !
       ! ncol assignment and checks
       !
       ncol_loc = linops(ilinop)%linop%ncol       
       if (  ncol_loc .le. 0 ) then
          write(msg,*) 'matrix ', ilinop, &
               'at block ', irow, icol,&
               'has ncol .le. 0 '       
          rc = IOerr(lun_err, err_val, 'init_block_linop', &
               etb(msg))
       end if
       
       if ( this%ncol_vec(icol) .eq. 0 ) then
          ! first assigment
          this%ncol_vec(icol) = ncol_loc
       else
          if (this%ncol_vec(icol) .ne. ncol_loc) then
             ! dimension match
             rc = IOerr(lun_err, err_val, 'init_block_linop', &
                  '  not consistent dimensions for rows value')
          end if
       end  if
    end do

    

    !
    ! set minimal operator properties
    !
    this%nrow         = sum(this%nrow_vec)
    this%ncol         = sum(this%ncol_vec)
    if(present (is_symmetric)) this%is_symmetric = is_symmetric

    !
    ! set real scaling factors
    !
    if (present(alphas)) then
       this%alphas = alphas
    else
       this%alphas = one
    end if
    
    
    allocate(&
         this%scr_nrow(this%nrow),&
         this%scr_ncol(this%ncol),&
         stat=res)
    if ( res .ne. 0) &
         rc = IOerr(lun_err, err_alloc, 'init_block_linop', &
         ' type block_linop member linops, nrow_vec, ncol_vec'//&
         ' block_structure scrt_nrow, scr_ncol')

    this%is_initialized =.True.


  end subroutine init_block_linop


  !>-----------------------------------------------------------------------
  !> Structure variable for definition of Saddle point Matrix
  !> M= ( A   B1^T )
  !>    ( B2 -C    )
  !> with A     := Real matrix of dimension n x n 
  !>      B1,B2 := Real matrix of dimension m x n (constraint matrix) 
  !>      C     := Real matrix of dimension m x m
  !> We follow the notation of 
  !> BMGHLJ05: Benzi, M. Golub H. and Liesen J.
  !> "Numerical solution of saddle point problems"
  !> Acta Numerica(2005) 
  !>------------------------------------------------------------------------


  !>-------------------------------------------------------------
  !> Static constructor.
  !> (procedure public for type stiffprec)
  !> Instantiate and initiliaze variable of type saddleprec
  !>
  !> usage:
  !>     call 'var'%init(lun_err,)
  !> where:
  !> \param[in] lun_err    -> integer. Error logical unit
  !> \param[in] matrix_A   -> class(abs_linop) Block 1,1
  !> \param[in] matrix_B1T -> class(abs_linop) Block 1,2
  !> \param[in] matrix_B2  -> class(abs_linop) Block 2,1
  !> \param[in] B1equalB2  -> logical. Flag to set if B1=B2
  !> \param[in,optional]
  !>            matrix_C   -> class(abs_linop) minus Block 2,2
  !<-------------------------------------------------------------
  subroutine init_saddle_point(this,&
       lun_err,&
       matrixA, matrixB1T,matrixB2,B1equalB2,&
       matrixC)
    use Globals
    implicit none
    class(block_linop), target, intent(inout) :: this
    integer,                  intent(in   ) :: lun_err
    class(abs_linop), target, intent(in   ) :: matrixA
    class(abs_linop), target, intent(in   ) :: matrixB1T
    class(abs_linop), target, intent(in   ) :: matrixB2
    logical,                  intent(in   ) :: B1equalB2
    class(abs_linop), target, optional,intent(in   ) :: matrixC

    ! List of matrices
    type(array_linop),allocatable :: components(:)
    ! Saddle point structure
    integer :: nnz
    integer :: block_structure(3,4)
    real(kind=double) :: alphas(4)
    logical :: is_sym
    
    if(this%is_initialized) call this%kill(lun_err)

    this%B1equalB2 =  B1equalB2

    allocate(components(4))

    block_structure(:,1) = (/1,1,1/)
    components(1)%linop => matrixA
    alphas(1)=one

    block_structure(:,2) = (/2,1,2/)
    components(2)%linop=> matrixB1T
    alphas(2)=one
    
    block_structure(:,3) = (/3,2,1/)
    components(3)%linop=> matrixB2
    alphas(3)=one

    if ( present (matrixC)) then
       nnz=4
       block_structure(:,4) = (/4,2,2/)
       components(4)%linop => matrixC
       alphas(4)=-one
    else
       nnz=3
    end if
    is_sym=B1equalB2.and.matrixA%is_symmetric .and. &
         components(4)%linop%is_symmetric


    !
    ! init block matrix at set it as saddle point matrix
    ! ( A  B1T )
    ! ( B2 -C  )
    !
    call this%init(lun_err, &
         nnz, components(1:nnz),&
         2, 2, &
         nnz, block_structure(1:3,1:nnz),&
         is_sym,&
         alphas= alphas)
    this%is_saddle_point=.True.

    this%is_initialized=.True.
    
    !
    ! free memory
    !
    deallocate(components)

  end subroutine init_saddle_point

 
  

  

  subroutine structure_transpose(this, transpose_block_structure)
    implicit none
    class(block_linop),          intent(in   ) :: this
    integer,                     intent(inout) :: transpose_block_structure(3,this%nnzblock)
    !local
    integer :: ilinop, iterm,irow,icol,i,j,indx,isgn,temp(3)

    do iterm = 1, this%nnzblock
       ilinop = this%block_structure(1,iterm)
       irow   = this%block_structure(2,iterm)
       icol   = this%block_structure(3,iterm)
       !
       ! swap column and row index
       !
       transpose_block_structure(1,iterm)= ilinop
       transpose_block_structure(2,iterm)= icol
       transpose_block_structure(3,iterm)= irow
    end do


    !
    ! sort in lexigraphic order
    !

    !  Initialize.
    i = 0
    indx = 0
    isgn = 0
    j = 0
    do 
       call global_heapsort(this%nnzblock, indx, i,j,isgn)
       if (indx .gt. 0 ) then
          ! SWAP ELEMENT 

          ! swap 
          temp(1:3)    = transpose_block_structure(1:3,i)
          transpose_block_structure(1:3,i) = &
               transpose_block_structure(1:3,j)
          transpose_block_structure(1:3,j) = temp(1:3)

       else if ( indx .lt. 0) then
          ! COMPARE
          isgn = 0
          if ( lexicographic_order( 2,&
               transpose_block_structure(2:3,i),&
               transpose_block_structure(2:3,j) ) ) then
             isgn = -1
          else
             isgn = 1
          end if
       else if ( indx .eq. 0 ) then
          exit
       end if
    end do
   

    
  end subroutine structure_transpose

  !>-------------------------------------------------------------
  !> Static constructor. Shorthand per init procedure
  !> for creating block diagonal matrix.
  !> Example
  !> M = ( A  0  0 )
  !>     ( 0 -B  0 )
  !>     ( 0  0  B )
  !>
  !> (procedure public for type block_linop)
  !> Instantiate variable of type block_linop
  !>
  !> usage:
  !>     call 'var'%init(lun_err, nlinops,linops [alphas])
  !>
  !> \param[in] lun_err           -> integer. Number of columns
  !> \param[in] nlinops           -> integer. Number of non-zero term
  !>                                     stored in the matrix
  !> \param[in] linops            -> type(array_linop). Dimensions = nlinops
  !>                                    List of operators on the diagonal
  !> \param[in,optinal] alphas     -> real.(dimensions=nnzblk)
  !>                                  Scalars premultiplynong each non
  !>                                  zeros lin. op. Default = 1.0
  !> Example of use for M:
  !>
  !> type(array_linop) :: linops(3)
  !> real(kind=double) :: alphas(3)
  !> type(block_linop) :: M
  !> linops(1)%linop => A
  !> linops(2)%linop => B
  !> linops(3)%linop => B
  !> alphas=(1,-1,1)
  !> call M%block_diagonal(0,3,linops,alphas)
  !>
  !<-------------------------------------------------------------
  subroutine init_block_diagonal(this, lun_err, &
       nlinops, linops,&
       alphas)
    use Globals
    implicit none
    !var
    class(block_linop),          intent(inout) :: this
    integer,                     intent(in   ) :: lun_err
    integer,                     intent(in   ) :: nlinops
    type(array_linop),           intent(in   ) :: linops(nlinops)
    real(kind=double), optional, intent(in   ) :: alphas(nlinops)
    !local
    logical :: rc
    integer :: i,res
    integer, allocatable :: block_structure(:,:)
    logical :: is_symmetric

    allocate(block_structure(3,nlinops),stat=res)
    if (res .ne. 0 ) &
         rc = IOerr(lun_err, err_alloc, 'block_diagonal', &
         'work array block_structure')

    !
    is_symmetric=.True.
    do i=1, nlinops
       is_symmetric = is_symmetric .and. linops(i)%linop%is_symmetric
       block_structure(1,i)=i
       block_structure(2,i)=i
       block_structure(3,i)=i
    end do
    
    call this%init(lun_err, &
       nlinops, linops,&
       nlinops, nlinops, &
       nlinops, block_structure, &
       is_symmetric,&
       alphas=alphas)

    deallocate(block_structure,stat=res)
    if (res .ne. 0 ) &
         rc = IOerr(lun_err, err_dealloc, 'block_diagonal', &
         ' work array block_structure')

    
  end subroutine init_block_diagonal


!!$  !>------------------------------------------------------
!!$  !> Subruotine to print linear operator properties
!!$  !> (public procedure for class adj_linop)
!!$  !> 
!!$  !> usage:
!!$  !>     call 'var'%info(lun)
!!$  !> where :
!!$  !>
!!$  !> \param[in   ] lun -> integer. I/O logical unit 
!!$  !<-------------------------------------------------------------
!!$  subroutine info_block_linop(this,lun)
!!$    implicit none
!!$    class(abs_linop),  intent(in   ) :: this
!!$    integer,           intent(in   ) :: lun
!!$    !local 
!!$    character(len=256) :: msg
!!$    logical :: all_named
!!$    integer :: ilinop,iterm,irow,col
!!$    character (len=256) :: out_row
!!$    character (len=256) :: format
!!$
!!$    all_named=.True.
!!$    do iterm=1,this%nnzblock 
!!$       ilinop = block_structure(1,iterm)
!!$       irow   = block_structure(2,iterm)
!!$       icol   = block_structure(3,iterm)
!!$
!!$       all_named = all_named .and. this%linop_list(ilinop)%name .ne. 'empty'
!!$    end do
!!$
!!$    if (all_named) then
!!$       max_length=0
!!$       do iterm=1,this%nnzblock
!!$          ilinop = block_structure(1,iterm)
!!$          irow   = block_structure(2,iterm)
!!$          icol   = block_structure(3,iterm)   
!!$          labels(iterm)=this%linop_list(ilinop)%name
!!$          max_length=max(max_length,len(etb(labels(iterm))))
!!$       end do
!!$       write(out_format,'(a,I1,a)')='(',max_length,',a)'
!!$
!!$       write(str0,out_format) '0',' '
!!$       do icol=1,this%block_ncol
!!$          row_out((icol-1)*max_length+1,icol*max_length) = str0
!!$       end do
!!$       
!!$       current_row=1
!!$       do iterm=1, this%nnzblock
!!$          irow   = block_structure(2,iterm)
!!$          icol   = block_structure(3,iterm)
!!$          if (irow.ne.current_row) then
!!$             write(str,'(a5,a1)') labels(iterm)),' '
!!$             row_out((icol-1)*max_length+1,icol*max_length) = str
!!$             write(lun,*) etb(row_out)
!!$             current_row=irow
!!$          else
!!$             
!!$          
!!$       do i=1, this%nblock_row-1
!!$          out_format=etb(etb(out_string)//'a,1x')
!!$       end do
!!$       out_format=etb(etb(out_format)//'a)')
!!$       
!!$    else
!!$       do iterm=1,this%nnzblock 
!!$          ilinop = block_structure(1,iterm)
!!$          irow   = block_structure(2,iterm)
!!$          icol   = block_structure(3,iterm)
!!$          write(labels(iterm),'(a,I0.3)')'M',ilinop
!!$       end do
!!$    end if
!!$
!!$    
!!$  end subroutine info_block_linop

    

  !>-------------------------------------------------------------
  !> Static destructor.
  !> (procedure public for type block_linop)
  !> Free memory and diassociate pointers
  !>
  !> usage:
  !>     call 'var'%init(lun_err)
  !>
  !> where:
  !> \param[in] lun_err               -> integer. Error logical unit
  !<-------------------------------------------------------------
  subroutine kill_block_linop(this, lun_err)
    use Globals
    implicit none
    !var
    class(block_linop),   intent(inout) :: this
    integer,           intent(in   ) :: lun_err
    ! local vars
    logical :: rc
    integer :: res

    this%nlinops      = 0
    this%nblock_row = 0
    this%nblock_col = 0
    this%nnzblock   = 0
    this%nrow = 0
    this%ncol = 0


    deallocate(&
         this%linop_list,&
         this%nrow_vec,&
         this%ncol_vec,&
         this%block_structure,&
         this%alphas,&
         stat=res)
    if (res .ne. 0) &
         rc = IOerr(lun_err, err_dealloc, 'kill_block_linop', &
         ' type block_linop member linops, nrow_vec, ncol_vec'//&
         ' block_structure scrt_nrow, scr_ncol')

    deallocate(&
         this%scr_nrow,&
         this%scr_ncol,&
         stat=res)
    if (res .ne. 0) &
    rc = IOerr(lun_err, err_dealloc, 'kill_blockslinop', &
    ' type block_linop member scr_nrow, scr_ncol')

    this%is_initialized = .False.

  end subroutine kill_block_linop

  
  

  !>-------------------------------------------------------------
  !> Concrete procedure procedure Mxv of parent type abs_linop.
  !> It performes the matrix-vector multiplication for type
  !> block_linop
  !> 
  !>         vec_out = (M) times (vec_in)
  !> (public procedure for class abs_linop)
  !> 
  !> usage:
  !>     call 'var'%matrix_times_vector(vec_in,vec_out,[info,lun_err])
  !>
  !> where 
  !> \param[in]          vec_in  -> real, dimension('var'%ncol)
  !>                                 vector to be multiplied
  !> \param[inout]       vec_out -> real, dimension('var'%nrow)
  !>                                 vector (M) times (vec_in) 
  !> \param[in,optional] info    -> integer. Info number
  !>                                 in case of error
  !> \param[in,optional] lun_err -> integer. Info number
  !>                                 in case of error  
  !<-------------------------------------------------------------
  recursive subroutine Mxv_block(this,vec_in,vec_out,info,lun_err)
       use Globals
       implicit none
       class(block_linop), intent(inout) :: this
       real(kind=double), intent(in   ) :: vec_in(this%ncol)
       real(kind=double), intent(inout) :: vec_out(this%nrow)
       integer,           intent(inout) :: info
       integer,           intent(in   ) :: lun_err

       !local
       logical :: rc
       integer :: iterm, ilinop, irow, icol
       character(len=1) :: itrans
       integer :: irowt, icolt
       integer :: in_begin,in_end, out_begin, out_end
       real(kind=double) :: dnrm2
       integer:: info_loc=0,lun_err_loc=0
       character(len=256) :: msg

       info = 0
       vec_out = zero
       this%scr_nrow = zero

       do iterm = 1, this%nnzblock
          ilinop   = this%block_structure(1,iterm)
          irow   = this%block_structure(2,iterm)
          icol   = this%block_structure(3,iterm)
                    
          in_begin  = sum(this%ncol_vec(1:icol)) - this%ncol_vec(icol) + 1
          in_end    = sum(this%ncol_vec(1:icol))

          out_begin = sum(this%nrow_vec(1:irow)) - this%nrow_vec(irow) + 1
          out_end   = sum(this%nrow_vec(1:irow))

          call this%linop_list(ilinop)%linop%matrix_times_vector(&
               vec_in(in_begin:in_end),&
               this%scr_nrow(out_begin:out_end),&
               info=info_loc,&
               lun_err=lun_err_loc)
          !
          ! check for errors
          !
          if (info_loc .ne. 0) then
             write(msg,*)  ' error in block (',irow,',',icol,')'
             rc = IOerr(lun_err, wrn_val, 'Mxv_block', &
                  etb(msg))
             info=info_loc
             return
          end if

          ! sum contribution
          call daxpy(&
               this%nrow_vec(irow),&
               this%alphas(iterm),&
               this%scr_nrow(out_begin:out_end),1,&
               vec_out(out_begin:out_end),1)
               
          
       end do
          
     end subroutine Mxv_block


  

end module LinearOperator
