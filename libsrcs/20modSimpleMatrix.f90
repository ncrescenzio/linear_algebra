module SimpleMatrix
  use Globals
  use LinearOperator
  use Matrix
  implicit none
  private
  public :: diag, vector2matrix

  !>-------------------------------------------------
  !> Structure variable for Real Diagonal Matrices
  !>------------------------------------------------
  type, extends(abs_matrix), public :: diagmat
     !> Vector y = prec_zero^{-1} uvec 
     !> Dimension(prec_zero%nequ)
     real(kind=double), allocatable :: diagonal(:)
   contains
     !> static constructor
     !> (procedure public for type diagmat)
     procedure, public, pass :: init => init_diagmat
     !> static constructor
     !> (procedure public for type diagmat)
     procedure, public, pass :: set => set_diagmat
     !> static destructor
     !> (procedure public for type diagmat)
     procedure, public, pass :: kill => kill_diagmat
     !> writing procedure
     !> (procedure public for type diagmat)
     procedure, public, pass :: write => write_diagmat
     !> Procedure to compute 
     !>         y = M * x 
     !> (public procedure for type blockmat)
     procedure, public, pass :: matrix_times_vector => Mxv_diagmat
     !> Procedure to compute 
     !>         y = M^T * x 
     !> with M^T the transposed of a matrix M
     !> (public procedure for type blockmat)
     procedure, public, pass :: matrix_transpose_times_vector  => MTxv_diagmat
     !> Fucntion to convert a vector into a diagonal matrix
     !> It act like init + set 
     !> (public procedure for type blockmat)
     procedure, public, nopass :: vector2matrix
  end type diagmat
  type, extends(abs_matrix), public :: scalmat
     !> Scalar on diagonal 
     !> Dimension(prec_zero%nequ)
     real(kind=double) :: scalar=zero
   contains
     !> static constructor
     !> (procedure public for type diagmat)
     procedure, public, pass :: init => init_scalmat
     !> static constructor
     !> (procedure public for type diagmat)
     procedure, public, pass :: eye => init_eye
     !> static constructor
     !> (procedure public for type diagmat)
     procedure, public, pass :: zeros => init_zeros
     !> static destructor
     !> (procedure public for type diagmat)
     procedure, public, pass :: kill => kill_scalmat
     !> Procedure to compute 
     !>         y = M * x 
     !> (public procedure for type blockmat)
     procedure, public, pass :: matrix_times_vector => Mxv_scalmat
     !> Procedure to compute 
     !>         y = M^T * x 
     !> with M^T the transposed of a matrix M
     !> (public procedure for type blockmat)
     procedure, public, pass :: matrix_transpose_times_vector => MTxv_scalmat
     !> Fucntion to convert a scalar into matrix
     !> It act like init + set 
     !> (public procedure for type blockmat)
     !procedure, public, nopass :: scalar2matrix
  end type scalmat

  type, extends(abs_linop), public :: rankk_mat
     integer :: nrankmax
     integer :: nrank
     !> Scalar on diagonal 
     !> Dimension(prec_zero%nequ)
     real(kind=double), allocatable :: Uvectors(:,:)
     !> Scalar on diagonal 
     !> Dimension(prec_zero%nequ)
     real(kind=double),allocatable :: Vvectors(:,:)
     
   contains
     !> static constructor
     !> (procedure public for type diagmat)
     procedure, public, pass :: init => init_rankk_mat
     !> static destructor
     !> (procedure public for type diagmat)
     procedure, public, pass :: kill => kill_rankk_mat
     !> static destructor
     !> (procedure public for type diagmat)
     procedure, public, pass :: set=> set_rankk_mat
     
     !> Procedure to compute 
     !>         y = M * x 
     !> (public procedure for type blockmat)
     procedure, public, pass :: matrix_times_vector => apply_rankk_mat
  end type rankk_mat

  !>-------------------------------------------------
  !> Structure variable for Identity matrix
  !>------------------------------------------------
  type, extends(abs_matrix), public :: permutation
     !> Permutation
     integer,allocatable :: perm(:)
     !> Inverse Permutaion
     integer,allocatable :: inverse_permutation(:)
   contains
     !> Static constructor 
     !> (procedure public for type permutation)
     procedure, public, pass :: init => init_permutation
     !> Static destructor
     !> (procedure public for type permutation)
     procedure, public, pass :: kill => kill_permutation
     !> Compute matrix times vector operatoration
     !> (public procedure for type permutation)
     procedure, public,  pass :: matrix_times_vector => apply_permutation
     !> Compute matrix transpose time vector operation
     !> (public procedure for type permutation)
     procedure, public,  pass :: matrix_transpose_times_vector => apply_permutation_transpose
  end type permutation

  interface diag
     module procedure vector2matrix
  end interface diag

contains
  !>-------------------------------------------------------------
  !> Static constructor.
  !> (procedure public for type scalmat)
  !>
  !> usage:
  !>     call 'var'%init( nrow, [ncol] )
  !>
  !> where:
  !> \param[in] nrow          -> integer. Number of rows
  !> \param[in] scalar        -> real. Scalar on diagonal
  !> \param[in,optional] ncol -> integer. Number of columns.
  !>                             If not passed ncol=nrow
  !<-------------------------------------------------------------
  subroutine init_scalmat(this,nrow,scalar,ncol)
     use Globals
     implicit none
     class(scalmat),     intent(inout) :: this
     integer,            intent(in   ) :: nrow
     real(kind=double),  intent(in   ) :: scalar
     integer, optional,  intent(in   ) :: ncol
     !local
     logical :: rc
     integer :: res

     this%nrow = nrow
     if ( present(ncol) )  then
        this%ncol = ncol
        this%is_symmetric = ncol.eq.nrow
     else
        this%ncol = nrow
        this%is_symmetric = .True.
     end if
     this%scalar = scalar
     if ( abs(scalar) < small )then
        this%name='0'
     else if ( abs(scalar-one) < small ) then
        this%name='I'
     else
        write(this%name,'(1pe8.2,a)') scalar,'I'
     end if
        
   end subroutine init_scalmat

   !>-------------------------------------------------------------
   !> Static constructor.
   !> (procedure public for type scalmat)
   !>
   !> usage:
   !>     call 'var'%eye( nrow, [ncol] )
   !>
   !> where:
   !> \param[in] nrow          -> integer. Number of rows
   !> \param[in,optional] ncol -> integer. Number of columns.
   !>                             If not passed ncol=nrow
   !<-------------------------------------------------------------
   subroutine init_eye(this,nrow,ncol)
     use Globals
     implicit none
     class(scalmat),     intent(inout) :: this
     integer,            intent(in   ) :: nrow
     integer, optional,  intent(in   ) :: ncol

     call this%init(nrow,one,ncol)

   end subroutine init_eye

    !>-------------------------------------------------------------
   !> Static constructor.
   !> (procedure public for type scalmat)
   !>
   !> usage:
   !>     call 'var'%eye( nrow, [ncol] )
   !>
   !> where:
   !> \param[in] nrow          -> integer. Number of rows
   !> \param[in,optional] ncol -> integer. Number of columns.
   !>                             If not passed ncol=nrow
   !<-------------------------------------------------------------
   subroutine init_zeros(this,nrow,ncol)
     use Globals
     implicit none
     class(scalmat),     intent(inout) :: this
     integer,            intent(in   ) :: nrow
     integer, optional,  intent(in   ) :: ncol

     call this%init(nrow,zero,ncol)

   end subroutine init_zeros

   !>-------------------------------------------------------------
   !> Static destructructor.
   !> (procedure public for type eye)
   !>
   !> usage:
   !>     call 'var'%kill()
   !<-------------------------------------------------------------
   subroutine kill_scalmat(this)
     use Globals
     implicit none
     class(scalmat),  intent(inout) :: this

     !
     ! reset to default
     !
     call this%to_default()
     
   end subroutine kill_scalmat

   !>-------------------------------------------------------------
   !> Matrix-vector multiplication procedure
   !>         vec_out = (M) times (vec_in)
   !> (public procedure for class abs_matrix)
   !> 
   !> usage:
   !>     call 'var'%Mxv(vec_in,vec_out,[info,lun_err])
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
   subroutine Mxv_scalmat(this,vec_in,vec_out,info,lun_err)
     use Globals
     implicit none
     class(scalmat), intent(inout) :: this
     real(kind=double), intent(in   ) :: vec_in(this%ncol)
     real(kind=double), intent(inout) :: vec_out(this%nrow)
     integer,           intent(inout) :: info
     integer,           intent(in   ) :: lun_err
     !local
     integer :: mindim

     mindim=min(this%ncol,this%nrow)
     
     vec_out = zero
     call daxpy(mindim,this%scalar,vec_in(1:mindim),1,vec_out(1:mindim),1)

   end subroutine Mxv_scalmat

   !>-------------------------------------------------------------
   !> Matrix transpose times vector multiplication procedure
   !>         vec_out = (M) times (vec_in)
   !> (public procedure for class abs_matrix)
   !> 
   !> usage:
   !>     call 'var'%MTxv(vec_in,vec_out,[info,lun_err])
   !>
   !> where 
   !> \param[in   ] vec_in  -> real, dimension('var'%ncol)
   !>                                   vector to be multiplied
   !> \param[inout] vec_out -> real, dimension('var'%nrow)
   !>                                    vector (M) times (vec_in) 
   !> \param[in   ] info    -> integer. Info number
   !>                                    in case of error 
   !> \param[in   ] lun_err -> integer. Info number
   !>                                    in case of error 
   !<-------------------------------------------------------------
   subroutine MTxv_scalmat(this,vec_in,vec_out,info,lun_err)
     use Globals
     implicit none
     class(scalmat),    intent(inout) :: this
     real(kind=double), intent(in   ) :: vec_in(this%ncol)
     real(kind=double), intent(inout) :: vec_out(this%nrow)
     integer,           intent(inout) :: info
     integer,           intent(in   ) :: lun_err
     !local
     integer :: mindim

     mindim=min(this%ncol,this%nrow)
     
     vec_out = zero
     call daxpy(mindim,this%scalar,vec_in(1:mindim),1,vec_out(1:mindim),1)
     info=0
     
   end subroutine MTxv_scalmat




  
  function vector2matrix(coeff) result (this)
    use Globals
    implicit none
    real(kind=double), intent(in  ) :: coeff(:)
    type(diagmat) :: this
    !local
    logical :: rc
    integer :: nequ

    nequ=size(coeff)
    if ( this%is_initialized ) then
       if ( this%nrow .ne. nequ) then
          call this%kill(6)
          call this%init(6, nequ)
       end if
    else
       call this%init(6, nequ)
    end if
    call this%set(6,coeff)

  end function vector2matrix


  !>-------------------------------------------------------------
  !> Static constructor.
  !> (procedure public for type diagmat)
  !>
  !> usage:
  !>     call 'var'%init( lun_err,nrow, [ncol] )
  !>
  !> where:
  !> \param[in] lun_err       -> integer. IO unit for error msg 
  !> \param[in] nrow          -> integer. Number of rows 
  !> \param[in,optional] ncol -> integer. Number of columns.
  !>                             If not passed ncol=nrow
  !<-------------------------------------------------------------
  subroutine init_diagmat(this,lun_err,nrow,ncol)
    use Globals
    implicit none
    class(diagmat),     intent(inout) :: this
    integer,            intent(in   ) :: lun_err
    integer,            intent(in   ) :: nrow
    integer,optional,   intent(in   ) :: ncol

    !local
    logical :: rc
    integer :: res,i
    integer, allocatable :: iwork(:)

    if ( this%is_initialized ) call this%kill(lun_err)
    
    ! set properties
    this%nrow = nrow
    if ( present(ncol) )  then
       this%ncol = ncol
       this%is_symmetric = ncol.eq.nrow
    else
       this%ncol = nrow
       this%is_symmetric = .True.
    end if

    ! copy coefficients
    allocate(this%diagonal(min(this%nrow,this%ncol)),stat=res)
    if(res .ne. 0)&
         rc = IOerr(lun_err, err_alloc, 'init_diagmat', &
         ' type diagmat member diagonal',res)

     this%is_initialized=.True.
    
  end subroutine init_diagmat


  !>-------------------------------------------------------------
  !> Procedure to set diagonal of diagmat
  !> (procedure public for type diagmat)
  !>
  !> usage:
  !>     call 'var'%set( diagonal )
  !>
  !> where:
  !> \param[in] diagonal -> real (dimension=min(nrow,ncol))
  !<-------------------------------------------------------------
  subroutine set_diagmat(this,lun_err,diagonal)
    use Globals
    implicit none
    class(diagmat),     intent(inout) :: this
    integer,            intent(in   ) :: lun_err
    real(kind=double),  intent(in   ) :: diagonal(min(this%nrow,this%ncol))
    !
    ! copy coefficients
    !
    this%diagonal(1:min(this%nrow,this%ncol))= diagonal(1:min(this%nrow,this%ncol))

  end subroutine set_diagmat

   !>--------------------------------------------------------
   !> Static destructructor.
   !> (procedure public for type diagmat)
   !>
   !> usage:
   !>     call 'var'%kill()
   !<--------------------------------------------------------
   subroutine kill_diagmat(this,lun_err)
     use Globals
     implicit none
     class(diagmat),  intent(inout) :: this
     integer,         intent(in   ) :: lun_err
     !local
     logical :: rc
     integer :: res

     !
     ! free memory
     !
     if ( this%is_initialized ) then
        deallocate(this%diagonal,stat=res)
        if(res .ne. 0) rc = IOerr(lun_err, &
             err_dealloc, 'kill_diagmat', &
             ' type diagmat member diagonal',res)
     end if

     !
     ! rest properties to default
     !
     call this%to_default()
   end subroutine kill_diagmat

   !>-------------------------------------------------------------
   !> Writing procedure.
   !> (procedure public for type diagmat)
   !>
   !> usage:
   !>     call 'var'%write(lun)
   !<-------------------------------------------------------------
   subroutine write_diagmat(this,lun)
     use Globals
     implicit none
     class(diagmat),  intent(in) :: this
     integer,         intent(in) :: lun
     !local
     integer :: inode

     write(lun,*) 'matrix coefficients' 
     do inode=1,this%nrow
        write(lun,'(i5,1pe15.6)') inode, this%diagonal(inode)
     end do
     
        
   end subroutine write_diagmat
   
   
   
   !>-------------------------------------------------------------
   !> Matrix-vector multiplication procedure
   !>         vec_out = (M) times (vec_in)
   !> (public procedure for class abs_matrix)
   !> 
   !> usage:
   !>     call 'var'%Mxv(vec_in,vec_out,[info])
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
   subroutine Mxv_diagmat(this,vec_in,vec_out,info,lun_err)
     use Globals
     implicit none
     class(diagmat), intent(inout) :: this
     real(kind=double), intent(in   ) :: vec_in(this%ncol)
     real(kind=double), intent(inout) :: vec_out(this%nrow)
     integer,           intent(inout) :: info
     integer,           intent(in   ) :: lun_err
     !local
     integer :: mindim

     info=0
     mindim=min(this%ncol,this%nrow)
     
     vec_out(1:mindim) = this%diagonal(1:mindim) * vec_in(1:mindim)

     if ( this%ncol < this%nrow ) vec_out(mindim+1:this%nrow) = zero
   end subroutine Mxv_diagmat

   !>-------------------------------------------------------------
   !> Matrix transpose times vector multiplication procedure
   !>         vec_out = (M) times (vec_in)
   !> (public procedure for class abs_matrix)
   !> 
   !> usage:
   !>     call 'var'%MTxv(vec_in,vec_out,[info,lun_err])
   !>
   !> where 
   !> \param[in   ] vec_in            -> real, dimension('var'%ncol)
   !>                                   vector to be multiplied
   !> \param[inout] vec_out           -> real, dimension('var'%nrow)
   !>                                    vector (M) times (vec_in) 
   !> \param[in   ] info    -> integer. Info number
   !>                                    in case of error 
   !> \param[in   ] lun_err -> integer. Info number
   !>                                    in case of error 
   !<-------------------------------------------------------------
   subroutine MTxv_diagmat(this,vec_in,vec_out,info,lun_err)
     use Globals
     implicit none
     class(diagmat), intent(inout) :: this
     real(kind=double), intent(in   ) :: vec_in(this%nrow)
     real(kind=double), intent(inout) :: vec_out(this%ncol)
     integer,           intent(inout) :: info
     integer,           intent(in   ) :: lun_err
     
     
     !local
     integer :: mindim

     mindim=min(this%ncol,this%nrow)
     
     vec_out(1:mindim) = this%diagonal(1:mindim) * vec_in(1:mindim)

     if ( this%ncol < this%nrow ) vec_out(mindim+1:this%nrow) = zero

     info=0

   end subroutine MTxv_diagmat

   
   

   subroutine init_rankk_mat(this,lun_err,len_Uvector,len_Vvector,nrankmax)
    use Globals
    implicit none
    class(rankk_mat),    intent(inout) :: this
    integer,         intent(in   ) :: lun_err
    integer,         intent(in   ) :: len_Uvector
    integer,         intent(in   ) :: len_Vvector
    integer,         intent(in   ) :: nrankmax
    !local
    logical :: rc
    integer :: res
    integer :: i


    this%nrow         = len_Uvector
    this%ncol         = len_Vvector
    this%nrankmax     = nrankmax
    
    allocate(&
         this%Uvectors(len_Uvector,this%nrankmax),&
         this%Vvectors(len_Vvector,this%nrankmax),&
         stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'init_rankk', &
         ' type rankk member updates work',res)
  end subroutine init_rankk_mat


  subroutine set_rankk_mat(this,lun_err,&
       nrank, Uvectors,Vvectors,is_symmetric)
    use Globals
    implicit none
    class(rankk_mat),        intent(inout) :: this
    integer,                 intent(in   ) :: lun_err
    integer,                 intent(in   ) :: nrank
    real(kind=double),       intent(in   ) :: Uvectors(this%nrow,nrank)
    real(kind=double),       intent(in   ) :: Vvectors(this%ncol,nrank)    
    logical,                 intent(in   ) :: is_symmetric
    !local
    logical :: rc
    integer :: res,i
    real(kind=double) :: ddot

    this%is_symmetric = is_symmetric
    
    this%nrank = nrank

    this%Uvectors(1:this%nrow, 1: nrank) = Uvectors  
    this%Vvectors(1:this%ncol, 1: nrank) = Vvectors 

  end subroutine set_rankk_mat


  subroutine kill_rankk_mat(this,lun_err)
    use Globals
    implicit none
    class(rankk_mat),     intent(inout) :: this
    integer,         intent(in   ) :: lun_err
    !local
    logical :: rc
    integer :: res 
    integer :: i

    deallocate(&
         this%Uvectors,&
         this%Vvectors,&
         stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_dealloc, 'init_rankk', &
         ' type rankk member updates work',res)

    !
    ! reset to default
    !
    call this%to_default()

  end subroutine kill_rankk_mat

  recursive subroutine apply_rankk_mat(this,vec_in,vec_out,info,lun_err)
    use Globals
    implicit none
    class(rankk_mat),  intent(inout) :: this
    real(kind=double), intent(in   ) :: vec_in(this%ncol)
    real(kind=double), intent(inout) :: vec_out(this%nrow)
    integer,           intent(inout) :: info
    integer,           intent(in   ) :: lun_err
    ! local
    integer :: i
    real(kind=double) :: ddot
    
    vec_out=zero
    do i = 1, this%nrank
       vec_out = vec_out + &
            ddot(this%ncol, &
            this%Vvectors(1:this%ncol,i),1,&
            vec_in(this%ncol),1) *&
            this%Uvectors(1:this%nrow,i)
    end do
    
  end subroutine apply_rankk_mat

  !>-------------------------------------------------------------
  !> Static constructor.
  !> (procedure public for type permuation)
  !>
  !> usage:
  !>     call 'var'%init( lun_err,nperm, perm )
  !>
  !> where:
  !> \param[in] nrow          -> integer. Number of rows 
  !> \param[in,optional] ncol -> integer. Number of columns.
  !>                             If not passed ncol=nrow
  !<-------------------------------------------------------------
  subroutine  init_permutation(this,lun_err, nperm, perm)
    implicit none
    class(permutation),  intent(inout) :: this
    integer,         intent(in   ) :: lun_err    
    integer,         intent(in   ) :: nperm
    integer,         intent(in   ) :: perm(nperm)
    !local
    logical :: rc
    integer :: res,i
    
    this%nrow = nperm
    this%ncol = nperm
    this%is_symmetric = .False.
    this%triangular = 'N'

    allocate(this%perm(nperm),this%inverse_permutation(nperm),stat=res)
    if(res .ne. 0)&
         rc = IOerr(lun_err, err_alloc, 'init_permutation', &
         ' type diagmat member diagonal',res)
    
    this%perm = perm
    do i = 1, nperm
       this%inverse_permutation(perm(i)) = i
    end do

    
    
    
  end subroutine init_permutation

  !>-------------------------------------------------------------
  !> Static destructructor.
  !> (procedure public for type eye)
  !>
  !> usage:
  !>     call 'var'%kill()
  !<-------------------------------------------------------------
  subroutine  kill_permutation(this,lun_err)
    implicit none
    class(permutation), intent(inout) :: this
    integer,            intent(in   ) :: lun_err
    !local
    logical :: rc
    integer :: res

    call this%to_default()
    deallocate(this%perm,this%inverse_permutation,stat=res)
    if(res .ne. 0)&
         rc = IOerr(lun_err, err_alloc, 'init_perm', &
         ' type diagmat member diagonal',res)

    !
    ! reset to default
    !
    call this%to_default()

    
  end subroutine kill_permutation

 
  !>-------------------------------------------------------------
  !> Matrix-vector multiplication procedure
  !>         vec_out = (M) times (vec_in)
  !> (public procedure for class abs_matrix)
  !> 
  !> usage:
  !>     call 'var'%Mxv(vec_in,vec_out,[info])
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
  subroutine apply_permutation(this,vec_in,vec_out,info,lun_err)
    use Globals
    class(permutation),   intent(inout) :: this
    real(kind=double), intent(in   ) :: vec_in(this%ncol)
    real(kind=double), intent(inout) :: vec_out(this%nrow)
    integer,           intent(inout) :: info
    integer,           intent(in   ) :: lun_err

    vec_out(this%perm) = vec_in 
    
    info=0
    
  end subroutine apply_permutation


  !>-------------------------------------------------------------
  !> Matrix transope-vector multiplication procedure
  !>         vec_out = (M) times (vec_in)
  !> (public procedure for class abs_matrix)
  !> 
  !> usage:
  !>     call 'var'%MTxv(vec_in,vec_out,[info])
  !>
  !> where 
  !> \param[in   ] vec_in  -> real, dimension('var'%ncol)
  !>                                   vector to be multiplied
  !> \param[inout] vec_out -> real, dimension('var'%nrow)
  !>                                    vector (M) times (vec_in) 
  !> \param[in   ] info    -> integer. Info number
  !>                                    in case of error 
  !> \param[in   ] lun_err -> integer. Info number
  !>                                    in case of error 
  !<-------------------------------------------------------------
  subroutine apply_permutation_transpose(this,vec_in,vec_out,info,lun_err)
    use Globals
    class(permutation), intent(inout) :: this
    real(kind=double),  intent(in   ) :: vec_in(this%nrow)
    real(kind=double),  intent(inout) :: vec_out(this%ncol)
    integer,            intent(inout) :: info
    integer,            intent(in   ) :: lun_err

    vec_out(this%inverse_permutation)  = vec_in
    info=0
    
  end subroutine apply_permutation_transpose
 end module SimpleMatrix
