module CombinedSparseMatrix
  use Globals
  use ScalableMatrix
  use SparseMatrix
  implicit none
  private 
  !> Structure variable containing the variables and the procedures
  !> to handle with a composed sparse matrix of the form
  !> 
  !>       M = \alpha A + \beta B^T D C
  !>
  !> where:
  !>   alpha,beta = real scalars
  !>   A          = sparse matrix
  !>   B,C        = with same sparsisty pattern
  !>   D          = diagonal matrix
  type, extends(scalable_mat), public :: combspmat
     !> Alpha coefficient
     real(kind=double) :: alpha=0.0d0
     !> Beta coefficient
     real(kind=double) :: beta=0.0d0
     !> Matrix A
     type(spmat), pointer :: A_matrix => null()
     !> Matrix B
     type(spmat), pointer :: B_matrix => null()
     ! Logical flag for the case B_amtrix=C_matrix
     logical :: BequalC=.false.
     !> Matrix C
     type(spmat), pointer :: C_matrix => null()
     !> Dimension = C_matrix%nrow 
     !> Diagonal Matrix    
     real(kind=double), pointer :: D_matrix(:) => null()
     !> Dimension (A_matrix%ncol).
     !> Working array
     real(kind=double), allocatable :: scr_ncol_A(:)
     !> Dimension (C_matrix%nrow).
     !> Working array
     real(kind=double), allocatable :: scr_ncol_D(:)
   contains
     !> static constructor
     !> (procedure public for type combspmat)
     procedure, public, pass :: init => init_combspmat
     !> static constructor
     !> (procedure public for type combspmat)
     procedure, public, pass :: set => set_combspmat
     !> static destructor
     !> (procedure public for type combspmat)
     procedure, public, pass :: kill => kill_combspmat
     !> Info procedure.
     !> (public procedure for type combspmat)
     procedure, public, pass :: info => info_combspmat
     !> Procedure to operate the diagonal scaling
     !> of matrix M by the diagonal matrix D
     !> out_matrix = D M D
     !> (public procedure for type combspmat)
     procedure, public , pass :: diagonal_scale
     !> Matrix vector moltiplication
     !>         y = (M) * x )
     !> (public procedure for type combspmat)
     procedure, public, pass :: matrix_times_vector => Mxv_combspmat
     !> Matrix vector moltiplication
     !>         y = (M^T) times (x)
     !> (public procedure for type combspmat)
     procedure, public, pass :: matrix_transpose_times_vector => MTxv_combspmat
     !> Get the diagonal of the matrix
     !> (public procedure for type combspmat)
     procedure, public, pass :: get_diagonal
  end type combspmat
contains
  !>-------------------------------------------------------------
  !> Static constructor.
  !> (procedure public for type spmat)
  !> Instantiate (allocate if necessary)
  !> and initilize (by also reading from input file)
  !> variable of type spmat
  !>
  !> usage:
  !>     call 'var'%init(lun_err, nrow, nterm )
  !>
  !> where:
  !> \param[in] lun_err            -> integer. Error logical unit
  !> \param[in] nrow               -> integer. Number of rows 
  !> \param[in] nrow               -> integer. Number of columns
  !> \param[in] nterm              -> integer. Number of non-zero term
  !>                                     stored in the matrix
  !> \param[in] (optional) dim_ker -> integer. Number of non-zero term
  !<-------------------------------------------------------------
  subroutine init_combspmat(this, &
       lun_err,&
       BequalC,&
       A_matrix,&
       B_matrix,&
       C_matrix,&
       D_matrix)
    use Globals
    use SparseMatrix
    implicit none
    !var
    class(combspmat),            intent(inout) :: this
    integer,                     intent(in   ) :: lun_err
    logical,                     intent(in   ) :: BequalC
    type(spmat),                 intent(in   ) :: A_matrix  
    type(spmat),                 intent(in   ) :: B_matrix
    type(spmat),                 intent(in   ) :: C_matrix
    real(kind=double),           intent(in   ) :: D_matrix(C_matrix%nrow) 
    
    ! local vars
    integer :: res
    logical :: rc
    integer :: i
    real(kind=double) :: min,max

    !  definition of abs_matrix quantities
    this%nrow          = A_matrix%nrow
    this%ncol          = A_matrix%ncol
    this%is_symmetric  = (A_matrix%is_symmetric) .and. (BequalC)
    this%triangular    = 'N'

    !
    ! allocation work array
    !
    if (.not. allocated(this%scr_ncol_A)  ) then
       allocate(this%scr_ncol_A(this%ncol),stat=res)
       if(res .ne. 0) rc = IOerr(lun_err, err_alloc, &
            'init_combspmat', &
            'type combspmat member scr_ncol_A',res)
    end if


    if (.not. allocated(this%scr_ncol_D)  ) then
       allocate(this%scr_ncol_D(C_matrix%nrow),stat=res)
       if(res .ne. 0) rc = IOerr(lun_err, err_alloc, &
            'init_combspmat', &
            'type combspmat member scr_ncol_D',res)
    end if

    
    
  end subroutine init_combspmat

  !>-------------------------------------------------------------
  !> Static constructor.
  !> (procedure public for type spmat)
  !> Instantiate (allocate if necessary)
  !> and initilize (by also reading from input file)
  !> variable of type spmat
  !>
  !> usage:
  !>     call 'var'%set(lun_err, nrow, nterm )
  !>
  !> where:
  !> \param[in] lun_err            -> integer. Error logical unit
  !> \param[in] nrow               -> integer. Number of rows 
  !> \param[in] nrow               -> integer. Number of columns
  !> \param[in] nterm              -> integer. Number of non-zero term
  !>                                     stored in the matrix
  !> \param[in] (optional) dim_ker -> integer. Number of non-zero term
  !<-------------------------------------------------------------
  subroutine set_combspmat(this, &
       lun_err,&
       BequalC,&
       alpha,&
       beta,&
       A_matrix,&
       B_matrix,&
       C_matrix,&
       D_matrix)
    use Globals
    use SparseMatrix
    implicit none
    !var
    class(combspmat),            intent(inout) :: this
    integer,                     intent(in   ) :: lun_err
    logical,                     intent(in   ) :: BequalC
    real(kind=double),           intent(in   ) :: alpha
    real(kind=double),           intent(in   ) :: beta
    type(spmat),       target,   intent(in   ) :: A_matrix  
    type(spmat),       target,   intent(in   ) :: B_matrix
    type(spmat),       target,   intent(in   ) :: C_matrix
    real(kind=double), target,   intent(in   ) :: D_matrix(C_matrix%nrow)
    ! local vars
    integer :: res
    logical :: rc
    integer :: i
    real(kind=double) :: min,max

    ! defition of type combspmat 
    this%alpha    = alpha
    this%beta     = beta
    this%BequalC  = BequalC

    this%A_matrix => A_matrix
    this%B_matrix => B_matrix
    this%C_matrix => C_matrix
    this%D_matrix => D_matrix(1:C_matrix%nrow)
    
  end subroutine set_combspmat


  !>-------------------------------------------------------------
  !> Static destructor.
  !> (procedure public for type spmat)
  !> deallocate all arrays for a var of type spmat
  !>
  !> usage:
  !>     call 'var'%kill(lun)
  !>
  !> where:
  !> \param[in] lun -> integer. I/O unit for error message output
  !<-----------------------------------------------------------
  subroutine kill_combspmat(this, lun_err)
    implicit none
    ! vars
    class(combspmat) :: this
    integer, intent(in) :: lun_err
    ! local vars
    integer :: res
    logical :: rc

    this%beta     = zero
    this%alpha    = zero
    this%A_matrix => null()
    this%B_matrix => null()
    this%C_matrix => null()
    this%D_matrix => null()

    !
    ! allocation work array
    !
    if ( allocated(this%scr_ncol_A)  ) then
       deallocate(this%scr_ncol_A,stat=res)
       if(res .ne. 0) rc = IOerr(lun_err, err_dealloc, &
            'init_combspmat', &
            'type combspmat member scr_ncol_A',res)
    end if

    if ( allocated(this%scr_ncol_D)  ) then
       deallocate(this%scr_ncol_D,stat=res)
       if(res .ne. 0) rc = IOerr(lun_err, err_dealloc, &
            'init_combspmat', &
            'type combspmat member scr_ncol_D',res)
    end if
    
  end subroutine kill_combspmat

  !>-------------------------------------------------------------
  !> Info procedure.
  !> (public procedure for type combspmat)
  !> Prints content of a variable of type combspmat
  !> 
  !> usage:
  !>     call 'var'%info(lun_out)
  !>
  !> where:
  !> \param[in] lun_out -> integer. I/O unit for output message
  !<-------------------------------------------------------------
  subroutine info_combspmat(this, lun)
    use Globals
    implicit none
    class(combspmat), intent(in) :: this
    integer, intent(in) :: lun

    write(lun,'(a,1e12.5)') 'alpha = ', this%alpha
    write(lun,'(a,1e12.5)') 'beta = ', this%beta
    write(lun,'(a)') 'Matrix_A '
    call this%A_matrix%info(lun)

    write(lun,'(a)') 'Matrix_B '
    call this%B_matrix%info(lun)

    write(lun,'(a)') 'Matrix_C '
    call this%C_matrix%info(lun)

  end subroutine info_combspmat



  !>-------------------------------------------------------------
  !> Procedure to compute Matrix vector product
  !>         vec_out = (matrix)  times (vec_in)
  !> (public procedure for type combspmat)
  !> 
  !> usage:
  !>     call 'var'%matrix_times_vector(vec_in,vec_out,info)
  !>
  !> where 
  !> \param[in   ] vec_in  -> real(dimension="var"%ncol). Array to be multiplied
  !> \param[inout] vec_out -> real(dimension="var"%ncol). Matrix times v
  !> \param[inout] info     -> integer. Info variable (useless in this case) 
  !>                                    It differs from  0 in case of error
  !<-------------------------------------------------------------
  subroutine Mxv_combspmat(this,vec_in, vec_out, info,lun_err)
    use Globals
    implicit none
    class(combspmat),      intent(inout) :: this
    real(kind=double), intent(in   ) :: vec_in(this%ncol)
    real(kind=double), intent(inout) :: vec_out(this%nrow)
    integer,           intent(inout) :: info
    integer,           intent(in   ) :: lun_err


    ! x  = vec_in
    ! y1 = A x
    ! y2 = B^T D C x
    ! y  = vec_out = alpha y1 + beta y2 
    vec_out=zero

    ! y1 = alpha A x
    call this%A_matrix%matrix_times_vector(vec_in,vec_out,info,lun_err)
    if (abs( this%alpha-one) .gt. small) then
       call dscal(this%A_matrix%nrow,this%alpha, vec_out,1)
    end  if
    
    ! compute y2
    ! y2 = B^T D C x
    call this%C_matrix%matrix_times_vector(vec_in,this%scr_ncol_D,info,lun_err)
    this%scr_ncol_D = this%scr_ncol_D * this%D_matrix
    call this%B_matrix%matrix_transpose_times_vector(this%scr_ncol_D,this%scr_ncol_A,info,lun_err)
    
    
    ! y = y1 + beta y2 
    call daxpy(this%ncol, this%beta, this%scr_ncol_A, 1 , vec_out, 1)  

  end subroutine Mxv_combspmat

  !>-------------------------------------------------------------
  !> Procedure to compute Matrix vector product
  !>         vec_out = (matrix)^T  times (vec_in)
  !> (public procedure for type combspmat)
  !> 
  !> usage:
  !>     call 'var'%matrix_times_vector(vec_in,vec_out,info)
  !>
  !> where 
  !> \param[in   ] vec_in  -> real(dimension="var"%ncol). Array to be multiplied
  !> \param[inout] vec_out -> real(dimension="var"%ncol). Matrix times v
  !> \param[inout] info     -> integer. Info variable (useless in this case) 
  !>                                    It differs from  0 in case of error
  !<-------------------------------------------------------------
  subroutine MTxv_combspmat(this,vec_in, vec_out, info,lun_err)
    use Globals
    implicit none
    class(combspmat),  intent(inout) :: this
    real(kind=double), intent(in   ) :: vec_in(this%nrow)
    real(kind=double), intent(inout) :: vec_out(this%ncol)
    integer,           intent(inout) :: info
    integer,           intent(in   ) :: lun_err
        

    ! x  = vec_in
    ! y1 = A^T x
    ! y2 = C^T D B x
    ! y  = vec_out = alpha y1 + beta y2 
    vec_out=zero

    ! y1 = alpha A^T x
    call this%A_matrix%matrix_transpose_times_vector(vec_in,vec_out,info,lun_err)
    if (abs( this%alpha-one) .gt. small) then
       call dscal(this%A_matrix%nrow,this%alpha, vec_out,1)
    end  if
    
    ! compute y2
    ! y2 = C^T D B x
    call this%B_matrix%matrix_times_vector(vec_in,this%scr_ncol_D,info,lun_err)
    this%scr_ncol_D = this%scr_ncol_D * this%D_matrix
    call this%C_matrix%matrix_transpose_times_vector(this%scr_ncol_D,this%scr_ncol_A,info,lun_err)


    ! y = y1 + beta y2 
    call daxpy(this%ncol, this%beta, this%scr_ncol_A, 1 , vec_out, 1)  

  end subroutine MTxv_combspmat



  !>-------------------------------------------------------
  !> Procedure to operate diagonal scaling of spmat
  !>     M => D M D 
  !> where D(M) is a diagonal matrix
  !> (pubblic procedure for type spmat, operate in place)
  !> 
  !> usage: var%MxD(diag_matrix)
  !>
  !> where: 
  !>\param[in] diag_matrix -> real(nrow). Diagonal mat. values
  !<-------------------------------------------------------
  subroutine diagonal_scale(this, lun_err, diagonal)
    use Globals
    implicit none
    class(combspmat),  intent(inout) :: this
    integer,           intent(in   ) :: lun_err
    real(kind=double), intent(in   ) :: diagonal(this%ncol)
    !local
    logical rc
    integer :: i
    
    if ( this%nrow .eq. this%ncol ) then
       rc = IOerr(lun_err, err_inp, &
            'diagonal_scale', &
            'combspmat is not squared ')
       return
    end if
    
    call this%A_matrix%diagonal_scale(lun_err,diagonal)
    call this%B_matrix%MxD(lun_err,diagonal )


    if ( .not. this%BequalC ) then
       call this%C_matrix%MxD(lun_err,diagonal )
    end if

  end subroutine diagonal_scale


  subroutine get_diagonal(this,diagonal)
    use Globals
    class(combspmat),  intent(in   ) :: this
    real(kind=double), intent(inout) :: diagonal(min(this%nrow,this%ncol))
    !local 
    integer :: i,j,start,finish,ind
    
    ! 
    ! diag=alpha diag(A)+beta diag(B^T D C)
    !


    ! 
    ! diag=alpha diag(A)
    !
    call this%A_matrix%get_diagonal(diagonal)
    if (abs( this%alpha-one) .gt. small) then
       call dscal(this%A_matrix%nrow,this%alpha, diagonal,1)
    end  if

    ! 
    ! diag= diag + beta diag(B^T D C)
    !
!!$    do i=1,this%C_matrix%nrow
!!$       start  = this%C_matrix%ia(i)
!!$       finish = this%C_matrix%ia(i+1)-1
!!$       do j = start, finish
!!$          ind = this%C_matrix%ja(j)
!!$          diagonal(ind) = diagonal(ind) + &
!!$               this%beta * (&
!!$               this%B_matrix%coeff(j) * &
!!$               this%D_matrix(i)     * &
!!$               this%C_matrix%coeff(j)   )
!!$       end do
!!$    end do     
    
  end subroutine get_diagonal

end module CombinedSparseMatrix


