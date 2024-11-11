module StiffnessMatrix
  use Globals
  use Matrix
  use SparseMatrix
  !use Graphs
  use ReadableMatrix
  use DenseMatrix
  implicit none
  private 
  !> Structure variable containing the variables and the procedures
  !> to handle with a composed sparse matrix of the form
  !> 
  !>       S(\Tdens) = A^T diag(\tdens) diag(inv_weigth) A + lasso* Id
  !>
  !> where:
  !>   A               = matrix
  !>   \tdens \weight  = positive vectors
  !>   lasso           = positive regularization
  type, extends(readable_mat), public :: stiffmat
     !> Sized of tdens-weight
     integer :: ntdens
     !> Matrix A
     class(abs_matrix), pointer :: matrix_A => null()
     !> Dimension = matrix_A%nrow 
     !> Diagonal Matrix    
     real(kind=double), allocatable :: tdens(:)
     !> Dimension = matrix_A%nrow 
     !> Diagonal Matrix    
     real(kind=double), allocatable :: sqrt_tdens(:)
     !> Dimension = matrix_A%nrow 
     !> Diagonal Matrix    
     real(kind=double), allocatable :: inv_weight(:)
     !> Dimension (matrix_A%ncol).
     !> Working array
     real(kind=double), allocatable :: scr_ncol_A(:)
     !> Dimension (matrix_A%nrow).
     !> Working array
     real(kind=double), allocatable :: scr_nrow_A(:)
     !> Dimension (matrix_A%nrow,matrix%ncol).
     !> Working matrix for assembly dense matrices
     real(kind=double), allocatable :: scr_tdens_A(:,:)
     !> Dimension (matrix_A%nrow,matrix%ncol).
     !> Working matrix for assembly dense matrices
     type(densemat) :: A_transpose
     !>-----------------------------------------------
     !> Varibles for explicit matrix A^T D A
     !>-----------------------------------------------
     !> Flag for assembled matrix
     logical :: is_assembled=.False.
     !> Flag for assembled matrix
     logical :: is_inexact=.False.
     !> Flag for assembled matrix
     class(abs_matrix), pointer :: assembled_matrix
     !> Working array
     type(densemat) :: dense
     !> Sparse matrix  case
     !> Working sparse matrix
     !> Use when matrix_A is of type
     !> 1-spmat 
     !> 2 adjmat
     type(spmat) :: sparse
     !> Assembler for sparse matrix
     !> Dimension()
     integer, allocatable :: assembler(:,:,:)
     !>-----------------------------------------------
     !> Column selection variables
     !> On/Off Flag for selection of rows
     logical :: select_tdens
     !> On/Off Flag for selection of rows
     integer :: nnz_tdens
     !> Indeces of rows of Matrix A to be computed
     integer, allocatable :: selection_tdens(:)
     !>---------------------------------------------
     real(kind=double) :: lasso=zero
   contains
     !> static constructor
     !> (procedure public for type stiffmat)
     procedure, public, pass :: init => init_stiffmat
      !> static constructor
     !> (procedure public for type stiffmat)
     procedure, public, pass :: set => set_stiffmat
     !> static destructor
     !> (procedure public for type stiffmat)
     procedure, public, pass :: kill => kill_stiffmat
     !> Info procedure.
     !> (public procedure for type stiffmat)
     procedure, public, pass :: info => info_stiffmat
     
     !> Matrix vector moltiplication
     !>         y = (M) * x )
     !> (public procedure for type stiffmat)
     procedure, public, pass :: matrix_times_vector => Mxv_stiffmat
     !> Matrix vector moltiplication
     !>         y = (M^T) times (x)
     !> (public procedure for type stiffmat)
     procedure, public, pass :: matrix_transpose_times_vector => MTxv_stiffmat
     !> Get the diagonal of the matrix
     !> (public procedure for type stiffmat)
     procedure, public, pass :: get_diagonal
     !> Assemble the matrix explitely 
     !> (public procedure for type stiffmat)
     procedure, public, pass :: assembly
!!$!> Procedure to operate the diagonal scaling
!!$!> Procedure to operate the diagonal scaling
!!$     !> of matrix M by the diagonal matrix D
!!$     !> out_matrix = D M D
!!$     !> (public procedure for type stiffmat)
!!$     procedure, public , pass :: diagonal_scale
  end type stiffmat
contains
  !>-------------------------------------------------------------
  !> Static constructor.
  !> (procedure public for type spmat)
  !> Instantiate (allocate if necessary)
  !> and initilize (by also reading from input file)
  !> variable of type spmat
  !>
  !> usage:
  !>     call 'var'%init(lun_err, matrix_A, tdens, inv_weight)
  !>
  !> where:
  !> \param[in] lun_err            -> integer. Error logical unit
  !> \param[in] nrow               -> integer. Number of rows 
  !> \param[in] nrow               -> integer. Number of columns
  !> \param[in] nterm              -> integer. Number of non-zero term
  !>                                     stored in the matrix
  !<-------------------------------------------------------------
  subroutine init_stiffmat(this, &
       lun_err,&
       matrix_A,&
       inv_weight)
    use Globals
    use Matrix
    implicit none
    !var
    class(stiffmat),             intent(inout) :: this
    integer,                     intent(in   ) :: lun_err
    class(abs_matrix), target,   intent(in   ) :: matrix_A  
    real(kind=double),           intent(in   ) :: inv_weight(matrix_A%nrow)
    
    ! local vars
    integer :: res,i,j,ntdens
    logical :: rc
    type(spmat) matrix_transpose

    !  definition of abs_matrix quantities
    this%nrow          = matrix_A%ncol
    this%ncol          = matrix_A%ncol
    this%is_symmetric  = .true.
    this%triangular    = 'N'
    this%ntdens        = matrix_A%nrow

    ! defition of type stiffmat 
    this%matrix_A  => matrix_A
    allocate(&
         this%tdens(matrix_A%nrow),&
         this%sqrt_tdens(matrix_A%nrow),&
         this%inv_weight(matrix_A%ncol),&
         stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, &
         'init_stiffmat', &
         'type stiffmat member tdens',res)
    this%tdens = zero

    this%inv_weight   = inv_weight

    !
    ! allocation work array
    !
    if (.not. allocated(this%scr_ncol_A)  ) then
       allocate(this%scr_ncol_A(this%ncol),stat=res)
       if(res .ne. 0) rc = IOerr(lun_err, err_alloc, &
            'init_stiffmat', &
            'type stiffmat member scr_ncol_A',res)
    end if

    if (.not. allocated(this%scr_nrow_A)  ) then
       allocate(this%scr_nrow_A(this%matrix_A%nrow),stat=res)
       if(res .ne. 0) rc = IOerr(lun_err, err_alloc, &
            'init_stiffmat', &
            'type stiffmat member scr_nrow_A',res)
    end if

    select type ( matrix => this%matrix_A )
    type is (densemat)
       write(*,*)' DENSE MATRIX'
       call this%dense%init(lun_err, &
            this%ncol, this%ncol,&
            ! optional arguments
            .True.,'N','U')
       
       if (.not. allocated(this%scr_tdens_A) ) then
          allocate(&
               this%scr_tdens_A(this%matrix_A%nrow,this%matrix_A%ncol),&
               stat=res)
          if(res .ne. 0) rc = IOerr(lun_err, err_alloc, &
               'init_stiffmat', & 
               'type stiffmat member scr_tdens_A',res) 
       end if
       
       if ( .not. this%A_transpose%is_initialized) then
          this%A_transpose = matrix
          call this%A_transpose%transpose(lun_err)
       end if
       
       
    type is (spmat)
       !
       ! use 
       !
       write(*,*)' SPARSE MATRIX'
       ntdens=this%ntdens
       matrix_transpose=matrix
       call matrix_transpose%transpose(lun_err)
       call this%sparse%MULT_MDN(lun_err, matrix_transpose , matrix, &
               ntdens, 100*ntdens ,this%tdens * this%inv_weight, matrix_transpose)
       call  matrix_transpose%kill(lun_err) 

       
!!$    type is (adjmat)
!!$       write(*,*)' INCIDENCE MATRIX GRAPH '
!!$       !
!!$       allocate(this%assembler(2,2,matrix%graph%ncell),stat=res)
!!$       if(res .ne. 0) rc = IOerr(lun_err, err_alloc, &
!!$            'init_stiffmat', & 
!!$            'type stiffmat member scr_tdens_A',res)
!!$       call matrix%preassembly_stiff(lun_err,&
!!$            this%sparse,this%assembler)
    end select

    !
    ! selection procedure
    !
    this%select_tdens = .False.
    allocate(this%selection_tdens(this%matrix_A%nrow),stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, &
         'init_stiffmat', &
         'type stiffmat member selection_tdens',res)

    
  end subroutine init_stiffmat

  !>-------------------------------------------------------------
  !> Static constructor.
  !> (procedure public for type spmat)
  !> Instantiate (allocate if necessary)
  !> and initilize (by also reading from input file)
  !> variable of type spmat
  !>
  !> usage:
  !>     call 'var'%init(lun_err, matrix_A, tdens, inv_weight)
  !>
  !> where:
  !> \param[in] lun_err            -> integer. Error logical unit
  !> \param[in] nrow               -> integer. Number of rows 
  !> \param[in] nrow               -> integer. Number of columns
  !> \param[in] nterm              -> integer. Number of non-zero term
  !>                                     stored in the matrix
  !<-------------------------------------------------------------
  subroutine set_stiffmat(this, tdens, &
       ! optional argument
       select_tdens, nnz, selection_tdens,&
       lasso )
    use Globals
    use Matrix
    implicit none
    !var
    class(stiffmat), intent(inout) :: this
    real(kind=double), intent(in)  :: tdens(this%matrix_A%nrow) 
    logical, optional, intent(in)  :: select_tdens
    integer, optional, intent(in)  :: nnz
    integer, optional, intent(in)  :: selection_tdens(this%matrix_A%nrow)
    real(kind=double), optional, intent(in) :: lasso
    
    this%tdens = tdens
    
    this%sqrt_tdens = sqrt(tdens*this%inv_weight)
    
    ! 
    !
    this%is_assembled = .False.

    if ( present( select_tdens ) ) then
       this%select_tdens    = select_tdens
       this%selection_tdens = selection_tdens
       this%nnz_tdens       = nnz
    end if

    ! lasso relazation
    if ( present(lasso) ) then
       this%lasso = lasso
    end if

  end subroutine set_stiffmat


  
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
  subroutine kill_stiffmat(this, lun_err)
    implicit none
    ! vars
    class(stiffmat) :: this
    integer, intent(in) :: lun_err
    ! local vars
    integer :: res
    logical :: rc

    this%matrix_A => null()
    
   
    if ( allocated(this%tdens)  ) then
       deallocate(this%tdens,stat=res)
       if(res .ne. 0) rc = IOerr(lun_err, err_dealloc, &
            'kill_stiffmat', &
            'type stiffmat member tdens',res)
    end if

    if ( allocated(this%sqrt_tdens)  ) then
       deallocate(this%sqrt_tdens,stat=res)
       if(res .ne. 0) rc = IOerr(lun_err, err_dealloc, &
            'kill_stiffmat', &
            'type stiffmat member sqrt_tdens',res)
    end if

    if ( allocated(this%scr_ncol_A)  ) then
       deallocate(this%scr_ncol_A,stat=res)
       if(res .ne. 0) rc = IOerr(lun_err, err_dealloc, &
            'kill_stiffmat', &
            'type stiffmat member scr_ncol_A',res)
    end if

    if ( allocated(this%scr_nrow_A)  ) then
       deallocate(this%scr_nrow_A,stat=res)
       if(res .ne. 0) rc = IOerr(lun_err, err_dealloc, &
            'kill_stiffmat', &
            'type stiffmat member scr_nrow_A',res)
    end if
    
    if ( allocated(this%scr_tdens_A)  ) then
       deallocate(this%scr_tdens_A,stat=res)
       if(res .ne. 0) rc = IOerr(lun_err, err_dealloc, &
            'kill_stiffmat', &
            'type stiffmat member scr_tdens_A',res)
    end if

    if ( this%dense%is_initialized) call this%dense%kill(lun_err) 
    
    if ( this%sparse%is_initialized) call this%sparse%kill(lun_err)

    if ( allocated(this%assembler)  ) then
       deallocate(this%assembler,stat=res)
       if(res .ne. 0) rc = IOerr(lun_err, err_dealloc, &
            'kill_stiffmat', &
            'type stiffmat member assembler',res)
    end if
    
  end subroutine kill_stiffmat

  !>-------------------------------------------------------------
  !> Info procedure.
  !> (public procedure for type stiffmat)
  !> Prints content of a variable of type stiffmat
  !> 
  !> usage:
  !>     call 'var'%info(lun_out)
  !>
  !> where:
  !> \param[in] lun_out -> integer. I/O unit for output message
  !<-------------------------------------------------------------
  subroutine info_stiffmat(this, lun)
    use Globals
    implicit none
    class(stiffmat), intent(in) :: this
    integer, intent(in) :: lun

    write(lun,'(a)') 'Nothing to declare'
    
  end subroutine info_stiffmat



  !>-------------------------------------------------------------
  !> Procedure to compute Matrix vector product
  !>         vec_out = (matrix)  times (vec_in)
  !> (public procedure for type stiffmat)
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
  subroutine Mxv_stiffmat(this,vec_in, vec_out, info,lun_err)
    use Globals
    implicit none
    class(stiffmat),   intent(inout) :: this
    real(kind=double), intent(in   ) :: vec_in(this%ncol)
    real(kind=double), intent(inout) :: vec_out(this%nrow)
    integer,           intent(inout) :: info
    integer,           intent(in   ) :: lun_err

    !
    real(kind=double) :: dnrm2
    integer :: iterm,irow
    
    info=0
    !if ( .not. this%is_assembled) then
    if ( this%select_tdens) then 
       select type ( matrix => this%matrix_A)
       type is (densemat)
          ! y1      = A vec_in
          ! y2      = diag(tdens) diag(inv_weight) y1
          ! vec_out = A^T y2

          ! y1 = select (A x)
          call matrix%select_Mxv(vec_in,this%scr_nrow_A, &
               this%nnz_tdens, this%selection_tdens)

          ! y2 = select(diag(tdens) diag(inv_weight) y1)
          do iterm=1,this%nnz_tdens
             irow=this%selection_tdens(iterm) 
             this%scr_nrow_A(irow) = &
                  this%tdens(irow) * &
                  this%inv_weight(irow) * &
                  this%scr_nrow_A(irow)
          end do
          ! vec_out = A^T select(y2)
          call matrix%MTxv_select(this%scr_nrow_A,vec_out,&
               this%nnz_tdens,this%selection_tdens)
       end select
    else
       ! y1      = A vec_in
       ! y2      = diag(tdens) diag(inv_weight) y1
       ! vec_out = A^T y2

       ! y1 = alpha A x
       call this%matrix_A%matrix_times_vector(vec_in,this%scr_nrow_A,info,lun_err)

       ! y2 = diag(tdens) diag(inv_weight) y1
       this%scr_nrow_A = this%tdens * this%inv_weight * this%scr_nrow_A

       ! vec_out = A^T y2
       call this%matrix_A%matrix_transpose_times_vector(this%scr_nrow_A,vec_out,info,lun_err)
       
    end if
       
    !
    ! add lasso contribution
    !
    if ( abs(this%lasso) .gt. small) then
       vec_out = vec_out + this%lasso * vec_in
    end if
       

  end subroutine Mxv_stiffmat

  !>-------------------------------------------------------------
  !> Procedure to compute Matrix vector product
  !>         vec_out = (matrix)^T  times (vec_in)
  !> (public procedure for type stiffmat)
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
  subroutine MTxv_stiffmat(this,vec_in, vec_out, info,lun_err)
    use Globals
    implicit none
    class(stiffmat),   intent(inout) :: this
    real(kind=double), intent(in   ) :: vec_in(this%nrow)
    real(kind=double), intent(inout) :: vec_out(this%ncol)
    integer,           intent(inout) :: info
    integer,           intent(in   ) :: lun_err


    call this%matrix_times_vector(vec_in,vec_out,info=info,lun_err=lun_err)

  end subroutine MTxv_stiffmat


  !>-------------------------------------------------------------
  !> Procedure to get the diagonal of a structure varible spmat
  !> of a matrix in ssr format
  !>
  !> usage : this%get_diagonal(diag)
  !> 
  !> \param[out] diag -> real. (dimension = this%nrow)
  !>                       Diagonal of spmat
  !>--------------------------------------------------------------
  subroutine get_diagonal(this,diagonal)
    use Globals
    use SparseMatrix
    use Densematrix
    implicit none
    class(stiffmat),   intent(in   ) :: this
    real(kind=double), intent(inout) :: diagonal(this%ncol)
    !local 
    integer i,j,ind,icol

    
    select type(matrix => this%matrix_A)
    class is (spmat)
       diagonal=zero
       do i=1,matrix%nrow
          do j=matrix%ia(i),matrix%ia(i+1)-1
             icol=matrix%ja(j)
             diagonal(icol) = diagonal(icol) + &
                  this%tdens(icol)*this%inv_weight(icol) * &
                  matrix%coeff(j)**2 
          end do
       end do
    class is (densemat)
       do j=1,matrix%ncol
          ! diag(s)_{j} = \sum_{k} ( A_{k,j} ) ^ 2 tdens(k) * inv_weight(k)
          diagonal(j) = zero
          do i=1,matrix%nrow
             diagonal(j) = diagonal(j) + &
                  this%tdens(i)*this%inv_weight(i) * &
                  matrix%coeff(i,j)**2 
          end do
       end do
    end select

    !
    ! add lasso contribution
    !
    if ( abs(this%lasso) .gt. 1.0d-12) then
       diagonal = diagonal + this%lasso 
    end if
          

  end subroutine get_diagonal

  !>-------------------------------------------------------------
  !> Procedure to get the diagonal of a structure varible spmat
  !> of a matrix in ssr format
  !>
  !> usage : this%get_diagonal(diag)
  !> 
  !> \param[out] diag -> real. (dimension = this%nrow)
  !>                       Diagonal of spmat
  !>--------------------------------------------------------------
  subroutine assembly(this,nnz,selection,upper_part,lasso)
    use Globals
    use SparseMatrix
    use Densematrix
    use Timing
    implicit none
    class(stiffmat), target,  intent(inout) :: this
    integer, intent(in ) :: nnz
    integer, intent(in ) :: selection(nnz)
    logical, optional, intent(in   ) :: upper_part
    real(kind=double), optional, intent(in) :: lasso
    
    !local 
    integer i,j,k,ind,icol,start,iterm,A_nrow,A_ncol,idiag
    integer :: ntdens
    type(Tim) :: zerot,copy, scale, mm
    real(kind=double) :: alpha
    type(spmat) :: matrix_at

    ntdens = this%ntdens

    this%is_assembled = .True.
    if ( present (lasso) ) this%is_inexact = .True.

    select type( matrix => this%matrix_A)
    type is (spmat)
       !
       ! replace this procedure with assembly_ATA 
       !
       matrix_at=matrix
       call matrix_at%transpose(0)
       call this%sparse%MULT_MDN(0, matrix_at , matrix, &
               ntdens, 100*ntdens ,this%tdens * this%inv_weight)
       call  matrix_at%kill(0)

       

       !
       ! add optional relaxation
       !
       if ( present (lasso) ) then
          this%is_inexact = .True.
          do i=1,this%sparse%nrow
             idiag = this%sparse%idiag(i)
             this%sparse%coeff(idiag) = this%sparse%coeff(idiag) + lasso
          end do
       end if
       
!!$    class is (adjmat)
!!$       !
!!$       ! create explicit matrix
!!$       !
!!$       this%assembled_matrix => this%sparse
!!$       call matrix%assembly_stiff( this%tdens * this%inv_weight, &
!!$            this%assembler,this%sparse)
!!$       
!!$       !
!!$       ! add optional relaxation
!!$       !
!!$       if ( present (lasso) ) then
!!$          this%is_inexact = .True.
!!$          do i=1,this%sparse%nrow
!!$             idiag = this%sparse%idiag(i)
!!$             this%sparse%coeff(idiag) = this%sparse%coeff(idiag) + lasso
!!$          end do
!!$       end if
    class is (densemat)
       ! fix label for assembled matrix
       this%assembled_matrix => this%dense
       
       !
       ! selection of rows 
       !
       if ( nnz .ne. this%ntdens ) then
          this%dense%coeff = zero
          ! out_{i,j} = \sum_{k} ( A_{k,i} ) tdens(k) * inv_weight(k)* A{k,j}
          this%dense%coeff = zero
          do iterm = 1,nnz
             k = selection(iterm)            
             alpha = this%tdens(k)*this%inv_weight(k)
             call dsyr('U',matrix%ncol,alpha,&
                     this%A_transpose%coeff(1,k),1,&
                     this%dense%coeff,matrix%ncol)
          end do
          do j=1,matrix%ncol
             ! diag(s)_{j} = \sum_{k} ( A_{k,j} ) ^ 2 tdens(k) * inv_weight(k)
             this%dense%coeff(j,j) = zero
             do i=1,matrix%nrow
                this%dense%coeff(j,j) = this%dense%coeff(j,j) + &
                     this%tdens(i)*this%inv_weight(i) * &
                     matrix%coeff(i,j)**2 
             end do
          end do
       else   
          ! copy matrix values
          this%scr_tdens_A = matrix%coeff

          ! compute matrix-matrix product : B =sqrt(diag(tdens)) A
          do j = 1, matrix%ncol
             do i = 1, matrix%nrow
                this%scr_tdens_A( i, j ) = this%scr_tdens_A( i, j ) * &
                     this%sqrt_tdens( i )
             end do
          end do
          ! compute matrix-matrix product : B^T B = A^T diag(tdens) A          
          call dsyrk('U','T', this%ncol, matrix%nrow,&
               one, this%scr_tdens_A,matrix%nrow,&
               zero,this%dense%coeff,&
               this%nrow)

          !
          ! found in
          ! https://stackoverflow.com/questions/47013581/blas-matrix-by-matrix-transpose-multiply
          !
          ! not efficient but if we increase efficiency can be used for selection
          !
!!$          this%dense%coeff = zero
!!$          do i = 1, matrix%nrow
!!$             call dsyr('U',matrix%ncol,this%tdens(i),&
!!$                  this%A_transpose(1,i),1,&
!!$                  this%dense%coeff,matrix%ncol)
!!$          end do         
      end if
      !
      ! add optional relaxation
      !
      if ( present (lasso) ) then
         if ( abs(lasso) .gt. small)  then
            this%is_inexact = .True.
            do i=1,this%dense%nrow
               this%dense%coeff(i,i) = this%dense%coeff(i,i) + lasso
            end do
         end if
      end if
    end select
  

  end subroutine assembly


end module StiffnessMatrix


