module DenseMatrix
  use Globals
  use ReadableMatrix
  implicit none
  private
  !>----------------------------------------------------------------------
  !> Structure variable containg member storing sparse real matrices
  !> in csr ( compress sparse row ) and ssr (symmetric sparse row format )
  !>----------------------------------------------------------------------
  type, extends(readable_mat),public :: densemat
     !> Storage system
     !> FULL,UPPER,LOWER
     character(len=1)  :: symmetric_storage='F'
     ! Dimension (nrow, ncol)
     ! Non zero elements of the sparse matrix
     real(kind=double), allocatable :: coeff(:,:)
   contains
     !> static constructor 
     !> (procedure public for type densemat)
     procedure, public, pass :: init => init_densemat
     !> static destructor
     !> (procedure public for type densemat)
     procedure, public, pass :: kill => kill_densemat
     !> Writing procedure.
     !> (public procedure for type densemat)
     procedure, public, pass :: write => write_densemat
     !> Reading procedure.
     !> (public procedure for type densemat)
     procedure, public, pass :: read => read_densemat
     !> Info procedure  sparse matrix
     !> (public procedure for type densemat)
     procedure, public, pass :: info => info_densemat
     !> Procedure to compute 
     !>         y = M * x 
     !> (public procedure for type densemat)
     procedure, public, pass :: matrix_times_vector => Mxv_densemat
     !> Procedure to compute 
     !>         y = M * x 
     !> (public procedure for type densemat)
     procedure, public, pass :: select_Mxv
     !> Procedure to compute 
     !>         y = M^T * x 
     !> with M^T the transposed of a matrix M
     !> (public procedure for type densemat)
     procedure, public, pass :: matrix_transpose_times_vector => MTxv_densemat
     !> Procedure to compute 
     !>         y = M * x 
     !> (public procedure for type densemat)
     procedure, public, pass :: MTxv_select
     !> Procedure to transpose matrix
     procedure, public, pass :: transpose
     !> Procedure to get the diagonal 
     procedure, public, pass :: get_diagonal
     !> Procedure for scaling matrix by a diagonal
     !procedure, public, pass :: diagonal_scaling
     !> Procedure to compute cholesky factor U
     !> with A=UU^T for A symmetric 
     procedure, public, pass :: cholesky
     !> Procedure to compute spectral decomposition
     !> of symmetric matrices
     procedure, public, pass :: spectral_decomposition
     !> Procedure to form explicitely 
     !> of symmetric matrices
     procedure, public, pass :: form => form_densemat
     !> Procedure to count elements above a
     !> given threshold
     procedure, public, pass :: count_above_threshold
     !> Procedure to form sprse matrix from those
     !> elements above a given threshold
     procedure, public, pass :: sparsify 
  end type densemat
contains
  !>-------------------------------------------------------------
  !> Static constructor.
  !> (procedure public for type densemat)
  !> Instantiate (allocate if necessary)
  !> and initilize (by also reading from input file)
  !> variable of type densemat
  !>
  !> usage:
  !>     call 'var'%init(lun_err, nrow, nterm )
  !>
  !> where:
  !> \param[in] lun_err               -> integer. Error logical unit
  !> \param[in] nrow                  -> integer. Number of rows 
  !> \param[in] ncol                  -> integer. Number of columns
  !> \param[in] (optional) is_sym     -> Logical. T/F flag for symmetric matrix
  !> \param[in] (optional) is_tri     -> Character. 'upper_triangular' and 
  !>                                     'lower_triangular' fix densemat to be  
  !>                                     upper or lower triangular
  !<-------------------------------------------------------------
  subroutine init_densemat(this, lun_err, &
       nrow, ncol,&
       ! optional arguments
       is_symmetric,triangular,symmetric_storage)
    use Globals
    implicit none
    !var
    class(densemat),              intent(inout) :: this
    integer,                     intent(in   ) :: lun_err
    integer,                     intent(in   ) :: nrow
    integer,                     intent(in   ) :: ncol
    logical,           optional, intent(in   ) :: is_symmetric
    character (len=1), optional, intent(in   ) :: triangular
    character (len=1), optional, intent(in   ) :: symmetric_storage
    ! local vars
    integer :: res
    logical :: rc

    if (this%is_initialized) call this%kill(lun_err)
    
    this%is_initialized  = .true.
    this%nrow            = nrow
    this%ncol            = ncol

    allocate(this%coeff(nrow,ncol),stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'init_densemat', &
         '  type densemat member coeff (array)',res)
    

    if ( present(is_symmetric) ) then
       this%is_symmetric   = is_symmetric
       if ( present(  symmetric_storage ) ) then
          this%symmetric_storage   =  symmetric_storage
       end if
    end if
       
    
    
    if ( present(triangular)) then
       if ( (triangular .eq. 'N') .or. &
            (triangular .eq. 'U') .or. &
            (triangular .eq. 'L') ) then
          this%triangular = triangular
       end if
    end if
          
  end subroutine init_densemat

  !>-------------------------------------------------------------
  !> Static destructor.
  !> (procedure public for type densemat)
  !> deallocate all arrays for a var of type densemat
  !>
  !> usage:
  !>     call 'var'%kill(lun)
  !>
  !> where:
  !> \param[in] lun -> integer. I/O unit for error message output
  !<-----------------------------------------------------------
  subroutine kill_densemat(this, lun_err)
    implicit none
    ! vars
    class(densemat),intent(inout) :: this
    integer,     intent(in   ) :: lun_err
    ! local vars
    integer :: res
    logical :: rc

    
    deallocate(this%coeff,stat=res)
    if (res.ne.0) rc=IOerr(lun_err, err_dealloc, 'kill_densemat', &
         'dealloc fail for densemat member coeff',res)

    call this%to_default()
    this%is_initialized = .false.

  end subroutine kill_densemat
  
  !>-------------------------------------------------------------
  !> Writing procedure.
  !> (public procedure for type densemat)
  !> Prints content of a variable of type spamat
  !> 
  !> usage:
  !>     call 'var'%write(lun)
  !>
  !> where:
  !> \param[in] lun -> integer. I/O unit for error message output
  !>
  !<-------------------------------------------------------------
  subroutine write_densemat(this, lun)
    use Globals
    implicit none
    class(densemat),    intent(in) :: this
    integer,           intent(in) :: lun
    !real(kind=double), optional, intent(in) :: rho
    ! loc. var
    integer i,j,m,n,ind

    write(lun,*) this%nrow, this%ncol,'! number of rows and columns     '
    write(lun,*) this%nrow * this%ncol, '! number of terms '
    ind=0
    do j=1,this%ncol
       do i=1,this%nrow
          ind=ind+1
          write(lun,1010) i,j,this%coeff(i,j)
       end do
    end do
1010 format(2i15,1pe24.16)

  end subroutine write_densemat


  !>-------------------------------------------------------------
  !> Writing procedure.
  !> (public procedure for type densemat)
  !> Prints content of a variable of type spamat
  !> 
  !> usage:
  !>     call 'var'%read(lun)
  !>
  !> where:
  !> \param[in] lun -> integer. I/O unit for error message output
  !>
  !<-------------------------------------------------------------
  subroutine read_densemat(this, lun_err,lun)
    use Globals
    implicit none
    class(densemat),    intent(inout) :: this
    integer,           intent(in   ) :: lun_err
    integer,           intent(in   ) :: lun
    
    ! loc. var
    logical :: rc
    integer :: res
    integer :: i,j,m,n,ind,nrow,ncol,nterm,l,iterm,icol,irow
    real(kind=double) :: rvalue

    read(lun,*,iostat=res) nrow, ncol
    if (res.ne.0) rc=IOerr(lun_err, err_inp, 'read_densemat', &
         'type densemat members nrow, ncol')
    read(lun,*,iostat=res) i
    if (res.ne.0) rc=IOerr(lun_err, err_inp, 'read_densemat', &
         'type densemat members nterm')
        
    call this%init(lun_err, nrow, ncol,&
         is_symmetric=(nrow .eq. ncol))
    
    ind=0
    do j=1,this%ncol
       do i=1,this%nrow
          ind=ind+1
          read(lun,*,iostat=res) irow, icol,rvalue
          if (res.ne.0) rc=IOerr(lun_err, err_inp, 'read_densemat', &
               'type densemat')
          this%coeff(i,j) = rvalue
       end do
    end do
    
  end subroutine read_densemat

  !>-------------------------------------------------------------
  !> Info procedure.
  !> (public procedure for type densemat)
  !> Prints content of a variable of type spamat
  !> 
  !> usage:
  !>     call 'var'%info(lun_out)
  !>
  !> where:
  !> \param[in] lun_out -> integer. I/O unit for output message
  !>
  !<-------------------------------------------------------------
  subroutine info_densemat(this, lun)
    use Globals
    implicit none
    class(densemat), intent(in) :: this
    integer , intent(in) :: lun
    ! loc. var
    character(len=256) :: state,mem,prop=' ' 
    character(len=3) :: sep=' | '

    if ( this%is_initialized ) then
       write(state,'(a)') 'Densemat allocated'
       if ( this%triangular .eq. 'U' ) write(prop,'(a)') 'upper-triangular'
       if ( this%triangular .eq. 'L' ) write(prop,'(a)') 'lower-triangular'
       if ( this%is_symmetric ) write(prop,'(a)') 'symmetric'
       
       write(lun,'(a)') etb(etb(state)//sep//etb(prop))

       write(lun,'(a,I8,a,a,I8)') &
            'nrows= ', this%nrow,sep,&      
            'ncol= ', this%ncol
    else
       write(lun,*) 'Densemat not initialized'
    end if

  end subroutine info_densemat

  !>-------------------------------------------------------------
  !> Procedure to cycle all element in matrix
  !> (public procedure for type densemat)
  !> 
  !> usage:
  !>     call 'var'%Mxv(vec_in,vec_out,[info])
  !>
  !> where 
  !> \param[in   ] vec_in          -> real. dimension('var'%ncol)
  !>                                  vector to be multiplied
  !> \param[inout] vec_out         -> real. dimension('var'%nrow)
  !>                                  vector (M) times (vec_in) 
  !> \param[in   ] (optional) info -> integer. Info number
  !>                                  in case of error   
  !<-------------------------------------------------------------
  subroutine get_coeff(this,iterm,irow,icol,rvalue)
    use Globals
    implicit none
    class(densemat),    intent(inout) :: this
    integer,           intent(in   ) :: iterm  
    integer,           intent(out  ) :: irow
    integer,           intent(out  ) :: icol
    real(kind=double), intent(out  ) :: rvalue

    irow   = (iterm+1) / this%ncol  
    icol   = (iterm+1) / this%nrow
    rvalue = this%coeff(irow,icol)
    
  end subroutine get_coeff

  !>-------------------------------------------------------------
  !> Procedure to compute Matrix vector product
  !>         vec_out = (M) times (vec_in)
  !> (public procedure for type densemat)
  !> 
  !> usage:
  !>     call 'var'%matrix_times_vector(vec_in,vec_out,[info])
  !>
  !> where 
  !> \param[in   ] vec_in          -> real. dimension('var'%ncol)
  !>                                  vector to be multiplied
  !> \param[inout] vec_out         -> real. dimension('var'%nrow)
  !>                                  vector (M) times (vec_in) 
  !> \param[in   ] (optional) info -> integer. Info number
  !>                                  in case of error   
  !<-------------------------------------------------------------
  subroutine Mxv_densemat(this,vec_in,vec_out, info,lun_err)
    use Globals
    implicit none
    class(densemat),   intent(inout) :: this
    real(kind=double), intent(in   ) :: vec_in(this%ncol)
    real(kind=double), intent(inout) :: vec_out(this%nrow)
    integer,           intent(inout) :: info
    integer,           intent(in   ) :: lun_err


    info=0

    if ( .not. this%is_symmetric ) then
       call dgemv ('N', this%nrow, this%ncol,&
            one, this%coeff, &
            this%nrow, vec_in, 1,&
            zero, vec_out, 1)
    else   
       if (this%symmetric_storage .eq. 'F') then
          call dsymv('U',&
               this%nrow,one,this%coeff,this%ncol,vec_in,1,zero,vec_out,1)
       else
          call dsymv(this%symmetric_storage,&
               this%nrow,one,this%coeff,this%ncol,vec_in,1,zero,vec_out,1)
       end if          
    end if
    
  end subroutine Mxv_densemat

  !>-------------------------------------------------------------
  !> Procedure to compute Matrix vector product
  !>         vec_out = (M) times (vec_in)
  !> (public procedure for type densemat)
  !> 
  !> usage:
  !>     call 'var'%Mxv(vec_in,vec_out,[info])
  !>
  !> where 
  !> \param[in   ] vec_in          -> real. dimension('var'%ncol)
  !>                                  vector to be multiplied
  !> \param[inout] vec_out         -> real. dimension('var'%nrow)
  !>                                  vector (M) times (vec_in) 
  !> \param[in   ] (optional) info -> integer. Info number
  !>                                  in case of error   
  !<-------------------------------------------------------------
  subroutine select_Mxv(this,vec_in,vec_out, nnz, selection_out,info)
    use Globals
    implicit none
    class(densemat),    intent(inout) :: this
    real(kind=double), intent(in   ) :: vec_in(this%ncol)
    real(kind=double), intent(inout) :: vec_out(this%nrow)
    integer,           intent(in   ) :: nnz
    integer,           intent(in   ) :: selection_out(nnz)
    integer, optional, intent(inout) :: info
    !local
    integer:: icol, irow, innz
    real(kind=double) :: vi
    
    
    vec_out = zero
    do icol = 1, this%ncol
       vi = vec_in(icol)
       do innz = 1, nnz
          irow = selection_out(innz)
          vec_out(irow) = vec_out(irow) + this%coeff(irow,icol) * vi
       end do
    end do
    
    if (present(info)) info=0

  end subroutine select_Mxv

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
  subroutine MTxv_densemat(this,vec_in,vec_out, info,lun_err)
    use Globals
    implicit none
    class(densemat),   intent(inout) :: this
    real(kind=double), intent(in   ) :: vec_in(this%nrow)
    real(kind=double), intent(inout) :: vec_out(this%ncol)
    integer,           intent(inout) :: info
    integer,           intent(in   ) :: lun_err

    info=0
    call dgemv ('T', this%nrow, this%ncol,&
         one, this%coeff, &
         this%nrow, vec_in, 1,&
         zero, vec_out, 1)
    
  end subroutine MTxv_densemat

  !>-------------------------------------------------------------
  !> Procedure to compute Matrix vector product
  !>         vec_out = (M) times (vec_in)
  !> (public procedure for type densemat)
  !> 
  !> usage:
  !>     call 'var'%MTxv(vec_in,vec_out,[info])
  !>
  !> where 
  !> \param[in   ] vec_in          -> real. dimension('var'%ncol)
  !>                                  vector to be multiplied
  !> \param[inout] vec_out         -> real. dimension('var'%nrow)
  !>                                  vector (M) times (vec_in) 
  !> \param[in   ] (optional) info -> integer. Info number
  !>                                  in case of error   
  !<-------------------------------------------------------------
  subroutine MTxv_select(this,vec_in,vec_out,nnz, selection_in, info)
    use Globals
    implicit none
    class(densemat),    intent(inout) :: this
    real(kind=double), intent(in   ) :: vec_in(this%nrow)
    real(kind=double), intent(inout) :: vec_out(this%ncol)
    integer,           intent(in   ) :: nnz
    integer,           intent(in   ) :: selection_in(this%nrow)
    integer, optional, intent(inout) :: info
    !local
    integer :: irow, icol, innz
    

    vec_out = zero
    do icol = 1,this%ncol
       do innz = 1, nnz
          irow = selection_in(innz)
          vec_out(icol) = vec_out(icol) + this%coeff(irow,icol) * vec_in(irow)
       end do
    end do

    if (present(info)) info=0
    
  end subroutine MTxv_select

  subroutine transpose(this,&
       lun_err)
    use Globals
    implicit none
    ! Input variables
    class(densemat), intent(inout) :: this
    integer,     intent(in   ) :: lun_err

    ! Local variables
    logical :: rc
    integer :: res
    integer :: i,j,ind
    character(len=1) :: new_triangular
    type(densemat) :: temp

    ! check if the matrix id symmetric
    if ( this%is_symmetric ) return
     
    ! set properties of transpose matrix
    if (this%triangular .eq. 'N') then
       new_triangular = 'N'
    else
       select case (this%triangular)
       case ('U')
          new_triangular = 'L'
       case ('L')
          new_triangular = 'U'
       end select
    end if

    ! initialized working matrix
    call temp%init(lun_err,&
         this%ncol,this%nrow,& 
         triangular=new_triangular) 
    
    do j=1,this%ncol
       do i=1,this%nrow
          temp%coeff(j,i)=this%coeff(i,j)
       end do
    end do

    call this%kill(lun_err)
       
    ! copy the working spmat into this
    select type(this)
    type is (densemat)
       this = temp
    end select
    
    call temp%kill(lun_err)


    

  end subroutine transpose

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
    implicit none
    class(densemat),    intent(in   ) :: this
    real(kind=double), intent(inout) :: diagonal(min(this%nrow,this%ncol))
    !local 
    integer i,j,ind,icol
    
    do i=1,min(this%nrow,this%ncol)
       diagonal(i)=this%coeff(i,i)
    end do
    
  end subroutine get_diagonal


  !>-------------------------------------------------------------
  !> Procedure to from matrix given from other format 
  !>
  !> usage : this%from_densemat(lun_err,matrixA)
  !> 
  !> \param[in] matrixA -> class(linop). Linear operator
  !>                       Supported type:
  !>                        scalmat,diagmat,
  !>                        spmat,linear_combination
  !>--------------------------------------------------------------
  subroutine form_densemat(this,lun_err,matrixA)
    use Globals
    use LinearOperator
    use SparseMatrix
    implicit none
    class(densemat),   intent(inout) :: this
    integer,          intent(in   ) :: lun_err
    class(abs_linop), intent(in   ) :: matrixA
    !local 
    integer i,j,ind,icol
    logical :: rc
    integer :: res

    call this%init(lun_err, &
         matrixA%nrow, matrixA%ncol,&
         ! optional arguments
         matrixA%is_symmetric, matrixA%triangular,'F')

    this%coeff=zero
    call form_coefficients(lun_err,matrixA,this%coeff,one)
    
  contains
    recursive subroutine form_coefficients(lun_err,matrixA,coeff,alpha)
      use Globals
      use LinearOperator
      use SparseMatrix
      use Matrix
      use SimpleMatrix
      implicit none
      integer,          intent(in   ) :: lun_err
      class(abs_linop), intent(in   ) :: matrixA
      real(kind=double), intent(inout) :: coeff(matrixA%nrow,matrixA%ncol)
      real(kind=double), intent(in  ) :: alpha
      
      !local 
      integer i,j,ind,irow,icol,ilinop
      integer :: in_begin,in_end, out_begin, out_end
      logical :: rc
      integer :: res
      select type( mat => matrixA)
      type is (spmat)
         do i=1,mat%nrow
            do j=mat%ia(i),mat%ia(i+1)-1
               coeff(i,mat%ja(j))=coeff(i,mat%ja(j))+alpha*mat%coeff(j)
            end do
         end do
      type is (scalmat)
         do i=1,min(mat%nrow,mat%ncol)
            coeff(i,i)=coeff(i,i)+alpha*mat%scalar
         end do         
      type is (diagmat)
         do i=1,min(mat%nrow,mat%ncol)
            coeff(i,i)=coeff(i,i)+alpha*mat%diagonal(i)
         end do
      type is (pair_linop)
         if ( mat%type .eq. 'LC')  then
            do ilinop = 1, mat%nlinops
               call form_coefficients(lun_err,  mat%linop_list(ilinop)%linop,coeff,mat%alphas(ilinop))
            end do
         else
            write(*,*) 'Assembly Linear product not defined yet'
         end if
      type is (block_linop)
         !
         ! to be tested
         !
         do i = 1, mat%nnzblock
            ilinop = mat%block_structure(1,i)
            irow   = mat%block_structure(2,i)
            icol   = mat%block_structure(3,i)

            in_end    = sum(mat%ncol_vec(1:icol))
            in_begin  = in_end - mat%ncol_vec(icol) + 1
            

            out_end   = sum(mat%nrow_vec(1:irow))
            out_begin = out_end - mat%nrow_vec(irow) + 1
            
            call form_coefficients(lun_err,  mat%linop_list(ilinop)%linop,&
                 coeff(in_begin:in_end,out_begin:out_end),mat%alphas(ilinop))
         end do
      class default
         if (res.ne.0) rc=IOerr(lun_err, err_inp, 'form_densemat', &
              'passe mLinearOperatoratrix type not supported')
      end select

    end subroutine form_coefficients

    
  end subroutine form_densemat

  !>-------------------------------------------------------------
  !> Procedure to get sparse matrix from
  !> a dense matrix selection those entries above a given
  !> threeshold 
  !>
  !> usage : this%sparsifty(lun_err,threshold,sparse)
  !> 
  !> \param[out] diag -> real. (dimension = this%nrow)
  !>                       Diagonal of spmat
  !>--------------------------------------------------------------
  subroutine sparsify(this,lun_err,threshold,sparse)
    use Globals
    use SparseMatrix
    class(densemat), intent(in) :: this
    integer, intent(in ) :: lun_err
    real(kind=double), intent(in) :: threshold
    type(spmat), intent(inout) :: sparse
    !local
    logical ::rc
    integer :: res
    integer :: irow,icol, ind,nnz
    integer, allocatable :: nnzinrow(:)

    !
    ! count nnz zeros in rows
    !
    allocate(nnzinrow(this%nrow),stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'sparsify', &
         ' temp array nnzinrow',res)   
    nnzinrow=0
    do icol=1,this%ncol
       do irow=1,this%nrow
          if ( abs( this%coeff(irow,icol) ) > threshold ) nnzinrow(irow)=nnzinrow(irow)+1
       end do
    end do
    nnz=sum(nnzinrow)
    
    !
    ! initialize the matrix
    !
    call sparse%init(lun_err,this%nrow,this%ncol,nnz,'csr',&
         is_symmetric=this%is_symmetric,&
         triangular=this%triangular)

    !
    ! set ia ja coeff of sparse matrix
    !
    sparse%ia(1)=1
    do irow=1,this%nrow
       sparse%ia(irow+1)=sparse%ia(irow)+nnzinrow(irow)
    end do
    nnzinrow=0
    do icol=1,this%ncol
       do irow=1,this%nrow
          if ( abs( this%coeff(irow,icol) ) > threshold ) then
             ind=sparse%ia(irow)+nnzinrow(irow)
             sparse%ja(ind)=icol
             sparse%coeff(ind)=this%coeff(irow,icol)
             nnzinrow(irow)=nnzinrow(irow)+1
          end if
       end do
    end do

    deallocate(nnzinrow,stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'sparsify', &
         ' temp array nnzinrow',res)   

      
  end subroutine sparsify

  !>-------------------------------------------------------------
  !> Procedure to count the number of elements above a given
  !> threeshold 
  !>
  !> usage : this%count_above_threshold(threshold)
  !> 
  !> \param[in] threshold -> real. Threshold
  !>
  !> \return    nnz       -> Number of elements above threshold 
  !>--------------------------------------------------------------
  function count_above_threshold(this,threshold) result (nnz)
    use Globals
    use SparseMatrix
    class(densemat), intent(in) :: this
    real(kind=double), intent(in) :: threshold
    integer :: nnz
    !local
    integer :: irow,icol

    nnz=0
    do icol=1,this%ncol
       do irow=1,this%nrow
          if ( abs( this%coeff(irow,icol) ) > threshold ) nnz=nnz+1
       end do
    end do
  end function count_above_threshold
  

  !>-------------------------------------------------------------
  !> Procedure to compute the cholesky factor of
  !> positive definitive matrix M.
  !> Only upper factor U , with M=U^T U, is computed
  !>
  !> usage : this%cholesky(lun_err,info,upper)
  !> 
  !> \param[in ] lun_err  -> integer. I/O unit for error msg
  !> \param[out] info     -> integer. Info variable 
  !>                           (info=0 == no errors)
  !> \param[inout] upper  -> type(densemat)
  !>                         Dense Matrix containg the
  !>                         U factor s.t. M=U^TU
  !>--------------------------------------------------------------
  subroutine cholesky(this,lun_err, info, upper_factor)
    use Globals
    implicit none
    class(densemat), intent(in) :: this
    integer, intent(in) :: lun_err
    integer, intent(out) :: info
    type(densemat), intent(inout) :: upper_factor
    !local 
    logical :: rc
    integer :: res,nequ
    
    
    !
    ! checks matrix
    !    
    if ( this%nrow .ne. this%ncol ) then
       if(res .ne. 0) rc = IOerr(lun_err, wrn_inp, 'cholesky', &
            ' not squared matrix passed')
       info = -1
       return
    end if
    nequ = this%nrow
    if ( this%nrow .le. 0 ) then
       if(res .ne. 0) rc = IOerr(lun_err, wrn_inp, 'cholesky', &
            ' Matrix dimensions less or equal than 0')
       info = -1
       return
    end if
    if ( .not. this%is_symmetric) then
       if(res .ne. 0) rc = IOerr(lun_err, wrn_inp, 'cholesky', &
            ' not symmetric matrix passed')
       info = -1
       return
    end if
    
    !
    ! check factor
    !
    if ( .not. upper_factor%is_initialized) then
       call upper_factor%init(lun_err, &
       nequ, nequ,&
       triangular='U',symmetric_storage='F')
    end if

    ! copy coefficient ( dpotrf works in place)
    upper_factor%coeff = this%coeff
    ! compute upper factor
    CALL DPOTRF( 'U', nequ, upper_factor%coeff, nequ, info )
        
  end subroutine cholesky

  !>-------------------------------------------------------------
  !> Procedure to compute eigenvalues and, optionally, eigenvectors
  !> for a dense matrix. Only Symmetric case is considered.
  !>
  !> usage : this%spectral_decomposition(lun_err,&
  !>                                     info,eigenvalues,&
  !>                                     [eigenvectors],&
  !>                                     [aux])
  !> 
  !> \param[in   ] lun_err  -> integer. I/O unit for error msg
  !> \param[inout] info     -> integer. Info variable 
  !>                           In output:
  !>                            (info =0 : no errors  )
  !>                            (info/=0 : check dsyev)
  !>                           In input: if info == -999
  !>                           we return in info the size of the
  !>                           work space required
  !> \param[inout]
  !>           eigenvalues  -> real(dimension=matrix%ncol)
  !>                           REAL eigenvalues
  !> \param[inout,optional]
  !>           eigenvectors -> real(dimension=matrix%ncol,matrix%ncol)
  !>                           Orthonormal eigenvectors
  !> \param[inout,optional]
  !>           aux          -> type(scrt) Scratch arrays containers
  !>--------------------------------------------------------------
  subroutine spectral_decomposition(this,&
       lun_err, info, eigenvalues, eigenvectors, aux)
    use Globals
    use Scratch
    implicit none
    class(densemat), intent(in) :: this
    integer, intent(in   ) :: lun_err
    integer, intent(inout) :: info
    real(kind=double),  intent(inout) :: eigenvalues(this%nrow)
    real(kind=double),  optional, intent(inout) :: eigenvectors(this%nrow,this%nrow)
    type(scrt), target, optional, intent(inout) :: aux
    !local 
    logical :: rc
    integer :: res,nequ,nraux
    type(scrt), target  :: aux_loc
    type(scrt), pointer :: aux_work
    integer :: n1,n2,n3,nb ,N,nwork
    integer :: ilaenv
    character(len=1) :: UPLO,JOB
    real(kind=double), pointer :: copy_coeff(:,:)

    
    !
    ! checks matrix
    !
    if ( this%nrow .ne. this%ncol ) then
       if(res .ne. 0) rc = IOerr(lun_err, wrn_inp, &
            ' spectral_decomposition', &
            ' not squared matrix passed')
       info = -1
       return
    end if
    N = this%nrow
    if ( .not. this%is_symmetric) then
       if(res .ne. 0) rc = IOerr(lun_err, wrn_inp, &
            ' spectral_decomposition', &
            ' not symmetric matrix passed')
       info = -1
       return
    end if
    !
    ! prepare work space 
    !
    UPLO=this%symmetric_storage
    if ( UPLO .eq. 'F') then
       UPLO='U'
    end if

    nb=ilaenv( 1, 'dsyev', UPLO, N1, N2, N3,this%nrow) ! n1,n2,n3 are ignored
    nwork=(nb+2)* this%nrow

    
    
    if ( present(eigenvectors) ) then
       JOB='V'
       nraux=nwork
    else
       !
       ! add space for copying coefficient of matrix
       !
       JOB='N'
       nraux=nwork+this%nrow**2
    end if
    !
    ! if (info is passed nothing is done and we return
    ! in info the dimension of the work space requiredx
    ! 
    if (info .eq. -999 ) then
       info=nraux
       !call aux%init(lun_err, 0 , nraux)
       return
    end if
    
    if ( present(aux)) then
       if (.not. aux%check(0,nraux) ) then
          info=-1
          return
       else
          aux_work => aux
       end if
    else
       call aux_loc%init(lun_err, 0 , nraux)
       aux_work => aux_loc
    end if

   
    if ( aux_work%nraux .eq. nwork ) then
       ! copy matrix ( dsyev works on place)
       eigenvectors = this%coeff 
       call dsyev( JOB, UPLO,&
            this%nrow, eigenvectors, this%nrow, eigenvalues,&
            aux_work%raux,  nwork , info )
    else
       copy_coeff (1:N,1:N) => aux_work%raux(nwork+1:aux_work%nraux)
       copy_coeff (1:N,1:N) = this%coeff 
       call dsyev( JOB, UPLO,&
            this%nrow,&
            copy_coeff,&
            this%nrow, eigenvalues,&
            aux_work%raux,  nwork, info )
    end if
       
       
    !
    ! free memory
    !
    if ( present(aux)) then
       aux_work => null()
    else
       aux_work => null()
       call aux_loc%kill(lun_err)
    end if

  end subroutine spectral_decomposition




  !>-------------------------------------------------------------
  !> Procedure to get the diagonal of a structure varible spmat
  !> of a matrix in ssr format
  !>
  !> usage : this%get_diagonal(diag)
  !> 
  !> \param[out] diag -> real. (dimension = this%nrow)
  !>                       Diagonal of spmat
  !>--------------------------------------------------------------
  subroutine apply_inverse(this,lun_err,info,rhs,sol,aux)
    use Globals
    use Scratch
    implicit none
    class(densemat),       intent(in   ) :: this
    integer,              intent(in   ) :: lun_err
    integer,              intent(inout) :: info
    real(kind=double),    intent(in   ) :: rhs(this%nrow)
    real(kind=double),    intent(inout) :: sol(this%ncol)
    type(scrt), optional,target, intent(inout) :: aux    
    !local 
    integer :: nrow,ncol
    character :: flag_unit ='N'
    type(scrt), pointer :: aux_loc
    type(scrt), target :: aux_temp
    real(kind=double),pointer :: mat_temp(:,:), rhs_temp(:,:)
    
    nrow = this%nrow
    ncol = this%ncol

    if (.not. present(aux) ) then
       call aux_temp%init(lun_err,0,ncol**2+ncol)
       aux_loc => aux_temp
    else
       aux_loc => aux
    end if

    
    if (this%triangular .ne. 'N') then
       mat_temp(1:ncol,1:ncol) => aux%raux(1:ncol**2)
       rhs_temp(1:1,1:ncol)      => aux%raux(ncol**2+1:ncol**2+ncol)
       
       mat_temp(:,:) = this%coeff
       rhs_temp(1,:) = rhs
       ! solve
       !   L: op( A )*X = alpha*B,   
       ! or
       !   R: X*op( A ) = alpha*B,
       ! 
       if (this%unitary_diag) flag_unit='U'
       call dtrsm('L', &      ! L or R
            this%triangular,& ! U or L, how triangualr matrix is stored
            'N',&             ! N or T , transpose or not
            flag_unit,&       ! N or U , unitary diag or not 
            nrow,&            ! number of rows of b
            1, &              ! number of column of b
            one, &            ! alpha
            mat_temp,&        ! A matrix
            nrow,&            ! number of rows of A to solve
            rhs_temp,&        ! B matrix
            nrow)             ! number of rows of B to solve
       
       sol = rhs_temp(1,:)
    end if
    

    aux_loc => null()
    if (.not. present(aux) ) then
       call aux_temp%kill(lun_err)
    end if
       
    
  end subroutine apply_inverse

    
end module DenseMatrix



