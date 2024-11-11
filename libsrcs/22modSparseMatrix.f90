module SparseMatrix
  use Globals
  use LinearOperator
  use Matrix
  use ScalableMatrix
  use ReadableMatrix
  use SimpleMatrix
  implicit none
  private
  type, public :: sppattern
     !> Number of rows
     integer :: nrow
     !> Nmber of columns
     integer :: ncol
     !> NUmber of non-zeros term 
     integer :: nterm
     ! Dimension (nrow+1)
     ! Pointers to the first term of sysamt in row
     integer, allocatable :: ia(:)
     ! Dimension (nterm)
     ! Pointers to the column index of coeff
     integer, allocatable :: ja(:)
!!$   contains
!!$     !> Procedure to initialize product matrix     
!!$     !> (procedure public for type sppattern)
!!$     procedure, public, pass :: init_spattern
!!$     !> Procedure to initialize product matrix     
!!$     !> (procedure public for type sppattern)
!!$     procedure, public, pass :: get_nnz_row
!!$     !> Procedure to initialize product matrix     
!!$     !> (procedure public for type sppattern)
!!$     procedure, public, pass :: get_row    
  end type sppattern
     
  
  public :: MtimesN,apply_perm, &
       spmat_constructor
  !>----------------------------------------------------------------------
  !> Structure variable containg member storing sparse real matrices
  !> in csr ( compress sparse row ) and ssr (symmetric sparse row format )
  !>----------------------------------------------------------------------
  type, extends(readable_mat),public :: spmat
     !> Number of non-zeros term
     integer :: nterm=0
     !> Flag for ordered coefficient
     logical :: is_sorted=.true.
     !> Flag for storage style
     !> 'csr'= Compressed Storage Row
     !> 'ssr'= Symmetric  Storege Row (use only
     !> for Symmetric Matrices, it stores only the
     !> upper triangular part)
     character(len=3) :: storage_system='csr'
     ! Dimension (nrow+1)
     ! Pointers to the first term of sysamt in row
     integer, allocatable  :: ia(:)
     ! Dimension (nrow)
     ! Pointers to the diagonal term in coeff
     ! May be not always be allocated
     integer, allocatable :: ind_diag(:)
     ! Dimension (nterm)
     ! Pointers to the column index of coeff
     integer, allocatable :: ja(:)
     ! Dimension (nterm)
     ! Non zero elements of the sparse matrix
     real(kind=double), allocatable :: coeff(:)
   contains
     !> static constructor 
     !> (procedure public for type spmat)
     procedure, public, pass :: init => init_spmat
     !> static destructor
     !> (procedure public for type spmat)
     procedure, public, pass :: kill => kill_spmat
     !> Reading procedure.
     !> (public procedure for type spmat)
     procedure, public, pass :: read => read_spmat
     !> Writing procedure.
     !> (public procedure for type spmat)
     procedure, public, pass :: write => write_spmat
     !> Info procedure  sparse matrix
     !> (public procedure for type spmat)
     procedure, public, pass :: info => info_spmat
     !> Build the pointers to the diagonal terms,
     !> for a specific row, given ia and ja
     !> (public procedure for type spmat)
     procedure, public, pass :: idiag
     !> Build the pointers to the diagonal terms,
     !> given ia and ja
     !> (public procedure for type spmat)
     procedure, public, pass :: build_ind_diag
     !> Build the Diagonal Matrix array
     !> (public procedure for type spmat)
     procedure, public, pass :: get_diagonal
     !> Procedure to read a given row 
     !> (public procedure for type spmat)
     procedure, public, pass :: get_row
     !> Convert CSR mat. into SSR
     !> It can initilize a new matrix
     !> (public procedure for type spmat)
     procedure, public , pass :: csr2ssr
     !> Convert SSR mat. into CSR
     !> It can initilize a new matrix
     !> (public procedure for type spmat)
     procedure, public , pass :: ssr2csr
     !> Permute a matrix
     !> (public procedure for type spmat)
     procedure, public , pass :: perm => perm_mat
     !> Procedure to set to 'value' all the non-diagonal
     !> elements of row 'index'
     !> (public procedure for type spmat)
     procedure, public , pass :: set_row
     !> Procedure to fix the irow-term of the solution to 
     !> given value, cmodifieng the matrix and the rhs
     !> (public procedure for type spmat)
     procedure, public , pass :: set_rowcol
     !> Procedure to operate the diagonal scaling
     !> of matrix M by the diagonal matrix D
     !> out_matrix = left_D M right_D
     !> (public procedure for type spmat)
     procedure, public , pass :: diagonal_scale
     !> Procedure to  multiply sparse matrix M 
     !> by the diagonal matrix D 
     !> out_matrix = D M D
     !> (public procedure for type spmat)
     procedure, public , pass :: MxD
     !> Procedure to multiply  diagonal matrix D 
     !> by the sparse matrix M 
     !> out_matrix = D M
     !> (public procedure for type spmat)
     procedure, public , pass :: DxM
     !> Procedure to compute 
     !>         y = M * x 
     !> (public procedure for type spmat)
     procedure, public, pass :: matrix_times_vector => mxv_spmat
     !> Procedure to compute 
     !>         y_i = M(i,:) * x
     !> for csr matrix only
     !> (public procedure for type spmat)
     procedure, public, pass :: Mixv
     !> Procedure to compute 
     !>         y = M^T * x 
     !> with M^T the transposed of a matrix M
     !> (public procedure for type spmat)
     procedure, public, pass :: matrix_transpose_times_vector => MTxv_spmat
     !> Procedure to compute 
     !>         y = y + alpha* M * e_i
     !> where e_i is the vector v with all vero and v(i) = 1.0 
     !> (public procedure for type spmat)
     procedure, public, pass :: aMxei
     !> Procedure to get the index in the member coeff(or ja) 
     !> of the entries (irow,icol) 
     !> return zero if not found
     !> (public procedure for type spmat)
     procedure, public, pass :: iterm
     !> Logical output for dimensions matching
     !> (public procedure for type spmat)
     procedure, public, pass :: check => check_spmat
     !> Procedure to sort arrays ja and coeff
     !> so that the colums indeces of ja of each row
     !> (j = ja(ia(irow),ja(ia(irow+1)-1) 
     !> are in increasing order 
     !> (public procedure for type spmat)
     procedure, public, pass :: sort => sort_spmat
     !> Procedure to sort arrays the passed portion
     !> ja and coeff so that the colums indeces of ja
     !> are sorted in increasing order 
     !> (public procedure for type spmat)
     procedure, private, nopass :: sort_row
     !> Procedure for selecting and/or prermuting
     !> rows of csr matrix
     !> (public procedure for type spmat)
     procedure, public, pass :: select_permute_rows
     !> Procedure for selecting and/or prermuting
     !> columns of csr matrix
     !> (public procedure for type spmat)
     procedure, public, pass :: select_permute_columns
     !> Procedure for trasposition (operate in place)
     !> (public procedure for type spmat)
     procedure, public, pass :: transpose => transpose_spmat 

     !> Procedure to solve system
     !>           M x = b
     !> with M triagular
     !> (public procedure for type spmat)
     procedure, public, pass :: solve_triangular
     !> Procedure to solve system
     !>           M x = b
     !> with M triagular
     !> (public procedure for type spmat)
     procedure, public, pass :: solve_triangular_unitdiag
     !> Procedure Build the approximate Cholesky factorization 
     !> of a symmetric sparse matrix M in the form
     !>          M = U^T U
     !> (public procedure for type spmat)
     procedure, public, pass :: incomplete_cholesky
     !> Procedure Build the approximate factorization 
     !> of a sparse matrix M in the form
     !>          M = L U
     !> with L and U lower and upper triangular matrix
     !> (public procedure for type spmat)
     procedure, public, pass :: incomplete_lu
     !> Procedure Build the approximate factorization 
     !> of a sparse matrix M in the form
     !>          M = L U
     !> with L and U lower and upper triangular matrix
     !> (public procedure for type spmat)
     procedure, public, nopass :: MtimesN
     !> Procedure that adds a row to a matrix
     !> ( A | 0 )
     !> ( w | 1 )
     !> where A is the matrix, w is a vector (dim=nrow)
     procedure, public, pass :: add_row
     !> Procedure Build the reverse Cuthill-Mckee ordering 
     !> for the spmat contiang the connection of a
     !> a genreal graph.
     !> (public procedure for type spmat)
     procedure, public, pass :: genrcm
     !> Procedure to compute the maximum bandwidth 
     !> for a general matrix.
     !> (public procedure for type spmat)
     procedure, public, pass :: bandwidth
     !> Procedure to create ps file containg the
     !> sparsity pattern of a sparse matrix.
     !> (public procedure for type spmat)
     procedure, public, pass :: plot2ps
     !> Procedure to create sparsity pattern
     !> of the normal matrix
     !> (procedure public for type spmat)
     procedure, public, pass :: pattern_ATA
     !> Procedure to create product matrix
     !> (procedure public for type spmat)
     procedure, public, pass :: mult_MDN
     !> Procedure to initialize product matrix
     !>    A = M x N
     !> It only builds the sparsity pattern
     !> A%coeff=zero
     !> (procedure public for type spmat)
     procedure, public, pass :: init_product_spmat
     !> Procedure to create product matrix
     !> (procedure public for type spmat)
     procedure, public, pass :: assembly_redirector
     !> Procedure to integer pointer assembler
     !> used in assembly_ATA procedure
     !> (procedure public for type spmat)
     procedure, public, pass :: assembly_assemblerATA
     !> Procedure to assembly matrix ATA using
     !> matrix A and assembler built in procedure 
     !> assembly_assemblerATA
     !> (procedure public for type spmat)
     procedure, public, pass :: assembly_ATA
     !> Procedure to compute 
     !>         N = N + alpha * M
     !> with M, N sparse matrices. Work in Place
     !> (public procedure for type spmat)
     procedure, public, pass :: aMpN
     !> Procedure convert sparse block matrix
     !> into one sparse matrix
     procedure, public, pass :: form_new_linop
     !> Procedure convert diagonal matrix into a
     !> into a sparse matrix
     procedure, public, pass :: form_diag
     !> Procedure convert sparse block matrix
     !> into one sparse matrix
     procedure, public, pass :: form_pair_linop
     !> Procedure convert sparse block matrix
     !> into one sparse matrix
     procedure, public, pass :: form_block
     !> Procedure convert sparse block matrix
     !> into one sparse matrix
     !procedure, public, pass :: form_linear_combination
     !> Procedure convert sparse block matrix
     !> into one sparse matrix
     !procedure, public, pass :: form_composed
  end type spmat

  

contains
  !>-------------------------------------------------------------
  !> Static constructor.
  !> (procedure public for type spmat)
  !> Instantiate (allocate if necessary)
  !> and initilize (by also reading from input file)
  !> variable of type spmat
  !>
  !> usage:
  !>     call 'var'%init(lun_err, nrow,ncol, nterm,&
  !>                    storage_system ,[is_symmetric,triangular] )
  !>
  !> where:
  !> \param[in] lun_err               -> integer. Error logical unit
  !> \param[in] nrow                  -> integer. Number of rows 
  !> \param[in] ncol                  -> integer. Number of columns
  !> \param[in] nterm                 -> integer. Number of non-zero term
  !>                                     stored in the matrix
  !> \param[in] storage_system        -> character(len=3) Storage system
  !>                                     csr (strongly recomended)
  !>                                     ssr 
  !> \param[in] (optional) is_symmetric -> Logical. T/F flag for symmetric matrix
  !> \param[in] (optional) triangular   -> Character. 
  !>                                       'N' = not triangular ( the default)
  !>                                       'U' = uppertriangular
  !>                                       'L' = lower_triangular' 
  !<-------------------------------------------------------------
  subroutine init_spmat(this, lun_err, &
       nrow, ncol, nterm,&
       storage_system,&
       ! optional arguments
       is_symmetric, &
       triangular,unitary_diag)
    use Globals
    implicit none
    !var
    class(spmat),                intent(inout) :: this
    integer,                     intent(in   ) :: lun_err
    integer,                     intent(in   ) :: nrow
    integer,                     intent(in   ) :: ncol
    integer,                     intent(in   ) :: nterm
    character(len=*),            intent(in   ) :: storage_system
    logical,           optional, intent(in   ) :: is_symmetric
    character (len=1), optional, intent(in   ) :: triangular
    logical,           optional, intent(in   ) :: unitary_diag

    ! local vars
    integer :: res
    logical :: rc


    if (this%is_initialized) call this%kill(lun_err)

    
    this%is_initialized  = .true.
    this%nrow            = nrow
    this%ncol            = ncol
    this%nterm           = nterm

    allocate(this%ia(nrow+1),stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'init_spmat', &
         '  type spmat member ia (array)',res)
    

    allocate(this%ja(nterm),stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'init_spmat', &
         '  type spmat member ja (array)',res)


    allocate(this%coeff(nterm),stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'init_spmat', &
         '  type spmat member coeff (array)',res)
    

    
    this%storage_system = storage_system

    if ( .not. ( (this%storage_system .eq. 'ssr') .or. &
         (this%storage_system .eq. 'csr') ) ) then
       rc = IOerr(lun_err, err_inp, 'init_spmat', &
            ' format '//storage_system//' nor supported. Use csr or ssr' )
       
    end if
       
    
    if (this%storage_system .eq. 'ssr') then
       this%is_symmetric=.true.
       this%triangular  = 'N'
       if ( present(is_symmetric) ) then
          if ( (.not.is_symmetric) .and. this%storage_system.eq.'ssr') &
               rc = IOerr(lun_err, err_inp, 'init_spmat', &
               '  conflict between storage_system and is_symmetric)')
          this%is_symmetric   = is_symmetric
       end if

    else
       
       if ( present(is_symmetric) ) then
          this%is_symmetric   = is_symmetric
       end if

       if ( present(triangular) ) then
          if ( ( triangular .eq. 'N') .or. &
               ( triangular .eq. 'U') .or. &
               ( triangular .eq. 'L') ) then
             this%triangular = triangular
          else
             rc = IOerr(lun_err, err_inp, 'init_spmat', &
                  '  wrong triangular passed ='//etb(triangular))
          end if
       end if
    end if


    if ( present(unitary_diag) ) this%unitary_diag = unitary_diag

    
  end subroutine init_spmat

  subroutine spmat_constructor(this, lun_err, &
       nrow, ncol, nterm,&
       storage_system,&
       ! optional arguments
       is_symmetric, &
       triangular,unitary_diag)
    use Globals
    implicit none
    !var
    type(spmat),                intent(inout) :: this
    integer,                     intent(in   ) :: lun_err
    integer,                     intent(in   ) :: nrow
    integer,                     intent(in   ) :: ncol
    integer,                     intent(in   ) :: nterm
    character(len=*),            intent(in   ) :: storage_system
    logical,           optional, intent(in   ) :: is_symmetric
    character (len=1), optional, intent(in   ) :: triangular
    logical,           optional, intent(in   ) :: unitary_diag

    call this%init( lun_err, &
       nrow, ncol, nterm,&
       storage_system,&
       ! optional arguments
       is_symmetric=is_symmetric, &
       triangular=triangular,&
       unitary_diag=unitary_diag)

  end subroutine spmat_constructor
  

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
  subroutine kill_spmat(this, lun_err)
    implicit none
    ! vars
    class(spmat),intent(inout) :: this
    integer,     intent(in   ) :: lun_err
    ! local vars
    integer :: res
    logical :: rc

    ! Matrix not initialized. nothing to do. return
    if ( this%nrow .le. 0 ) return 

    deallocate(this%ia,stat=res)
    if (res.ne.0) rc=IOerr(lun_err, err_dealloc, 'kill_spmat', &
         'dealloc fail for spmat members ia',res)
    deallocate(this%ja,stat=res)
    if (res.ne.0) rc=IOerr(lun_err, err_dealloc, 'kill_spmat', &
         'dealloc fail for spmat members ja',res)

    deallocate(this%coeff,stat=res)
    if (res.ne.0) rc=IOerr(lun_err, err_dealloc, 'kill_spmat', &
         'dealloc fail for spmat member coeff',res)


    if ( allocated(this%ind_diag) ) then
       deallocate(this%ind_diag, stat=res)
       if(res .ne. 0) rc = IOerr(lun_err, err_dealloc, 'kill_spmat', &
            '  type spmat member ia (array)',res)

    end if

   
    !
    ! reset to default
    !
    call this%to_default()
    this%nterm       = 0
    

  end subroutine kill_spmat

  !>-------------------------------------------------------------
  !> Reading procedure.
  !> (public procedure for type spmat)
  !> Read content of a variable of type spamat
  !> 
  !> usage:
  !>     call 'var'%read(lun)
  !>
  !> where:
  !> \param[in] lun -> integer. I/O unit for error message output
  !>
  !<-------------------------------------------------------------
  subroutine read_spmat(this,lun_err,input_file,format)
    use Globals
    implicit none
    class(spmat),      intent(inout) :: this
    integer,           intent(in   ) :: lun_err
    type(file),        intent(in   ) :: input_file
    character(len=*), optional, intent(in) :: format
    ! loc. var
    logical :: rc
    integer :: res, lun
    integer:: iterm,icol,irow,irow_tmp,nrow_max,nterm_max,nterm,nrow,ncol
    character(len=256) :: inputs
    integer , allocatable :: ia(:),ja(:)
    real(kind=double) , allocatable :: coeff(:) 
    
    lun = input_file%lun

    if (this%is_initialized) call this%kill(lun_err)
    
    if (.not. present(format)) then
       read(lun,*) this%nrow  ! '! number of rows '
       read(lun,*) this%ncol  ! '! number of colums '
       read(lun,*) this%nterm      !' ! number of non-zero terms '
       allocate(this%ia(this%nrow+1),this%ja(this%nterm),&
            this%coeff(this%nterm), stat=res)
       if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'read_spmat', &
            '  type spmat member ia, ja, coeff (array)',res)
       irow = 1
       this%ia(1) = 1
       do iterm=1,this%nterm
          read(lun,*) irow_tmp,this%ja(iterm),this%coeff(iterm)
          if(irow_tmp.gt.irow) then
             irow = irow+1
             this%ia(irow) = iterm
          end if
       end do
       this%ia(this%nrow+1) = this%nterm+1

       this%is_initialized = .true.
    else
       if ( format .eq. 'matlab') then
          nrow_max=10000000
          nterm_max=100000000
          allocate(ia(nrow_max+1),ja(nterm_max),&
               coeff(nterm_max), stat=res)
          if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'read_spmat', &
               '  type spmat member ia, ja, coeff (array)',res)

          nterm=0
          nrow=0
          ncol=0
          do iterm=1,nterm_max
             read(lun,*,iostat=res) irow_tmp,ja(iterm),coeff(iterm)
             if (res.ne.0) exit
             nrow=max(nrow,irow_tmp)
             ncol=max(ncol,ja(iterm))
             nterm=nterm+1
             if(irow_tmp.gt.irow) then
                irow = irow+1
                ia(irow) = iterm
             end if
          end do

          call this%init(lun_err, &
               nrow, ncol, nterm,&
               'csr')
          this%ia(1:nrow+1)   = ia(1:nrow+1)
          this%ja(1:nterm)    = ja(1:nterm)
          this%coeff(1:nterm) = coeff(1:nterm)
                    
          
          deallocate(ia,ja, coeff, stat=res)
          if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'read_spmat', &
               '  type spmat member ia, ja, coeff (array)',res)


          
       end if
    end if
      
  end subroutine read_spmat

  !>-------------------------------------------------------------
  !> Writing procedure.
  !> (public procedure for type spmat)
  !> Prints content of a variable of type spamat
  !> 
  !> usage:
  !>     call 'var'%write(lun)
  !>
  !> where:
  !> \param[in] lun -> integer. I/O unit for error message output
  !>
  !<-------------------------------------------------------------
  subroutine write_spmat(this, lun, format)
    use Globals
    implicit none
    class(spmat),      intent(in) :: this
    integer,           intent(in) :: lun
    character(len=*), optional, intent(in) :: format
    ! loc. var
    integer i,j,m,n,ind
    character(len=256) :: out_format='fortran'
    
    if ( present(format) ) out_format =etb(format)
    
    select case (out_format)
    case ('fortran')
       write(lun,*) this%nrow, '! number of rows     '
       write(lun,*) this%ncol, '! number of columns     '
       write(lun,*) this%nterm,'! number of non-zero terms '
       do i=1,this%nrow
          m=this%ia(i)
          n=this%ia(i+1)-1
          do j=m,n
             write(lun,1010) i,this%ja(j),this%coeff(j)
          end do
       end do
    case ('matlab')
       do i=1,this%nrow
          m=this%ia(i)
          n=this%ia(i+1)-1
          write(lun,1010) (i,this%ja(ind),this%coeff(ind),ind=m,n)
       end do

    end select
1010   format(2i15,1pe24.16)       

  end subroutine write_spmat

  !>-------------------------------------------------------------
  !> Info procedure.
  !> (public procedure for type spmat)
  !> Prints content of a variable of type spamat
  !> 
  !> usage:
  !>     call 'var'%info(lun_out)
  !>
  !> where:
  !> \param[in] lun_out -> integer. I/O unit for output message
  !>
  !<-------------------------------------------------------------
  subroutine info_spmat(this, lun)
    use Globals
    implicit none
    class(spmat), intent(in) :: this
    integer, intent(in) :: lun
    ! loc. var
    character(len=256) :: state,mem,prop
    character(len=3) :: sep=' | '

    if ( this%is_initialized ) then
       if ( this%name .eq. 'empty') then
          write(state,'(a)') 'Sparse matrix allocated'
       else
          write(state,'(a,a,a)') 'Sparse matrix ',etb(this%name),' allocated'
       end if
       mem=this%storage_system
       if ( this%triangular .ne. 'N') then
          if ( this%triangular .eq. 'U' ) then
             write(prop,'(a)') 'upper-triangular'
          else 
             write(prop,'(a)') 'lower-triangular'
          end if
       else
          if ( this%is_symmetric ) then  
             write(prop,'(a)') 'symmetric'
          else
             write(prop,'(a)') 'non-symmetric'
          end if
       end if

       write(lun,'(a)') etb(etb(state)//sep//etb(mem)//sep//etb(prop))

       write(lun,'(a,I8,a,a,I8,a,a,I8)') 'nrows= ', this%nrow,sep,&      
            'ncol= ', this%ncol,sep,&      
            'non-zero terms= ', this%nterm     
    else
       write(lun,*) 'Spmat not initialized'
    end if

  end subroutine info_spmat


  !>-------------------------------------------------------------
  !> check procedure.
  !> (public procedure for type spmat)
  !> verifies if a matrix structure sastisfy the array bounds
  !> for given, nrow, ncol, nterm
  !> 
  !> usage:
  !>     sym = 'var'%check(nrow,ncol,nterm)
  !>
  !> \param[in] (optional) nrow          -> integer. Nmb. of rows required
  !> \param[in] (optional) ncol          -> integer. Nmb. of columnss required
  !> \param[in] (optional) nterm         -> integer. Nmb. of non-zero 
  !>                                       terms required
  !> \param[in] (optional)            -> integer. Storage system
  !> \param[in] (optional) is_symmetric  -> Logical. pass ".true." to  
  !>                                        ask if spmat is symmetric
  !> \param[in] (optional) triangular         -> Character(len=1). pass "U" or "L"
  !>                                        to ask if spmat is upper or lower 
  !>                                        triangular
  !> \return    check_spmat -> logical. true/false flag if Spmat respects
  !>                                   dimension bound   e               
  !<------------------------------------------------------------
  function check_spmat(this, &
       nrow, ncol, nterm,&
       storage_system,is_symmetric,triangular)
    use Globals
    implicit none
    class(spmat),                intent(in) :: this
    integer,           optional, intent(in) :: nrow
    integer,           optional, intent(in) :: ncol
    integer,           optional, intent(in) :: nterm
    character(len=3),  optional, intent(in) :: storage_system
    logical          , optional, intent(in) :: is_symmetric
    character (len=1), optional, intent(in) :: triangular


    logical                       :: check_spmat

    check_spmat = .true.


    if ( present(nrow) .and. (nrow .ne. this%nrow ) ) then
       check_spmat=.false. 
    end if

    if ( present(ncol) .and. (ncol .ne. this%ncol ) ) then
       check_spmat=.false. 
    end if

    if ( present(nterm) .and. (nterm .ne. this%nterm) ) then
       check_spmat = .false.
    end if


    if ( present(storage_system) .and. &
         (storage_system .ne. this%storage_system) ) then
       check_spmat = .false.
    end if

  
    if ( ( present(is_symmetric)  ) .and. &
         ( .not. this%is_symmetric ) ) then
       check_spmat = .false.
    end if

    if ( ( present (triangular)  ) .and. &
         ( this%triangular .ne. triangular  ) ) then
       check_spmat = .false.
    end if




  end function check_spmat


  !>-------------------------------------------------------------
  !> Procedure for identify the pointer the diagonal 
  !> term(irow) of a matrix 
  !>
  !> usage : this%idiag(irow)
  !> 
  !> where:
  !> \param[in] irow  -> integer. Index of the row
  !> \return    idiag -> integer. Index of the diagonal term of the
  !>                              irow^th row
  !>                              idiag=0 if the diagonal term is null
  !>--------------------------------------------------------------
  function idiag( this, irow )
    use Globals
    implicit none
    class(spmat), intent(in) :: this
    integer,      intent(in) :: irow
    integer :: idiag
    !local 
    integer :: j

    idiag=0
    if (allocated(this%ind_diag) ) then
       idiag = this%ind_diag(irow)
    else
       if ( this%storage_system .eq. 'ssr' ) then
          !
          ! use ssr structure
          !
          idiag = this%ia(irow)

          if ( this%ja(idiag) .ne. irow ) idiag = 0 
       else if ( this%storage_system == 'csr' ) then 
          if ( this%triangular .ne. 'N' ) then
             ! use triangular structure
             select case( this%triangular ) 
             case ('L') 
                idiag = this%ia(irow+1)-1
             case ('U') 
                idiag = this%ia(irow)
             end select
             if ( this%ja(idiag) .ne. irow ) idiag = 0
          else
             ! search in row
             j = this%ia(irow)
             idiag = 0
             do j = this%ia(irow),this%ia(irow+1) - 1
                if (this%ja(j) .eq. irow) then
                   idiag=j
                   exit
                end if
             end do
          end if
       end if
    end if
  end function idiag


  !>-------------------------------------------------------------
  !> Procedure to initilized and build the member ind_diag
  !> of structure spamt     
  !>
  !> usage : this%build_diag(lun_err)
  !>
  !> \param[in] lun_err -> integer. I/O unit for err. msg. 
  !>--------------------------------------------------------------
  subroutine build_ind_diag(this,lun_err)
    use Globals
    implicit none
    class(spmat),      intent(inout) :: this
    integer,           intent(in   ) :: lun_err
    !local 
    logical rc
    integer res
    integer i
    integer :: diag_loc(this%nrow)
    
    do i = 1,this%nrow
       diag_loc(i) = this%idiag(i)
    end do

    if ( .not. allocated(this%ind_diag)) then
       allocate (this%ind_diag(this%nrow),stat=res) 
       if(res .ne. 0) rc = IOerr(lun_err, err_alloc, &
            'build_ind_diag', &
            'type spmat member ind_diag',res)
    end if
    this%ind_diag = diag_loc

  end subroutine build_ind_diag

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
    class(spmat),      intent(in   ) :: this
    real(kind=double), intent(inout) :: diagonal(min(this%nrow,this%ncol))
    !local 
    integer i, ind
    
    diagonal=zero   
    do i=1,min(this%nrow,this%ncol)
       ind = this%idiag(i)
       if (ind .ne. 0 ) diagonal(i) = this%coeff(ind)
    end do

  end subroutine get_diagonal


  !>------------------------------------------------------------------
  !> Procedure to transform matrix structure from ssr into csr format.
  !> Reinitialized the matrix.
  !>
  !> usage: call var%ssr2csr(lun_err)
  !>
  !> where:
  !> \param[in ] lun_err -> integer. I\O unit for error message
  !>------------------------------------------------------------------   
  subroutine ssr2csr(this,lun_err)
    use Globals
    implicit none
    class(spmat), intent(inout) :: this
    integer,      intent(in   ) :: lun_err

    !local
    integer i, j ,m
    integer nrow, nterm, nterm_out
    type(spmat) :: upper, lower

    nrow =this%nrow
    nterm=this%nterm
    nterm_out=2*nterm-nrow

    ! create a copy and change to csr 
    select type (this)
    type is (spmat) 
       upper = this
    end select
    upper%storage_system = 'csr'
    upper%triangular     = 'U'

    ! create workin matrix
    ! We used the following trick: change the storage system
    ! to upper to build a lower triangular matrix scr_transpose
    lower =  upper
    call lower%transpose(lun_err)

    ! renitialized 
    call this%kill(lun_err)
    call this%init(lun_err, &
         nrow, nrow, nterm_out,&
         'csr',&
         is_symmetric=.true.)

    ! the matrix is build this=upper+lower-diag
    m=0
    this%ia(1)=1
    do i = 1,nrow
       do j = lower%ia(i),lower%ia(i+1)-2              
          m=m+1
          this%ja(m) = lower%ja(j)
          this%coeff(m) = lower%coeff(j)              
       end do
       do j = upper%ia(i),upper%ia(i+1)-1
          m=m+1
          this%ja(m) = upper%ja(j)
          this%coeff(m) = upper%coeff(j)
       end do
       this%ia(i+1) = m+1
    end do


    ! free memory 
    call upper%kill(lun_err)
    call lower%kill(lun_err)


  end subroutine ssr2csr

  !>---------------------------------------------------------------
  !> Procedure to transform a symmetric matrix from csr into ssr
  !> format. Initialize the varaible mat_out if it is not
  !> 
  !> usage: call var%ssr2csr(lun_err,mat_out)
  !>
  !!> where:
  !> \param[in ] lun_err            -> integer. I\O unit for error message
  !> \param[out] (optional) mat_out -> type(spmat). Matrix in csr format
  !<---------------------------------------------------------------------
  subroutine csr2ssr(this, lun_err)
    use Globals
    implicit none
    class(spmat), intent(inout) :: this
    integer,      intent(in   ) :: lun_err

    !local
    integer :: i, j , m
    integer :: nrow, ncol, nterm_out
    type(spmat) :: ssr        

    if ( .not. this%is_symmetric  ) then
       write(lun_err,*) ' Not symmetric matrix'
       write(lun_err,*) ' Convertion not possible'
       stop
    end if

    nrow      = this%nrow
    ncol      = this%ncol
    nterm_out = (this%nterm + nrow) / 2

    ! initialized and set properteis of working spmat
    call ssr%init(lun_err,nrow,ncol, nterm_out, &
         'ssr')
    
    m=0
    ssr%ia(1)=1
    do i = 1, nrow       
       do j = this%idiag(i), this%ia(i+1)-1
          m=m+1
          ssr%ja(m)    = this%ja(j)
          ssr%coeff(m) = this%coeff(j)
       end do
       ssr%ia(i+1) = m + 1
    end do

    ! copy the working spmat into this
    select type (this)
    type is (spmat) 
       this = ssr
    end select

    ! free memory
    call ssr%kill(lun_err)

  end subroutine csr2ssr

    !>-----------------------------------------------------------
    !> Procedure to build the matrix PAP' for given sparse matrix 
    !> A and a permutation perm. It initializes spmat mat_out.
    !>
    !> usage : call var%perm(lun_err,nrow,perm,iperm)
    !>
    !> where:
    !> \param[in ] lun_err -> integer. I\O unit for error message
    !> \param[in ] nrow    -> integer. Nmb of rows and columns
    !> \param[in ] perm    -> integer(nrow). Permutation
    !> \param[in ] iperm   -> integer(nrow). Inverse Permutation
    !>-----------------------------------------------------------
    subroutine perm_mat(this, lun_err,  perm , inv_perm )
      use Globals
      implicit none
      class(spmat), intent(inout)  :: this
      integer,      intent(in   )  :: lun_err
      integer,      intent(in   )  :: perm(this%nrow)
      integer,      intent(in   )  :: inv_perm(this%nrow)

      !local
      logical :: rc
      integer :: res,info
      integer :: i
      type(spmat) :: full, sorted

      if ( this%storage_system .eq. 'ssr' ) then
         ! local copy in ssr format
         select type (this)
         type is (spmat)
            full = this
            call full%ssr2csr(lun_err)
         end select

         ! create a sorted the local copy
         call sorted%init(lun_err,&
              full%nrow,full%nrow,full%nterm,'csr')
         
         call apply_perm(lun_err,&
              full%nrow,full%nterm,&
              perm,inv_perm,&
              full%ia,full%ja,full%coeff,&
              sorted%ia,sorted%ja,sorted%coeff,&
              info)

         
         
         call full%kill(lun_err)
         ! convert sorted mat into ssr format and copy into this
         call sorted%csr2ssr(lun_err) 
         select  type (this)
         type is (spmat)
            this = sorted
         end select
         call sorted%kill(lun_err)

         

      else         
         call sorted%init(&
              lun_err,this%nrow,this%nrow,this%nterm,'csr')
         
         call apply_perm(lun_err,&
              this%nrow,this%nterm,&
              perm,inv_perm,&
              this%ia,this%ja,this%coeff,&
              sorted%ia,sorted%ja,sorted%coeff,&
              info)    
         select type (this)
         type is (spmat)
            this = sorted
         end select
         call sorted%kill(lun_err)
         
      end if
      if (info .ne. 0) then
         write(lun_err,*) 'Error in procedure perm_mat'
         stop
      end if
      
    

    end subroutine perm_mat
    
    subroutine apply_perm(lun_err,&
           nn,nt,&
           perm,iperm,&
           iat_in,ja_in,coef_in,&
           iat_out,ja_out,coef_out,&
           info)
        !--------------------------------------
        !
        ! Creates a permuted matrix
        !
        !--------------------------------------

        use Globals

        implicit none

        ! Interfaces
        integer                        :: lun_err
        integer, intent(in)            :: nn,nt
        integer, intent(in)            :: perm(nn),iperm(nn)
        integer, intent(in)            :: iat_in(nn+1),ja_in(nt)
        real(kind=double), intent(in)  :: coef_in(nt)

        integer, intent(out)           :: info
        integer, intent(out)           :: iat_out(nn+1),ja_out(nt)
        real(kind=double), intent(out) :: coef_out(nt)

        integer                        :: i,j,k,m1,m2
        integer, allocatable           :: irow_tmp(:)
        !
        integer, allocatable           :: iat_tmp(:),ja_tmp(:)
        real(kind=double), allocatable :: coef_tmp(:)
        !
        !local
        logical :: rc
        integer :: res
        info = 0

        allocate(irow_tmp(nt),iat_tmp(nn+1),ja_tmp(nt),coef_tmp(nt),stat=res)
        if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'init_spmat', &
             '  type spmat member ia (array)',res)

        k = 0
        do i = 1,nn
           m1 = iat_in(iperm(i))
           m2 = iat_in(iperm(i)+1)-1
           do j = m1,m2
              k = k + 1
              irow_tmp(k) = i
              ja_out(k) = perm(ja_in(j))
              coef_out(k) = coef_in(j)
           end do
        end do
        call irow2iat(nn,nt,irow_tmp,iat_out)
        call TranspMat(nn,nt,&
             iat_out,ja_out,coef_out,&
             nn,k,iat_tmp,ja_tmp,coef_tmp,&
             irow_tmp)
        call TranspMat(nn,nt,&
             iat_tmp,ja_tmp,coef_tmp,&
             nn,k,iat_out,ja_out,coef_out,&
             irow_tmp)

        deallocate(irow_tmp,iat_tmp,ja_tmp,coef_tmp,stat=info)
        if (info .ne. 0) return

      end subroutine apply_perm
      subroutine irow2iat(n,nterm,irow,iat)
        !-----------------------------------------------------------------
        !
        !  Subroutine: irow2iat
        !
        !  Coded by Carlo Janna
        !  April 2012
        !
        !  Purpose: Build the topology vector iat 
        !  from a row indices list irow
        !
        !  Variables:
        !
        !  n     : # of rows of a given matrix stored in CSR format
        !  nterm : # of non-zeroes of a given matrix stored in CSR format
        !  irow  : row indices of non-zeroes of a matrix in CSR format
        !  iat   : integer array of the pointers to the begin of each row
        !
        !----------------------------------------------------------------

        implicit none

        !Input variables
        integer, intent(in) :: n,nterm
        integer, intent(in) :: irow(nterm)

        !Output variables
        integer, intent(out) :: iat(n+1)

        !Local variables
        integer             :: j,k,irow_old,irow_new

        !------------------------------------------------------------

        irow_old = 0
        do k=1,nterm
           irow_new = irow(k)
           if ( irow_new .gt. irow_old ) then
              do j = irow_old+1,irow_new
                 iat(j)=k
              enddo
              irow_old = irow_new
           end if
        end do
        k = nterm+1
        do j = irow_old+1,n+1
           iat(j) = k
        enddo

      end subroutine irow2iat
      subroutine TranspMat(n,nterm,&
         iat,ja,mat,&
         m,nterm_out,&
         iat_T,ja_T,mat_T,&
         IW)
      !----------------------------------------------------------------
      !
      !  Subroutine: TranspMat
      !
      !  Coded by Carlo Janna
      !  September 2011
      !
      !  Purpose: Transpose a (n x m) matrix.
      !
      !---------------------------------------------------------------
      use Globals
      
      implicit none
      ! Input variables
      integer, intent(in)            :: n,m,nterm
      integer, intent(in)            :: iat(n+1),ja(nterm)

      ! Output variables
      integer, intent(out)           :: nterm_out
      integer, intent(out)           :: iat_T(m+1),ja_T(nterm)
      real(kind=double), intent(in)  :: mat(nterm)
      real(kind=double), intent(out) :: mat_T(nterm)

      ! Work array (given as input)
      integer, intent(out)           :: IW(m+1)

      ! Local variables
      integer                        :: i,j,ind

      ! Initialize pointers
      do i = 1,m+1
         iat_T(i) = 0
      enddo

      ! Count non-zeroes for each column of the input matrix
      do i = 1,n
         do j = iat(i),iat(i+1)-1
            iat_T(ja(j)) = iat_T(ja(j)) + 1
         enddo
      enddo

      ! Reset pointers
      IW(1) = 1
      do i = 2,m+1
         IW(i) = IW(i-1) + iat_T(i-1)
      enddo
      call SCOPY(m+1,IW,1,iat_T,1)


      ! Transpose coeffiecients
      do i = 1,n
         do j = iat(i),iat(i+1)-1
            ind = IW(ja(j))
            ja_T(ind) = i
            mat_T(ind) = mat(j)
            IW(ja(j)) = ind+1
         enddo
      enddo

      nterm_out = nterm

    end subroutine TranspMat

  !>-------------------------------------------------------------
  !> Procedure for reset a row and optionally compute
  !> the correction of rhs in order to fix the solution
  !> at irow with value sol_irow
  !> (public procedure for type spmat, works in place)
  !> 
  !> usage call var%fix_sol_irow(irow, nrow, sol_irow,rhs)
  !> 
  !> where:
  !> \param[in ] irow              -> integer. Index of 
  !>                                     row where fix the diag.
  !> \param[in ] nrow              -> integer. Nmb of eqs
  !> \param[in ] (optin.) sol_irow -> real. Sol in irow
  !> \param[out] (optin.) rhs      -> real(nrow). rhs to correct
  !>----------------------------------------------------------------
  function iterm(this,irow,icol)  result (ind)
    use Globals
    implicit none   

    class(spmat),       intent(inout) :: this
    integer,            intent(in   ) :: irow
    integer,            intent(in   ) :: icol
    integer :: ind
    !local
    integer :: i1,i2,nel,pos,irow_loc,icol_loc,j
    
    !
    ! set the default output to zero
    !
    ind = 0
    
    select case (this%storage_system)
    case ('ssr')
       !
       ! switch indeces
       !
       if ( icol > irow ) then
          irow_loc = icol
          icol_loc = irow
       else
          irow_loc = irow
          icol_loc = icol
       end if
       !
       ! set bound and lenght of row indeces
       !
       i1=this%ia(irow_loc)
       i2=this%ia(irow_loc+1)-1
       nel=i2-i1+1
       
       !
       ! find icol-element position
       !
       pos=ifind(nel,this%ja(i1:i2),icol_loc)
       ind=i1+pos-1
    case ('csr')
       !
       ! set bound and lenght of row indeces
       !
       i1=this%ia(irow)
       i2=this%ia(irow+1)-1
       nel=i2-i1+1

       !
       ! find icol-element position
       !
       pos=ifind(nel,this%ja(i1:i2),icol)
       ind=i1+pos-1
       !write(0,*) 'i,j',irow,icol,'pos,nel',pos, nel,'term', ind, 'check',this%ja(ind)-icol
    end select
  end function iterm

  !>-------------------------------------------------------------
  !> Procedure to set a 'value' to all the non-diagonal
  !> elements of row 'index'
  !> (public procedure for type spmat, works in place)
  !> 
  !> usage call var%set_row(index, value)
  !> 
  !> where:
  !> \param[in ] index  -> integer. Index of row to be fixed
  !> \param[in ] value  -> real. value of the constant 
  !>----------------------------------------------------------------
  subroutine set_row(this,index,alpha)
    use Globals
    implicit none   
    class(spmat),       intent(inout) :: this
    integer,            intent(in   ) :: index
    real(kind=double),  intent(in   ) :: alpha

    !local 
    integer :: i,j,ind

    do j = this%ia(index), this%ia(index+1)-1
       this%coeff(j)=alpha
    end do

  end subroutine set_row
  
  !>-------------------------------------------------------------
  !> Procedure to set a 'value' to all the non-diagonal
  !> elements of row and column 'index'
  !> (public procedure for type spmat, works in place)
  !> 
  !> usage call var%set_rowcol(index, value)
  !> 
  !> where:
  !> \param[in ] index  -> integer. Index of row to be fixed
  !> \param[in ] value  -> real. value of the constant 
  !>----------------------------------------------------------------
  subroutine set_rowcol(this,index,alpha)
    use Globals
    implicit none   

    class(spmat),       intent(inout) :: this
    integer,            intent(in   ) :: index
    real(kind=double),  intent(in   ) :: alpha

    !local 
    integer :: i,j,ind
    
    if ( this%is_symmetric ) then
       select case ( this%storage_system) 
       case ('ssr') 
          !
          ! brute force researchs for upper column terms
          !
          do i=1, index
             do j = this%ia(i), this%ia(i+1)-1
                if (this%ja(j) .eq. index) this%coeff(j)=alpha
             end do
          end do
          !
          ! set to alpha current coefficient A(index,:)  
          !
          do j=this%ia(index), this%ia(index+1)-1
             this%coeff(j) = alpha
          end do
       case ('csr' )
          !
          ! set to alpha current coefficient A(index,:)  
          !
          do j=this%ia(index), this%ia(index+1)-1
             this%coeff(j) = alpha
          end do

          !
          ! find the index of the term A(j,index) 
          !
          do j=this%ia(index), this%ia(index+1)-1
             ind=this%iterm(this%ja(j),index)
             this%coeff(ind) = alpha
          end do
       end select
    else  
       !
       ! brute force researchs
       !
       do i=1, this%nrow
          do j = this%ia(i), this%ia(i+1)-1
             if (this%ja(j) .eq. index) this%coeff(j)=alpha
          end do
       end do
    end if
    
  end subroutine set_rowcol

  !>-------------------------------------------------------------
  !> Procedure to compute Matrix vector product
  !>         vec_out = vec_out+ alpha * (M) times (e_i)
  !> where is the vector
  !>   e_i=(0,.....,1,.....0)
  !>                ^
  !>                i
  !> (public procedure for type spmat)
  !> 
  !> usage:
  !>     call 'var'%Mxv(index,alpha,vec_out)
  !>
  !> where 
  !> \param[in   ] lun_err -> integer. I/O unit for err. msg. 
  !> \param[in   ] index   -> integer. index of the vector e_i
  !> \param[in   ] alpha   -> real. scalar multiplinge vector e_i
  !> \param[inout] vec_out -> real. dimension('var'%nrow)
  !>                            vector (M) times (vec_in) 
  !<-------------------------------------------------------------
  subroutine aMxei(this,lun_err,index,alpha,vec_out)
    use Globals
    implicit none
    class(spmat),      intent(inout) :: this
    integer,           intent(in   ) :: lun_err
    integer,           intent(in   ) :: index
    real(kind=double), intent(in   ) :: alpha
    real(kind=double), intent(inout) :: vec_out(this%nrow)

    ! local 
    logical :: rc
    integer :: stat,res
    integer :: j, m, mm,info
    real(kind=double), allocatable:: scr(:)
    
    if ( this%is_symmetric) then    
       select case(this%storage_system) 
       case('ssr') 
          m=this%ia(index)
          mm=this%ia(index+1)-1
          vec_out(index)=vec_out(index)+this%coeff(m)*alpha
          do j=m+1,mm
             vec_out(index)=vec_out(index)+this%coeff(j)*alpha
             vec_out(this%ja(j))=vec_out(this%ja(j))+this%coeff(j)*alpha
          end do
       case('csr') 
          m=this%ia(index)
          mm=this%ia(index+1)-1
          do j=m,mm
             vec_out(this%ja(j))=vec_out(this%ja(j))+this%coeff(j)*alpha
          end do
       end select
    else
       allocate(scr(this%nrow+this%ncol),stat=res)
       if (res.ne. 0) rc = IOerr(lun_err, err_alloc, 'aMxei', &
            'scratch array',res)
       !
       ! define vector alpha * e_i
       !
       scr=zero
       scr(index)=alpha
       call this%Mxv(&
            scr(1:this%ncol),&
            scr(1+this%ncol:this%nrow+this%ncol))
       vec_out=vec_out+scr(1+this%ncol:this%nrow+this%ncol)
       deallocate(scr,stat=res)
       if (res.ne. 0) rc = IOerr(lun_err, err_dealloc, 'aMxei', &
            'scratch array',res)
    end if

  end subroutine aMxei
  

  !>-------------------------------------------------------
  !> Procedure to operate diagonal scaling of spmat
  !>     M => left_D M right_D
  !> where D(M) is a diagonal matrix
  !> (pubblic procedure for type spmat, operate in place)
  !> 
  !> usage: var%diag_scale(lun_err,diag_matrix)
  !>
  !> where: 
  !>\param[in] lun_err     -> integer. I/O unit for err. msg. 
  !>\param[in] diag_matrix -> real(nrow). Diagonal mat. values
  !<-------------------------------------------------------
  subroutine diagonal_scale(this, lun_err, diagonal)
    use Globals
    implicit none
    class(spmat),      intent(inout) :: this
    integer,           intent(in   ) :: lun_err
    real(kind=double), intent(in   ) :: diagonal(this%ncol)
    !local
    logical rc
    integer res,i,j
    real(kind=double) ::  di

    if ( this%ncol .ne. this%nrow ) then 
       rc = IOerr(lun_err, err_val, 'diagonal_scale', &
            'spmat is not squared',res)
    end if


    do i = 1, this%nrow
       di = diagonal(i)
       do j = this%ia(i), this%ia(i+1)-1
          this%coeff(j) = &
               di * this%coeff(j) * diagonal( this%ja(j) )
       end do
    end do


  end subroutine diagonal_scale


  !>-------------------------------------------------------
  !> Procedure to operate matrix / matrix product.
  !>   B = D * A
  !> where D is a diagonal matrix and A is a sparse matrix
  !> (pubblic procedure for type spmat, operate in place)
  !> It converts the storage scheme in case of 
  !> A=symmetric matrix in ssr
  !> 
  !> 
  !> usage: var%DxM(lun_err,diag)
  !>
  !> where: 
  !>\param[in] lun_err -> integer. I/O unit for err. msg. 
  !>\param[in] diag    -> real. dimension(nequ). 
  !>                         Diagonal mat. values
  !<-------------------------------------------------------
  subroutine DxM( this, lun_err, diag)
    use Globals
    implicit none
    class(spmat),       intent(inout) :: this
    integer,            intent(in   ) :: lun_err
    real(kind= double), intent(in   ) :: diag(this%nrow)
    !local 
    integer i
    real(kind=double) :: diag_scal

    if (this%storage_system .eq. 'ssr' ) then
       call this%ssr2csr(lun_err)
    end if
    do i=1, this%nrow
       diag_scal = diag(i)
       this%coeff( this%ia(i): this%ia(i+1)-1) =  &
            diag_scal * &
            this%coeff( this%ia(i):this%ia(i+1)-1 )
    end do

    !
    ! in general this is not true
    ! 
    this%is_symmetric = .false.
    this%unitary_diag = .false.
  end subroutine DxM

  !>-------------------------------------------------------
  !> Procedure to operate matrix / matrix product.
  !>   B = A * D
  !> where D is a diagonal matrix and A is a sparse matrix
  !> (pubblic procedure for type spmat, operate in place)
  !> It converts the storage scheme in case of 
  !> A=symmetric matrix in ssr
  !> 
  !> 
  !> usage: var%DxM(lun_err,diag)
  !>
  !> where: 
  !>\param[in] lun_err -> integer. I/O unit for err. msg. 
  !>\param[in] diag    -> real. dimension(nequ). 
  !>                         Diagonal mat. values
  !<-------------------------------------------------------
  subroutine MxD( this, lun_err, diag)
    use Globals
    implicit none
    class(spmat),       intent(inout) :: this
    integer,            intent(in   ) :: lun_err
    real(kind= double), intent(in   ) :: diag(this%ncol)
    !local 
    integer i, j

    if (this%storage_system .eq. 'ssr' ) then
       call this%ssr2csr(lun_err)
    end if

    do i=1, this%nrow
       do j = this%ia(i), this%ia(i+1)-1
          this%coeff(j) =  &
               diag(this%ja(j)) * this%coeff( j )
       end do
    end do
    

    !
    ! in general this is not true
    ! 
    this%is_symmetric = .false.
    this%unitary_diag = .false.
  end subroutine MxD



  !>-------------------------------------------------------------
  !> Procedure to compute Matrix vector product
  !>         vec_out = (M) times (vec_in)
  !> (public procedure for type spmat)
  !> 
  !> usage:
  !>     call 'var'%matrix_times_vector(vec_in,vec_out,[info])
  !>
  !> where 
  !> \param[in   ] vec_in          -> real. dimension('var'%ncol)
  !>                                  vector to be multiplied
  !> \param[inout] vec_out         -> real. dimension('var'%nrow)
  !>                                  vector (M) times (vec_in) 
  !> \param[in   ]  info -> integer. Info number
  !>                                  in case of error   
  !<-------------------------------------------------------------
  subroutine Mxv_spmat(this,vec_in,vec_out, info,lun_err)
    use Globals
    implicit none
    class(spmat),      intent(inout) :: this
    real(kind=double), intent(in   ) :: vec_in(this%ncol)
    real(kind=double), intent(inout) :: vec_out(this%nrow)
    integer,           intent(inout) :: info
    integer,           intent(in   ) :: lun_err


    info=0

    select case(this%storage_system) 
    case('ssr') 
       call axbsym(this%nrow,this%ncol,this%nterm,this%ia,this%ja,&
            this%coeff,vec_in,vec_out)
    case('csr') 
       call axbnsym(this%nrow,this%ncol,this%nterm,this%ia,this%ja,&
            this%coeff,vec_in,vec_out)
    end select
  contains

    ! vec_out = M * vec_in 
    subroutine axbsym(nrow,ncol,nterm,ia,ja,coef1,xvec,bvec)
      use Globals
      implicit none
      integer :: nrow,ncol,nterm
      integer :: ja(nterm),ia(nrow+1)
      real(kind=double) ::coef1(nterm),xvec(ncol),bvec(nrow)
      integer :: i,k,m,mm

      bvec=zero
      do k=1,nrow
         m=ia(k)
         mm=ia(k+1)-1
         bvec(k)=bvec(k)+coef1(m)*xvec(ja(m))
         do i=m+1,mm
            bvec(k)=bvec(k)+coef1(i)*xvec(ja(i))
            bvec(ja(i))=bvec(ja(i))+coef1(i)*xvec(k)
         end do
      end do

    end subroutine axbsym

    ! vec_out = M * vec_in 
    subroutine axbnsym(nrow,ncol,nterm,ia,ja,coef1,xvec,bvec)
      use GLobals
      implicit none
      integer nrow,ncol,nterm
      integer i,k
      real(kind=double)  :: coef1(nterm),xvec(ncol),bvec(nrow)
      integer ja(nterm),ia(nrow+1)

      bvec = zero
      do k=1,nrow
         do i=ia(k),ia(k+1)-1
            bvec(k)=bvec(k)+coef1(i)*xvec(ja(i))
         end do
      end do
    end subroutine axbnsym

  end subroutine Mxv_spmat

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
  !> \param[in   ] (optional) info -> integer. Info number
  !>                                  in case of error   
  !<-------------------------------------------------------------
  subroutine Mixv(this,irow,vec_in,y_out, info,lun_err)
    use Globals
    implicit none
    class(spmat),      intent(inout) :: this
    integer,           intent(in   ) :: irow
    real(kind=double), intent(in   ) :: vec_in(this%ncol)   
    real(kind=double), intent(inout) :: y_out
    integer, optional, intent(inout) :: info
    integer, optional, intent(in   ) :: lun_err
    !local
    integer :: i

    if (present(info)) info=0

    select case(this%storage_system) 
    case('ssr') 
       stop
    case('csr') 
       y_out = zero
       do i=this%ia(irow),this%ia(irow+1)-1
          y_out=y_out+this%coeff(i)*vec_in(this%ja(i))
       end do
    end select

  end subroutine Mixv


  !>-------------------------------------------------------------
  !> Matrix transope-vector multiplication procedure
  !>         vec_out = (M) times (vec_in)
  !> (public procedure for class abs_matrix)
  !> 
  !> usage:
  !>     call 'var'%MTxv(vec_in,vec_out,[info])
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
  subroutine MTxv_spmat(this,vec_in,vec_out, info,lun_err)
    use Globals
    implicit none
    class(spmat),      intent(inout) :: this
    real(kind=double), intent(in   ) :: vec_in(this%nrow)
    real(kind=double), intent(inout) :: vec_out(this%ncol)
    integer,           intent(inout) :: info
    integer,           intent(in   ) :: lun_err

    select case ( this%storage_system ) 
    case ('ssr')
       call axbsym(this%nrow,this%ncol,this%nterm,this%ia,this%ja,&
            this%coeff,vec_in,vec_out)
    case('csr')
       call axbnsym_transpose(this%nrow,this%ncol,&
            this%nterm,this%ia,this%ja,&
            this%coeff,vec_in,vec_out)
    end select
  contains

    ! vec_out = M * vec_in 
    subroutine axbsym(nrow,ncol,nterm,ia,ja,coef1,xvec,bvec)
      use Globals
      implicit none
      integer :: nrow,ncol,nterm
      integer :: ja(nterm),ia(nrow+1)
      real(kind=double) ::coef1(nterm),xvec(ncol),bvec(nrow)
      integer :: i,k,m,mm

      bvec=zero
      do k=1,nrow
         m=ia(k)
         mm=ia(k+1)-1
         bvec(k)=bvec(k)+coef1(m)*xvec(ja(m))
         do i=m+1,mm
            bvec(k)=bvec(k)+coef1(i)*xvec(ja(i))
            bvec(ja(i))=bvec(ja(i))+coef1(i)*xvec(k)
         end do
      end do

    end subroutine axbsym

    ! vec_out = M^t * vec_in 
    subroutine axbnsym_transpose(nrow,ncol,nterm,ia,ja,coef1,xvec,bvec)
      use GLobals
      implicit none
      integer nrow,ncol,nterm
      integer i,k,m,mm
      real(kind=double)  coef1(nterm),xvec(nrow),bvec(ncol)
      integer ja(nterm),ia(nrow+1)

      bvec = zero
      do k=1,nrow
         m=ia(k)
         mm=ia(k+1)-1
         do i=m,mm
            bvec(ja(i))=bvec(ja(i))+coef1(i)*xvec(k)
         end do
      end do
    end subroutine axbnsym_transpose

  end subroutine MTxv_spmat


  subroutine transpose_spmat(this, lun_err,perm)
    use Globals
    implicit none
    ! Input variables
    class(spmat),     intent(inout) :: this
    integer,          intent(in   ) :: lun_err
    integer,optional, intent(inout) :: perm(this%nterm)

    ! Local variables
    logical :: rc
    integer :: res
    integer :: i,j,ind
    character(len=1) :: new_triangular
    integer, allocatable :: IW(:),permutation(:)
    type(spmat) :: transpose
    
    ! check if the matrix id symmetric
    if ( this%is_symmetric ) return
    ! check if the matrix uses csr storage
    if ( this%storage_system  .ne. 'csr' ) then
       write(lun_err,*) 'Matrix not in csr format'
       stop
    end if
    ! wok array
    allocate (&
         IW(this%ncol+1),&
         permutation(this%nterm),&
         stat=res) 
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'transpose matrix', &
         'work array IW',res)


    ! set properties of transpose matrix
    select case (this%triangular)
    case ('N')
       new_triangular = 'N'
    case ('U')
       new_triangular = 'L'
    case ('L')
       new_triangular = 'U'
    end select

    ! initialized working matrix
    call transpose%init(lun_err,&
         this%ncol,this%nrow,& 
         this%nterm,&
         'csr',&
         triangular=new_triangular)  

    ! Initialize pointers
    transpose%ia = 0

    ! Count non-zeroes for each column of the input matrix
    do i = 1,this%nrow
       do j = this%ia(i),this%ia(i+1)-1
          transpose%ia(this%ja(j)) = transpose%ia(this%ja(j)) + 1
       enddo
    enddo

    ! Reset pointers
    IW(1) = 1
    do i = 2,transpose%nrow+1
       IW(i) = IW(i-1) + transpose%ia(i-1)
    enddo
    call SCOPY(transpose%nrow+1,IW,1,transpose%ia,1)


    ! Transpose coeffiecients
    do i = 1,this%nrow
       do j = this%ia(i),this%ia(i+1)-1
          ind = IW(this%ja(j))
          transpose%ja(ind) = i
          transpose%coeff(ind) = this%coeff(j)
          permutation(j) = ind 
          IW(this%ja(j)) = ind+1
       end do
    end do

    ! copy the working spmat into this
    select type(this)
    type is (spmat)
       this = transpose
    end select
    ! copy if present
    if (present(perm) ) perm=permutation


    ! free memory
    call transpose%kill(lun_err)
    deallocate (IW,permutation,stat=res) 
    if(res .ne. 0) rc = IOerr(lun_err, err_dealloc, 'transpose matrix', &
         'work array IW perm',res)

  end subroutine transpose_spmat


  !>-----------------------------------------------------------------
  !> Subroutine to solve linear system in the form
  !>      M x = b    or     M^T x = b
  !> with M triangular
  !>-------------------------------------------------------
  subroutine solve_triangular(this,lun_err,rhs,sol,&
                                ! optional arguments
       transpose, &
       inverse_diagonal)
    use Globals
    !use Timing
    implicit none
    class(spmat),                intent(in   ) :: this
    integer,                     intent(in   ) :: lun_err
    real(kind=double),           intent(in   ) :: rhs(this%nrow)
    real(kind=double),           intent(inout) :: sol(this%ncol)
    character(len=1),  optional, intent(in   ) :: transpose    
    real(kind=double), optional, intent(in   ) :: inverse_diagonal(this%nrow)
    ! local
    logical :: rc
    integer :: i,j,k,m,nequ
    integer :: idiag,ifirst,ilast,irow
    real(kind=double) :: scr
    !type(Tim) :: trisol1,trisol2



    nequ = this%nrow
    !
    ! checks
    !
    if ( (this%triangular .eq. 'N') ) then
       rc = IOerr(lun_err, err_inp, 'solve_triangular', &
            'non-triangualr matrix passed')
    end if

    if ( this%ncol .ne. this%nrow ) then
       rc = IOerr(lun_err, err_inp, 'solve_triangular', &
            'non-squared matrix passed')
    end if


    !
    ! case without inverse of the diagonal
    !
    if ( .not. present (inverse_diagonal) ) then 
       if ( present(transpose) .and. (transpose .eq. 'T') ) then
          select case (this%triangular)
          case('U')
             !
             ! solve system U^t x = rhs
             !
             !call trisol1%init()
             !call trisol1%set('start')
             sol = zero
             do k = 1, nequ
                idiag = this%ia(k)
                ilast = this%ia(k+1) - 1
                sol(k) = ( rhs(k) - sol(k) ) / this%coeff(this%ia(k))
                do m = idiag+1, ilast
                   sol(this%ja(m)) = sol(this%ja(m)) + &
                        this%coeff(m) * sol(k)
                end do
             end do
             !call  trisol1%set('stop')
             !call  trisol1%info(lun_err,'      no  diag inv U^T')
             !call  trisol1%kill()
          case('L')
             !
             rc = IOerr(lun_err, err_inp, 'solve_triangular', &
                  ' not implemated for L^t x = b')
          end select
       else
          !
          ! solve with given matrix 
          !
          select case (this%triangular)
          case('L')
             !
             ! solve system L x = rhs
             !
             sol=zero
             do k = 1, nequ
                ! select range
                ifirst = this%ia(k)
                idiag  = this%ia(k+1) - 1  
                sol(k) = rhs(k)
                do m = ifirst,idiag-1
                   sol(k) = sol(k) - this%coeff(m) * sol(this%ja(m))
                end do
                sol(k) = sol(k) / this%coeff(idiag)
             end do
          case('U')
             !
             ! solve system U x = rhs
             !
             !call trisol1%init()
             !call trisol1%set('start')
             do k = nequ,1,-1
                idiag  = this%ia(k  )
                ilast  = this%ia(k+1) - 1
                ! sol(k) =  a_{k,k}^{-1} (
                !           rhs(k) - 
                !           \sum_{j} a_{irow,j} sol_{j} )
                scr = zero
                !
                ! EF change the order of the loop ??
                ! use do m = ilast, idiag+1,-1
                !
                do m = ilast, idiag+1,-1
                   scr = scr + this%coeff(m)*sol(this%ja(m))
                end do
                sol(k) = ( rhs(k) - scr ) / this%coeff(idiag)
             end do
             !call  trisol1%set('stop')
             !call  trisol1%info(lun_err,'     no  diag inv U')
             !call  trisol1%kill()
          end select
       end if
    else
       !
       ! use given inverse_diagonal 
       ! 
       if ( present(transpose) .and. (transpose.eq.'T') ) then
          select case (this%triangular)
          case('U')
             !
             ! solve system U^t x = rhs
             !
             sol = zero
             !
             do k = 1, nequ
                ! line range
                idiag = this%ia(k)
                ilast = this%ia(k+1) - 1
                sol(k) = ( rhs(k) - sol(k) ) * inverse_diagonal(k)
                do m = idiag+1, ilast
                   sol(this%ja(m)) = sol(this%ja(m)) + &
                        this%coeff(m) * sol(k) 
                end do
             end do
!!$             do k = 1, nequ
!!$                sol(k) = sol(k) * inverse_diagonal(k)
!!$             end do
          case('L')
             !
             rc = IOerr(lun_err, err_inp, 'solve_triangular', &
                  ' not implemated for L^t x = b')
          end select
       else
          !
          ! solve with given matrix 
          !
          select case (this%triangular)
          case('L')
             !
             ! solve system L x = rhs
             !
             sol=zero

             do irow = 1, nequ
                ! select line range
                ifirst = this%ia(irow)
                idiag  = this%ia(irow+1) - 1  

                ! commpute  
                !   sol(k) =  a_{k,k}^{-1} (
                !               rhs(k) - 
                !               \sum_{j \neq irow } a_{irow,j} sol_{j} )
                sol(irow)=rhs(irow)
                do m = ifirst,idiag-1
                   sol(irow) = sol(irow) - this%coeff(m) * sol(this%ja(m))
                end do
                sol(irow) = sol(irow) * inverse_diagonal(irow)
             end do
          case('U')
             !
             ! solve system U x = rhs
             !
             do irow = nequ,1,-1
                idiag  = this%ia(irow  )
                ilast  = this%ia(irow+1) - 1
                ! sol(k) =  a_{k,k}^{-1} (
                !           rhs(k) - 
                !           \sum_{j} a_{irow,j} sol_{j} )
                sol(irow) = rhs(irow)
                do m = ilast, idiag+1, -1
                   sol(irow) = sol(irow) - this%coeff(m) * sol(this%ja(m))
                end do
                sol(irow) = sol(irow) * inverse_diagonal(irow) 
             end do
          end select
       end if
    end if

  end subroutine solve_triangular


  !>-----------------------------------------------------------------
  !> Subroutine to solve triangular linear system in the form
  !>      M x = b    or     M^T x = b
  !> with M triangular with diagonal equal one
  !>-------------------------------------------------------
  subroutine solve_triangular_unitdiag(this,lun_err,rhs,sol,transpose)
    use Globals
    !use Timing
    implicit none
    class(spmat),                intent(in   ) :: this
    integer,                     intent(in   ) :: lun_err
    real(kind=double),           intent(in   ) :: rhs(this%nrow)
    real(kind=double),           intent(inout) :: sol(this%ncol)
    character(len=1),  optional, intent(in   ) :: transpose    
    ! local
    logical :: rc
    integer :: i,j,k,m,nequ
    integer :: idiag,ifirst,ilast,irow
    real(kind=double) :: scr
    !type(Tim) :: trisol1,trisol2
        

    nequ = this%nrow
    !
    ! checks
    !
    if ( (this%triangular .eq. 'N') ) then
       rc = IOerr(lun_err, err_inp, 'solve_triangular', &
            'non-triangualr matrix passed')
    end if

    if ( this%nrow .ne. this%ncol ) then
       rc = IOerr(lun_err, err_inp, 'solve_triangular', &
            'non-squared matrix passed')
    end if

    if ( .not. this%unitary_diag ) then
       rc = IOerr(lun_err, err_inp, 'solve_triangular', &
            'not unitary diagonal matrix passed')
    end if

    if ( present(transpose) .and. (transpose .eq. 'T') ) then
       select case (this%triangular)
       case('U')
          !call trisol1%init()
          !call trisol1%set('start')
          !
          ! solve system U^t x = rhs
          !
          sol = zero
          do k = 1, nequ
             idiag = this%ia(k)
             ilast = this%ia(k+1) - 1
             ! removed inversion
             !sol(k) = ( rhs(k) - sol(k) ) / this%coeff(idiag)
             sol(k) = rhs(k) - sol(k) 
             do m = idiag+1, ilast
                sol(this%ja(m)) = sol(this%ja(m)) + &
                     this%coeff(m) * sol(k)
             end do
          end do
          !call  trisol1%set('stop')
          !call  trisol1%info(lun_err,'    yes diag inv U^T')
          !call  trisol1%kill()
          return
       case('L')
          !
          rc = IOerr(lun_err, err_inp, 'solve_triangular', &
               ' not implemated for L^t x = b')
       end select
    else
       !
       ! solve with given matrix 
       !
       select case (this%triangular)
       case('L')
          !
          ! solve system L x = rhs
          !
          sol=zero
          do k = 1, nequ
             ! select range
             ifirst = this%ia(k)
             idiag  = this%ia(k+1) - 1  
             sol(k)=rhs(k)
             do m = ifirst,idiag-1
                sol(k) = sol(k) - this%coeff(m) * sol(this%ja(m))
             end do
             ! removed inversion{
             ! sol(k) = sol(k) / this%coeff(idiag)
             ! }
          end do
       case('U')
          !
          ! solve system U x = rhs
          !
          do k = nequ,1,-1
             ! sol(k) =  a_{k,k}^{-1} (
             !           rhs(k) - 
             !           \sum_{j} a_{k,j} sol_{j} )
             !
             idiag  = this%ia(k  )
             ilast  = this%ia(k+1) - 1

             sol(k) = rhs(k)
             do m = ilast, idiag+1, -1
                sol(k) = sol(k) - this%coeff(m)*sol(this%ja(m))
             end do
             ! removed{
             !sol(k) = sol(k) / this%coeff(idiag)
             !} 
          end do
       end select
    end if

  end subroutine solve_triangular_unitdiag

  !>-------------------------------------------------------
  !> Subroutine for incomplete cholesky factorization
  !> of a symmetric matrix A. 
  !> It computes the approximated factor U giving
  !>                  A ~ U^t U
  !> with U a upper triangualr matrix in csr format.
  !> If the optional array "diagonal" is passed the
  !> the output are
  !> 1- an upper triangular matrix with one on the diagonal
  !>            upper    = diag^{-1}(U) U
  !> 2 -the diagonal array containg the diagonal of U 
  !             diagonal = diag(U)
  !> usage:
  !>     call 'var'%incomplete_cholesky(vec_in,vec_out,[info])
  !>
  !> where
  !> \param[in   ] lun_err    -> integer. Info number
  !>                              in case of error 
  !> \param[in   ] nfillin    -> integer. Number of non-zeros
  !>                              per row
  !> \param[in   ] tol_fillin -> real. Drop terms below
  !>                              this tolerance
  !> \param[in   ] job        -> integer. Job to be performed
  !>                              in incoplete factorization
  !> \param[inout] info       -> integer. Info err number
  !> \param[inout] upper      -> type(spmat). Upper factor
  !>                              such that A~=U^t U. If argument
  !>                              diagonal upper=diag^{-1}(U) U
  !>                               
  !> \param[inout] diagonal        -> real(dimension=matrix%nrow)
  !>                              diagonal of U
  !>-------------------------------------------------------
  subroutine incomplete_cholesky (this,&
       lun_err,n_fillin,tol_fillin,job,&
       info,&
       upper,&
       ! optional arguments
       diagonal)
    use Globals
    implicit none
    class(spmat), target,        intent(in   ) :: this
    integer,                     intent(in   ) :: lun_err
    integer,                     intent(in   ) :: n_fillin
    real(kind=double),           intent(in   ) :: tol_fillin
    integer,                     intent(in   ) :: job
    integer,                     intent(inout) :: info
    type(spmat),                 intent(inout) :: upper
    real(kind=double), optional, intent(inout) :: diagonal(this%nrow)
    !local
    logical :: rc
    integer :: res
    integer :: nterm_max, nterm ,nterm_out
    integer :: irow,nrow
    integer, allocatable :: jlu(:),jw(:),iscr(:),iend(:)
    integer, allocatable :: iw1(:),jw1(:),jw2(:),jw3(:)
    integer, target,   allocatable :: iwork(:)
    real(kind=double), allocatable :: alu(:),w(:)
    integer, pointer :: idiag(:)
    type(spmat) :: ssr_temp

    info = 0
    !
    ! checks
    !
    if ( this%nrow .ne. this%ncol ) &
         rc = IOerr(lun_err, err_inp, 'incomplete_cholesky ', &
         'not square matrix passed')

    if ( .not. (this%is_symmetric) ) &
         rc = IOerr(lun_err, err_inp, 'incomplete_cholesky ', &
         'not symmetric  matrix passed')
    
    !
    ! dimensions
    !
    nrow      = this%nrow
    nterm_max = this%nterm + (n_fillin+1) * nrow
    !
    ! work arrays for both storage system
    !
    allocate(&
         alu(nterm_max), &
         jlu(nterm_max),&
         jw1(nrow),&
         jw2(nrow),&
         jw3(nrow+1),&
         w(nrow+1),&
         jw(2*nrow),&
         iscr(nrow+1),&
         stat=res)
    if ( res .ne. 0) &
         rc = IOerr(lun_err, err_alloc, 'incomplete_cholesky ', &
         'work arrays ',res)

    
    if ( n_fillin .gt. 0) then      
       if ( this%storage_system .eq. 'ssr' ) then
          !
          ! 1 : assign integer pointer idiag, ia, ja, iw1
          ! 

          !
          ! 1.1 : allocate extra work arrays 
          !
          allocate(&
               iw1(this%nterm),&     ! ja_transpose
               iwork(nrow+1),&  ! ia_work to cycle ja_transpose
               stat=res)
          if ( res .ne. 0) &
               rc = IOerr(lun_err, err_alloc, 'incomplete_cholesky ', &
               'work arrays iw1, iwork',res)

          !
          ! 1.2 : build iw1, iwork
          ! iw1 shuold act like ja for the lower part of matrix
          ! iwork contains the incides of the startin indices for each row
          !
          
          call my_unpackja(nrow,this%nterm, this%ia, this%ja, iscr(1:nrow),&
               iw1, iwork)

          !
          ! 2 : compute factorization
          !
          
          info = 0
          call my_ssr_incomplete_cheolesky(&
               this%nrow,&
               this%nterm,&
               this%coeff,&
               this%ja,&
               this%ia,&
               iwork,&
               n_fillin,tol_fillin,&
               alu,&
               jlu,&
               nterm_max,&
               w,&
               iw1,& ! already computed 
               jw1,&
               jw2,&
               jw3,&
               info,lun_err)

          !
          ! free local memory
          !
          deallocate(iw1, iwork, stat=res)
          if ( res .ne. 0) &
               rc = IOerr(lun_err, err_dealloc, 'incomplete_cholesky ', &
               'work arrays iw1, iwork',res)

       else if ( this%storage_system .eq. 'csr') then
          !
          ! 1 : assign integer pointer idiag, ia, ja, iw1
          ! 
          
          !
          ! diagonal pointer
          !
          if ( .not. allocated(this%ind_diag) ) then
             !
             ! 1.1 : allocate extra work arrays 
             !
             allocate(&
                  iwork(nrow),&  ! idiag
                  stat=res)
             if ( res .ne. 0) &
                  rc = IOerr(lun_err, err_alloc, 'incomplete_cholesky ', &
                  'work arrays ',res)
             !
             ! build diagonal
             !
             do irow = 1, nrow
                iwork(irow) = this%idiag(irow)
                if ( iwork(irow) .eq. 0) &
                     rc = IOerr(lun_err, err_inp, 'incomplete_cholesky ', &
                     'not zero diagonal term on line', irow)
             end do
             idiag => iwork
          else
             idiag => this%ind_diag
          end if

          !
          ! compute factorization
          !
          info = 0
          call my_csr_incomplete_cheolesky(&
               lun_err, job, &
               this%nrow,&
               this%nterm,&
               this%coeff,&
               this%ja,&
               this%ia,&
               idiag,&
               n_fillin,tol_fillin,&
               alu,&
               jlu,&
               nterm_max,&
               w,&
               jw1,&
               jw2,&
               jw3,&
               info,lun_err)
          
          !
          ! free local memory
          !
          idiag => null()
          if ( allocated(iwork) ) then
             deallocate(iwork,stat=res)
              if ( res .ne. 0) &
                  rc = IOerr(lun_err, err_dealloc, 'incomplete_cholesky ', &
                  'work array iwork',res)
          end if
       end if
                        
       if (info.eq.0) then
          !
          ! convert mssr to csr
          !
          call mss2scaled_upper(lun_err,&
               nrow,&
               nterm_max,&            
               jlu,&
               alu,&
               upper)
                    

          !
          ! aluU(1:nrow) contains the square of the diagonal term
          ! of the upper factor of the choelsky factorization
          !

          if ( present (diagonal) ) then
             !
             ! case for D=diag(U),  \tilde{U}=diag(U)^{-1} U 
             !
             diagonal(1:nrow) = sqrt( alu(1:nrow) )
          else
             !
             ! case for U factorization
             !
             ! use alu as scratch array to compute roots
             alu(1:nrow) = sqrt( alu(1:nrow) )
             call upper%DxM( lun_err,alu(1:nrow) )
             upper%unitary_diag = .false.
          end if


       else if ( info .lt. 0) then
          !
          ! info in case of lu factorization error 
          !
          write(lun_err,*) ' Error in IC construction' 
          select case (info) 
          case (-1)
             rc = IOerr(lun_err, wrn_out, 'build_ic_fillin', &
                  ' The elimination process has generated'//&
                  ' a row in U whose length is .gt.  n)')
          case (-2)  
             rc = IOerr(lun_err, wrn_out, 'build_ic_fillin', &
                  'The matrix U overflows the array U')
          case (-3)  
             rc = IOerr(lun_err, wrn_out, 'build_ic_fillin', &
                  ' Illegal value for n_fillin')
          case(-4)   
             rc = IOerr(lun_err, wrn_out, 'build_ic_fillin', &
                  'zero row encountered')
          case(-5) 
             rc = IOerr(lun_err, wrn_out, 'build_ic_fillin', &
                  'Non Positive Diagonal Element')
          case(-6) 
             rc = IOerr(lun_err, wrn_out, 'build_ic_fillin', &
                  ' Non Positive Pivot Element')
          end select
       end if

    else
       !
       ! No-fillin algorithm 
       !
       if ( this%storage_system .eq. 'ssr' ) then
          !
          ! init upper factor
          !
          if ( upper%is_initialized ) then !!!!!
             call upper%kill(lun_err)
          end if
          call upper%init(lun_err, &
               nrow, nrow, this%nterm,&
               'csr',&
               is_symmetric=.false.,&
               triangular='U')

          !
          ! build U factor 
          !
          upper%ia = this%ia
          upper%ja = this%ja
          call kersh_loc(lun_err,nrow,this%nterm,job,&
               this%ia,&
               this%ja,&
               this%coeff,&
               upper%coeff,&
               info)
       else
          !
          ! assign pointer to diagonal terms of this
          !
          select type( this )
          type is (spmat)
             ssr_temp = this
          end select          
          call ssr_temp%csr2ssr(lun_err)

          !
          ! init upper factor
          !
          if ( upper%is_initialized ) then!!!!!
             call upper%kill(lun_err)
          end if
          call upper%init(lun_err, &
               nrow, nrow, ssr_temp%nterm,&
               'csr',&
               is_symmetric=.false.,&
               triangular='U')
          

          !
          ! build U factor 
          !
          call kersh_loc(lun_err,nrow,ssr_temp%nterm,job,&
               ssr_temp%ia,&
               ssr_temp%ja,&
               ssr_temp%coeff,&
               upper%coeff,&
               info)
          upper%ia = ssr_temp%ia
          upper%ja = ssr_temp%ja


          !
          ! free memory
          ! 
          call ssr_temp%kill(lun_err)

       end if
       if ( info .ne. 0 ) return

       !
       ! U => diag(U) , diag(U)^{-1} U 
       !
       if ( present (diagonal) ) then
          !
          ! case for D=diag(U),  \tilde{U}=diag(U)^{-1} U 
          !
          call upper%get_diagonal(diagonal)
          diagonal = one / diagonal
          call upper%DxM(lun_err, diagonal)
          upper%unitary_diag = .true. 
          diagonal = one / diagonal
       end if
    end if

    !
    ! free memory
    !
    
     deallocate(&
         alu, &
         jlu,&
         jw1,&
         jw2,&
         jw3,&
         w,&
         jw,&
         iscr,&
         stat=res)
     if ( res .ne. 0) &
          rc = IOerr(lun_err, err_dealloc, 'incomplete_cholesky ', &
          'work arrays',res)

 
  contains
    !--------------------------------------------------------------
    ! Subroutine computing incomplete upper triangular factor
    ! of the Cholesky decomposition 
    !            A ~ U^T U
    ! Output matrix is stored in MSS( modified storage stystem) 
    ! containg 
    !         diag(U)**2     [ in array U(1     :nequ) ]
    !         diag(U)^{-1} U [ in array U(nequ+2:iwk ) ]
    !---------------------------------------------------------------
    subroutine my_csr_incomplete_cheolesky(lun_err,job,&
         n,nterm,&
         a,ja,ia,idiag,&
         lfil,droptol,&
         U,ju,iwk,&
         w,&
         jw1,jw2,jw3,ierr,iout)
      !------------------------------------------------------------
      !      
      implicit none 
      !     
      integer lun_err,job
      integer n, nterm 
      integer lfil, iwk, ierr, iout
      integer ja(nterm), ia(n+1),idiag(n),iend(n)
      integer ju(iwk), jw1(n), jw2(n), jw3(n+1)   
      !      
      real*8  a(nterm), U(iwk), w(n+1), droptol
      logical :: rc
      !      
      !------------------------------------------------------------
      !                    *** SYMILUT preconditioner ***         
      !      incomplete U^tU factorization with dual truncation mechanism 
      !-------------------------------------------------------------------
      !     Author: Carlo Janna *December, 12, 2005                       
      !-------------------------------------------------------------------
      ! PARAMETERS                                                        
      !
      ! on entry:
      !========== 
      ! n       = integer. The row dimension of the matrix A.  
      !
      ! nterm   = integer. The number of non-zero of the matrix A. 
      !
      ! a       = ceofficient of matrix 
      ! ja      = column index 
      ! first_row = indecx of the first nonzero element of a for each row
      ! idaig     = 
      !
      ! lfil    = integer. Fill-in parameter. Each row of L and each row
      !           of U will have a maximum of lfil elements (excluding the 
      !           diagonal element). lfil must be .ge. 0.
      !           ** WARNING: THE MEANING OF LFIL HAS CHANGED 
      !           WITH RESPECT T0 EARLIER VERSIONS. 
      !
      ! droptol = real*8. Sets the threshold for dropping small terms in the
      !           factorization. See below for details on dropping strategy.
      !
      !  
      ! iwk     = integer. The lengths of arrays U and ju. If the arrays
      !           are not big enough to store the U^tU factorizations, symilut
      !           will stop with an error message. Last term in ju is used 
      !           to deactivate jw3 work array. The space avaliable for the
      !           factorization is (iwk-1).
      !
      ! On return:
      !===========
      !
      ! U       = matrix stored in Modified Symmetric Sparse Row (MSSR) format
      !           containing the U factor. 
      !           U(1:n) contains the squared of the diagonal of the factor U.
      !           Extra-Diagonal elements are stored in
      !           U(n+2:nterm) and they scaled by the inverse of the diagonal of U
      !           M = D+T with D=diag(U)^2 T=D(U)^{-1} U
      !           U(n+1) is unused.
      !
      ! ju      = integer array. The first n+1 position contain the pointers to
      !           the beginning of each row of U after the diagonal. ju(n+2:*)
      !           are the column indeces of the matrix factor U (excluding the
      !           diagonal entry). [ju(n+1)-1] is the dimension of vectors U and
      !           ju.
      !
      ! ierr    = integer. Error message with the following meaning.
      !           ierr  = 0    --> successful return.
      !           ierr  = -1   --> Error. input matrix may be wrong.
      !                            (The elimination process has generated a
      !                            row in U whose length is .gt.  n).
      !           ierr  = -2   --> The matrix U overflows the array U.
      !           ierr  = -3   --> Illegal value for lfil.
      !           ierr  = -4   --> zero row encountered.
      !           ierr  = -5   --> Non Positive Diagonal Element
      !           ierr  = -6   --> Non Positive Pivot Element
      !
      ! work arrays:
      !=============
      ! iw1     = integer work array of lenght nterm+1
      ! jw1     = integer work array of length n (max number of non-zero).
      ! jw2     = integer work array of length n.
      ! jw3     = integer work array of length n.
      ! w       = real work array of length n+1 (max number of non-zero).
      !  
      !----------------------------------------------------------------------
      ! w, jw1      store the working array [1:ii-1 = L-part, ii:n = u] 
      ! jw2         stores nonzero indicators
      ! jw3         stores pointers to the "interesting" part of row of U 
      !             factor
      ! iw1         stores lines to be eliminated for each line.
      !             first n positions are lenght of each line of iw1
      !                                 "POINTERS???"
      ! 
      ! Notes:
      ! -----
      ! The diagonal elements of the input matrix must be  nonzero (at least
      ! 'structurally'). 
      !
      !----------------------------------------------------------------------* 
      !---- Dual drop strategy works as follows.                             *
      !                                                                      *
      !     1) Theresholding in L and U as set by droptol. Any element whose *
      !        magnitude is less than some tolerance (relative to the abs    *
      !        value of diagonal element in u) is dropped.                   *
      !                                                                      *
      !     2) Keeping only the largest lfil elements in the i-th row of U   * 
      !        (excluding diagonal elements).                                *
      !                                                                      *
      ! Flexibility: one  can use  droptol=0  to get  a strategy  based on   *
      ! keeping  the largest  elements in  each row  of U. Taking            *
      ! droptol .ne.  0 but lfil=n will give  the usual threshold strategy   *
      ! (however, fill-in is then umpredictible).                            *
      !----------------------------------------------------------------------*
      !
      !  Locals
      integer  ju0, rowHead, rowMid , rowEnd, jrow , jcol, jpos
      integer  lenu , lenu0 , lenl , lenght, ncut
      integer  i,  ii, j,  jj, k, kk, j1, j2 , j3 , j4
      character(len=256) :: msg
      !      
      real*8   tnorm, abstol, s, fact, zero 
      !      
      !
      !-----------------------------------------------------------------------
      !  Initialize ju0 (points to next element to be added to U,ju)
      !  and pointer array.
      !-----------------------------------------------------------------------
      ju0 = n+2 
      ju(1) = ju0
      jw3(1) = ju0 
      !
      !  Initialize nonzero indicator array. 
      !
      do j=1,n 
         jw2(j)  = 0
      enddo
      !
      !  Set ju(iwk) to deactivate jw3 pointers
      !
      ju(iwk) = -1      
      !--------------------------------------------------------------
      !  Beginning of main loop.
      !-------------------------------------------------------------
      jw1=0

      do ii = 1, n
         j1 = idiag(ii)
         j2 = ia(ii+1) - 1
         !j2 = iend(ii)
         j3 = ia(ii)
         j4 = idiag(ii)-1
         !  Calculate the norm of the ii-th row         
         tnorm = 0.0d0
         do  k=j1,j2
            tnorm = tnorm+abs(a(k))
         end do
         !         
         if (tnorm .eq. 0.0) goto 999
         tnorm = tnorm/real(j2-j1+1)
         !     
         !  Unpack Upper and Lower part of row of A in array w 
         !   
         lenu  = j2 - j1 + 1
         lenu0 = lenu
         lenl  = j4 - j3 + 1         
         k = ii-1          
         do j = j1,j2         
            w(j-j1+ii)   = a(j)
            jw1(j-j1+ii) = ja(j)
            k = k + 1         
            jw2(ja(j))  = k
         end do          

         k = 0         
         do j = j3,j4   
            ! in case of csr iw1=ja
            jw1(j-j3+1)  = ja(j)
            k = k + 1         
            jw2(ja(j))  = k
         end do
!!$         write(*,'(10(1x,e8.1))') (w(k),k=1,n)  
!!$         write(*,'(10I6)') (jw1(k),k=1,n) 
!!$         write(*,'(10I6)') (jw2(k),k=1,n) 
!!$         write(*,*) ' ' 
!
         !----BLAS-Method-----------------------------------------------------
         !  Unpack Upper and Lower part of row of A in array w 
         !         lenu  = j2 - j1 + 1
         !         lenu0 = lenu
         !         lenl  = j4 - j3 + 1
         !         call DCOPY(lenu,a(j1),1,w(ii),1)
         !         call SCOPY(lenu,ja(j1),1,jw1(ii),1)
         !         call SCOPY(lenl,iw1(j3),1,jw1(1),1)
         !   
         !  Set non-zero indicators
         !         
         !         k = ii-1
         !         do j = j1,j2
         !            k = k + 1
         !            jw2(ja(j)) = k
         !         enddo
         !         k = 0
         !         do j = j3,j4
         !            k = k + 1
         !            jw2(iw1(j)) = k
         !         enddo
         !----End of BLAS-Method----------------------------------------------
         !
         !  Eliminate previous rows
         !
         jj = 0 
150      jj = jj+1
         if (jj .gt. lenl) goto 160
         !
         !  Determine the smallest column index
         !         
         jrow = jw1(jj)
         k = jj
         do j = jj+1,lenl
            if (jw1(j) .lt. jrow) then
               jrow = jw1(j)
               k = j
            endif
         enddo

         if (k .ne. jj) then
            ! Exchange in jw1
            j = jw1(jj)
            jw1(jj) = jw1(k)
            jw1(k) = j
            ! Exchange in jw1
            jw2(jrow) = jj
            jw2(j) = k
         endif
         !write(*,*) 'ii, jrow',ii, jrow

         !
         !  Zero out element in row by setting jw2(jrow) to zero
         !
         jw2(jrow) = 0
         !         
         !  Get the leading term of the row if it exists a(jj,ii)
         !  
         rowHead = ju(jrow)
         rowMid  = jw3(jrow)
         if (ju(rowMid) .eq. ii) then  
            !------------------------------------------------------------------
            !  Mettere la possibilità di saltare 
            !  questo ciclo se il termine pivotale
            !  è troppo piccolo: fact = U(rowMid)
            !-----------------------------------------------------------------
            !
            ! Exists, mark fill-in terms of the current row before the
            ! diagonal            
            !
            do k = rowHead,rowMid-1
               jcol = ju(k)
               jpos = jw2(jcol)
               if (jpos .eq. 0) then
                  ! This is a fill-in term
                  lenl = lenl+1
                  if (lenl .gt. n) goto 996
                  jw1(lenl) = jcol
                  jw2(jcol) = lenl
               endif
            enddo
            !                    
            ! Combine current row ii with row jrow
            !

            fact = U(rowMid)*U(jrow)
            ! fact = U(rowMid) * U(jrow) * sqrt(U(jrow))
            rowEnd = ju(jrow+1) - 1
            do k = rowMid,rowEnd
               !
               !write(*,*) ii, rowMid,ju(rowMid), jrow,  ju(k), k
               s = fact * U(k)
               jcol = ju(k)
               jpos = jw2(jcol)
               if (jpos .eq. 0) then
                  !                          
                  ! This is a fill-in element
                  !
                  lenu = lenu + 1
                  if (lenu .gt. n) goto 996
                  i = ii+lenu-1
                  jw1(i) = jcol
                  jw2(jcol) = i
                  w(i) = -s
                  !                     
               else
                  !
                  ! This is not a fill-in element
                  !
                  w(jpos) = w(jpos) - s
                  !                     
               endif
               !                  
            enddo
            !               
            ! Update pointer jw3
            !
            jw3(jrow) = jw3(jrow) + 1
            if (jw3(jrow) .gt. rowEnd) then
               ! Deactivate pointer
               jw3(jrow) = iwk
            endif
            !               
         endif

         !           
         goto 150
160      continue
         !
         !  Reset non-zero indicators
         !
         do k = 1,lenu
            jw2(jw1(ii+k-1)) = 0
         enddo
         !
         !  Store the diagonal term
         !
         s = w(ii)
         if (s .le. zero) then
            !write(iout,100) ii,s
            select case (job)
            case (0)
               ! breaks factorization   
               write(msg,*) ' non-positive pivot ', s , ' found at line',ii 
               rc =IOerr(lun_err, wrn_val, 'incomplete_cholesky', &
                 etb(msg) )
               ierr=ii
               return
            case (1)
               !
               ! s is replacede as in ilut by Saad
               !
               ! breaks factorization  
               !write(msg,*) ' non-positive pivot ', s , &
               !     ' found at line',ii ,&
               !     ' repalced with=  ', (0.0001 + droptol)*tnorm

               !rc = IOerr(lun_err, wrn_val, 'incomplete_cholesky', &
               !  etb(msg))
               ierr=-6              
               s =(0.0001 + droptol)*tnorm
               w(ii)= s
               ierr = 0
            case (2)
               !
               ! s is replacede as in ilut by Saad
               !  
               if ( abs(s) < large ) then
                  s =(0.0001 + droptol)*tnorm
                  w(ii)= s
                  ierr = 0
               end if
            end select
         endif
         U(ii) = s
         s = 1.d0/s
         !
         !  Apply dropping strategy
         !
         abstol = sqrt(s) * tnorm * droptol  ! Set value of true tolerance
         lenght = 0 
         do k = ii+1,ii+lenu-1
            if (abs(w(k)) .gt. abstol) then
               lenght = lenght+1
               w(ii+lenght)   = w(k)
               jw1(ii+lenght) = jw1(k)
            endif
         enddo
         !
         !  Set lenu = number of elements bigger than tolerance        
         !
         lenu = lenght
         ncut = min0(lenu,lenu0+lfil-1)
         !----------------------------------------------------------------
         !  May be interesting to save a fixed number of non-zero 
         !  for each line
         ! ncut = min0(lenu,lfil)         
         !----------------------------------------------------------------
         !
         !  Select the ncut biggest elements only if lenu gt zero
         !
         if (lenu .gt. 0) then
            jpos = ii+1
            call qsplit(w(jpos),jw1(jpos),lenu,ncut)
            !          
            !  Order them in increasing jw1
            !
            !call HPSORT(jw1(jpos),w(jpos),ncut)
            call sort_row(ncut,&
                 jw1(jpos:jpos+ncut-1),&
                 w(jpos:jpos+ncut-1))
            !
            !  Scale and store in the factor
            !         
            if ((ju0+ncut-1) .gt. iwk-1) goto 997
            do k = 0,ncut-1
               ! EF
               U(ju0+k)  = s*w(jpos+k)
               ju(ju0+k) = jw1(jpos+k) 
            enddo
            !            
            !----BLAS-Method-------------------------------------------------
            !  Scale by the pivotal term and store
            !            if ((ju0+ncut-1) .gt. iwk-1) goto 997
            !            call DCOPY(ncut,w(jpos),1,U(ju0),1)
            !            call DSCAL(ncut,s,w(jpos),1)
            !            call SCOPY(ncut,jw1(jpos),1,ju(ju0),1)
            !----End of BLAS-Method------------------------------------------
            !            
         endif
         !
         !  Update pointer to beginning of next row of U    
         !
         ju0 = ju0 + ncut
         ju(ii+1) = ju0 
         jw3(ii+1) = ju0 
         !
      enddo


      !-----------------------------------------------------------------------
      !     end main loop
      !-----------------------------------------------------------------------
      !
      !     Successful run
      !
      ierr = 0
      !     
      return
      !
      !     incomprehensible error. Matrix must be wrong.
      !     
996   ierr = -1
      return
      !    
      !     insufficient storage in U.
      !     
997   ierr = -2
      return
      !     
      !     illegal lfil entered.
      !     
998   ierr = -3
      return
      !     
      !     zero row encountered
      !     
999   ierr = -4
      return
      !----------------end-of-ilut--------------------------------------------
100   format(' NULL OR NEGATIVE DIAGONAL ELEMENT: I,J =',I8,2X,E16.5)
101   format(' NORM OF ROW =  ',E16.8)
      !-----------------------------------------------------------------------
    end subroutine my_csr_incomplete_cheolesky


    !--------------------------------------------------------------
    ! Subroutine computing incomplete upper triangular factor
    ! of the Cholesky decomposition 
    !            A ~ U^T U
    ! Output matrix is stored in MSS( modified storage stystem) 
    ! containg 
    !         diag(U)**2     [ in array U(1     :nequ) ]
    !         diag(U)^{-1} U [ in array U(nequ+2:iwk ) ]
    !---------------------------------------------------------------
!!$    subroutine general_csr_incomplete_cheolesky(iout,ierr,&
!!$         nrow , irow ,nnz_row,&
!!$         ja_row,coeff_row,&
!!$         lfil,droptol,&
!!$         iwk,U,ju,&
!!$         ju0,jw1,jw2,jw3)
!!$      !------------------------------------------------------------
!!$      !      
!!$      implicit none 
!!$      ! comunication varaibles
!!$      integer, intent(in)              :: iout
!!$      integer, intent(inout)           :: ierr
!!$      ! input rows
!!$      integer, intent(in)              :: nrow
!!$      integer, intent(in)              :: irow
!!$      integer, intent(in)              :: nnz_row
!!$      ! inputs used as work arrays
!!$      integer, intent(inout)           :: ja_row(nrow)
!!$      real(kind=double), intent(inout) :: coeff_row(nrow)
!!$      ! factorization controls
!!$      integer, intent(in)              :: lfil
!!$      real(kind=double), intent(in   ) :: droptol
!!$      ! Factors
!!$      integer, intent(inout)           :: iwk
!!$      integer, intent(inout)           :: ju(iwk)
!!$      real(kind=double), intent(inout) :: U(iwk)
!!$      ! work arrays
!!$      integer, intent(inout)           :: ju0
!!$      integer, intent(inout)           :: jw1(nrow), jw2(nrow), jw3(nrow+1)
!!$      !      
!!$      !      
!!$      !------------------------------------------------------------
!!$      !                    *** SYMILUT preconditioner ***         
!!$      !      incomplete U^tU factorization with dual truncation mechanism 
!!$      !-------------------------------------------------------------------
!!$      !     Author: Carlo Janna *December, 12, 2005                       
!!$      !-------------------------------------------------------------------
!!$      ! PARAMETERS                                                        
!!$      !
!!$      ! on entry:
!!$      !========== 
!!$      ! n       = integer. The row dimension of the matrix A.  
!!$      !
!!$      ! nterm   = integer. The number of non-zero of the matrix A. 
!!$      !
!!$      ! a       = ceofficient of matrix 
!!$      ! ja      = column index 
!!$      ! first_row = indecx of the first nonzero element of a for each row
!!$      ! idaig     = 
!!$      !
!!$      ! lfil    = integer. Fill-in parameter. Each row of L and each row
!!$      !           of U will have a maximum of lfil elements (excluding the 
!!$      !           diagonal element). lfil must be .ge. 0.
!!$      !           ** WARNING: THE MEANING OF LFIL HAS CHANGED 
!!$      !           WITH RESPECT T0 EARLIER VERSIONS. 
!!$      !
!!$      ! droptol = real*8. Sets the threshold for dropping small terms in the
!!$      !           factorization. See below for details on dropping strategy.
!!$      !
!!$      !  
!!$      ! iwk     = integer. The lengths of arrays U and ju. If the arrays
!!$      !           are not big enough to store the U^tU factorizations, symilut
!!$      !           will stop with an error message. Last term in ju is used 
!!$      !           to deactivate jw3 work array. The space avaliable for the
!!$      !           factorization is (iwk-1).
!!$      !
!!$      ! On return:
!!$      !===========
!!$      !
!!$      ! U       = matrix stored in Modified Symmetric Sparse Row (MSSR) format
!!$      !           containing the U factor. 
!!$      !           U(1:n) contains the squared of the diagonal of the factor U.
!!$      !           Extra-Diagonal elements are stored in
!!$      !           U(n+2:nterm) and they scaled by the inverse of the diagonal of U
!!$      !           M = D+T with D=diag(U)^2 T=D(U)^{-1} U
!!$      !           U(n+1) is unused.
!!$      !
!!$      ! ju      = integer array. The first n+1 position contain the pointers to
!!$      !           the beginning of each row of U after the diagonal. ju(n+2:*)
!!$      !           are the column indeces of the matrix factor U (excluding the
!!$      !           diagonal entry). [ju(n+1)-1] is the dimension of vectors U and
!!$      !           ju.
!!$      !
!!$      ! ierr    = integer. Error message with the following meaning.
!!$      !           ierr  = 0    --> successful return.
!!$      !           ierr  = -1   --> Error. input matrix may be wrong.
!!$      !                            (The elimination process has generated a
!!$      !                            row in U whose length is .gt.  n).
!!$      !           ierr  = -2   --> The matrix U overflows the array U.
!!$      !           ierr  = -3   --> Illegal value for lfil.
!!$      !           ierr  = -4   --> zero row encountered.
!!$      !           ierr  = -5   --> Non Positive Diagonal Element
!!$      !
!!$      ! work arrays:
!!$      !=============
!!$      ! iw1     = integer work array of lenght nterm+1
!!$      ! jw1     = integer work array of length n (max number of non-zero).
!!$      ! jw2     = integer work array of length n.
!!$      ! jw3     = integer work array of length n.
!!$      ! w       = real work array of length n+1 (max number of non-zero).
!!$      !  
!!$      !----------------------------------------------------------------------
!!$      ! w, jw1      store the working array [1:ii-1 = L-part, ii:n = u] 
!!$      ! jw2         stores nonzero indicators
!!$      ! jw3         stores pointers to the "interesting" part of row of U 
!!$      !             factor
!!$      ! iw1         stores lines to be eliminated for each line.
!!$      !             first n positions are lenght of each line of iw1
!!$      !                                 "POINTERS???"
!!$      ! 
!!$      ! Notes:
!!$      ! -----
!!$      ! The diagonal elements of the input matrix must be  nonzero (at least
!!$      ! 'structurally'). 
!!$      !
!!$      !----------------------------------------------------------------------* 
!!$      !---- Dual drop strategy works as follows.                             *
!!$      !                                                                      *
!!$      !     1) Theresholding in L and U as set by droptol. Any element whose *
!!$      !        magnitude is less than some tolerance (relative to the abs    *
!!$      !        value of diagonal element in u) is dropped.                   *
!!$      !                                                                      *
!!$      !     2) Keeping only the largest lfil elements in the i-th row of U   * 
!!$      !        (excluding diagonal elements).                                *
!!$      !                                                                      *
!!$      ! Flexibility: one  can use  droptol=0  to get  a strategy  based on   *
!!$      ! keeping  the largest  elements in  each row  of U. Taking            *
!!$      ! droptol .ne.  0 but lfil=n will give  the usual threshold strategy   *
!!$      ! (however, fill-in is then umpredictible).                            *
!!$      !----------------------------------------------------------------------*
!!$      !
!!$      !  Locals
!!$      integer  rowHead, rowMid , rowEnd, jrow , jcol, jpos
!!$      integer  lenu , lenu0 , lenl , lenght, ncut
!!$      integer  i,  ii, j,  jj, k, kk, j1, j2 , j3 , j4, n, idiag
!!$      !      
!!$      real*8   tnorm, abstol, s, fact
!!$      !      
!!$      !
!!$      
!!$      ! TODO move outside
!!$      ! !-----------------------------------------------------------------------
!!$      ! !  Initialize ju0 (points to next element to be added to U,ju)
!!$      ! !  and pointer array.
!!$      ! !-----------------------------------------------------------------------
!!$      ! ju0 = n+2 
!!$      ! ju(1) = ju0
!!$      ! jw3(1) = ju0 
!!$      ! !
!!$      ! !  Initialize nonzero indicator array. 
!!$      ! !
!!$      ! do j=1,n 
!!$      !    jw2(j)  = 0
!!$      ! enddo
!!$      ! !
!!$      ! !  Set ju(iwk) to deactivate jw3 pointers
!!$      ! !
!!$      ! ju(iwk) = -1      
!!$      ! !--------------------------------------------------------------
!!$      ! !  Beginning of main loop.
!!$      ! !-------------------------------------------------------------
!!$      ! jw1=0
!!$
!!$
!!$      ! lenu  = nnzterm  in row after  diagonal (included)
!!$      ! lenl  = nnzterm  in row before diagonal 
!!$      ! 
!!$      ! work = ( 0...0.. a_ii (other nnz term packed after diagonal...0 ) size = nequ 
!!$      ! jw1  = ( column index packed around diagonal                     ) size = nequ
!!$      ! jw2  = ( ___1____2____3___4___5__) 
!!$      !             |    |    |   |   |
!!$      !             (nonzeros elements)
!!$      
!!$      !
!!$      ! assembly work array jw2
!!$      !
!!$      k=0
!!$      do i=1,nnz_row
!!$         k=k+1
!!$         jw2(ja_row(i))=k
!!$         if ( ja_row(i) .eq. irow ) then 
!!$            idiag = k
!!$            lenl  = k - 1
!!$         end if
!!$      end do
!!$     
!!$      lenu = nnz_row - lenl
!!$      
!!$      !  Calculate the norm of the ii-th row         
!!$      tnorm = zero
!!$      do  k = idiag, idiag + lenu - 1 
!!$         tnorm = tnorm+abs(coeff_row(k))
!!$      end do              
!!$      if (tnorm .eq. 0.0) goto 999
!!$      tnorm = tnorm/real(j2-j1+1)
!!$
!!$      ii = irow
!!$      n  = nrow 
!!$      ! do ii = 1, n
!!$      !    j1 = idiag(ii)
!!$      !    j2 = ia(ii+1) - 1
!!$      !    !j2 = iend(ii)
!!$      !    j3 = ia(ii)
!!$      !    j4 = idiag(ii)-1
!!$         
!!$      !    !  Calculate the norm of the ii-th row         
!!$      !    tnorm = 0.0d0
!!$      !    do  k=j1,j2
!!$      !       tnorm = tnorm+abs(a(k))
!!$      !    end do
!!$      !    !         
!!$      !    if (tnorm .eq. 0.0) goto 999
!!$      !    tnorm = tnorm/real(j2-j1+1)
!!$      !    !     
!!$      !    !  Unpack Upper and Lower part of row of A in array w 
!!$      !    !   
!!$      !    lenu  = j2 - j1 + 1
!!$      !    lenu0 = lenu
!!$      !    lenl  = j4 - j3 + 1         
!!$      !    k = ii-1          
!!$      !    do j = j1,j2         
!!$      !       w(j-j1+ii)   = a(j)
!!$      !       jw1(j-j1+ii) = ja(j)
!!$      !       k = k + 1         
!!$      !       jw2(ja(j))  = k
!!$      !    end do
!!$          
!!$
!!$      !    k = 0         
!!$      !    do j = j3,j4   
!!$      !       ! in case of csr iw1=ja
!!$      !       jw1(j-j3+1)  = ja(j)
!!$      !       k = k + 1         
!!$      !       jw2(ja(j))  = k
!!$      !    end do
!!$      !
!!$      ! Core algorithm for incomplete cholesky facorization
!!$      !
!!$
!!$
!!$      !
!!$      !  Eliminate previous rows
!!$      !
!!$      jj = 0 
!!$150   jj = jj+1
!!$      if (jj .gt. lenl) goto 160
!!$      !
!!$      !  Determine the smallest column index
!!$      !         
!!$      jrow = ja_row(jj)
!!$      k = jj
!!$      do j = jj+1,lenl
!!$         if (ja_row(j) .lt. jrow) then
!!$            jrow = ja_row(j)
!!$            k = j
!!$         endif
!!$      enddo
!!$      
!!$      if (k .ne. jj) then
!!$         ! Exchange in ja_row
!!$         j = ja_row(jj)
!!$         ja_row(jj) = ja_row(k)
!!$         ja_row(k) = j
!!$         ! Exchange in ja_row
!!$         jw2(jrow) = jj
!!$         jw2(j) = k
!!$      endif
!!$
!!$      !
!!$      !  Zero out element in row by setting jw2(jrow) to zero
!!$      !
!!$      jw2(jrow) = 0
!!$      !         
!!$      !  Get the leading term of the row if it exists a(jj,ii)
!!$      !  
!!$      rowHead = ju(jrow)
!!$      rowMid  = jw3(jrow)
!!$      if (ju(rowMid) .eq. ii) then  
!!$         !------------------------------------------------------------------
!!$         !  Mettere la possibilità di saltare 
!!$         !  questo ciclo se il termine pivotale
!!$         !  è troppo piccolo: fact = U(rowMid)
!!$         !-----------------------------------------------------------------
!!$         !
!!$         ! Exists, mark fill-in terms of the current row before the
!!$         ! diagonal            
!!$         !
!!$         do k = rowHead,rowMid-1
!!$            jcol = ju(k)
!!$            jpos = jw2(jcol)
!!$            if (jpos .eq. 0) then
!!$               ! This is a fill-in term
!!$               lenl = lenl+1
!!$               if (lenl .gt. n) goto 996
!!$               ja_row(lenl) = jcol
!!$               jw2(jcol) = lenl
!!$            endif
!!$         enddo
!!$         !                    
!!$         ! Combine current row ii with row jrow
!!$         !
!!$
!!$         fact = U(rowMid)*U(jrow)
!!$         ! fact = U(rowMid) * U(jrow) * sqrt(U(jrow))
!!$         rowEnd = ju(jrow+1) - 1
!!$         do k = rowMid,rowEnd
!!$            !
!!$            !write(*,*) ii, rowMid,ju(rowMid), jrow,  ju(k), k
!!$            s = fact * U(k)
!!$            jcol = ju(k)
!!$            jpos = jw2(jcol)
!!$            if (jpos .eq. 0) then
!!$               !                          
!!$               ! This is a fill-in element
!!$               !
!!$               lenu = lenu + 1
!!$               if (lenu .gt. n) goto 996
!!$               i = ii+lenu-1
!!$               ja_row(i) = jcol
!!$               jw2(jcol) = i
!!$               coeff_row(i) = -s
!!$               !                     
!!$            else
!!$               !
!!$               ! This is not a fill-in element
!!$               !
!!$               coeff_row(jpos) = coeff_row(jpos) - s
!!$               !                     
!!$            endif
!!$            !                  
!!$         enddo
!!$         !               
!!$         ! Update pointer jw3
!!$         !
!!$         jw3(jrow) = jw3(jrow) + 1
!!$         if (jw3(jrow) .gt. rowEnd) then
!!$            ! Deactivate pointer
!!$            jw3(jrow) = iwk
!!$         endif
!!$         !               
!!$      endif
!!$      !           
!!$      goto 150
!!$160   continue
!!$      !
!!$      !  Reset non-zero indicators
!!$      !
!!$      do k = 1,lenu
!!$         jw2(ja_row(ii+k-1)) = 0
!!$      enddo
!!$      !
!!$      !  Store the diagonal term
!!$      !
!!$      s = coeff_row(ii)         
!!$      if (s .le. zero) then
!!$         write(iout,100) ii,s
!!$         write(6,100) ii,s
!!$         write(6,*) 'ilu error'
!!$         !stop
!!$         ierr = -5  
!!$         write(iout,101) tnorm
!!$         return
!!$      endif
!!$      U(ii) = s
!!$      s = 1.d0/s
!!$      !
!!$      !  Apply dropping strategy
!!$      !
!!$      abstol = s * tnorm * droptol  ! Set value of true tolerance
!!$      lenght = 0 
!!$      do k = ii+1,ii+lenu-1
!!$         if (abs(coeff_row(k)) .gt. abstol) then
!!$            lenght = lenght+1
!!$            coeff_row(ii+lenght)   = coeff_row(k)
!!$            ja_row(ii+lenght) = ja_row(k)
!!$         endif
!!$      enddo
!!$      !
!!$      !  Set lenu = number of elements bigger than tolerance        
!!$      !
!!$      lenu = lenght
!!$      ncut = min0(lenu,lenu0+lfil-1)
!!$      !----------------------------------------------------------------
!!$      !  May be interesting to save a fixed number of non-zero 
!!$      !  for each line
!!$      ! ncut = min0(lenu,lfil)         
!!$      !----------------------------------------------------------------
!!$      !
!!$      !  Select the ncut biggest elements only if lenu gt zero
!!$      !
!!$      if (lenu .gt. 0) then
!!$         jpos = ii+1
!!$         call qsplit(coeff_row(jpos),ja_row(jpos),lenu,ncut)
!!$         !          
!!$         !  Order them in increasing ja_row
!!$         !
!!$         !call HPSORT(ja_row(jpos),w(jpos),ncut)
!!$         call sort_row(ncut,&
!!$              ja_row(jpos:jpos+ncut-1),&
!!$              coeff_row(jpos:jpos+ncut-1))
!!$         !
!!$         !  Scale and store in the factor
!!$         !         
!!$         if ((ju0+ncut-1) .gt. iwk-1) goto 997
!!$         do k = 0,ncut-1
!!$            ! EF
!!$            U(ju0+k)  = s*coeff_row(jpos+k)
!!$            ju(ju0+k) = ja_row(jpos+k) 
!!$         enddo
!!$         !            
!!$         !----BLAS-Method-------------------------------------------------
!!$         !  Scale by the pivotal term and store
!!$         !            if ((ju0+ncut-1) .gt. iwk-1) goto 997
!!$         !            call DCOPY(ncut,w(jpos),1,U(ju0),1)
!!$         !            call DSCAL(ncut,s,w(jpos),1)
!!$         !            call SCOPY(ncut,ja_row(jpos),1,ju(ju0),1)
!!$         !----End of BLAS-Method------------------------------------------
!!$         !            
!!$      endif
!!$      !
!!$      !  Update pointer to beginning of next row of U    
!!$      !
!!$      ju0 = ju0 + ncut
!!$      ju(ii+1) = ju0 
!!$      jw3(ii+1) = ju0 
!!$      
!!$      !
!!$      !
!!$      !     Successful run
!!$      !
!!$      ierr = 0
!!$      !     
!!$      return
!!$      !
!!$      !     incomprehensible error. Matrix must be wrong.
!!$      !     
!!$996   ierr = -1
!!$      return
!!$      !    
!!$      !     insufficient storage in U.
!!$      !     
!!$997   ierr = -2
!!$      return
!!$      !     
!!$      !     illegal lfil entered.
!!$      !     
!!$998   ierr = -3
!!$      return
!!$      !     
!!$      !     zero row encountered
!!$      !     
!!$999   ierr = -4
!!$      return
!!$      !----------------end-of-ilut--------------------------------------------
!!$100   format(' NULL OR NEGATIVE DIAGONAL ELEMENT: I,J =',I8,2X,E16.5)
!!$101   format(' NORM OF ROW =  ',E16.8)
!!$      !-----------------------------------------------------------------------
!!$    end subroutine general_csr_incomplete_cheolesky

    

    !--------------------------------------------------------------
    ! Subroutine computing incomplete upper triangular factor
    ! of the Cholesky decomposition 
    !            A ~ U^T U
    ! Output matrix is stored in MSS( modified storage stystem) 
    ! containg 
    !         diag(U)**2     [ in array U(1     :nequ) ]
    !         diag(U)^{-1} U [ in array U(nequ+2:iwk ) ]
    !---------------------------------------------------------------
    subroutine my_ssr_incomplete_cheolesky(n,nterm,&
         a,ja,ia,ia_work,&
         lfil,droptol,&
         U,ju,iwk,&
         w,&
         iw1,jw1,jw2,jw3,ierr,iout)
      !------------------------------------------------------------
      !      
      implicit none 
      !      
      integer n, nterm 
      integer lfil, iwk, ierr, iout
      integer ja(nterm), ia(n+1), ia_work(n+1)
      integer ju(iwk), iw1(nterm), jw1(n), jw2(n), jw3(n+1)   
      !      
      real*8  a(nterm), U(iwk), w(n+1), droptol
      !      
      !------------------------------------------------------------
      !                    *** SYMILUT preconditioner ***         
      !      incomplete U^tU factorization with dual truncation mechanism 
      !-------------------------------------------------------------------
      !     Author: Carlo Janna *December, 12, 2005                       
      !-------------------------------------------------------------------
      ! PARAMETERS                                                        
      !
      ! on entry:
      !========== 
      ! n       = integer. The row dimension of the matrix A.  
      !
      ! nterm   = integer. The number of non-zero of the matrix A. 
      !
      ! a       = ceofficient of matrix 
      ! ja      = column index 
      ! first_row = indecx of the first nonzero element of a for each row
      ! idaig     = 
      !
      ! lfil    = integer. Fill-in parameter. Each row of L and each row
      !           of U will have a maximum of lfil elements (excluding the 
      !           diagonal element). lfil must be .ge. 0.
      !           ** WARNING: THE MEANING OF LFIL HAS CHANGED 
      !           WITH RESPECT T0 EARLIER VERSIONS. 
      !
      ! droptol = real*8. Sets the threshold for dropping small terms in the
      !           factorization. See below for details on dropping strategy.
      !
      !  
      ! iwk     = integer. The lengths of arrays U and ju. If the arrays
      !           are not big enough to store the U^tU factorizations, symilut
      !           will stop with an error message. Last term in ju is used 
      !           to deactivate jw3 work array. The space avaliable for the
      !           factorization is (iwk-1).
      !
      ! On return:
      !===========
      !
      ! U       = matrix stored in Modified Symmetric Sparse Row (MSSR) format
      !           containing the U factor. 
      !           U(1:n) contains the squared of the diagonal of the factor U.
      !           Extra-Diagonal elements are stored in
      !           U(n+2:nterm) and they scaled by the inverse of the diagonal of U
      !           M = D+T with D=diag(U)^2 T=D(U)^{-1} U
      !           U(n+1) is unused.
      !
      ! ju      = integer array. The first n+1 position contain the pointers to
      !           the beginning of each row of U after the diagonal. ju(n+2:*)
      !           are the column indeces of the matrix factor U (excluding the
      !           diagonal entry). [ju(n+1)-1] is the dimension of vectors U and
      !           ju.
      !
      ! ierr    = integer. Error message with the following meaning.
      !           ierr  = 0    --> successful return.
      !           ierr  = -1   --> Error. input matrix may be wrong.
      !                            (The elimination process has generated a
      !                            row in U whose length is .gt.  n).
      !           ierr  = -2   --> The matrix U overflows the array U.
      !           ierr  = -3   --> Illegal value for lfil.
      !           ierr  = -4   --> zero row encountered.
      !           ierr  = -5   --> Non Positive Diagonal Element
      !
      ! work arrays:
      !=============
      ! iw1     = integer work array of lenght nterm+1
      ! jw1     = integer work array of length n (max number of non-zero).
      ! jw2     = integer work array of length n.
      ! jw3     = integer work array of length n.
      ! w       = real work array of length n+1 (max number of non-zero).
      !  
      !----------------------------------------------------------------------
      ! w, jw1      store the working array [1:ii-1 = L-part, ii:n = u] 
      ! jw2         stores nonzero indicators
      ! jw3         stores pointers to the "interesting" part of row of U 
      !             factor
      ! iw1         stores lines to be eliminated for each line.
      !             first n positions are lenght of each line of iw1
      !                                 "POINTERS???"
      ! 
      ! Notes:
      ! -----
      ! The diagonal elements of the input matrix must be  nonzero (at least
      ! 'structurally'). 
      !
      !----------------------------------------------------------------------* 
      !---- Dual drop strategy works as follows.                             *
      !                                                                      *
      !     1) Theresholding in L and U as set by droptol. Any element whose *
      !        magnitude is less than some tolerance (relative to the abs    *
      !        value of diagonal element in u) is dropped.                   *
      !                                                                      *
      !     2) Keeping only the largest lfil elements in the i-th row of U   * 
      !        (excluding diagonal elements).                                *
      !                                                                      *
      ! Flexibility: one  can use  droptol=0  to get  a strategy  based on   *
      ! keeping  the largest  elements in  each row  of U. Taking            *
      ! droptol .ne.  0 but lfil=n will give  the usual threshold strategy   *
      ! (however, fill-in is then umpredictible).                            *
      !----------------------------------------------------------------------*
      !
      !  Locals
      integer  ju0, rowHead, rowMid , rowEnd, jrow , jcol, jpos
      integer  lenu , lenu0 , lenl , lenght, ncut
      integer  i,  ii, j,  jj, k, kk, j1, j2 , j3 , j4
      !      
      real*8   tnorm, abstol, s, fact, zero 
      !      
      parameter(zero = 0.0)
      !
      !-----------------------------------------------------------------------
      !  Initialize ju0 (points to next element to be added to U,ju)
      !  and pointer array.
      !-----------------------------------------------------------------------
      ju0 = n+2 
      ju(1) = ju0
      jw3(1) = ju0 
      !
      !  Initialize nonzero indicator array. 
      !
      do j=1,n 
         jw2(j)  = 0
      enddo
      !
      !  Set ju(iwk) to deactivate jw3 pointers
      !
      ju(iwk) = -1      
      !--------------------------------------------------------------
      !  Beginning of main loop.
      !-------------------------------------------------------------
      jw1=0

      do ii = 1, n
         j1 = ia(ii)
         j2 = ia(ii+1) - 1
         j3 = ia_work(ii)
         j4 = ia_work(ii+1)-1

         !  Calculate the norm of the ii-th row         
         tnorm = 0.0d0
         do  k=j1,j2
            tnorm = tnorm+abs(a(k))
         end do
         !         
         if (tnorm .eq. 0.0) goto 999
         tnorm = tnorm/real(j2-j1+1)
         !     
         !  Unpack Upper and Lower part of row of A in array w 
         !   
         lenu  = j2 - j1 + 1
         lenu0 = lenu
         lenl  = j4 - j3 + 1         
         k = ii-1          
         do j = j1,j2         
            w(j-j1+ii)   = a(j)
            jw1(j-j1+ii) = ja(j)
            k = k + 1         
            jw2(ja(j))  = k
         end do
          

         k = 0         
         do j = j3,j4   
            ! in case of csr iw1=ja
            jw1(j-j3+1)  = iw1(j)
            k = k + 1         
            jw2(iw1(j))  = k
         end do
!!$         write(*,*) ii
!!$         write(*,*) 'range',j1,j2, j3,j4
!!$         write(*,'(10(1x,e8.1))') (w(k),k=1,n)  
!!$         write(*,'(10I6)') (jw1(k),k=1,n) 
!!$         write(*,'(10I6)') (jw2(k),k=1,n) 
!!$         write(*,*) ' ' 
!
         !----BLAS-Method-----------------------------------------------------
         !  Unpack Upper and Lower part of row of A in array w 
         !         lenu  = j2 - j1 + 1
         !         lenu0 = lenu
         !         lenl  = j4 - j3 + 1
         !         call DCOPY(lenu,a(j1),1,w(ii),1)
         !         call SCOPY(lenu,ja(j1),1,jw1(ii),1)
         !         call SCOPY(lenl,iw1(j3),1,jw1(1),1)
         !   
         !  Set non-zero indicators
         !         
         !         k = ii-1
         !         do j = j1,j2
         !            k = k + 1
         !            jw2(ja(j)) = k
         !         enddo
         !         k = 0
         !         do j = j3,j4
         !            k = k + 1
         !            jw2(iw1(j)) = k
         !         enddo
         !----End of BLAS-Method----------------------------------------------
         !
         !  Eliminate previous rows
         !
         jj = 0 
150      jj = jj+1
         if (jj .gt. lenl) goto 160
         !
         !  Determine the smallest column index
         !         
         jrow = jw1(jj)
         k = jj
         do j = jj+1,lenl
            if (jw1(j) .lt. jrow) then
               jrow = jw1(j)
               k = j
            endif
         enddo

         if (k .ne. jj) then
            ! Exchange in jw1
            j = jw1(jj)
            jw1(jj) = jw1(k)
            jw1(k) = j
            ! Exchange in jw1
            jw2(jrow) = jj
            jw2(j) = k
         endif
         !write(*,*) 'ii, jrow',ii, jrow

         !
         !  Zero out element in row by setting jw2(jrow) to zero
         !
         jw2(jrow) = 0
         !         
         !  Get the leading term of the row if it exists a(jj,ii)
         !  
         rowHead = ju(jrow)
         rowMid  = jw3(jrow)
         if (ju(rowMid) .eq. ii) then  
            !------------------------------------------------------------------
            !  Mettere la possibilità di saltare 
            !  questo ciclo se il termine pivotale
            !  è troppo piccolo: fact = U(rowMid)
            !-----------------------------------------------------------------
            !
            ! Exists, mark fill-in terms of the current row before the
            ! diagonal            
            !
            do k = rowHead,rowMid-1
               jcol = ju(k)
               jpos = jw2(jcol)
               if (jpos .eq. 0) then
                  ! This is a fill-in term
                  lenl = lenl+1
                  if (lenl .gt. n) goto 996
                  jw1(lenl) = jcol
                  jw2(jcol) = lenl
               endif
            enddo
            !                    
            ! Combine current row ii with row jrow
            !

            fact = U(rowMid)*U(jrow)
            ! fact = U(rowMid) * U(jrow) * sqrt(U(jrow))
            rowEnd = ju(jrow+1) - 1
            do k = rowMid,rowEnd
               !
               !write(*,*) ii, rowMid,ju(rowMid), jrow,  ju(k), k
               s = fact * U(k)
               jcol = ju(k)
               jpos = jw2(jcol)
               if (jpos .eq. 0) then
                  !                          
                  ! This is a fill-in element
                  !
                  lenu = lenu + 1
                  if (lenu .gt. n) goto 996
                  i = ii+lenu-1
                  jw1(i) = jcol
                  jw2(jcol) = i
                  w(i) = -s
                  !                     
               else
                  !
                  ! This is not a fill-in element
                  !
                  w(jpos) = w(jpos) - s
                  !                     
               endif
               !                  
            enddo
            !               
            ! Update pointer jw3
            !
            jw3(jrow) = jw3(jrow) + 1
            if (jw3(jrow) .gt. rowEnd) then
               ! Deactivate pointer
               jw3(jrow) = iwk
            endif
            !               
         endif
         !           
         goto 150
160      continue
         !
         !  Reset non-zero indicators
         !
         do k = 1,lenu
            jw2(jw1(ii+k-1)) = 0
         enddo
         !
         !  Store the diagonal term
         !
         s = w(ii)         
         if (s .le. zero) then
            write(iout,100) ii,s
            write(6,100) ii,s
            write(6,*) 'ilu error'
            !stop
            ierr = -5  
            write(iout,101) tnorm
            return
         endif
         U(ii) = s
         s = 1.d0/s
         !
         !  Apply dropping strategy
         !
         abstol = s * tnorm * droptol  ! Set value of true tolerance
         lenght = 0 
         do k = ii+1,ii+lenu-1
            if (abs(w(k)) .gt. abstol) then
               lenght = lenght+1
               w(ii+lenght)   = w(k)
               jw1(ii+lenght) = jw1(k)
            endif
         enddo
         !
         !  Set lenu = number of elements bigger than tolerance        
         !
         lenu = lenght
         ncut = min0(lenu,lenu0+lfil-1)
         !----------------------------------------------------------------
         !  May be interesting to save a fixed number of non-zero 
         !  for each line
         ! ncut = min0(lenu,lfil)         
         !----------------------------------------------------------------
         !
         !  Select the ncut biggest elements only if lenu gt zero
         !
         if (lenu .gt. 0) then
            jpos = ii+1
            call qsplit(w(jpos),jw1(jpos),lenu,ncut)
            !          
            !  Order them in increasing jw1
            !
            !call HPSORT(jw1(jpos),w(jpos),ncut)
            call sort_row(ncut,&
                 jw1(jpos:jpos+ncut-1),&
                 w(jpos:jpos+ncut-1))
            !
            !  Scale and store in the factor
            !         
            if ((ju0+ncut-1) .gt. iwk-1) goto 997
            do k = 0,ncut-1
               ! EF
               U(ju0+k)  = s*w(jpos+k)
               ju(ju0+k) = jw1(jpos+k) 
            enddo
            !            
            !----BLAS-Method-------------------------------------------------
            !  Scale by the pivotal term and store
            !            if ((ju0+ncut-1) .gt. iwk-1) goto 997
            !            call DCOPY(ncut,w(jpos),1,U(ju0),1)
            !            call DSCAL(ncut,s,w(jpos),1)
            !            call SCOPY(ncut,jw1(jpos),1,ju(ju0),1)
            !----End of BLAS-Method------------------------------------------
            !            
         endif
         !
         !  Update pointer to beginning of next row of U    
         !
         ju0 = ju0 + ncut
         ju(ii+1) = ju0 
         jw3(ii+1) = ju0 
         !
      enddo


      !-----------------------------------------------------------------------
      !     end main loop
      !-----------------------------------------------------------------------
      !
      !     Successful run
      !
      ierr = 0
      !     
      return
      !
      !     incomprehensible error. Matrix must be wrong.
      !     
996   ierr = -1
      return
      !    
      !     insufficient storage in U.
      !     
997   ierr = -2
      return
      !     
      !     illegal lfil entered.
      !     
998   ierr = -3
      return
      !     
      !     zero row encountered
      !     
999   ierr = -4
      return
      !----------------end-of-ilut--------------------------------------------
100   format(' NULL OR NEGATIVE DIAGONAL ELEMENT: I,J =',I8,2X,E16.5)
101   format(' NORM OF ROW =  ',E16.8)
      !-----------------------------------------------------------------------
    end subroutine my_ssr_incomplete_cheolesky

    subroutine my_unpackja(n,nterm,ia,ja,point,iw,first_row)
      !----------------------------------------------------------------------  
      !  Subroutine to unpack array ja and create array iw = "transposed"
      !  of array ja      
      !-----------------------------------------------------------------------
      implicit none
      !Input variables
      integer n,nterm
      integer ia(n+1),ja(nterm),point(n)
      !Output variables
      integer iw(nterm+1)      
      integer first_row(n+1) 
      !Local variables
      integer i,j,k,mm,nn   ,m   
      !
      !
      !Initialise pointers
      !
      do i = 1,n
         point(i) = 0
      enddo
      !
      ! Evaluate lenght of pre-diagonal terms for
      ! each line of array iw
      !
      do i = 1,n
         mm = ia(i)+1
         nn = ia(i+1)-1
         do j = mm,nn
            point(ja(j)) = point(ja(j)) + 1
         enddo
      enddo
      
      !
      ! Set pointers first_row
      !    
      first_row(1) = 1
      do i = 2,n+1
         !write(*,*) i-1, point(i-1)
         first_row(i) = first_row(i-1) + point(i-1)
         point(i-1)   = first_row(i-1)-1  
         !write(*,*) i, first_row(i)         
      enddo

      
      !
      !  Set iw1
      !
      m=0
      iw=0
      do i = 1,n
         mm = ia(i) + 1
         nn = ia(i+1)-1
         do j = mm,nn
            k = ja(j)
            point(k) = point(k) + 1            
            iw( point(k) ) = i
         enddo
      enddo
      
!!$      iw(1)=iw(2)-1
!!$      do i = 1,n
!!$         mm = first_row(i)
!!$         nn = first_row(i+1)-1
!!$         !write(*,*) 'range,',   mm, nn 
!!$         do j = mm,nn
!!$            write(*,*) i,iw(j)                    
!!$         enddo
!!$      enddo

      !
      return
      !
    end subroutine my_unpackja

    !----------------------------------------------------------------------- 
    subroutine qsplit(a,ind,n,ncut)
      use Globals
      integer,           intent(in   ):: n
      integer,           intent(inout):: ind(n)
      integer,           intent(in   ):: ncut
      real(kind=double), intent(inout) :: a(n)
      !-----------------------------------------------------------------------
      !     does a quick-sort split of a real array.
      !     on input a(1:n). is a real array
      !     on output a(1:n) is permuted such that its elements satisfy:
      !
      !     abs(a(i)) .ge. abs(a(ncut)) for i .lt. ncut and
      !     abs(a(i)) .le. abs(a(ncut)) for i .gt. ncut
      !
      !     ind(1:n) is an integer array which permuted in the same way as a(*).
      !-----------------------------------------------------------------------
      !local
      real(kind=double) :: tmp, abskey
      integer itmp, first, last
      integer i,j,mid
      logical keep_cycling
      !-----
      first = 1
      last = n
      if (ncut .lt. first .or. ncut .gt. last) return
      !
      !     outer loop -- while mid .ne. ncut do
      !

      keep_cycling = .true.
      do while ( keep_cycling ) 
         mid = first
         abskey = abs(a(mid))
         do j=first+1, last
            if (abs(a(j)) .gt. abskey) then
               mid = mid+1
               !     interchange
               tmp = a(mid)
               itmp = ind(mid)
               a(mid) = a(j)
               ind(mid) = ind(j)
               a(j)  = tmp
               ind(j) = itmp
            endif
         end do
         !
         !     interchange
         !
         tmp = a(mid)
         a(mid) = a(first)
         a(first)  = tmp
         !
         itmp = ind(mid)
         ind(mid) = ind(first)
         ind(first) = itmp
         !
         !     test for while loop
         !
         if (mid .eq. ncut) keep_cycling = .false.
         if (mid .gt. ncut) then
            last = mid-1
         else
            first = mid+1
         end if
      end do
    end subroutine qsplit



    !
    ! incomplete Choleski decompostion with no-fillin
    !
    subroutine kersh_loc(iout,nequ,nterm,job,ia,ja,sysmat,prec,info)
      use Globals
      implicit none
      integer  iout,nequ,nterm,job
      integer  ia(nequ+1),ja(nterm)
      integer  i,j,k,kk,k1,i1,j1,k2
      integer :: nterm_prec
      real(kind=double) :: prec(nterm),sysmat(nterm)
      real(kind=double) :: a
      integer, intent(inout) :: info

      info=0
      !
      do k=1,nterm
         prec(k) = zero
      end do
      !
      do kk=1,nequ-1
         k = ia(kk)
         a = sysmat(k) - prec(k)
         if (a.le.zero) then
            info = info+1
            !write(iout,100) kk,a
            !write(iout,101) prec(ia(kk-1))
            !a = (prec(ia(kk-1)))**2
            
            select case (job) 
            case (0) 
               return 
            case (1)
               write(iout,100) kk ,a
               write(iout,101) prec(ia(kk-1))
               a = (prec(ia(kk-1)))**2
            case (2)              
               write(iout,100) nequ ,a
               if ( abs(a) < large ) then 
                  write(iout,101) prec(ia(nequ-1))
                  a = (prec(ia(kk-1)))**2
               else
                  return
               end if
            end select
         end if
         prec(k) = sqrt(a)
         !
         i = ia(kk) + 1
         j = ia(kk+1) - 1
         !
         do k1 = i,j
            prec(k1) = (sysmat(k1) - prec(k1)) / prec(k)
         end do
         !
         do k2 = i,j-1
            j1 = ia(ja(k2))
            prec(j1) = prec(j1) + prec(k2)**2
            i1 = k2 + 1
            j1 = j1 + 1
            do while (j1.lt.ia(ja(k2)+1).and.i1.le.j)
               if (ja(j1).eq.ja(i1)) then
                  prec(j1) = prec(j1) + prec(k2) * prec(i1)
                  i1 = i1 + 1
                  j1 = j1 + 1
               else if (ja(j1).lt.ja(i1)) then
                  j1 = j1 + 1
               else if (ja(j1).gt.ja(i1)) then
                  i1 = i1 + 1
               end if
            end do
         end do
         !
         if (j.ge.i) prec(ia(ja(j))) = prec(ia(ja(j))) + prec(k2)**2
      end do
      !
      k = ia(nequ) 
      a = sysmat(k) - prec(k)
      if (a.le.zero) then
         select case (job) 
         case (0) 
            info = nequ
            return 
         case (1)
            write(iout,100) nequ ,a
            write(iout,101) prec(ia(nequ-1))
            a = (prec(ia(nequ-1)))**2
         case (2)              
            write(iout,100) nequ ,a
            if ( abs(a) < large ) then 
               write(iout,101) prec(ia(nequ-1))
               a = (prec(ia(nequ-1)))**2
            else
               info = nequ
               return
            end if
         end select
      end if
      prec(k) = sqrt(a)
      !
      return
100   format(' ELEMENTO DIAGONALE DI L NULLO,I,J =',I5,2X,E16.5)
101   format(' ELEMENTO DIAGONALE PRECEDENTE =  ',E16.8)
    end subroutine kersh_loc


  end subroutine incomplete_cholesky


      subroutine  mss2scaled_upper(lun_err,&
         nrow,&
         iwk,&
         ju,&
         U,&
         upper)
      use Globals
      !-------------------------------------------------------------------
      ! Change the format of a given matrix from Modified Symmetric Sparse
      ! Row format to CSR. The upper matrix is 
      !             upper = diag(U)^{-1} U
      ! thus its has one on the diagonal
      !-------------------------------------------------------------------
      implicit none
      ! Input variables      
      integer,           intent(in   ) :: lun_err
      integer,           intent(in   ) :: nrow
      integer,           intent(in   ) :: iwk
      integer,           intent(in   ) :: ju(iwk)
      real(kind=double), intent(in   ) :: U(iwk)
      ! Output variables
      type(spmat),       intent(inout) :: upper
      !  Local variables
      integer i,j,k
      integer ind,lenght,nterm
      real(kind=double) :: d_term


      !
      ! init upper factor
      !
      if ( upper%is_initialized ) then
         call upper%kill(lun_err)
      end if
      nterm = ju(nrow+1)-1
      call upper%init(lun_err, &
           nrow, nrow, nterm,&
           'csr',&
           is_symmetric=.false.,&
           unitary_diag=.true.,&
           triangular='U')

      !
      ! copy into the csr the matrix diag(U)^{-1} U
      ! 
      upper%ia(1) = 1
      do i = 1,nrow
         ! diagonal
         k      = ju(i)
         lenght = ju(i+1)-k
         ind           = upper%ia(i)
         upper%ia(i+1) = ind+lenght+1
         upper%coeff(ind) = one
         upper%ja(ind)    = i
         ! extra diagonal
         ind = ind+1
         do j = 0,lenght-1
            upper%coeff(ind+j)  = U(k+j)
            upper%ja(ind+j)     = ju(k+j)
         end do
      end do

    end subroutine mss2scaled_upper

  

  !-------------------------------------------------------
  ! Subroutine for incomplete LU factorization
  !-------------------------------------------------------
  subroutine incomplete_lu(this,&
       lun_err,n_fillin,tol_fillin,job,&
       info,&
       lower, upper)
    use Globals
    implicit none
    class(spmat),               intent(in   ) :: this
    integer,                    intent(in   ) :: lun_err
    integer,                    intent(in   ) :: n_fillin
    real(kind=double),          intent(in   ) :: tol_fillin
    integer,                    intent(in   ) :: job
    integer,                    intent(inout) :: info
    type(spmat),                intent(inout) :: lower
    type(spmat),                intent(inout) :: upper

    !local
    logical :: rc
    integer :: res, i, j
    integer :: nterm_max, nterm 
    integer :: nrow
    integer,           allocatable :: jlu(:),ju(:),jw(:),iscr(:)
    integer,           allocatable :: iw1(:),jw1(:),jw2(:),jw3(:)
    integer,           allocatable :: idiag(:),iwork(:)
    real(kind=double), allocatable :: alu(:),w(:),rwork(:)
    type(spmat)  :: lower_loc


    nrow = this%nrow
    nterm = this%nterm
    
    if ( this%storage_system .eq. 'ssr') then
       rc = IOerr(lun_err, err_inp, 'incomplete_lu', &
            'storage system ssr not supported yet')
    end if
    

    if ( n_fillin .eq. 0 ) then
       !
       ! use ILU0 from Saad
       !

       !
       ! prepare work arrays and diagonal pointer
       !
       allocate (iwork(nrow), idiag(nrow),rwork(nterm),stat=res)
       if( res .ne. 0) rc = IOerr(lun_err, err_alloc, 'incomplete_lu', &
            'work arrays idiag iwork rwork' )

       !
       ! build LU decomposition in Modified Compressed Storage
       !
       call loc_ilu0(lun_err,&
            nrow,nterm,info,this%ia,this%ja,idiag,iwork,&
            this%coeff,rwork)
!!$       do i=1,nrow
!!$          write(60,*) i,this%ja(idiag(i))
!!$          do j=this%ia(i),this%ia(i+1)-1
!!$             write(60,*) i,this%ja(j),rwork(j)
!!$          end do
!!$       end do
       
       !
       ! convert to scr format
       ! init. lower and upper matrix
       !
       call ilu02lu(lun_err,nrow,nterm,&
            rwork,this%ia,this%ja,idiag,&
            lower, upper)

       !
       ! prepare work arrays and diagonal pointer
       !
       deallocate (iwork, idiag,rwork,stat=res)
       if( res .ne. 0)  rc = IOerr(lun_err, err_dealloc, 'ilu0', &
            'work arrays idiag iwork rwork' )
    else
       !
       ! use ILUT from Saad
       !
       ! prepare work arrays 
       nrow      = this%nrow
       nterm_max = 2 * (n_fillin+1) * nrow

              
       allocate(&
            alu(nterm_max), &
            jlu(nterm_max),&
            ju(nrow),&
            w(nrow+1),&
            jw(2*nrow),&
            iscr(nrow+1),&
            stat=res)
       if (res.ne.0) rc = IOerr(lun_err, err_alloc, 'incomplete_lu', &
            'work arrays alu, jlu ju w jw',res)

       !
       ! non-symmetric case 
       !
       info =0
       call loc_ilut(nrow,&
            this%coeff,&
            this%ja,&
            this%ia,&
            n_fillin,&
            tol_fillin,&
            job,&
            alu,jlu,&
            ju,&
            nterm_max,&
            w,jw,info)

       if ( info .eq. 0 ) then
          !
          ! convert msr into csr format
          ! init. lower and upper matrix
          !    
          call ilut2lu(lun_err,nrow,nterm_max,&
               alu,jlu,ju,&! jw, &
               lower, upper)
       else
          !
          ! info in case of lu factorization error 
          !
          write(lun_err,*) ' Error in ILUT construction' 
          select case (info) 
          case (-1)
             rc = IOerr(lun_err, err_val, 'build_ic_fillin', &
                  ' The elimination process has generated'//&
                  ' a row in U whose length is .gt.  n)')
          case (-2)  
             rc = IOerr(lun_err, err_val, 'build_ic_fillin', &
                  'The matrix U overflows the array U')
          case (-3)  
             rc = IOerr(lun_err, err_val, 'build_ic_fillin', &
                  ' Illegal value for n_fillin')
          case(-4)   
             rc = IOerr(lun_err, err_val, 'build_ic_fillin', &
                  'zero row encountered')
          case(-5) 
             rc = IOerr(lun_err, err_val, 'build_ic_fillin', &
                  'Non Positive Diagonal Element')
          end select
          if (info .gt. 0) then
             rc = IOerr(lun_err, err_val, 'build_ic_fillin', &
                  'Zero pivot encountered at step number=',info)
          end if
          return
       end if


       !
       ! free memory
       !
       if ( lower_loc%is_initialized ) call lower_loc%kill(lun_err)

       deallocate(&
            alu, &
            jlu,&
            ju,&
            w,&
            jw,&
            stat=res)
       if (res.ne.0) rc = IOerr(lun_err, err_dealloc, 'incomplete_lu', &
            'work arrays alu, jlu ju w jw',res)  
    end if

  contains
    !----------------------------------------------------------------------- 
    subroutine qsplit(a,ind,n,ncut)
      use Globals
      integer,           intent(in   ):: n
      integer,           intent(inout):: ind(n)
      integer,           intent(in   ):: ncut
      real(kind=double), intent(inout) :: a(n)
      !-----------------------------------------------------------------------
      !     does a quick-sort split of a real array.
      !     on input a(1:n). is a real array
      !     on output a(1:n) is permuted such that its elements satisfy:
      !
      !     abs(a(i)) .ge. abs(a(ncut)) for i .lt. ncut and
      !     abs(a(i)) .le. abs(a(ncut)) for i .gt. ncut
      !
      !     ind(1:n) is an integer array which permuted in the same way as a(*).
      !-----------------------------------------------------------------------
      !local
      real(kind=double) :: tmp, abskey
      integer itmp, first, last
      integer i,j,mid
      logical keep_cycling
      !-----
      first = 1
      last = n
      if (ncut .lt. first .or. ncut .gt. last) return
      !
      !     outer loop -- while mid .ne. ncut do
      !

      keep_cycling = .true.
      do while ( keep_cycling ) 
         mid = first
         abskey = abs(a(mid))
         do j=first+1, last
            if (abs(a(j)) .gt. abskey) then
               mid = mid+1
               !     interchange
               tmp = a(mid)
               itmp = ind(mid)
               a(mid) = a(j)
               ind(mid) = ind(j)
               a(j)  = tmp
               ind(j) = itmp
            endif
         end do
         !
         !     interchange
         !
         tmp = a(mid)
         a(mid) = a(first)
         a(first)  = tmp
         !
         itmp = ind(mid)
         ind(mid) = ind(first)
         ind(first) = itmp
         !
         !     test for while loop
         !
         if (mid .eq. ncut) keep_cycling = .false.
         if (mid .gt. ncut) then
            last = mid-1
         else
            first = mid+1
         end if
      end do
    end subroutine qsplit

    !
    !************************* ILU0 *************************************
    !
    !>--------------------------------------------------------------
    !> Subroutine for the computation of the lu decompostion
    !> in the mss format
    !> DIAGONAL IS NOT INVERTED
    !>--------------------------------------------------------------
    subroutine loc_ilu0(lun_err,&
         nequ,nterm,info,ia,ja,idiag,iwork,&
         sysmat,prec)
      !
      !  calculates the ILU(0) factorization of matrix SYSMAT
      !  (algorithm 10.4 page 277 in Saad 1996)
      !
      use Globals
      implicit none
      integer :: lun_err,nequ,nterm,info   
      integer :: irow,jrow,ind1,ind2,j,jnd,jw
      integer :: ia(nequ+1),ja(nterm),idiag(nequ),iwork(nequ)
      real(kind=double) :: sysmat(nterm),prec(nterm)
      real(kind=double) :: mult


      info = 0 
      !
      !  initializes prec to sysmat
      !
      call dcopy(nterm,sysmat,1,prec,1)

      do j=1,nequ
         iwork(j)=0
      end do

      do irow=1,nequ

         ind1 = ia(irow)
         ind2 = ia(irow+1)-1

         do j=ind1,ind2
            iwork(ja(j)) = j
         end do

         jnd=ind1
         jrow = ja(jnd)
         do while (jnd.le.ind2 .and. jrow.lt.irow)
            
            mult = prec(jnd)*prec(idiag(jrow))
            prec(jnd) = mult

            do j=idiag(jrow)+1,ia(jrow+1)-1
               jw=iwork(ja(j))
               if(jw.ne.0) prec(jw)=prec(jw)-mult*prec(j)
            end do

            jnd=jnd+1
            jrow = ja(jnd)

         end do

         !
         !  set pointer to diagonal element
         !
         idiag(irow)=jnd
         !
         !   check if diagonal element exists and MULT is nonzero 
         !
         if(jrow.ne.irow .or. prec(jnd).eq.zero) then
            write(lun_err,*)' ILU(0) error at row',irow
            info=irow
            return
         end if

         prec(jnd)=one/prec(jnd)

         do j=ind1,ind2
            iwork(ja(j))=0
         end do

      end do

      !
      ! invert diagonal back
      !
      do j=1,nequ
         prec(idiag(j)) = one / prec(idiag(j))
      end do

    end subroutine loc_ilu0

    !>---------------------------------------------------
    !> Converts the lud factors from the mss format to 
    !> two matrix (upper, lower)  in csr format 
    !>       A ~ L U
    !> L having one on the diagonal
    !>---------------------------------------------------
    subroutine ilu02lu(lun_err,nrow,nterm,&
         a_lu,ia_lu,ja_lu,idiag,&
         lower, upper)
      use Globals
      implicit none

      integer,           intent(in   ) :: lun_err
      integer,           intent(in   ) :: nrow
      integer,           intent(in   ) :: nterm
      real(kind=double), intent(in   ) :: a_lu(nterm)
      integer,           intent(in   ) :: ia_lu(nrow+1)
      integer,           intent(in   ) :: ja_lu(nterm)
      integer,           intent(in   ) :: idiag(nrow)
      type(spmat),       intent(inout) :: lower, upper
      !local
      logical :: rc
      integer :: res 
      integer :: irow, j,start,i
      integer :: low_ind,up_ind
      integer, allocatable :: lower_ia(:),upper_ia(:)
      integer, allocatable :: lower_ja(:),upper_ja(:)
      real(kind=double), allocatable :: lower_coeff(:),upper_coeff(:)


      !
      ! allocation of work arrays 
      !
      allocate(&
           lower_ia(nrow+1),&
           upper_ia(nrow+1),&
           lower_ja(nterm),&
           upper_ja(nterm),&
           lower_coeff(nterm),&
           upper_coeff(nterm),&
           stat=res)
      if (res.ne.0) rc = IOerr(lun_err, err_alloc, 'msr2ldu', &
           'work arrays lower_ia, lower_ja, lower_coeff'//&
           ' upper_ia, upper_ja, upper_coeff',res) 

      !
      ! init counter
      !
      low_ind= 0
      up_ind = 0
      lower_ia(1) = 1
      upper_ia(1) = 1

      do irow=1,nrow     
         !
         ! copy lower matrix
         !
         do j = ia_lu(irow), idiag(irow)-1 
            low_ind              = low_ind+1
            lower_ja(low_ind)    = ja_lu(j)
            lower_coeff(low_ind) = a_lu(j)
         end do
         ! add one on the diagonal
         low_ind              = low_ind+1
         lower_ja(low_ind)    = irow
         lower_coeff(low_ind) = one
         ! set next index of first element 
         lower_ia(irow+1) = low_ind+1
         !write(*,*) 'irow, lower_ia(irow+1)',irow, lower_ia(irow+1)
         !write(*,*) lower_coeff(lower_ia(irow):lower_ia(irow+1)-1)

         !
         ! copy upper matrix
         !
         ! the diagonal term
         up_ind              = up_ind+1
         upper_ja(up_ind)    = irow
         upper_coeff(up_ind) = a_lu(idiag(irow))
         ! extra diagonal terms
         do j = idiag(irow)+1, ia_lu(irow+1)-1
            up_ind              = up_ind+1
            upper_ja(up_ind)    = ja_lu(j)
            upper_coeff(up_ind) = a_lu(j)
         end do
         ! set next index of first element
         upper_ia(irow+1) = up_ind+1
      end  do

      !
      ! init and copy lower sparse matrix
      !
      if (lower%is_initialized) call lower%kill(lun_err)
      call lower%init(lun_err, &
           nrow, nrow, low_ind,&
           'csr',&
           is_symmetric=.false.,&
           unitary_diag = .true.,&
           triangular='L')

      lower%ia               = lower_ia
      lower%ja   (1:low_ind) = lower_ja   (1:low_ind)
      lower%coeff(1:low_ind) = lower_coeff(1:low_ind)
      lower%is_sorted = .true.
      call lower%sort()


      !
      ! init and copy upper sparse matrix
      !
      if (upper%is_initialized) call  upper%kill(lun_err)
      !
      call upper%init(lun_err, &
           nrow, nrow, up_ind,&
           'csr',&
           is_symmetric=.false.,&
           triangular='U')

      upper%ia              = upper_ia
      upper%ja(1:up_ind)    = upper_ja(1:up_ind)
      upper%coeff(1:up_ind) = upper_coeff(1:up_ind)
      upper%is_sorted = .true.
      call upper%sort()

      !
      ! free memory 
      !
      deallocate(&
           lower_ia,&
           upper_ia,&
           lower_ja,&
           upper_ja,&
           lower_coeff,&
           upper_coeff,&
           stat=res)
      if (res.ne.0) rc = IOerr(lun_err, err_dealloc, 'msr2ldu', &
           'work arrays lower_ia, lower_ja, lower_coeff'//&
           ' upper_ia, upper_ja, upper_coeff',res) 

    end subroutine ilu02lu


    
    !---------------------------------------------------------
    !                          S P A R S K I T  
    !---------------------------------------------------------
    !                   ITERATIVE SOLVERS MODULE 
    !---------------------------------------------------------
    ! This Version Dated: August 13, 1996. 
    !  Warning: meaning of some  arguments have changed 
    ! w.r.t. earlier versions. 
    ! Some Calling sequences may also have changed
    !--------------------------------------------------------
    subroutine loc_ilut(n,a,ja,ia,lfil,droptol,job,alu,jlu,ju,iwk,w,jw,ierr)
      !--------------------------------------------------------------------
      implicit none 
      integer n ,job
      real*8 a(*),alu(*),w(n+1),droptol
      integer ja(*),ia(n+1),jlu(*),ju(n),jw(2*n),lfil,iwk,ierr
      !--------------------------------------------------------------------*
      !                      *** ILUT preconditioner ***                   *
      !      incomplete LU factorization with dual truncation mechanism    *
      !--------------------------------------------------------------------*
      !     Author: Yousef Saad *May, 5, 1990, Latest revision, August 1996*
      !--------------------------------------------------------------------*
      ! PARAMETERS                                                           
      !-----------                                                           
      !
      ! on entry:
      !========== 
      ! n       = integer. The row dimension of the matrix A. The matrix 
      !
      ! a,ja,ia = matrix stored in Compressed Sparse Row format.           
      !
      ! lfil    = integer. The fill-in parameter. Each row of L and each row
      !           of U will have a maximum of lfil elements (excluding the 
      !           diagonal element). lfil must be .ge. 0.
      !           ** WARNING: THE MEANING OF LFIL HAS CHANGED WITH RESPECT TO
      !           EARLIER VERSIONS. 
      !
      ! droptol = real*8. Sets the threshold for dropping small terms in the
      !           factorization. See below for details on dropping strategy.
      !
      !  
      ! iwk     = integer. The lengths of arrays alu and jlu. If the arrays
      !           are not big enough to store the ILU factorizations,i lut
      !           will stop with an error message. 
      !
      ! On return:
      !===========
      !
      ! alu,jlu = matrix stored in Modified Sparse Row (MSR)
      !           format containing the L and U factors together.
      !           The diagonal (stored in alu(1:n) ) is not inverted.
      !           Each i-th row of the alu,jlu matrix
      !           contains the i-th row of L (excluding the diagonal entry=1)
      !           followed by the i-th row of U.
      !
      ! ju      = integer array of length n containing the pointers to
      !           the beginning of each row of U in the matrix alu,jlu.
      !
      ! ierr    = integer. Error message with the following meaning.
      !           ierr  = 0   --> successful return.
      !           ierr .gt. 0 --> zero pivot encountered at step number ierr.
      !           ierr  = -1  --> Error. input matrix may be wrong.
      !                            (The elimination process has generated a
      !                            row in L or U whose length is .gt.  n.)
      !           ierr  = -2  --> The matrix L overflows the array al.
      !           ierr  = -3  --> The matrix U overflows the array alu.
      !           ierr  = -4  --> Illegal value for lfil.
      !           ierr  = -5  --> zero row encountered.
      !
      ! work arrays:
      !=============
      ! jw      = integer work array of length 2*n.
      ! w       = real work array of length n 
      !  
      !----------------------------------------------------------------------
      ! w, ju (1:n) store the working array [1:ii-1 = L-part, ii:n = u] 
      ! jw(n+1:2n)  stores nonzero indicators
      ! 
      ! Notes:
      ! ------
      ! The diagonal elements of the input matrix must be  nonzero (at least
      ! 'structurally'). 
      !
      !--------------------------------------------------------------------* 
      !---- Dual drop strategy works as follows.                           *
      !                                                                    *
      !     1) Theresholding in L and U as set by droptol. Any element     *
      !        whose magnitude is less than some tolerance (relative to    *
      !        the abs value of diagonal element in u) is dropped.         *
      !                                                                    *
      !     2) Keeping only the largest lfil elements in the i-th row of L * 
      !        and the largest lfil elements in the i-th row of U          *
      !        (excluding diagonal elements).                              *
      !                                                                    *
      ! Flexibility: one  can use  droptol=0  to get  a strategy  based on *
      ! keeping  the largest  elements in  each row  of L  and U.   Taking *
      ! droptol .ne.  0 but lfil=n will give  the usual threshold strategy *
      ! (however, fill-in is then unpredictible).                          *
      !--------------------------------------------------------------------*
      !     locals
      integer ju0,k,j1,j2,j,ii,i,lenl,lenu,jj,jrow,jpos,len 
      real*8 tnorm, t, abs, s, fact 
      if (lfil .lt. 0) goto 998
      !---------------------------------------------------------------------
      !     initialize ju0 (points to next element to be added to alu,jlu)
      !     and pointer array.
      !---------------------------------------------------------------------
      ju0 = n+2
      jlu(1) = ju0
      !
      !     initialize nonzero indicator array. 
      !
      do 1 j=1,n
         jw(n+j)  = 0
1     end do
      !---------------------------------------------------------------------
      !     beginning of main loop.
      !---------------------------------------------------------------------
      do 500 ii = 1, n
         j1 = ia(ii)
         j2 = ia(ii+1) - 1
         !write(0,*) j1,j2
         tnorm = 0.0d0
         do 501 k=j1,j2
            tnorm = tnorm+abs(a(k))
501      end do
         if (tnorm .eq. 0.0) goto 999
         tnorm = tnorm/real(j2-j1+1)
         !     
         ! unpack L-part and U-part of row of A in arrays w 
         !     
         lenu = 1
         lenl = 0
         jw(ii) = ii
         w(ii) = 0.0
         jw(n+ii) = ii
         !
         do 170  j = j1, j2
            k = ja(j)
            t = a(j)
            if (k .lt. ii) then
               lenl = lenl+1
               jw(lenl) = k
               w(lenl) = t
               jw(n+k) = lenl
            else if (k .eq. ii) then
               w(ii) = t
            else
               lenu = lenu+1
               jpos = ii+lenu-1 
               jw(jpos) = k
               w(jpos) = t
               jw(n+k) = jpos
            endif
170      end do
         jj = 0
         len = 0 
         !     
         !     eliminate previous rows
         !     
150      jj = jj+1
         if (jj .gt. lenl) goto 160
         !---------------------------------------------------
         ! in order to do the elimination in the correct 
         ! order we must select  the smallest 
         !  column index among jw(k), k=jj+1, ..., lenl.
         !---------------------------------------------------
         jrow = jw(jj)
         k = jj
         !     
         ! determine smallest column index
         !     
         do 151 j=jj+1,lenl
            if (jw(j) .lt. jrow) then
               jrow = jw(j)
               k = j
            end if
151      end do
         !
         if (k .ne. jj) then
            ! exchange in jw
            j = jw(jj)
            jw(jj) = jw(k)
            jw(k) = j
            ! exchange in jr
            jw(n+jrow) = jj
            jw(n+j) = k
            ! exchange in w
            s = w(jj)
            w(jj) = w(k)
            w(k) = s
         endif
         !
         ! zero out element in row by setting 
         ! jw(n+jrow) to zero.
         !     
         jw(n+jrow) = 0
         !
         ! get the multiplier for row to be eliminated (jrow).
         !     
         fact = w(jj)*alu(jrow)
         if (abs(fact) .le. droptol) goto 150
         !     
         ! combine current row and row jrow
         !
         do 203 k = ju(jrow), jlu(jrow+1)-1
            s = fact*alu(k)
            j = jlu(k)
            jpos = jw(n+j)
            if (j .ge. ii) then
               !     
               !     dealing with upper part.
               !     
               if (jpos .eq. 0) then
                  !
                  ! this is a fill-in element
                  !     
                  lenu = lenu+1
                  if (lenu .gt. n) goto 995
                  i = ii+lenu-1
                  jw(i) = j
                  jw(n+j) = i
                  w(i) = - s
               else
                  !
                  ! this is not a fill-in element 
                  !
                  w(jpos) = w(jpos) - s
               endif
            else
               !     
               ! dealing  with lower part.
               !     
               if (jpos .eq. 0) then
                  !
                  !     this is a fill-in element
                  !     
                  lenl = lenl+1
                  if (lenl .gt. n) goto 995
                  jw(lenl) = j
                  jw(n+j) = lenl
                  w(lenl) = - s
               else
                  !     
                  !     this is not a fill-in element 
                  !     
                  w(jpos) = w(jpos) - s
               endif
            endif
203      end do
         !     
         ! store this pivot element -- 
         ! (from left to right -- no danger of
         ! overlap with the working elements in L (pivots). 
         !     
         len = len+1 
         w(len) = fact
         jw(len)  = jrow
         goto 150

         !
         ! reset double-pointer to zero (U-part)
         !     
160      do k=1, lenu
            jw(n+jw(ii+k-1)) = 0
         end do
         !     
         ! update L-matrix
         !     
         lenl = len 
         len = min0(lenl,lfil)
         !     
         ! sort by quick-split
         !
         call qsplit (w,jw,lenl,len)
         !
         ! store L-part
         ! 
         do 204 k=1, len 
            if (ju0 .gt. iwk) goto 996
            alu(ju0) =  w(k)
            jlu(ju0) =  jw(k)
            ju0 = ju0+1
204      end do
         !     
         !save pointer to beginning of row ii of U
         !     
         ju(ii) = ju0
         !
         ! update U-matrix -- first apply dropping strategy 
         !
         len = 0
         do k=1, lenu-1
            if (abs(w(ii+k)) .gt. droptol*tnorm) then 
               len = len+1
               w(ii+len) = w(ii+k) 
               jw(ii+len) = jw(ii+k) 
            endif
         end do
         lenu = len+1
         len = min0(lenu,lfil)
         !
         call qsplit (w(ii+1), jw(ii+1), lenu-1,len)
         !
         !     copy
         ! 
         t = abs(w(ii))
         if (len + ju0 .gt. iwk) goto 997
         do 302 k=ii+1,ii+len-1 
            jlu(ju0) = jw(k)
            alu(ju0) = w(k)
            t = t + abs(w(k) )
            ju0 = ju0+1
302      end do
         !     
         !     store inverse of diagonal element of u
         !   
         if (w(ii) .eq. 0.0) then
            select case (job) 
            case (0)
               info=ii
               return
            case (1)
               w(ii) = (0.0001 + droptol)*tnorm
            end select
         end if
         alu(ii) = 1.0d0/w(ii) 
         !     
         !     update pointer to beginning of next row of U.
         !     
         jlu(ii+1) = ju0
         !--------------------------------------------------------------
         !     end main loop
         !--------------------------------------------------------------
500   end do
      ! EF
      ! get the not inverted diagonal
      do ii=1,n
         alu(ii)=one/alu(ii)
      end do

      ierr = 0
      return
      !
      !     incomprehensible error. Matrix must be wrong.
      !     
995   ierr = -1
      return
      !     
      !     insufficient storage in L.
      !     
996   ierr = -ju0
      return
      !     
      !     insufficient storage in U.
      !     
997   ierr = -len-ju0
      return
      !     
      !     illegal lfil entered.
      !     
998   ierr = -4
      return
      !     
      !     zero row encountered
      !     
999   ierr = -5
      return
      !----------------end-of-ilut----------------------
      !-------------------------------------------------

    end subroutine loc_ilut

    subroutine ilut2lu(lun_err,nrow,nterm,&
         a_lu,ja_lu,ia_u,&
         lower, upper)
      use Globals
      implicit none

      integer,           intent(in   ) :: lun_err
      integer,           intent(in   ) :: nrow
      integer,           intent(in   ) :: nterm
      real(kind=double), intent(in   ) :: a_lu(nterm)
      integer,           intent(in   ) :: ja_lu(nterm)
      integer,           intent(in   ) :: ia_u(nrow)
      type(spmat),  intent(inout) :: lower, upper
      !local
      logical :: rc
      integer :: res 
      integer :: irow, j,start,i
      integer :: low_ind,up_ind
      integer, allocatable :: lower_ia(:),upper_ia(:)
      integer, allocatable :: lower_ja(:),upper_ja(:)
      real(kind=double), allocatable :: lower_coeff(:),upper_coeff(:)


      !
      ! allocation of work arrays 
      !
      allocate(&
           lower_ia(nrow+1),&
           upper_ia(nrow+1),&
           lower_ja(nterm),&
           upper_ja(nterm),&
           lower_coeff(nterm),&
           upper_coeff(nterm),&
           stat=res)
      if (res.ne.0) rc = IOerr(lun_err, err_alloc, 'msr2ldu', &
           'work arrays lower_ia, lower_ja, lower_coeff'//&
           ' upper_ia, upper_ja, upper_coeff',res) 

      !
      ! init counter
      !
      low_ind=0
      up_ind=0
      lower_ia(1)=1
      upper_ia(1)=1

!!$    do i=1,nrow
!!$       write(*,*) i,a_lu(i),ja_lu(i)
!!$    end do

      do irow=1,nrow     
         !
         ! copy lower matrix
         !
         !write(*,*) 'lu',irow
         do j = ja_lu(irow), ia_u(irow)-1 
            low_ind              = low_ind+1
            lower_ja(low_ind)    = ja_lu(j)
            lower_coeff(low_ind) = a_lu(j)
            !write(*,*) irow, ja_lu(j), a_lu(j) 
            !write(*,*) irow, lower_ja(low_ind),lower_coeff(low_ind)
         end do
         ! add one on the diagonal
         low_ind              = low_ind+1
         lower_ja(low_ind)    = irow
         lower_coeff(low_ind) = one
         ! set next index of first elemen 
         lower_ia(irow+1) = low_ind+1
         !write(*,*) 'irow, lower_ia(irow+1)',irow, lower_ia(irow+1)
         !write(*,*) lower_coeff(lower_ia(irow):lower_ia(irow+1)-1)

         !
         ! copy upper matrix
         !
         ! add one on the diagonal
         up_ind              = up_ind+1
         upper_ja(up_ind)    = irow
         upper_coeff(up_ind) = a_lu(irow)
         do j = ia_u(irow), ja_lu(irow+1)-1
            up_ind              = up_ind+1
            upper_ja(up_ind)    = ja_lu(j)
            upper_coeff(up_ind) = a_lu(j)
         end do
         upper_ia(irow+1) = up_ind+1
      end  do

      !
      ! init and copy lower sparse matrix
      !
      if (lower%is_initialized) call lower%kill(lun_err)
      call lower%init(lun_err, &
           nrow, nrow, low_ind,&
           'csr',&
           is_symmetric=.false.,&
           unitary_diag=.true.,&
           triangular='L')
      lower%ia               = lower_ia
      lower%ja   (1:low_ind) = lower_ja   (1:low_ind)
      lower%coeff(1:low_ind) = lower_coeff(1:low_ind)
      lower%is_sorted = .false.

      call lower%sort()


      if (upper%is_initialized) call  upper%kill(lun_err)
      !
      ! init and copy upper sparse matrix
      !
      call upper%init(lun_err, &
           nrow, nrow, up_ind,&
           'csr',&
           is_symmetric=.false.,&
           triangular='U')

      upper%ia              = upper_ia
      upper%ja(1:up_ind)    = upper_ja(1:up_ind)
      upper%coeff(1:up_ind) = upper_coeff(1:up_ind)

      upper%is_sorted = .false.
      call upper%sort()

      !
      ! free memory 
      !
      deallocate(&
           lower_ia,&
           upper_ia,&
           lower_ja,&
           upper_ja,&
           lower_coeff,&
           upper_coeff,&
           stat=res)
      if (res.ne.0) rc = IOerr(lun_err, err_dealloc, 'msr2ldu', &
           'work arrays lower_ia, lower_ja, lower_coeff'//&
           ' upper_ia, upper_ja, upper_coeff',res) 

    end subroutine ilut2lu
  end subroutine incomplete_lu

  !>--------------------------------------------------------------
  !> Subroutine for sorting the arrays ja and coeff
  !> to be ordered according to column index
  !>--------------------------------------------------------------
  subroutine sort_spmat(this)
    use Globals
    implicit none
    class(spmat), intent(inout) :: this
    !local
    integer :: irow,start,finish,nnz

    if (this%is_sorted) return

    do irow = 1, this%nrow
       ! select matrix line
       ! works both for ssr an csr format
       start  = this%ia(irow)
       finish = this%ia(irow+1)-1
       nnz    = finish-start+1
       ! order 
       if (nnz .gt. 1 ) then 
          call sort_row(nnz,&
               this%ja(start:finish),&
               this%coeff(start:finish))
       end if
    end do
    this%is_sorted = .true.
  contains
    subroutine sort_matrix_line(nnz,ja,coeff)
      use Globals
      implicit none
      integer,           intent(in   ) :: nnz
      integer,           intent(inout) :: ja(nnz)
      real(kind=double), intent(inout) :: coeff(nnz)
      !local 
      integer :: i,j, indx,isgn, itemp
      real(kind=double) :: rtemp

      !  Initialize.
      i = 0
      indx = 0
      isgn = 0
      j = 0
      do 
         call global_heapsort(nnz, indx, i,j,isgn)
         if (indx .gt. 0 ) then
            ! SWAP ELEMENT 

            ! swap column indeces
            itemp = ja(i)
            ja(i) = ja(j)
            ja(j) = itemp

            ! swap real nnzero coeff
            rtemp    = coeff(i)
            coeff(i) = coeff(j)
            coeff(j) = rtemp
         else if ( indx .lt. 0) then
            ! COMPARE
            isgn = 0
            if ( ja(i) .lt.  ja(j) ) isgn = -1
            if ( ja(i) .gt.  ja(j) ) isgn = 1
         else if ( indx .eq. 0 ) then
            exit
         end if
      end do
    end subroutine sort_matrix_line
  end subroutine sort_spmat


  subroutine MtimesN(matrix_m,matrix_n, out_matrix)
    use Globals
    implicit none
    type(spmat), intent(in) :: matrix_m
    type(spmat), intent(in) :: matrix_n
    real(kind=double), intent(inout) :: out_matrix(matrix_m%nrow,matrix_n%ncol)
    !local
    integer irow,icol, k, ind,ind2
    real(kind=double) :: b,c

    out_matrix = zero
    do irow=1,matrix_m%nrow
       do ind=matrix_m%ia(irow),matrix_m%ia(irow+1)-1
          b = matrix_m%coeff(ind)
          k = matrix_m%ja(ind)
          do ind2 = matrix_n%ia(k),matrix_n%ia(k+1)-1
             c    = matrix_n%coeff(ind2)
             icol = matrix_n%ja(ind2)
             out_matrix(irow,icol) = out_matrix(irow,icol) + &
                  b * c
          end do
       end do
    end do

  end subroutine MtimesN
  
  subroutine add_row(this,lun_err,vec_w)
    use Globals
    implicit none
    ! vars
    class(spmat),      intent(inout) :: this
    integer,           intent(in   ) :: lun_err
    real(kind=double), intent(in   ) :: vec_w(this%ncol)
    ! local vars
    integer :: icol
    integer :: nrow, ncol, nterm
    type(spmat) :: matrix_out
    
    call matrix_out%init(lun_err,this%nrow+1,this%ncol+1,&
         this%nterm+this%ncol+1,this%storage_system,&
         is_symmetric=.false.,&
         unitary_diag=.false.)

    if(this%triangular.eq.'L') matrix_out%triangular='L'

    nrow = this%nrow
    ncol = this%ncol
    nterm = this%nterm

    matrix_out%ia(1:nrow+1) = this%ia
    matrix_out%ia(nrow+2)   = nterm+ncol+2
    
    matrix_out%ja(1:nterm) = this%ja
    do icol = 1,ncol+1
       matrix_out%ja(icol+nterm) = icol
    end do

    matrix_out%coeff(1:nterm) = this%coeff
    do icol = 1,ncol
       matrix_out%coeff(icol+nterm) = vec_w(icol)
    end do
    matrix_out%coeff(nterm+ncol+1) = sum(vec_w)    

    select type (this)
    type is (spmat)
       this = matrix_out
    end select
    
    call matrix_out%kill(lun_err)   
    
  end subroutine add_row

    function bandwidth ( this , perm, iperm) 
    !************************************************************************
    !
    !! ADJ_BANDWIDTH computes the bandwidth of an adjacency matrix.
    !  Enrico Facca 2018-09-09
    !  Al inputs included the spamt format
    !
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    11 March 2005
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Alan George and Joseph Liu,
    !    Computer Solution of Large Sparse Positive Definite Systems,
    !    Prentice Hall, 1981.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
    !
    !    Input, integer ( kind = 4 ) ADJ_NUM, the number of adjacency entries.
    !
    !    Input, integer ( kind = 4 ) ADJ_ROW(NODE_NUM+1).  Information about row I
    !    is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
    !
    !    Input, integer ( kind = 4 ) ADJ(ADJ_NUM), the adjacency structure.
    !    For each row, it contains the column indices of the nonzero entries.
    !
    !    Output, integer ( kind = 4 ) ADJ_BANDWIDTH, the bandwidth of the adjacency
    !    matrix.
    !
    use Globals
    implicit none
    class(spmat),      intent(in   ) :: this
    integer, optional, intent(in   ) :: perm(this%nrow),iperm(this%nrow)
    integer ( kind = 4 ) :: bandwidth

    ! local 
    integer ( kind = 4 ) band_hi
    integer ( kind = 4 ) band_lo
    integer ( kind = 4 ) col
    integer ( kind = 4 ) i
    integer ( kind = 4 ) j

    band_lo = 0
    band_hi = 0


    if ( present(perm) .and. present(iperm) ) then
       do i = 1, this%nrow
          do j = this%ia(perm(i)), this%ia(perm(i)+1)-1
             col = iperm(this%ja(j))
             band_lo = max ( band_lo, i - col )
             band_hi = max ( band_hi, col - i )
          end do

       end do
    else
       do i = 1, this%nrow
          do j = this%ia(i), this%ia(i+1)-1
             col = this%ja(j)
             band_lo = max ( band_lo, i - col )
             band_hi = max ( band_hi, col - i )
          end do
       end do
    end if
          

    bandwidth = band_lo + 1 + band_hi

    return
  end function bandwidth

  
  subroutine genrcm ( this, lun_err,perm, inv_perm)
    !*************************************************************************
    !
    !! GENRCM finds the reverse Cuthill-Mckee ordering for a general graph.
    !  Enrico Facca :
    !   2018-09-09 The inputs of the algorithm are now included in the 
    !   csr format. Each irow-line of the matrix contians the index of the 
    !   graph nodes that are connnected to irow. All the indeces are listed in 
    !   thsi%ja
    !  
    !
    !  Discussion:
    !
    !    For each connected component in the graph, the routine obtains
    !    an ordering by calling RCM.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    04 January 2003
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Alan George, Joseph Liu.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Alan George and Joseph Liu,
    !    Computer Solution of Large Sparse Positive Definite Systems,
    !    Prentice Hall, 1981.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
    !
    !    Input, integer ( kind = 4 ) ADJ_NUM, the number of adjacency entries.
    !
    !    Input, integer ( kind = 4 ) ADJ_ROW(NODE_NUM+1).  Information about row
    !    I is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
    !
    !    Input, integer ( kind = 4 ) ADJ(ADJ_NUM), the adjacency structure.
    !    For each row, it contains the column indices of the nonzero entries.
    !
    !    Output, integer ( kind = 4 ) PERM(NODE_NUM), the RCM ordering.
    !
    !  Local Parameters:
    !
    !    Local, integer LEVEL_ROW(NODE_NUM+1), the index vector for a level
    !    structure.  The level structure is stored in the currently unused
    !    spaces in the permutation vector PERM.
    !
    !    Local, integer MASK(NODE_NUM), marks variables that have been numbered.
    !
    use GLobals
    implicit none

    class(spmat), intent(in   ) :: this
    integer,      intent(in   ) :: lun_err    
    integer,      intent(inout) :: perm(this%nrow)
    integer,      intent(inout) :: inv_perm(this%nrow)
        
    ! local
    logical :: rc
    integer :: res
    integer ( kind = 4 ) adj_num
    integer ( kind = 4 ) node_num
    integer ( kind = 4 ) i
    integer ( kind = 4 ) iccsze
    integer ( kind = 4 ), allocatable :: mask(:)
    integer ( kind = 4 ) level_num
    integer ( kind = 4 ), allocatable ::level_row(:)
    integer ( kind = 4 ) num
    integer ( kind = 4 ) root

    node_num = this%nrow
    adj_num  = this%nterm

    ! allocate work arrays
    allocate(&
         mask(node_num),&
         level_row(node_num+1),&
         stat=res)
    if (res .ne. 0) &
         rc = IOerr(lun_err, err_alloc, 'genrcm', &
         'work arrays mask, level_row',res) 
    
    !build permutation 
    
    mask(1:node_num) = 1

    num = 1
    
    do i = 1, node_num
       !
       !  For each masked connected component...
       !
       if ( mask(i) /= 0 ) then

          root = i
          !
          !  Find a pseudo-peripheral node ROOT.  The level structure found by
          !  ROOT_FIND is stored starting at PERM(NUM).
          !
          call root_find ( root, adj_num, this%ia, this%ja,&
               mask, level_num, &
               level_row, perm(num), node_num )
          !
          !  RCM orders the component using ROOT as the starting node.
          !
          call rcm ( root, adj_num, this%ia, this%ja,&

               mask, perm(num), iccsze, &
               node_num )

          num = num + iccsze
          !
          !  We can stop once every node is in one of the connected components.
          !
          if ( node_num < num ) then
             exit
          end if

       end if

    end do

    do i = 1, node_num
       inv_perm(perm(i)) = i
    end do

    
    ! clean memory
    deallocate(&
         mask,&
         level_row,&
         stat=res)

    if (res .ne. 0) &
         rc = IOerr(lun_err, err_dealloc, 'genrcm', &
         'work arrays mask, level_row',res) 

  contains
    subroutine root_find ( root, adj_num, adj_row, adj, mask, level_num, &
         level_row, level, node_num )

      !************************************************************************0
      !
      !! ROOT_FIND finds a pseudo-peripheral node.
      !
      !  Discussion:
      !
      !    The diameter of a graph is the maximum distance (number of edges)
      !    between any two nodes of the graph.
      !
      !    The eccentricity of a node is the maximum distance between that
      !    node and any other node of the graph.
      !
      !    A peripheral node is a node whose eccentricity equals the
      !    diameter of the graph.
      !
      !    A pseudo-peripheral node is an approximation to a peripheral node;
      !    it may be a peripheral node, but all we know is that we tried our
      !    best.
      !
      !    The routine is given a graph, and seeks pseudo-peripheral nodes,
      !    using a modified version of the scheme of Gibbs, Poole and
      !    Stockmeyer.  It determines such a node for the section subgraph
      !    specified by MASK and ROOT.
      !
      !    The routine also determines the level structure associated with
      !    the given pseudo-peripheral node; that is, how far each node
      !    is from the pseudo-peripheral node.  The level structure is
      !    returned as a list of nodes LS, and pointers to the beginning
      !    of the list of nodes that are at a distance of 0, 1, 2, ...,
      !    NODE_NUM-1 from the pseudo-peripheral node.
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license.
      !
      !  Modified:
      !
      !    28 October 2003
      !
      !  Author:
      !
      !    Original FORTRAN77 version by Alan George, Joseph Liu.
      !    FORTRAN90 version by John Burkardt.
      !
      !  Reference:
      !
      !    Alan George and Joseph Liu,
      !    Computer Solution of Large Sparse Positive Definite Systems,
      !    Prentice Hall, 1981.
      !
      !    Norman Gibbs, William Poole, Paul Stockmeyer,
      !    An Algorithm for Reducing the Bandwidth and Profile of a Sparse Matrix,
      !    SIAM Journal on Numerical Analysis,
      !    Volume 13, pages 236-250, 1976.
      !
      !    Norman Gibbs,
      !    Algorithm 509: A Hybrid Profile Reduction Algorithm,
      !    ACM Transactions on Mathematical Software,
      !    Volume 2, pages 378-387, 1976.
      !
      !  Parameters:
      !
      !    Input/output, integer ( kind = 4 ) ROOT.  On input, ROOT is a node in the
      !    the component of the graph for which a pseudo-peripheral node is
      !    sought.  On output, ROOT is the pseudo-peripheral node obtained.
      !
      !    Input, integer ( kind = 4 ) ADJ_NUM, the number of adjacency entries.
      !
      !    Input, integer ( kind = 4 ) ADJ_ROW(NODE_NUM+1).  Information about row I
      !    is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
      !
      !    Input, integer ( kind = 4 ) ADJ(ADJ_NUM), the adjacency structure.
      !    For each row, it contains the column indices of the nonzero entries.
      !
      !    Input, integer ( kind = 4 ) MASK(NODE_NUM), specifies a section subgraph.  
      !    Nodes for which MASK is zero are ignored by FNROOT.
      !
      !    Output, integer ( kind = 4 ) LEVEL_NUM, is the number of levels in the 
      !    level structure rooted at the node ROOT.
      !
      !    Output, integer ( kind = 4 ) LEVEL_ROW(NODE_NUM+1), 
      !    integer LEVEL(NODE_NUM), the level structure array pair containing the 
      !    level structure found.
      !
      !    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
      !
      implicit none

      integer ( kind = 4 ) adj_num
      integer ( kind = 4 ) node_num

      integer ( kind = 4 ) adj(adj_num)
      integer ( kind = 4 ) adj_row(node_num+1)
      integer ( kind = 4 ) iccsze
      integer ( kind = 4 ) j
      integer ( kind = 4 ) jstrt
      integer ( kind = 4 ) k
      integer ( kind = 4 ) kstop
      integer ( kind = 4 ) kstrt
      integer ( kind = 4 ) level(node_num)
      integer ( kind = 4 ) level_num
      integer ( kind = 4 ) level_num2
      integer ( kind = 4 ) level_row(node_num+1)
      integer ( kind = 4 ) mask(node_num)
      integer ( kind = 4 ) mindeg
      integer ( kind = 4 ) nabor
      integer ( kind = 4 ) ndeg
      integer ( kind = 4 ) node
      integer ( kind = 4 ) root
      !
      !  Determine the level structure rooted at ROOT.
      !
      call level_set ( root, adj_num, adj_row, adj, mask, level_num, &
           level_row, level, node_num )
      !
      !  Count the number of nodes in this level structure.
      !
      iccsze = level_row(level_num+1) - 1
      !
      !  Extreme case:
      !    A complete graph has a level set of only a single level.
      !    Every node is equally good (or bad).
      !
      if ( level_num == 1 ) then
         return
      end if
      !
      !  Extreme case:
      !    A "line graph" 0--0--0--0--0 has every node in its only level.
      !    By chance, we've stumbled on the ideal root.
      !
      if ( level_num == iccsze ) then
         return
      end if
      !
      !  Pick any node from the last level that has minimum degree
      !  as the starting point to generate a new level set.
      !
      do

         mindeg = iccsze

         jstrt = level_row(level_num)
         root = level(jstrt)

         if ( jstrt < iccsze ) then

            do j = jstrt, iccsze

               node = level(j)
               ndeg = 0
               kstrt = adj_row(node)
               kstop = adj_row(node+1)-1

               do k = kstrt, kstop
                  nabor = adj(k)
                  if ( 0 < mask(nabor) ) then
                     ndeg = ndeg+1
                  end if
               end do

               if ( ndeg < mindeg ) then
                  root = node
                  mindeg = ndeg
               end if

            end do

         end if
         !
         !  Generate the rooted level structure associated with this node.
         !
         call level_set ( root, adj_num, adj_row, adj, mask, level_num2, &
              level_row, level, node_num )
         !
         !  If the number of levels did not increase, accept the new ROOT.
         !
         if ( level_num2 <= level_num ) then
            exit
         end if

         level_num = level_num2
         !
         !  In the unlikely case that ROOT is one endpoint of a line graph,
         !  we can exit now.
         !
         if ( iccsze <= level_num ) then
            exit
         end if

      end do

      return
    end subroutine root_find

    subroutine level_set ( root, adj_num, adj_row, adj, mask, level_num, &
         level_row, level, node_num )

      !*****************************************************************************80
      !
      !! LEVEL_SET generates the connected level structure rooted at a given node.
      !
      !  Discussion:
      !
      !    Only nodes for which MASK is nonzero will be considered.
      !
      !    The root node chosen by the user is assigned level 1, and masked.
      !    All (unmasked) nodes reachable from a node in level 1 are
      !    assigned level 2 and masked.  The process continues until there
      !    are no unmasked nodes adjacent to any node in the current level.
      !    The number of levels may vary between 2 and NODE_NUM.
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license.
      !
      !  Modified:
      !
      !    28 October 2003
      !
      !  Author:
      !
      !    Original FORTRAN77 version by Alan George, Joseph Liu.
      !    FORTRAN90 version by John Burkardt.
      !
      !  Reference:
      !
      !    Alan George and Joseph Liu,
      !    Computer Solution of Large Sparse Positive Definite Systems,
      !    Prentice Hall, 1981.
      !
      !  Parameters:
      !
      !    Input, integer ( kind = 4 ) ROOT, the node at which the level structure
      !    is to be rooted.
      !
      !    Input, integer ( kind = 4 ) ADJ_NUM, the number of adjacency entries.
      !
      !    Input, integer ( kind = 4 ) ADJ_ROW(NODE_NUM+1).  Information about row I
      !    is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
      !
      !    Input, integer ( kind = 4 ) ADJ(ADJ_NUM), the adjacency structure.
      !    For each row, it contains the column indices of the nonzero entries.
      !
      !    Input/output, integer ( kind = 4 ) MASK(NODE_NUM).  On input, only nodes 
      !    with nonzero MASK are to be processed.  On output, those nodes which were
      !    included in the level set have MASK set to 1.
      !
      !    Output, integer ( kind = 4 ) LEVEL_NUM, the number of levels in the level
      !    structure.  ROOT is in level 1.  The neighbors of ROOT
      !    are in level 2, and so on.
      !
      !    Output, integer ( kind = 4 ) LEVEL_ROW(NODE_NUM+1), LEVEL(NODE_NUM), the 
      !    rooted level structure.
      !
      !    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
      !
      implicit none

      integer ( kind = 4 ) adj_num
      integer ( kind = 4 ) node_num

      integer ( kind = 4 ) adj(adj_num)
      integer ( kind = 4 ) adj_row(node_num+1)
      integer ( kind = 4 ) i
      integer ( kind = 4 ) iccsze
      integer ( kind = 4 ) j
      integer ( kind = 4 ) jstop
      integer ( kind = 4 ) jstrt
      integer ( kind = 4 ) lbegin
      integer ( kind = 4 ) level_num
      integer ( kind = 4 ) level_row(node_num+1)
      integer ( kind = 4 ) level(node_num)
      integer ( kind = 4 ) lvlend
      integer ( kind = 4 ) lvsize
      integer ( kind = 4 ) mask(node_num)
      integer ( kind = 4 ) nbr
      integer ( kind = 4 ) node
      integer ( kind = 4 ) root

      mask(root) = 0
      level(1) = root
      level_num = 0
      lvlend = 0
      iccsze = 1
      !
      !  LBEGIN is the pointer to the beginning of the current level, and
      !  LVLEND points to the end of this level.
      !
      do

         lbegin = lvlend + 1
         lvlend = iccsze
         level_num = level_num + 1
         level_row(level_num) = lbegin
         !
         !  Generate the next level by finding all the masked neighbors of nodes
         !  in the current level.
         !
         do i = lbegin, lvlend

            node = level(i)
            jstrt = adj_row(node)
            jstop = adj_row(node+1)-1

            do j = jstrt, jstop

               nbr = adj(j)

               if ( mask(nbr) /= 0 ) then
                  iccsze = iccsze + 1
                  level(iccsze) = nbr
                  mask(nbr) = 0
               end if

            end do

         end do
         !
         !  Compute the current level width (the number of nodes encountered.)
         !  If it is positive, generate the next level.
         !
         lvsize = iccsze - lvlend

         if ( lvsize <= 0 ) then
            exit
         end if

      end do

      level_row(level_num+1) = lvlend + 1
      !
      !  Reset MASK to 1 for the nodes in the level structure.
      !
      mask(level(1:iccsze)) = 1

      return
    end subroutine level_set

    subroutine rcm ( root, adj_num, adj_row, adj, mask, perm, iccsze, node_num )

      !*****************************************************************************80
      !
      !! RCM renumbers a connected component by the reverse Cuthill McKee algorithm.
      !
      !  Discussion:
      !
      !    The connected component is specified by a node ROOT and a mask.
      !    The numbering starts at the root node.
      !
      !    An outline of the algorithm is as follows:
      !
      !    X(1) = ROOT.
      !
      !    for ( I = 1 to N-1)
      !      Find all unlabeled neighbors of X(I),
      !      assign them the next available labels, in order of increasing degree.
      !
      !    When done, reverse the ordering.
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license.
      !
      !  Modified:
      !
      !    05 December 2008
      !
      !  Author:
      !
      !    Original FORTRAN77 version by Alan George, Joseph Liu.
      !    FORTRAN90 version by John Burkardt.
      !
      !  Reference:
      !
      !    Alan George and Joseph Liu,
      !    Computer Solution of Large Sparse Positive Definite Systems,
      !    Prentice Hall, 1981.
      !
      !  Parameters:
      !
      !    Input, integer ( kind = 4 ) ROOT, the node that defines the connected 
      !    component.  It is used as the starting point for the RCM ordering.
      !
      !    Input, integer ( kind = 4 ) ADJ_NUM, the number of adjacency entries.
      !
      !    Input, integer ( kind = 4 ) ADJ_ROW(NODE_NUM+1).  Information about row I
      !    is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
      !
      !    Input, integer ( kind = 4 ) ADJ(ADJ_NUM), the adjacency structure.
      !    For each row, it contains the column indices of the nonzero entries.
      !
      !    Input/output, integer ( kind = 4 ) MASK(NODE_NUM), a mask for the nodes. 
      !    Only those nodes with nonzero input mask values are considered by the
      !    routine.  The nodes numbered by RCM will have their mask values
      !    set to zero.
      !
      !    Output, integer ( kind = 4 ) PERM(NODE_NUM), the RCM ordering.
      !
      !    Output, integer ( kind = 4 ) ICCSZE, the size of the connected component
      !    that has been numbered.
      !
      !    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
      !
      !  Local parameters:
      !
      !    Workspace, integer DEG(NODE_NUM), a temporary vector used to hold
      !    the degree of the nodes in the section graph specified by mask and root.
      !
      implicit none

      integer ( kind = 4 ) adj_num
      integer ( kind = 4 ) node_num

      integer ( kind = 4 ) adj(adj_num)
      integer ( kind = 4 ) adj_row(node_num+1)
      integer ( kind = 4 ) deg(node_num)
      integer ( kind = 4 ) fnbr
      integer ( kind = 4 ) i
      integer ( kind = 4 ) iccsze
      integer ( kind = 4 ) j
      integer ( kind = 4 ) jstop
      integer ( kind = 4 ) jstrt
      integer ( kind = 4 ) k
      integer ( kind = 4 ) l
      integer ( kind = 4 ) lbegin
      integer ( kind = 4 ) lnbr
      integer ( kind = 4 ) lperm
      integer ( kind = 4 ) lvlend
      integer ( kind = 4 ) mask(node_num)
      integer ( kind = 4 ) nbr
      integer ( kind = 4 ) node
      integer ( kind = 4 ) perm(node_num)
      integer ( kind = 4 ) root
      !
      !  Find the degrees of the nodes in the component specified by MASK and ROOT.
      !
      call degree ( root, adj_num, adj_row, adj, mask, deg, iccsze, perm, node_num )

      mask(root) = 0

      if ( iccsze <= 1 ) then
         return
      end if

      lvlend = 0
      lnbr = 1
      !
      !  LBEGIN and LVLEND point to the beginning and
      !  the end of the current level respectively.
      !
      do while ( lvlend < lnbr )

         lbegin = lvlend + 1
         lvlend = lnbr

         do i = lbegin, lvlend
            !
            !  For each node in the current level...
            !
            node = perm(i)
            jstrt = adj_row(node)
            jstop = adj_row(node+1) - 1
            !
            !  Find the unnumbered neighbors of NODE.
            !
            !  FNBR and LNBR point to the first and last neighbors
            !  of the current node in PERM.
            !
            fnbr = lnbr + 1

            do j = jstrt, jstop

               nbr = adj(j)

               if ( mask(nbr) /= 0 ) then
                  lnbr = lnbr + 1
                  mask(nbr) = 0
                  perm(lnbr) = nbr
               end if

            end do
            !
            !  If no neighbors, skip to next node in this level.
            !
            if ( lnbr <= fnbr ) then
               cycle
            end if
            !
            !  Sort the neighbors of NODE in increasing order by degree.
            !  Linear insertion is used.
            !
            k = fnbr

            do while ( k < lnbr )

               l = k
               k = k + 1
               nbr = perm(k)

               do while ( fnbr < l )

                  lperm = perm(l)

                  if ( deg(lperm) <= deg(nbr) ) then
                     exit
                  end if

                  perm(l+1) = lperm
                  l = l-1

               end do

               perm(l+1) = nbr

            end do

         end do

      end do
      !
      !  We now have the Cuthill-McKee ordering.  Reverse it.
      !
      call i4vec_reverse ( iccsze, perm )

      return
    end subroutine rcm

    subroutine i4vec_reverse ( n, a )

      !*****************************************************************************80
      !
      !! I4VEC_REVERSE reverses the elements of an integer vector.
      !
      !  Example:
      !
      !    Input:
      !
      !      N = 5,
      !      A = ( 11, 12, 13, 14, 15 ).
      !
      !    Output:
      !
      !      A = ( 15, 14, 13, 12, 11 ).
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license.
      !
      !  Modified:
      !
      !    26 July 1999
      !
      !  Author:
      !
      !    John Burkardt
      !
      !  Parameters:
      !
      !    Input, integer ( kind = 4 ) N, the number of entries in the array.
      !
      !    Input/output, integer ( kind = 4 ) A(N), the array to be reversed.
      !
      implicit none

      integer ( kind = 4 ) n

      integer ( kind = 4 ) a(n)
      integer ( kind = 4 ) i

      do i = 1, n/2
         call i4_swap ( a(i), a(n+1-i) )
      end do

      return
    end subroutine i4vec_reverse
    subroutine i4_swap ( i, j )

      !*****************************************************************************80
      !
      !! I4_SWAP swaps two integer values.
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license.
      !
      !  Modified:
      !
      !    30 November 1998
      !
      !  Author:
      !
      !    John Burkardt
      !
      !  Parameters:
      !
      !    Input/output, integer ( kind = 4 ) I, J.  On output, the values of I and
      !    J have been interchanged.
      !
      implicit none

      integer ( kind = 4 ) i
      integer ( kind = 4 ) j
      integer ( kind = 4 ) k

      k = i
      i = j
      j = k

      return
    end subroutine i4_swap

    subroutine degree ( root, adj_num, adj_row, adj, mask, deg, iccsze, ls, &
         node_num )

      !*****************************************************************************80
      !
      !! DEGREE computes the degrees of the nodes in the connected component.
      !
      !  Discussion:
      !
      !    The connected component is specified by MASK and ROOT.
      !    Nodes for which MASK is zero are ignored.
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license.
      !
      !  Modified:
      !
      !   05 January 2003
      !
      !  Author:
      !
      !    Original FORTRAN77 version by Alan George, Joseph Liu.
      !    FORTRAN90 version by John Burkardt.
      !
      !  Reference:
      !
      !    Alan George and Joseph Liu,
      !    Computer Solution of Large Sparse Positive Definite Systems,
      !    Prentice Hall, 1981.
      !
      !  Parameters:
      !
      !    Input, integer ( kind = 4 ) ROOT, the node that defines the 
      !    connected component.
      !
      !    Input, integer ( kind = 4 ) ADJ_NUM, the number of adjacency entries.
      !
      !    Input, integer ( kind = 4 ) ADJ_ROW(NODE_NUM+1).  Information about row I 
      !    is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
      !
      !    Input, integer ( kind = 4 ) ADJ(ADJ_NUM), the adjacency structure.
      !    For each row, it contains the column indices of the nonzero entries.
      !
      !    Input, integer ( kind = 4 ) MASK(NODE_NUM), is nonzero for those nodes 
      !    which are to be considered.
      !
      !    Output, integer ( kind = 4 ) DEG(NODE_NUM), contains, for each  node in 
      !    the connected component, its degree.
      !
      !    Output, integer ( kind = 4 ) ICCSIZE, the number of nodes in the connected
      !    component.
      !
      !    Output, integer ( kind = 4 ) LS(NODE_NUM), stores in entries 1 through 
      !    ICCSIZE the nodes in the connected component, starting with ROOT, and 
      !    proceeding by levels.
      !
      !    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
      !
      implicit none

      integer ( kind = 4 ) adj_num
      integer ( kind = 4 ) node_num

      integer ( kind = 4 ) adj(adj_num)
      integer ( kind = 4 ) adj_row(node_num+1)
      integer ( kind = 4 ) deg(node_num)
      integer ( kind = 4 ) i
      integer ( kind = 4 ) iccsze
      integer ( kind = 4 ) ideg
      integer ( kind = 4 ) j
      integer ( kind = 4 ) jstop
      integer ( kind = 4 ) jstrt
      integer ( kind = 4 ) lbegin
      integer ( kind = 4 ) ls(node_num)
      integer ( kind = 4 ) lvlend
      integer ( kind = 4 ) lvsize
      integer ( kind = 4 ) mask(node_num)
      integer ( kind = 4 ) nbr
      integer ( kind = 4 ) node
      integer ( kind = 4 ) root
      !
      !  The sign of ADJ_ROW(I) is used to indicate if node I has been considered.
      !
      ls(1) = root
      adj_row(root) = -adj_row(root)
      lvlend = 0
      iccsze = 1
      !
      !  LBEGIN is the pointer to the beginning of the current level, and
      !  LVLEND points to the end of this level.
      !
      do

         lbegin = lvlend + 1
         lvlend = iccsze
         !
         !  Find the degrees of nodes in the current level,
         !  and at the same time, generate the next level.
         !
         do i = lbegin, lvlend

            node = ls(i)
            jstrt = -adj_row(node)
            jstop = abs ( adj_row(node+1) ) - 1
            ideg = 0

            do j = jstrt, jstop

               nbr = adj(j)

               if ( mask(nbr) /= 0 ) then

                  ideg = ideg + 1

                  if ( 0 <= adj_row(nbr) ) then
                     adj_row(nbr) = -adj_row(nbr)
                     iccsze = iccsze + 1
                     ls(iccsze) = nbr
                  end if

               end if

            end do

            deg(node) = ideg

         end do
         !
         !  Compute the current level width.
         !
         lvsize = iccsze - lvlend
         !
         !  If the current level width is nonzero, generate another level.
         !
         if ( lvsize == 0 ) then
            exit
         end if

      end do
      !
      !  Reset ADJ_ROW to its correct sign and return.
      !
      do i = 1, iccsze
         node = ls(i)
         adj_row(node) = -adj_row(node)
      end do

      return
    end subroutine degree

  end subroutine genrcm
  
  
  subroutine sort_row(nnz,ja,coeff)
      use Globals
      implicit none
      integer,           intent(in   ) :: nnz
      integer,           intent(inout) :: ja(nnz)
      real(kind=double), intent(inout) :: coeff(nnz)
      !local 
      integer :: i,j, indx,isgn, itemp
      real(kind=double) :: rtemp

      if (nnz.lt.2) return


      !  Initialize.
      i = 0
      indx = 0
      isgn = 0
      j = 0
      do 
         call global_heapsort(nnz, indx, i,j,isgn)
         if (indx .gt. 0 ) then
            ! SWAP ELEMENT 

            ! swap column indeces
            itemp = ja(i)
            ja(i) = ja(j)
            ja(j) = itemp

            ! swap real nnzero coeff
            rtemp    = coeff(i)
            coeff(i) = coeff(j)
            coeff(j) = rtemp
         else if ( indx .lt. 0) then
            ! COMPARE
            isgn = 0
            if ( ja(i) .lt.  ja(j) ) isgn = -1
            if ( ja(i) .gt.  ja(j) ) isgn = 1
         else if ( indx .eq. 0 ) then
            exit
         end if
      end do
    end subroutine sort_row


    subroutine plot2ps ( this, job, iounit )

      !*****************************************************************************80
      !
      !! PLTMTPS creates a PostScript plot of a sparse matrix.
      !
      !  Discussion:
      !
      !    This routine creates a 'PS' file for plotting the pattern of
      !    a sparse matrix stored in general sparse format. It can be used
      !    for inserting matrix plots in a text. The size of the plot can be
      !    7in x 7in or 5 in x 5in ..
      !
      !    1) Plots square as well as rectangular matrices.
      !    2) Does not writer a caption yet.
      !    3) No bounding box put in yet
      !
      !  Modified:
      !
      !    07 January 2004
      !
      !  Author:
      !
      !    Paul Frederickson
      !
      !  Parameters:
      !
      !    Input, integer ( kind = 4 ) NROW, the row dimension of the matrix.
      !
      !    Input, integer ( kind = 4 ) NCOL, the column dimension of the matrix.
      !
      !    Input, integer ( kind = 4 ) MODE, indicates the matrix storage mode:
      !    0, by rows;
      !    1, by columns.
      !
      ! ja     = column indices of nonzero elements when matrix is
      !         stored rowise. Row indices if stores column-wise.
      ! ia     = integer ( kind = 4 ) array of containing the pointers to the
      !         beginning of the columns in arrays a, ja.
      !
      ! title  = character*72 = title of matrix test ( character a*72 ).
      ! key    = character*8  = key of matrix
      ! type   = character*3  = type of matrix.
      !
      ! job, integer ( kind = 4 ). tells pltmt whether or not to reduce
      ! the plot.
      !           if enabled then the standard size of 7in will be
      !           replaced by a 5in plot.
      !          job = 0 : do not reduce
      !          job = 1 : reduce plot to 5 inches.
      !
      ! iounit = logical unit number where to write the matrix into.
      !
      use Globals
      implicit none
      
      class(spmat),     intent(in ) :: this
      integer,          intent(in ) :: job
      ! local
      real ( kind = 8 ) delta
      integer ( kind = 4 ) ii
      integer ( kind = 4 ) ilast
      integer ( kind = 4 ) iounit
      integer ( kind = 4 ) istart


      integer ( kind = 4 ) k
      character ( len = 8 ) key
      integer ( kind = 4 ) m
      integer ( kind = 4 ) maxdim
      integer ( kind = 4 ) mode
      integer ( kind = 4 ) n
      integer ( kind = 4 ) ncol
      integer ( kind = 4 ) nrow
      integer ( kind = 4 ) nnz
      character ( len = 3 ) type

      nrow = this%nrow
      ncol = this%ncol
      
      mode = 1
      type = this%storage_system

      if ( mode == 0 ) then
         n = nrow
      else
         n = ncol
      end if

      nnz = this%ia(n+1) - this%ia(1)
      maxdim = max ( nrow, ncol )
      m = 1 + maxdim
      !
      !  Keep this test as in old pltmt (for future changes).
      !
      if ( mod ( job, 10 ) == 1 ) then
         delta = 72.0D+00 * 5.0D+00 / ( 2.0D+00 + maxdim )
      else
         delta = 72.0D+00 * 7.0D+00 / (2.0D+00 + maxdim )
      end if

      write(iounit,*)'%!PS'
      write(iounit,*)' gsave 50 50 translate'
      write(iounit,*) delta, delta, ' scale'
      write(iounit,*) ' 0.25 setlinewidth'

      if ( mod ( job, 10 ) == 1 ) then
         write (iounit,*) ' 23 55 translate'
      else
         write (iounit,*) ' 2 35 translate'
      end if

      write(iounit,*) ' newpath'
      write(iounit,*) 0,0,' moveto'
      write(iounit,*) m,0,' lineto'
      write(iounit,*) m,m,' lineto'
      write(iounit,*) 0,m,' lineto'
      write(iounit,*) ' closepath stroke'
      write(iounit,*) ' 1 1 translate'
      write(iounit,*) ' 0.5 setlinewidth'
      write(iounit,*) ' /p {moveto 0 -.25 rmoveto '
      write(iounit,*) '            0  .50 rlineto stroke} def'
      !
      !  Plotting loop
      !
      do ii = 1, n

         istart = this%ia(ii)
         ilast  = this%ia(ii+1)-1

         if ( mode /= 0 ) then

            do k = istart, ilast
               write(iounit,*) ii-1, nrow-this%ja(k), ' p'
            end do

         else
            !             y = xnrow - real ( ii, kind = 8 )
            do k = istart, ilast
               !               x = real ( this%ja(k) - 1, kind = 8 )
               write(iounit,*) this%ja(k)-1, nrow-ii, ' p'
            end do

         end if

      end do

      write(iounit,*)' showpage grestore'
130   format('Dimension: ',i4,' x ',i4',  Nonzero elements: ',i5)
      return
    end subroutine plot2ps


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
    subroutine pattern_ATA(matrix_A, &
         lun_err,ncol_max,nterm_max, &
         matrix_ATA,&
         matrix_AT)
      implicit none
      ! vars
      class(spmat), intent(in   ) :: matrix_A
      integer,      intent(in   ) :: lun_err
      integer,      intent(in   ) :: ncol_max
      integer,      intent(in   ) :: nterm_max
      type(spmat), intent(inout) :: matrix_ATA
      type(spmat), optional, intent(in   ) :: matrix_AT
 
      ! local vars
      integer :: res
      logical :: rc
      integer :: n1,n2,ncol ,nrow, nterm,nmerged
      integer :: i,i1,i2
      integer :: j,j1,j2,k
      integer :: nwork, nel, nb
      integer,  allocatable :: ia(:)
      integer,  allocatable :: ja(:)
      integer,  allocatable :: array_b(:)
      integer , allocatable :: merged(:)
      integer , allocatable :: count(:)

      ncol = matrix_A%ncol
      nrow = matrix_A%nrow 
      
      allocate(&
              ia(ncol+1),&
              ja(ncol_max*ncol),&
              count(ncol),&
              merged(2*ncol_max),&
              array_b(2*ncol_max),&
              stat=res)
      if ( .not. present( matrix_AT) ) then
         ja=0
         count=0
         ia(1)=1
         nterm=0
         do i=1,matrix_A%nrow
            !
            ! Add column indeces of rank one matrix A(i,:) A(i,:)^T 
            !
            ia(i+1) = ia(i) + ncol_max
         end do


         do i=1,matrix_A%nrow
            !
            i1=matrix_A%ia(i)
            i2=matrix_A%ia(i+1)-1
            nb=i2-i1+1
            array_b(1:nb) = matrix_A%ja(i1:i2)

            do k=i1,i2
               j=matrix_A%ja(k)            
               call union_sorted(&
                    count(j), ja(ia(j):ia(j)+count(j)-1),&
                    nb, array_b(1:nb), &
                    nmerged, merged)


               nterm = nterm + nmerged - count(j)


               !
               ! update row(j) with merged list of columns if limits
               ! are satistied
               !
               if ( (nmerged > ncol_max) ) then
                  rc = IOerr(lun_err, err_val , 'pattern_ATA', &
                       'exceeded the max number of element for row =',j)
               else if ( nterm > nterm_max  ) then
                  rc = IOerr(lun_err, err_val , 'pattern_ATA', &
                       'exceeded the max number of nonzeros elements')
               else
                  count(j) = nmerged
                  ja(ia(j):ia(j) + count(j) -1 ) = merged(1:nmerged)
               end if
            end do
         end do

      else
!!$         allocate(&
!!$              rows_done(matrix_AT%nrow),&
!!$              stat=res)
!!$         rows_done=.False.
!!$         
!!$         do i = 1, matrix_AT%nrow
!!$            do j=matrix_AT%ia(i),matrix_AT%ia(i+1)-1
!!$               irow = matrix_AT%ja(j)
!!$               
!!$               !
!!$               i1=matrix_A%ia(i)
!!$               i2=matrix_A%ia(i+1)-1
!!$               nb=i2-i1+1
!!$               array_b(1:nb) = matrix_A%ja(i1:i2)
!!$               
!!$               do k=i1,i2
!!$               j=matrix_A%ja(k)            
!!$               call diff_inter_sorted(&
!!$                    count(j), ja(ia(j):ia(j)+count(j)-1),&
!!$                    nb, array_b(1:nb), &
!!$                    ndiff,  diff,&
!!$                    ninter, inter)
!!$               
!!$               do k=1,ndiff
!!$                  
!!$               
!!$
!!$
!!$               nterm = nterm + nmerged - count(j)
               
               

      end if

      
      

      call matrix_ATA%init(lun_err,&
           ncol, ncol, nterm,&
           'csr',&
           is_symmetric = .True. )

      
      matrix_ATA%ia(1)=1
      do i=1,ncol
         matrix_ATA%ia(i+1)=matrix_ATA%ia(i)+count(i)
         matrix_ATA%ja(&
              matrix_ATA%ia(i):matrix_ATA%ia(i+1)-1) = &
              ja( ia(i):ia(i)+count(i)-1)
      end do

    contains
      subroutine union_sorted(&
           na, array_a,&
           nb, array_b,&
           nmerged, merged)
        implicit none
        integer, intent(in   ) :: na
        integer, intent(in   ) :: array_a(na)
        integer, intent(in   ) :: nb
        integer, intent(in   ) :: array_b(nb)
        integer, intent(inout) :: nmerged
        integer, intent(inout) :: merged(na+nb)
        ! local 
        integer :: iii,jjj,counter

        !
        ! if one array has zero length copy the other into merged 
        !
        if ( na .eq. 0 ) then 
           nmerged = nb
           merged(1:nb)  = array_b(1:nb)
           return
        end if

        if ( nb .eq. 0 ) then 
           nmerged = na
           merged(1:na)  = array_a(1:na)
           return
        end if
        !
        ! moved along the two arrays
        !
        counter = 0
        iii = 1
        jjj = 1
        do while ( ( iii .le. na) .and. ( jjj .le. nb ) )
           !
           ! if the smallest is in the first array
           ! copy and advance in the first array
           ! Otherwise advance copy from the second and 
           ! andvance
           !
           if ( array_a(iii)  .lt. array_b(jjj) ) then
              counter=counter+1
              merged(counter) = array_a(iii)
              iii = iii+1
           else if ( array_a(iii)  .eq. array_b(jjj) ) then
              counter=counter+1
              merged(counter) = array_b(jjj)
              jjj = jjj+1
              iii = iii+1
           else
              counter=counter+1
              merged(counter) = array_b(jjj)
              jjj = jjj+1
           end if
           !write(*,*) 'i,j',counter, merged(counter) 
        end do
!!$        if (i .lt. na) then
!!$           merged(counter:counter + na -i + 1 ) = array_a(i:na) 
!!$           counter = counter + na - i  
!!$        end if
!!$        if (j .lt. nb) then
!!$           merged(counter:counter + nb -j + 1 ) = array_b(j:nb) 
!!$           counter = counter + nb - j  
!!$        end if

        do while( iii .le. na) 
           counter=counter+1
           !write(*,*) 'i,',counter, array_a(iii)
           merged(counter) = array_a(iii)
           iii = iii+1
        end do
        do while( jjj .le. nb) 
           counter=counter+1
           !write(*,*) 'j,',counter, array_b(jjj)
           merged(counter) = array_b(jjj)
           jjj = jjj+1
        end do
        nmerged=counter



      end subroutine union_sorted

      subroutine diff_inter_sorted(&
           na, array_a,&
           nb, array_b,&
           ndiff,diff,&
           ninter,inter)
        implicit none
        integer, intent(in   ) :: na
        integer, intent(in   ) :: array_a(na)
        integer, intent(in   ) :: nb
        integer, intent(in   ) :: array_b(nb)
        integer, intent(inout) :: ndiff
        integer, intent(inout) :: diff(na+nb)
        integer, intent(inout) :: ninter
        integer, intent(inout) :: inter(na+nb)
        ! local 
        integer :: iii,jjj,counter_inter,counter_diff,irow

        !
        ! if one array has zero length copy the other into merged 
        !
        if ( na .eq. 0 ) then 
           ndiff = nb
           diff(1:nb)  = array_b(1:nb)
           ninter = 0
           return
        end if

        if ( nb .eq. 0 ) then 
           ndiff = na
           diff(1:na)  = array_a(1:na)
           ninter = 0
           return
        end if
        !
        ! moved along the two arrays
        !
        counter_inter = 0
        counter_diff  = 0
        iii = 1
        jjj = 1
        do while ( ( iii .le. na) .and. ( jjj .le. nb ) )
           !
           ! if the smallest is in the first array
           ! copy and advance in the first array
           ! Otherwise advance copy from the second and 
           ! andvance
           !
           if ( array_a(iii)  .lt. array_b(jjj) ) then
              counter_diff=counter_diff+1
              merged(counter_diff) = array_a(iii)
              iii = iii+1
           else if ( array_a(iii)  .eq. array_b(jjj) ) then
              counter_inter=counter_diff+1
              merged(counter_inter) = array_b(jjj)
              jjj = jjj+1
              iii = iii+1
           else 
              counter_diff=counter_diff+1
              merged(counter_diff) = array_b(jjj)
              jjj = jjj+1
           end if
           !write(*,*) 'i,j',counter, merged(counter) 
        end do

        do while( iii .le. na) 
           counter_diff=counter_diff+1
           merged(counter_diff) = array_a(iii)
           iii = iii+1
        end do
        do while( jjj .le. nb) 
           counter_diff=counter_diff+1
           !write(*,*) 'j,',counter, array_b(jjj)
           merged(counter_diff) = array_b(jjj)
           jjj = jjj+1
        end do
        ndiff  = counter_diff
        ninter = counter_inter

      end subroutine diff_inter_sorted

    end subroutine pattern_ATA

    !>-------------------------------------------------------------
    !> Procedure to initialize and build sparse matrix 
    !> from matrix-matrix multiplication of sparse matrices
    !> (procedure public for type spmat)
    !> Initialize and compute all members of 
    !>            MDN = M x D x N
    !> usage:
    !>     call 'var'%mult_MDN(lun_err,&
    !>                         matrix_M, matrix_N,&
    !>                         max_nnz_row, max_nnz_total,&
    !>                         diagonal, &
    !>                         [mat_NT])
    !>
    !> where:
    !> \param[in] lun_err           -> integer. I/O unit for error message
    !> \param[in] matrix_M          -> type(spmat) First  (left) factor of 
    !>                                   multiplication MDM = M D N
    !> \param[in] matrix_M          -> type(spmat) Last  (right) factor
    !>                                   multiplication MDM = M D N
    !> \param[in] max_nnz_row       -> integer. Max number of non-zeros for row
    !> \param[in] max_nnz_total     -> integer. Max number of non-zeros of MDN
    !> \param[in] daigonal          -> real(dim=matrix_N%nrow). Coefficent of
    !>                                 diagonal matrix D
    !> \param[in] [optional] mat_NT -> type(spamt) 
    !<-----------------------------------------------------------
    subroutine mult_MDN( matrix_MDN, &
         lun_err, &
         matrix_M, matrix_N , &
         max_nnz_row, max_nnz_total, &
         ! optional
         diagonal,&
         mat_NT)
     implicit none
     ! vars
     class(spmat), intent(inout) ::  matrix_MDN
     integer,     intent(in   ) :: lun_err
     type(spmat), intent(in   ) ::matrix_M
     type(spmat), target, intent(in   ) :: matrix_N
     !integer,     intent(in   ) :: max_connection
     integer,     intent(in   ) :: max_nnz_row
     integer,     intent(in   ) :: max_nnz_total
     real(kind=double), optional, intent(in) :: diagonal(matrix_N%nrow)
     type(spmat), target, optional :: mat_NT

     ! local
     logical :: rc
     integer :: res,i,j
     integer :: n1, n, m, ntermax ,n2,info,nterm_MDN
     ! matrix used to 
     type(spmat), target :: loc_matrix_NT
     class(spmat), pointer :: matrix_N_csc
     
    
     integer, allocatable ::  WN1(:)
     integer, allocatable ::  JWN1(:)
     integer, allocatable ::  WN2(:)
     integer, allocatable ::  JWN2(:)
     real(double) , allocatable ::  WR1(:)

     integer , allocatable :: MDN_ia(:), MDN_ja(:)
     real(double) , allocatable :: MDN_coeff(:)
     
     !
     ! defined matrix_N in CSC (compressed store column) format 
     ! which is the transpose of matrix_N
     !
     if ( .not. present(mat_NT) ) then
        loc_matrix_NT = matrix_N
        call loc_matrix_NT%transpose(lun_err)
        matrix_N_csc => loc_matrix_NT
     else
        matrix_N_csc => mat_NT
     end if
         
     !
     ! define matrices dimensions and work arrays
     !
     n1= matrix_M%nrow
     m=matrix_N%nrow
     n2=matrix_N%ncol

     allocate( &
          WN1(m),&
          JWN1(m),&
          WN2(n2),&
          JWN2(n2),&
          WR1(m),&
          stat=res)
     if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'mult_MDN', &
          ' work arrays WN1,JWN1,WN2,JWN2 ')
     
     allocate( MDN_ia(n1+1),&
          MDN_ja(max_nnz_total),&
          MDN_coeff(max_nnz_total),&
          stat=res)
     if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'mult_MDN', &
          ' work MDN_ia, MDN_ja, MDN_coeff ')
     
     !
     ! assembly
     !
     if (present(diagonal)) then
        call Sparse_PMDM( matrix_M%nrow, matrix_M%ncol,matrix_N%ncol,&
             matrix_M%nterm,matrix_N%nterm,&
             max_nnz_total,nterm_MDN,&       
             matrix_M%ia, matrix_M%ja,&
             matrix_N%ia,matrix_N%ja,&
             matrix_N_csc%ia,matrix_N_csc%ja,&
             diagonal,&
             MDN_ia,MDN_ja,&
             WN1,WN2,JWN1,JWN2,&
             matrix_M%coeff,matrix_N_csc%coeff,&
             MDN_coeff,WR1,info)
     else
        call Sparse_PMM( matrix_M%nrow, matrix_M%ncol,matrix_N%ncol,&
             matrix_M%nterm,matrix_N%nterm,&
             max_nnz_total,nterm_MDN,&       
             matrix_M%ia, matrix_M%ja,&
             matrix_N%ia,matrix_N%ja,&
             matrix_N_csc%ia,matrix_N_csc%ja,&
             MDN_ia,MDN_ja,&
             WN1,WN2,JWN1,JWN2,&
             matrix_M%coeff,matrix_N_csc%coeff,&
             MDN_coeff,WR1,info)
     end if
     
     !
     ! init matrix and copy matrix
     !
     if ( matrix_MDN%is_initialized) call matrix_MDN%kill(lun_err)
     call matrix_MDN%init(lun_err, &
          n1,n2, nterm_MDN,&
          'csr',&
          is_symmetric=.True.)  
     matrix_MDN%ia(1:n1+1) = MDN_ia(1:n1+1)
     matrix_MDN%ja(1:nterm_MDN) = MDN_ja(1:nterm_MDN)
     matrix_MDN%coeff(1:nterm_MDN) = MDN_coeff(1:nterm_MDN)
     
     !
     ! fre memory
     !
     deallocate( &
          WN1,&
          JWN1,&
          WN2,&
          JWN2,&
          WR1,&
          stat=res)
     if(res .ne. 0) rc = IOerr(lun_err, err_dealloc, 'mult_MDN', &
          ' work arrays WN1,JWN1,WN2,JWN2 ')

     deallocate( &
          MDN_ia,&
          MDN_ja,&
          MDN_coeff,&
          stat=res)
     if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'mult_MDN', &
          ' work MDN_ia, MDN_ja, MDNcoeff ')

     if ( loc_matrix_NT%is_initialized) call loc_matrix_NT%kill(lun_err)
     matrix_N_csc => null()


   contains

     subroutine Sparse_PMDM(n1,m,n2,&
          nterm_F1,nterm_F2,&
          ntermax_P,nterm_P,&       
          ptr_F1,jcol_F1,&
          ptr_F2,jcol_F2,&
          ptc_F2,irow_F2,&
          diag,&
          ptr_P,jcol_P,&
          WN1,WN2,JWN1,JWN2,&
          mat_F1,mat_F2,mat_P,WR1,info)
       !-----------------------------------------------------------------------
       !  Esegue il prodotto matrice-matrice in modalit sparsa.
       !  La prima matrice "fattore" memorizzata in formato CSR
       !  La seconda matrice "fattore" memorizzata in formato CSC, ma  anche
       !  data la struttura del suo formato CSR
       !
       !  Variabili in Input:
       !
       !  F1                    : primo fattore n1*m
       !  ptr_F1,jcol_F1,mat_F1 : struttura e coefficienti di F1 per righe
       !  F2                    : secondo fattore m*n2
       !  ptr_F2,jcol_F2        : struttura di F2 per righe
       !  ptc_F2,irow_F2,mat_F2 : struttura e coefficienti di F2 per colonne
       !
       !  Variabili in Output
       !
       !  P                     : matrice prodotto
       !  ptr_P,jcol_P,mat_P    : struttura e coefficienti di P per righe
       !
       !  work arrays
       !
       !  WN1     : viene usato come indicatore di non-zero per la riga i di F1
       !  WN2     : viene usato come indicatore di colonna intercettata (pu�
       !            essere una variabile logica (testare efficienza)
       !  JWN1    : indici di colonna della riga i di F1
       !  WR1     : coefficienti della riga i di F1
       !  JWN2    : lista delle colonne intercettate
       !
       !  NB: JWN1,JWN2,WR1 possono essere dimensionati "piccoli"
       !
       !-----------------------------------------------------------------------
       use Globals
       implicit none
       !  Varibili in Input
       integer, intent(in) ::   n1,m,n2
       integer, intent(in) ::   nterm_F1,nterm_F2,ntermax_P
       integer, intent(in) ::   ptr_F1(n1+1),jcol_F1(nterm_F1)
       integer, intent(in) ::   ptr_F2(m+1) ,jcol_F2(nterm_F2)
       integer, intent(in) ::   ptc_F2(n2+1),irow_F2(nterm_F2)
       real(kind=double), intent(in) :: mat_F1(nterm_F1)
       real(kind=double), intent(in) ::  mat_F2(nterm_F2)
       real(kind=double), intent(in) :: diag(m)
       !  Variabili in Output
       integer , intent(inout) ::  info
       integer , intent(inout) :: nterm_P
       integer , intent(inout) :: ptr_P(n1+1),jcol_P(ntermax_P)
       real(kind=double), intent(inout) ::  mat_P(ntermax_P)
       !  Work arrays
       integer , intent(inout) :: WN1(m),JWN1(m),WN2(n2),JWN2(n2)
       real(kind=double), intent(inout) :: WR1(m)
       !  Variabili locali
       integer  i,j,k,ii,jj,kk,ISTRT,ISTOP,JSTRT,JSTOP
       integer  ind_P,ind_Psav,ind1,ind2
       !
       !  Inizializza a zero i work arrays
       !

       WN1=0
       WN2=0
       ptr_P(1) = 1
       ind_P = 0
       info = 0
       !
!  Ciclo sulle righe di F1
!
      do i = 1,n1
!
!  Esplora la riga i-esima
!
         ISTRT = ptr_F1(i)
         ISTOP = ptr_F1(i+1)-1
         ind1 = 0
         ind2 = 0
         do j = ISTRT,ISTOP
!
            jj = jcol_F1(j)
!           Carica in WN1,JWN1 e WR1 il termine j-esimo della riga
            ind1 = ind1 + 1
            JWN1(ind1) = jj        ! indice di colonna
            WN1(jj) = ind1         ! indicatore di non-zero che punta
                                   ! alla posizione del coef. in WR1
            WR1(ind1) = mat_F1(j)  ! coefficiente
!
!           Estrai le colonne di F2 "intercettate" dal termine jj
!           ---> Esplora la riga jj di F2
!
            JSTRT = ptr_F2(jj)
            JSTOP = ptr_F2(jj+1)-1
            do k = JSTRT,JSTOP
               ii = jcol_F2(k)  ! indice di colonna di F2
               if (WN2(ii) .eq. 0) then
!                 Se tale colonna non compare nella lista di colonne
!                 intercettate, aggiungila alla lista
                  WN2(ii) = 1
                  ind2 = ind2 + 1
                  JWN2(ind2) = ii
               endif
            enddo  ! ciclo su k
!
         enddo  ! ciclo su j
!
!  Creata la lista delle colonne intercettate, adesso esegui i prodotti
!  scalari
!
         ind_Psav = ind_P+1
         do j = 1,ind2
            ind_P = ind_P+1
            mat_P(ind_P) = 0.d0      ! Inizializza il coef. di mat. prod.
            jj = JWN2(j)             ! Colonna di F2 intercettata
            jcol_P(ind_P) = jj       ! Indice di colonna della mat. prod.
            JSTRT = ptc_F2(jj)       ! Puntatore alla testa della colonna
            JSTOP = ptc_F2(jj+1)-1   ! Puntatore alla fine della colonna
            do k = JSTRT,JSTOP
               kk = WN1(irow_F2(k))  ! Termine della riga di F1
               if (kk .gt. 0) then
                  mat_P(ind_P) = mat_P(ind_P) + mat_F2(k)*WR1(kk)*diag(irow_F2(k))
               endif
            enddo  ! ciclo su k
         enddo  ! ciclo su j
!
!  Controllo di smatriciamento
!
         if (ind_P .gt. ntermax_P) then
            info = i  ! Indice della riga prima dello smatriciamento
            return
         endif
!
!  Ordina la riga
!
         call sort_row(ind_P+1-ind_Psav,jcol_P(ind_Psav),mat_P(ind_Psav))

!
!  Assegna il punatore
!
         ptr_P(i+1) = ind_P+1
!
!  Fine del ciclo sull'i-esima riga di F1
!  Azzeramento (sparso) delle liste di indicatori WN1 e WN2
!
         do j = 1,ind1
            WN1(JWN1(j)) = 0
         enddo
         do j = 1,ind2
            WN2(JWN2(j)) = 0
         enddo
!
      enddo  ! ciclo su i
!
!  Calcola il numero di non-zeri della matrice prodotto
!
      nterm_P = ptr_P(n1+1)-1
!
      return
!
    end subroutine Sparse_PMDM

    subroutine Sparse_PMM(n1,m,n2,&
          nterm_F1,nterm_F2,&
          ntermax_P,nterm_P,&       
          ptr_F1,jcol_F1,&
          ptr_F2,jcol_F2,&
          ptc_F2,irow_F2,&
          ptr_P,jcol_P,&
          WN1,WN2,JWN1,JWN2,&
          mat_F1,mat_F2,mat_P,WR1,info)
       !-----------------------------------------------------------------------
       !  Esegue il prodotto matrice-matrice in modalit sparsa.
       !  La prima matrice "fattore" memorizzata in formato CSR
       !  La seconda matrice "fattore" memorizzata in formato CSC, ma  anche
       !  data la struttura del suo formato CSR
       !
       !  Variabili in Input:
       !
       !  F1                    : primo fattore n1*m
       !  ptr_F1,jcol_F1,mat_F1 : struttura e coefficienti di F1 per righe
       !  F2                    : secondo fattore m*n2
       !  ptr_F2,jcol_F2        : struttura di F2 per righe
       !  ptc_F2,irow_F2,mat_F2 : struttura e coefficienti di F2 per colonne
       !
       !  Variabili in Output
       !
       !  P                     : matrice prodotto
       !  ptr_P,jcol_P,mat_P    : struttura e coefficienti di P per righe
       !
       !  work arrays
       !
       !  WN1     : viene usato come indicatore di non-zero per la riga i di F1
       !  WN2     : viene usato come indicatore di colonna intercettata (pu�
       !            essere una variabile logica (testare efficienza)
       !  JWN1    : indici di colonna della riga i di F1
       !  WR1     : coefficienti della riga i di F1
       !  JWN2    : lista delle colonne intercettate
       !
       !  NB: JWN1,JWN2,WR1 possono essere dimensionati "piccoli"
       !
       !-----------------------------------------------------------------------
       use Globals
       implicit none
       !  Varibili in Input
       integer, intent(in) ::   n1,m,n2
       integer, intent(in) ::   nterm_F1,nterm_F2,ntermax_P
       integer, intent(in) ::   ptr_F1(n1+1),jcol_F1(nterm_F1)
       integer, intent(in) ::   ptr_F2(m+1) ,jcol_F2(nterm_F2)
       integer, intent(in) ::   ptc_F2(n2+1),irow_F2(nterm_F2)
       real(kind=double), intent(in) :: mat_F1(nterm_F1)
       real(kind=double), intent(in) ::  mat_F2(nterm_F2)
       !  Variabili in Output
       integer , intent(inout) ::  info
       integer , intent(inout) :: nterm_P
       integer , intent(inout) :: ptr_P(n1+1),jcol_P(ntermax_P)
       real(kind=double), intent(inout) ::  mat_P(ntermax_P)
       !  Work arrays
       integer , intent(inout) :: WN1(m),JWN1(m),WN2(n2),JWN2(n2)
       real(kind=double), intent(inout) :: WR1(m)
       !  Variabili locali
       integer  i,j,k,ii,jj,kk,ISTRT,ISTOP,JSTRT,JSTOP
       integer  ind_P,ind_Psav,ind1,ind2
       !
       !  Inizializza a zero i work arrays
       !

       WN1=0
       WN2=0
       ptr_P(1) = 1
       ind_P = 0
       info = 0
       !
!  Ciclo sulle righe di F1
!
      do i = 1,n1
!
!  Esplora la riga i-esima
!
         ISTRT = ptr_F1(i)
         ISTOP = ptr_F1(i+1)-1
         ind1 = 0
         ind2 = 0
         do j = ISTRT,ISTOP
!
            jj = jcol_F1(j)
!           Carica in WN1,JWN1 e WR1 il termine j-esimo della riga
            ind1 = ind1 + 1
            JWN1(ind1) = jj        ! indice di colonna
            WN1(jj) = ind1         ! indicatore di non-zero che punta
                                   ! alla posizione del coef. in WR1
            WR1(ind1) = mat_F1(j)  ! coefficiente
!
!           Estrai le colonne di F2 "intercettate" dal termine jj
!           ---> Esplora la riga jj di F2
!
            JSTRT = ptr_F2(jj)
            JSTOP = ptr_F2(jj+1)-1
            do k = JSTRT,JSTOP
               ii = jcol_F2(k)  ! indice di colonna di F2
               if (WN2(ii) .eq. 0) then
!                 Se tale colonna non compare nella lista di colonne
!                 intercettate, aggiungila alla lista
                  WN2(ii) = 1
                  ind2 = ind2 + 1
                  JWN2(ind2) = ii
               endif
            enddo  ! ciclo su k
!
         enddo  ! ciclo su j
!
!  Creata la lista delle colonne intercettate, adesso esegui i prodotti
!  scalari
!
         ind_Psav = ind_P+1
         do j = 1,ind2
            ind_P = ind_P+1
            mat_P(ind_P) = 0.d0      ! Inizializza il coef. di mat. prod.
            jj = JWN2(j)             ! Colonna di F2 intercettata
            jcol_P(ind_P) = jj       ! Indice di colonna della mat. prod.
            JSTRT = ptc_F2(jj)       ! Puntatore alla testa della colonna
            JSTOP = ptc_F2(jj+1)-1   ! Puntatore alla fine della colonna
            do k = JSTRT,JSTOP
               kk = WN1(irow_F2(k))  ! Termine della riga di F1
               if (kk .gt. 0) then
                  mat_P(ind_P) = mat_P(ind_P) + mat_F2(k)*WR1(kk)
               endif
            enddo  ! ciclo su k
         enddo  ! ciclo su j
!
!  Controllo di smatriciamento
!
         if (ind_P .gt. ntermax_P) then
            info = i  ! Indice della riga prima dello smatriciamento
            return
         endif
!
!  Ordina la riga
!
         call sort_row(ind_P+1-ind_Psav,jcol_P(ind_Psav),mat_P(ind_Psav))

!
!  Assegna il punatore
!
         ptr_P(i+1) = ind_P+1
!
!  Fine del ciclo sull'i-esima riga di F1
!  Azzeramento (sparso) delle liste di indicatori WN1 e WN2
!
         do j = 1,ind1
            WN1(JWN1(j)) = 0
         enddo
         do j = 1,ind2
            WN2(JWN2(j)) = 0
         enddo
!
      enddo  ! ciclo su i
!
!  Calcola il numero di non-zeri della matrice prodotto
!
      nterm_P = ptr_P(n1+1)-1
!
      return
!
    end subroutine Sparse_PMM
    
  end subroutine MULT_MDN


  !>-------------------------------------------------------------
    !> Procedure to initialize and build sparse matrix 
    !> from matrix-matrix multiplication of sparse matrices
    !> (procedure public for type spmat)
    !> Initialize and compute all members of 
    !>            MDN = M x D x N
    !> usage:
    !>     call 'var'%mult_MDN(lun_err,&
    !>                         matrix_M, matrix_N,&
    !>                         max_nnz_row, max_nnz_total,&
    !>                         diagonal, &
    !>                         [mat_NT])
    !>
    !> where:
    !> \param[in] lun_err           -> integer. I/O unit for error message
    !> \param[in] matrix_M          -> type(spmat) First  (left) factor of 
    !>                                   multiplication MDM = M D N
    !> \param[in] matrix_M          -> type(spmat) Last  (right) factor
    !>                                   multiplication MDM = M D N
    !> \param[in] max_nnz_row       -> integer. Max number of non-zeros for row
    !> \param[in] max_nnz_total     -> integer. Max number of non-zeros of MDN
    !<-----------------------------------------------------------
    subroutine init_product_spmat( matrix_MDN, &
         lun_err, &
         matrix_M, matrix_N , &
         max_nnz_row, max_nnz_total)
      implicit none
     ! vars
     class(spmat), intent(inout) ::  matrix_MDN
     integer,     intent(in   ) :: lun_err
     type(spmat), intent(in   ) ::matrix_M
     type(spmat), target, intent(in   ) :: matrix_N
     !integer,     intent(in   ) :: max_connection
     integer,     intent(in   ) :: max_nnz_row
     integer,     intent(in   ) :: max_nnz_total
   
     !  local 
     integer  i,j,k,ii,jj,kk,ISTRT,ISTOP,JSTRT,JSTOP
     integer  ind_P,ind_Psav,ind1,ind2
     logical :: rc
     integer :: res
     integer :: n1, n, m, ntermax ,n2,info,nterm_MDN
     integer, allocatable ::  WN1(:)
     integer, allocatable ::  JWN1(:)
     integer, allocatable ::  WN2(:)
     integer, allocatable ::  JWN2(:)
     integer , allocatable :: MDN_ia(:), MDN_ja(:)
       
     if ( matrix_M%ncol .ne. matrix_N%nrow) then
        rc = IOerr(lun_err, wrn_val, 'init_product_spmat', &
          ' matrix dimension mismatch')
        return
     end if
     
     n1= matrix_M%nrow
     m=matrix_N%nrow
     n2=matrix_N%ncol
     
     !
     ! define matrices dimensions and work arrays
     !
     n1= matrix_M%nrow
     m=matrix_N%nrow
     n2=matrix_N%ncol

     allocate( &
          WN1(m),&
          JWN1(m),&
          WN2(n2),&
          JWN2(n2),&
          stat=res)
     if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'mult_MDN', &
          ' work arrays WN1,JWN1,WN2,JWN2 ')
     
     allocate( &
          MDN_ia(n1+1),&
          MDN_ja(max_nnz_total),&
          stat=res)
     if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'mult_MDN', &
          ' work MDN_ia, MDN_ja, MDN_coeff ')
     
     !
     !  Inizializza a zero i work arrays
     !
     WN1=0
     WN2=0
     MDN_ia(1) = 1
     ind_P = 0
     info = 0
     !
     !  Ciclo sulle righe di F1
     !
     do i = 1,n1
        !
        !  Esplora la riga i-esima
        !
        ISTRT = matrix_M%ia(i)
        ISTOP = matrix_M%ia(i+1)-1
        ind1 = 0
        ind2 = 0
        do j = ISTRT,ISTOP
           !
           jj = matrix_M%ja(j)
           ! Carica in WN1,JWN1  il termine j-esimo della riga
           ind1 = ind1 + 1
           JWN1(ind1) = jj        ! indice di colonna
           WN1(jj) = ind1         ! indicatore di non-zero che punta
           !
           ! Estrai le colonne di F2 "intercettate" dal termine jj
           ! ---> Esplora la riga jj di F2
           !
           JSTRT = matrix_N%ia(jj)
           JSTOP = matrix_N%ia(jj+1)-1
           do k = JSTRT,JSTOP
              ii = matrix_N%ja(k)  ! indice di colonna di F2
              if (WN2(ii) .eq. 0) then
                 ! Se tale colonna non compare nella lista di colonne
                 ! intercettate, aggiungila alla lista
                 WN2(ii) = 1
                 ind2 = ind2 + 1
                 JWN2(ind2) = ii
              endif
           enddo  ! ciclo su k
           !
        enddo  ! ciclo su j
        !
        !
        !  Ordina la riga
        !
        call isort(ind_P+1-ind_Psav,MDN_ja(ind_Psav))

        !
        !  Assegna il punatore
        !
        MDN_ia(i+1) = ind_P+1
        !
        !  Fine del ciclo sull'i-esima riga di F1
        !  Azzeramento (sparso) delle liste di indicatori WN1 e WN2
        !
        do j = 1,ind1
           WN1(JWN1(j)) = 0
        enddo
        do j = 1,ind2
           WN2(JWN2(j)) = 0
        enddo
        !
     enddo  ! ciclo su i
     !
     !  Calcola il numero di non-zeri della matrice prodotto
     !
     nterm_MDN = MDN_ia(n1+1)-1

     !
     ! init matrix and copy matrix
     !
     call matrix_MDN%init(lun_err, &
          n1,n2, nterm_MDN,&
          'csr',&
          is_symmetric=.True.)  
     matrix_MDN%ia(1:n1+1) = MDN_ia(1:n1+1)
     matrix_MDN%ja(1:nterm_MDN) = MDN_ja(1:nterm_MDN)
     matrix_MDN%coeff(1:nterm_MDN) = zero

     !
     ! fre memory
     !
     deallocate( &
          WN1,&
          JWN1,&
          WN2,&
          JWN2,&
          stat=res)
     if(res .ne. 0) rc = IOerr(lun_err, err_dealloc, 'mult_MDN', &
          ' work arrays WN1,JWN1,WN2,JWN2 ')

     deallocate( &
          MDN_ia,&
          MDN_ja,&
          stat=res)
     if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'mult_MDN', &
          ' work MDN_ia, MDN_ja, MDNcoeff ')

   end subroutine init_product_spmat


  !>-------------------------------------------------------------
  !> Procedure to build integer arrays that redirect the non zeros
  !> term of matrix A into the non-zeros term of matrix B.
  !> The sparsity patterno of matrix B must contain the sparsity 
  !> pattern of A.
  !> (procedure public for type spmat)
  !>
  !> usage:
  !>     call 'var'%assembly_redirector(lun_err,&
  !>                         dst_matrix,&
  !>                         redirector)
  !>
  !> where:
  !> \param[in] lun_err            -> integer. I/O unit for error message
  !> \param[in] dst_matrix -> type(spmat) Matrix having a greater or equal
  !>                                  sparsity pattern of matrix var 
  !> \param[in] redirector         -> integer. (dim=this%nterm)
  !>                                  Integer pointer redirectiong nonzero term of 
  !>                                  sparse matric this those of dst_matrix
  !>                                   dest%ja(redicter(ind)) = this%ja(ind)
  !<-----------------------------------------------------------
  subroutine assembly_redirector( this, &
       lun_err, &
       dst_matrix, &
       redirector)
    use Globals
    implicit none
    class(spmat), intent(in   ) :: this
    integer,      intent(in   ) :: lun_err
    type(spmat),  intent(in   ) :: dst_matrix
    integer,      intent(inout) :: redirector(this%nterm)
    ! local
    integer :: irow, j, jcol, pos, start, finish,nel

    do irow=1,this%nrow
       !
       ! defined destionation column indees bounds
       !
       start  = dst_matrix%ia(irow)
       finish = dst_matrix%ia(irow+1)-1
       nel = finish-start+1
       pos = 1
       
       !
       ! search column indeces position
       ! 
       do j=this%ia(irow),this%ia(irow+1)-1
          jcol = this%ja(j)
          nel = finish-start+1 
          pos = ifind(nel , dst_matrix%ja(start:finish),jcol)
          redirector(j) = start + pos - 1
          start = start + pos
       end do
    end do
  end subroutine assembly_redirector

  !>-------------------------------------------------------------
  !> Procedure create the pointer to assembly matrix
  !> for assembly matrix ATA given pattern ATA (computed via MmultN)
  !> and matrix A. The pointer contains the indeces of the 
  !> element of ATA for a given element of A. (see assembly_ATA for usage)
  !> The pointer assembly can be used for assembly
  !> ATA with D diagonal matrix
  !>
  !> usage:
  !>     call 'var'%mult_MDN(lun_err,&
  !>                         matrix_M, matrix_N,&
  !>                         max_nnz_row, max_nnz_total,&
  !>                         diagonal, &
  !>                         [mat_NT])
  !>
  !> where:
  !> \param[in] lun_err           -> integer. I/O unit for error message
  !> \param[in] matrixATA         -> type(spmat) Assembled matrix ATA
  !> \param[in] matrixA           -> type(spmat) Matrix A
  !> \param[in] nassembler        -> integer. Number od element in 
  !> \param[in] assembler         -> integer. (dimension=nassembler)
  !>                                 Pointer allocated in this soubotuine
  !>                                 
  !<-----------------------------------------------------------
  subroutine assembly_assemblerATA( matrixATA, &
       lun_err, &
       matrix_A, &
       nassembler, assembler)
    implicit none
    class(spmat),         intent(in   ) :: matrixATA
    integer,              intent(in   ) :: lun_err
    type(spmat),          intent(in   ) :: matrix_A
    integer,              intent(inout) :: nassembler
    integer, allocatable, intent(inout) :: assembler(:)
    ! loc
    logical :: rc
    integer :: res
    integer :: i,j,k,start,finish
    integer :: irow, k1,k2, icol, jcol
    integer :: nel, nel_ATA
    nassembler=0
    do irow=1,matrix_A%nrow
       nel=matrix_A%ia(irow+1)-matrix_A%ia(irow) 
       nassembler =  nassembler +  nel **2
    end do

    if ( allocated(assembler) ) then
       deallocate (assembler,stat=res)
       if(res .ne. 0) rc = IOerr(lun_err, err_dealloc, 'assembly_assemblerATA', &
          ' assembler')
    end if
    allocate (assembler(nassembler),stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_dealloc, 'assembly_assemblerATA', &
         ' assembler')

    do irow=1,matrix_A%nrow
       nel = matrix_A%ia(irow+1)-matrix_A%ia(irow)
       do k1=1,nel
          icol = matrix_A%ja(k1)
          do k2=1,nel
             jcol = matrix_A%ja(k2)
             k=k+1
             nel_ATA = matrixATA%ia(icol+1)-matrixATA%ia(icol)
             start  = matrixATA%ia(icol)
             finish = matrixATA%ia(icol)
             assembler(k)=ifind(nel_ATA,matrixATA%ja(start:finish),jcol )
          end do
       end do
    end do



  end subroutine assembly_assemblerATA

  !>-------------------------------------------------------------
  !> Procedure to assembly matrix ATA given pattern ATA (computed via assembly assembler)
  !> and matrix A. The pointer contains the indeces of the 
  !> element of ATA for a given element of A. (see assembly_ATA for usage)
  !> The pointer assembly can be used for assembly
  !> ATA with D diagonal matrix
  !>
  !> usage:
  !>     call 'var'%mult_MDN(lun_err,&
  !>                         matrix_M, matrix_N,&
  !>                         max_nnz_row, max_nnz_total,&
  !>                         diagonal, &
  !>                         [mat_NT])
  !>
  !> where:
  !> \param[in] lun_err           -> integer. I/O unit for error message
  !> \param[in] matrixATA         -> type(spmat) Assembled matrix ATA
  !> \param[in] matrixA           -> type(spmat) Matrix A
  !> \param[in] nassembler        -> integer. Number od element in 
  !> \param[in] assembler         -> integer. (dimension=nassembler)
  !>                                 Pointer allocated in this soubotuine
  !>                                 
  !<-----------------------------------------------------------
  subroutine assembly_ATA( matrixATA, &
       lun_err, &
       matrix_A, &
       nassembler, assembler,&
       diagonal)
    implicit none
    class(spmat),               intent(inout) :: matrixATA
    integer,                    intent(in   ) :: lun_err
    type(spmat),                intent(in   ) :: matrix_A
    integer,                    intent(in   ) :: nassembler
    integer,                    intent(in   ) :: assembler(nassembler)
    real(kind=double), optional,intent(in   ) :: diagonal(matrix_A%nrow)
    
    ! loc
    logical :: rc
    integer :: res
    integer :: i,j,k,ind
    integer :: irow, k1,k2, icol, jcol
    integer :: nel, nel_ATA
    
    if  ( present(diagonal) ) then
       do irow=1,matrix_A%nrow
          nel=matrix_A%ia(irow+1)-matrix_A%ia(irow)
          do k1=1,nel
             icol = matrix_A%ja(k1)
             do k2=1,nel
                jcol = matrix_A%ja(k2)
                k=k+1
                ind = assembler(k)
                matrixATA%coeff(k) = matrixATA%coeff(k) + matrix_A%coeff(icol) * matrix_A%coeff(jcol)
             end do
          end do
       end do
    else
       do irow=1,matrix_A%nrow
          nel=matrix_A%ia(irow+1)-matrix_A%ia(irow)
          do k1=1,nel
             icol = matrix_A%ja(k1)
             do k2=1,nel
                jcol = matrix_A%ja(k2)
                k=k+1
                ind = assembler(k)
                matrixATA%coeff(k) = matrixATA%coeff(k) + &
                     matrix_A%coeff(icol) * diagonal(irow) * matrix_A%coeff(jcol)
             end do
          end do
       end do
    end if




  end subroutine assembly_ATA

    

    
    !>-------------------------------------------------------------
    !> Reading procedure.
    !> (public procedure for type spmat)
    !> Read content of a variable of type spamat
    !> 
    !> usage:
    !>     call 'var'%read(lun)
    !>
    !> where:
    !> \param[in] lun -> integer. I/O unit for error message output
    !>
    !<-------------------------------------------------------------
    subroutine get_row(this,lun_err,&
         irow,&
         nnz,coeff_row,ja_row,&
         idiag)
      use Globals
      implicit none
      class(spmat), intent(in)       :: this
      integer, intent(in)            :: lun_err
      integer, intent(in)            :: irow
      integer, intent(out)           :: nnz
      real(kind=double), intent(out) :: coeff_row(this%ncol)
      integer, optional, intent(out) :: ja_row(this%ncol)
      integer, optional, intent(out) :: idiag


      ! loc. var
      logical :: rc
      integer :: res, lun
      integer:: iterm,icol

      if ( this%storage_system .eq. 'ssr' ) then
         if(res .ne. 0) rc = IOerr(lun_err, err_inp, 'get_row', &
              ' ssr format not supported')
      end if

      ! get number of non-zero terms
      nnz   = this%ia(irow+1) - this%ia(irow) 

      ! copy column indeces  
      if ( present(ja_row)) &
           ja_row(1:nnz) = this%ja(this%ia(irow):this%ia(irow+1)-1)

      ! copy non-zero terms
      coeff_row(1:nnz) = this%coeff(this%ia(irow):this%ia(irow+1)-1)

      ! optional evaluation of  iadig
      if ( present(idiag) )  idiag = this%idiag(irow)

    end subroutine get_row

!!$  !>-------------------------------------------------------
!!$  !> Subroutine for incomplete cholesky factorization
!!$  !> of a symmetric matrix A. 
!!$  !> It computes the approximated factor U giving
!!$  !>                  A ~ U^t U
!!$  !> with U a upper triangualr matrix in csr format.
!!$  !> If the optional array "diagonal" is passed the
!!$  !> the output are
!!$  !> 1- an upper triangular matrix with one on the diagonal
!!$  !>            upper    = diag^{-1}(U) U
!!$  !> 2 -the diagonal array containg the diagonal of U 
!!$  !             diagonal = diag(U) 
!!$  !>-------------------------------------------------------
!!$  subroutine new_incomplete_cholesky (this,&
!!$       lun_err,n_fillin,tol_fillin,info,&
!!$       upper,&
!!$       ! optional arguements
!!$       diagonal,aux_saved)
!!$    use Globals
!!$    use Scratch
!!$    implicit none
!!$    class(spmat), target,        intent(in   ) :: this
!!$    integer,                     intent(in   ) :: lun_err
!!$    integer,                     intent(in   ) :: n_fillin
!!$    real(kind=double),           intent(in   ) :: tol_fillin
!!$    integer,                     intent(inout) :: info
!!$    type(spmat),                 intent(inout) :: upper
!!$    real(kind=double), optional, intent(inout) :: diagonal(this%nrow)
!!$    type(scrt), optional, target,intent(inout) :: aux_saved
!!$
!!$    !local
!!$    logical :: rc
!!$    integer :: res
!!$    integer :: nterm_max, nterm ,nterm_out,niaux,nraux
!!$    integer :: irow,nrow
!!$    integer :: ibegin, iend, idiag, nnz
!!$    integer :: ju0
!!$    integer,           pointer :: ja_row(:), ja_centered(:)
!!$    integer,           pointer :: jlu(:)
!!$    integer,           pointer :: jw1(:),jw2(:),jw3(:)
!!$    real(kind=double), pointer :: alu(:)
!!$    real(kind=double), pointer :: coeff_row(:)
!!$    real(kind=double), pointer :: coeff_centered(:)
!!$    type(scrt), pointer :: aux
!!$    type(scrt), target  :: aux_loc 
!!$
!!$    !
!!$    ! checks
!!$    !
!!$    if ( this%nrow .ne. this%ncol  ) &
!!$         rc = IOerr(lun_err, err_inp, 'incomplete_cholesky ', &
!!$         'not square matrix passed')
!!$
!!$    if ( .not. (this%is_symmetric) ) &
!!$         rc = IOerr(lun_err, err_inp, 'incomplete_cholesky ', &
!!$         'not symmmetric matrix passed')
!!$    
!!$    !
!!$    ! dimensions
!!$    !
!!$    nrow      = this%nrow
!!$    nterm_max = this%nterm + (n_fillin+1) * nrow
!!$
!!$    !----------------------------------------------
!!$    ! scratch setting
!!$    !
!!$    niaux = nterm_max + 5 * nrow + 1
!!$    nraux = nterm_max + 3 * nrow
!!$    !
!!$    ! work arrays for both storage system
!!$    !
!!$    if ( present( aux_saved ) ) then
!!$       if ( .not. aux%check(niaux,nraux) ) then
!!$          call aux_loc%init(lun_err, niaux,nraux)
!!$          aux => aux_loc
!!$       else
!!$          aux => aux_saved
!!$       end if
!!$    else
!!$       call aux_loc%init(lun_err, niaux,nraux)
!!$       aux => aux_loc
!!$    end if
!!$    !
!!$    ! redirect real arrays
!!$    !
!!$    iend = 0
!!$    call aux%range(nterm_max,ibegin,iend)
!!$    alu             => aux%raux(ibegin:iend)
!!$    call aux%range(nterm_max,ibegin,iend)
!!$    coeff_row       => aux%raux(ibegin:iend)
!!$    call aux%range(nrow,ibegin,iend)
!!$    coeff_centered  => aux%raux(ibegin:iend)
!!$
!!$    !
!!$    ! redirect integer arrays
!!$    !
!!$    iend = 0
!!$    call aux%range(nterm_max,ibegin,iend)
!!$    jlu => aux%iaux(ibegin:iend)
!!$    call aux%range(nrow,ibegin,iend)
!!$    ja_row => aux%iaux(ibegin:iend)
!!$    call aux%range(nrow,ibegin,iend)
!!$    ja_centered  => aux%iaux(ibegin:iend)
!!$    call aux%range(nrow,ibegin,iend)
!!$    jw1 => aux%iaux(ibegin:iend)
!!$    call aux%range(nrow,ibegin,iend)
!!$    jw2 => aux%iaux(ibegin:iend)
!!$    call aux%range(nrow+1,ibegin,iend)
!!$    jw3 => aux%iaux(ibegin:iend)
!!$    !----------------------------------------------
!!$    
!!$    !----------------------------------------------
!!$    ! cycle all rows and generate factorization
!!$    !
!!$    info = 0  
!!$    do irow = 1, this%nrow
!!$       ! get irow^th row of the matrix
!!$       call this%get_row(lun_err, irow,nnz, coeff_row, ja_row,idiag)
!!$       
!!$       !
!!$       ! convert ja_row and coeff_row in the format required
!!$       ! by general_csr_incomplete_cheolesky algorithm
!!$       !
!!$       ! set portion of work array
!!$       ibegin = irow - ( nnz - idiag ) 
!!$       iend   = irow + ( nnz - idiag )
!!$       ! create work arrays for incomplete cholesky
!!$       ja_centered   (ibegin:iend) = ja_row(:)
!!$       coeff_centered(ibegin:iend) = coeff_row(:)
!!$       
!!$       !
!!$       ! compute factorization
!!$       !
!!$       call general_csr_incomplete_cheolesky(&
!!$            lun_err,info,&
!!$            this%nrow,&
!!$            irow,&
!!$            nnz,&
!!$            ja_centered,&
!!$            coeff_centered,&
!!$            n_fillin,tol_fillin,&
!!$            nterm_max,&
!!$            alu,&
!!$            jlu,&
!!$            ju0,&
!!$            jw1,&
!!$            jw2,&
!!$            jw3)
!!$
!!$       !
!!$       ! info in case of error
!!$       ! 
!!$       if (info .ne. 0 ) then
!!$          select case (info) 
!!$          case (-1)
!!$             rc = IOerr(lun_err, wrn_out, 'build_ic_fillin', &
!!$                  ' The elimination process has generated'//&
!!$                  ' a row in U whose length is .gt.  n)')
!!$          case (-2)  
!!$             rc = IOerr(lun_err, wrn_out, 'build_ic_fillin', &
!!$                  'The matrix U overflows the array U')
!!$          case (-3)  
!!$             rc = IOerr(lun_err, wrn_out, 'build_ic_fillin', &
!!$                  ' Illegal value for n_fillin')
!!$          case(-4)   
!!$             rc = IOerr(lun_err, wrn_out, 'build_ic_fillin', &
!!$                  'zero row encountered')
!!$          case(-5) 
!!$             rc = IOerr(lun_err, wrn_out, 'build_ic_fillin', &
!!$                  'Non Positive Diagonal Element')
!!$          end select
!!$          !
!!$          ! return for failure
!!$          !
!!$          return 
!!$       end if
!!$    end do
!!$    !----------------------------------------------
!!$    ! successfull factorization
!!$    !----------------------------------------------
!!$
!!$
!!$    !----------------------------------------------
!!$    ! convert mssr to csr : 
!!$    ! upper will store the matrix D(U)^{-1} U
!!$    call mss2scaled_upper(lun_err,&
!!$         nrow,&
!!$         nterm_max,&            
!!$         jlu,&
!!$         alu,&
!!$         upper)
!!$    !----------------------------------------------
!!$
!!$    !-----------------------------------------------
!!$    ! Prepare matrix 
!!$    !    U 
!!$    ! or matrices 
!!$    !    D(U), D(U)^{-1} U
!!$    !
!!$    ! Array alu(1:nrow) contains the square 
!!$    ! of the diagonal term  of the upper 
!!$    ! factor of the choelsky factorization
!!$    !
!!$    if ( present (diagonal) ) then
!!$       !
!!$       ! case for D=diag(U),  \tilde{U}=diag(U)^{-1} U 
!!$       !
!!$       diagonal(1:nrow) = sqrt( alu(1:nrow) )
!!$    else
!!$       !
!!$       ! case for U factorization
!!$       !
!!$       ! use alu as scratch array to compute roots
!!$       alu(1:nrow) = sqrt( alu(1:nrow) )
!!$       call upper%DxM( lun_err,alu(1:nrow) )
!!$       upper%unitary_diag = .false.
!!$    end if
!!$    
!!$    !
!!$    ! free memory
!!$    !
!!$    aux => null()
!!$    if ( aux_loc%is_initialized) call aux_loc%kill(lun_err)
!!$    
!!$  end subroutine new_incomplete_cholesky

  !*******************************

! Code to compute the pattern of lower triangle of C= A^T * A, without forming C
! C held in CSC format

! A is m x n, held in CSC format.
! There is no checking of the user's data.

   subroutine ata_pattern(job,m,n,ptr,row,rowptr,col,limit,nzc,ptrc,rowc,mask,&
              flag,stat,count)

   integer, intent(in) :: job ! 
      ! 1 to true if only lower triangular part of C is required
      ! 2 to true if only lower triangular part of C is required but rowc
      !   allocated to be large enough for upper and lower parts
   integer, intent(in) :: m  ! rows in A
   integer, intent(in) :: n  ! cols in A
   integer, intent(in) ::  ptr(n+1) ! column pointers for A 
   integer, intent(in) ::  row(:) ! row indices of A
   integer, intent(out) :: rowptr(m+1) ! row pointers for A
   integer, intent(out) :: col(:) ! col indices of A 
   integer, intent(in) :: limit(2) ! used to terminate the routine if C
     ! found to be too dense. 
     ! If a col of C has more than limit(1) entries, computation terminates
     ! If C has more than limit(2) entries, computation terminates.
     ! If negative, no limit
   integer, intent(out) :: nzc ! no. of entries in (lower triangle of) C
   integer, intent(out) ::  ptrc(n+1) ! column pointers for C 
   integer, intent(out), allocatable :: rowc(:)
   integer, intent(out) :: mask(m) ! work array
   integer, intent(out) :: flag ! error flag
     ! -1 : allocation error
     ! -2 : deallocation error
     ! -3 : a column of C has more than limit(1) entries (if lower = .true. this
     !      is the lower triangular part of C only)
     ! -4 : number of entries in C exceeds limit(2). 
   integer, intent(out) :: stat ! stat parameter
   integer, intent(out), optional :: count(n) ! only used if job = 1,2. In this
     ! case, accumulates counts in rows so that we can get a count for
     ! the number of entries in a col of C (upper and lower parts)

   integer :: i,i1,i2,ii,j,j1,j2,jcol,jj,k
   integer :: lenrowc,lenc

   flag = 0
   stat = 0

   ! Set rowptr, col to hold A in CSR format
   mask(1:m) = 0

   do j = 1, n
      do i = ptr(j), ptr(j+1)-1
         k = row(i)
         mask(k) = mask(k) + 1
      end do
   end do

   rowptr(1) = 1
   do i = 1, m
      rowptr(i+1) = rowptr(i) + mask(i)
      mask(i) = rowptr(i)
   end do

   do j = 1, n
      do i = ptr(j), ptr(j+1)-1
         k = row(i)
         jj = mask(k)
         col(jj) = j
         mask(k) = jj + 1
      end do
   end do

   ! First compute how many entries in (lower triangular part of) C.
   mask(1:m) = 0
   if (present(count)) count(1:n) = 0

   ! loop over cols of C
   nzc = 0
   ptrc(1) = 1
      do j = 1,n
         i1 = ptr(j)
         i2 = ptr(j+1) - 1
         do i = i1,i2
            ii = row(i)
            j1 = rowptr(ii)
            j2 = rowptr(ii+1) - 1
            do jj = j1,j2
               jcol = col(jj)
                if (jcol.le.j .and. mask(jcol).lt.j) then
                   mask(jcol) = j
                   nzc = nzc + 1
                   ! by symmetry, there is an entry in col jcol, row j
                   if (jcol.ne.j) count(jcol) = count(jcol) + 1
                end if
             end do
          end do
          ptrC(j+1) = nzc + 1
          lenc = ptrC(j+1) - ptrC(j)
          lenc = lenc + count(j) ! this allows us to include upper triangular 
                                 ! part (by symmetry, equal to entries in row j)
          if (lenc.gt.limit(1)) then
             flag = -3
             return
          end if

         ! check size of C (upper and lower triangular parts)
         if (nzc.gt.(limit(2)-n)/2) then
           flag = -4
           return
         end if

       end do

   ! reinitialise mask
   mask(1:m) = 0

   lenrowc = nzc
   if (job.eq.2) lenrowc = 2*nzc - n

   ! allocate space for the pattern of C
     allocate(rowc(lenrowc),stat=stat)
     if (stat.ne.0) then
       flag = -1
       return
     end if

   !call matmulC_pattern(n, ptr, row, rowptr, col, ptrc, rowc, count)

   end subroutine ata_pattern

!****************************************************************************
   subroutine matmulC_pattern(n,ia,ja,ib,jb,ic,jc,mask)

    ! compute the product of two sparse matrices A=(ia,ja,a) and
    ! B = (ib,jb,b) in CSR format. Output C in CSR format.
    ! Here we compute only the upper triangle of C
    ! since we are interested in A = B^T so that C = B^T * B is symmetric.
    ! Of course, upper triangle of C in CSR format = lower triangle in CSC
    ! format
    !!! note: as we only want the case where C is square, we have 
    ! taken row size of A = col size of B

       ! n: row size of A = column size of B
    integer, intent (in) :: n
       ! ia: row pointer of A
       ! ja: column indices of A
       ! ib: row pointer of B
       ! jb: column indices of B
    integer, intent(in) :: ia(*),ib(*),ja(*),jb(*)
       ! ic: row pointer of C
       ! jc: column indices of C
    integer, intent(out) :: ic(*),jc(*)
       ! a: entry values of A
       ! b: entry values of B
       ! c: entry values of C
       ! working array
    integer, intent(out) :: mask(*) ! length n

    integer :: nz,i,j,k,icol,icol_add,neigh

    ! initialise the mask array which is an array that
    ! has non-zero value if the
    ! column index already exist, in which case
    ! the value is the index of that column

    ! code computes one row of C at a time
    mask(1:n) = 0
    nz = 0
    ic(1) = 1
    do i = 1,n
       do j = ia(i),ia(i+1)-1
          neigh = ja(j)
          do k = ib(neigh),ib(neigh+1)-1
             icol_add = jb(k)

             ! cycle if entry is in lower triangle
             if (icol_add.lt.i) cycle

             icol = mask(icol_add)
             if (icol.eq.0) then
                nz = nz + 1
                jc(nz) = icol_add
                mask(icol_add) = nz     ! add mask
             end if
          end do
       end do
       ! done with row i, so set mask to zero again
       do j = ic(i),nz
          mask(jc(j)) = 0
       end do
       ic(i+1) = nz + 1
    end do

  end subroutine matmulC_pattern 

  !>-------------------------------------------------------------
  !> Procedure to compute select and/or permute
  !> rows from csr matrix. Given a sparse matrix "A"
  !> the new matrix will contain only the rows
  !> listed in "row_list".
  !> Example
  !> A=( 1 0 0 0 3 ) P = (3,1) new=( 5 4 0 6 0 )
  !>   ( 4 0 5 5 6 )               ( 1 0 0 0 3 )
  !>   ( 5 4 0 6 0 )
  !> 
  !> (public procedure for type spmat)
  !> 
  !> usage:
  !>     call 'var'%select_permute_row(vec_in,vec_out,[info])
  !>
  !> where 
  !> \param[in   ] newrows         -> integer. Number of rows
  !>                                  to be selected or permuted
  !>                                  vector to be multiplied
  !> \param[in   ] row_list        -> integer(dimesion=nrows)
  !>                                  List of the rows of the original
  !>                                  matrix
  !> \param[inout] new             -> type(spmat) New Matrix
  !> \param[in   ] (optional) lun_err -> integer. Info number
  !>                                  in case of error   
  !<-------------------------------------------------------------
  subroutine select_permute_rows(this, new_nrow,row_list ,new,lun_err)
    use Globals
    implicit none
    class(spmat),  intent(in   ) :: this
    integer,       intent(in   ) :: new_nrow
    integer,       intent(in   ) :: row_list(new_nrow)
    class(spmat),  intent(inout) :: new
    integer,  optional, intent(in   ) :: lun_err
    
    !local
    logical :: rc
    integer :: i,j,res,lun
    integer :: new_start,new_finish, start, finish,new_nterm
    integer, allocatable :: count(:)

    if (present(lun_err))then
       lun=lun_err
    else
       lun=0
    end if
    
    ! count new non-zeros
    new_nterm=0
    do i=1,new_nrow
       new_nterm=new_nterm+this%ia(row_list(i)+1)-this%ia(row_list(i))
    end do

    !
    ! init new matrix
    !
    if ( new%is_initialized) call new%kill(lun)
    call new%init(lun,new_nrow,this%ncol,new_nterm,'csr')

    ! new ia
    new%ia(1)=1
    do i=1,new_nrow
       new%ia(i+1)=new%ia(i)+this%ia(row_list(i)+1)-this%ia(row_list(i))
    end do

    ! new ja and coeff
    do i=1,new_nrow
       new_start=new%ia(i)
       new_finish=new%ia(i+1)-1

       start=this%ia(row_list(i))
       finish=this%ia(row_list(i)+1)-1
    
       new%ja   (new_start:new_finish)=this%ja(start:finish)
       new%coeff(new_start:new_finish)=this%coeff(start:finish)
    end do
    
  end subroutine SELECT_PERMUTE_ROWS

  !>-------------------------------------------------------------
  !> Procedure to compute select and/or permute
  !> rows from csr matrix. Given a sparse matrix "A"
  !> the new matrix will contain only the rows
  !> listed in "row_list".
  !> Example
  !> A=( 1 0 0 0 3 ) P = (3,1) new=( 5 4 0 6 0 )
  !>   ( 4 0 5 5 6 )               ( 1 0 0 0 3 )
  !>   ( 5 4 0 6 0 )
  !> 
  !> (public procedure for type spmat)
  !> 
  !> usage:
  !>     call 'var'%select_permute_row(vec_in,vec_out,[info])
  !>
  !> where 
  !> \param[in   ] newrows         -> integer. Number of rows
  !>                                  to be selected or permuted
  !>                                  vector to be multiplied
  !> \param[in   ] row_list        -> integer(dimesion=nrows)
  !>                                  List of the rows of the original
  !>                                  matrix
  !> \param[inout] new             -> type(spmat) New Matrix
  !> \param[in   ] (optional) lun_err -> integer. Info number
  !>                                  in case of error   
  !<-------------------------------------------------------------
  subroutine remove_non_zeros(this, nterm2remove,irows,jcols,lun_err)
    use Globals
    implicit none
    class(spmat),    intent(in   ) :: this
    integer,         intent(in   ) :: nterm2remove
    integer,         intent(in   ) :: irows(nterm2remove)
    class(spmat),    intent(inout) :: jcols(nterm2remove)
    integer,optional,intent(in   ) :: lun_err
    
    !local
    logical :: rc
    integer :: i,j,res,lun
    integer :: new_start,new_finish, start, finish,new_nterm
    integer, allocatable :: count(:)

    if (present(lun_err))then
       lun=lun_err
    else
       lun=0
    end if
    
    
  end subroutine remove_non_zeros


  !>-------------------------------------------------------------
  !> Procedure to compute select and/or permute
  !> rows from csr matrix. Given a sparse matrix "A"
  !> the new matrix will contain only the rows
  !> listed in "row_list".
  !> Example
  !> A=( 1 0 0 0 3 ) column list = (3,1) new=( 0 1 )
  !>   ( 4 0 5 5 6 )                         ( 5 4 )
  !>   ( 5 4 0 6 0 )                         ( 0 5 )
  !> 
  !> (public procedure for type spmat)
  !> 
  !> usage:
  !>     call 'var'%select_permute_columns(vec_in,vec_out,[info])
  !>
  !> where 
  !> \param[in   ] newrows         -> integer. Number of rows
  !>                                  to be selected or permuted
  !>                                  vector to be multiplied
  !> \param[in   ] column_list      -> integer(dimesion=nrows)
  !>                                  List of the columns of the original
  !>                                  matrixxs
  !> \param[inout] new             -> type(spmat) New Matrix
  !> \param[in   ] (optional) lun_err -> integer. Info number
  !>                                  in case of error   
  !<-------------------------------------------------------------
  subroutine select_permute_columns(this, new_ncol,column_list ,new,lun_err)
    use Globals
    implicit none
    class(spmat),  intent(in   ) :: this
    integer,       intent(in   ) :: new_ncol
    integer,       intent(in   ) :: column_list(new_ncol)
    class(spmat),  intent(inout) :: new
    integer,  optional, intent(in   ) :: lun_err
    
    !local
    logical :: rc
    integer :: i,j,res,lun
    integer :: new_start,new_finish, start, finish,new_nterm
    integer, allocatable :: count(:)
    type(spmat) :: temp
    
    if (present(lun_err))then
       lun=lun_err
    else
       lun=0
    end if

    select type (this)
    type is (spmat)
       temp=this
    end select

    ! in csr format is easier to
    ! 1-transpose the matrix
    ! 2-select/permute by rows
    ! 3-transpose back
    call temp%transpose(lun)
    call select_permute_rows(temp,new_ncol,column_list,new)    
    call new%transpose(lun)

    ! free memory
    call temp%kill(lun)
    
  end subroutine SELECT_PERMUTE_COLUMNS
  
  

  !>-------------------------------------------------------------
  !> Procedure to sum to compute
  !>
  !> N=alpha*M+N
  !>
  !> Supports different format for M
  !> (like axpy in blas library)
  !> (public procedure for type spmat)
  !> 
  !> usage:
  !>     call 'var'%(alpha,matrix_M,matrix_N,[Mcoeff2Ncoeff])
  !>
  !> where 
  !> \param[in   ] alpha           -> real. Scaling factor for M
  !> \param[in   ] matrix_M        -> type(identity,diagmat,spmat)
  !>                                  Matrix to add to N
  !> \param[in   ] (optional) info -> integer. dimension=matrix_M%nterm
  !>                                  Used only for M=spmat with different
  !>                                  non-zeros pattern
  !>                                  Map from non-zeros of M into
  !>                                  non-zeros of N  
  !<-------------------------------------------------------------
  subroutine aMpN(matrix_N,alpha, matrix_M, Mcoeff2Ncoeff,same_pattern)
    use Globals
    implicit none
    class(spmat),      intent(inout) :: matrix_N
    real(kind=double), intent(in   ) :: alpha
    class(abs_linop),  intent(in   ) :: matrix_M
    integer, optional, intent(in   ) :: Mcoeff2Ncoeff(:)
    logical, optional, intent(in   ) :: same_pattern
    !local
    logical :: rc
    integer :: i,j
    integer :: irow, jcol, pos, start, finish,nel,ind


    if ( ( matrix_M%nrow .ne. matrix_N%nrow) .or. &
         ( matrix_M%ncol .ne. matrix_N%ncol)) then
       call matrix_M%info(0)
       call matrix_N%info(0)
       rc = IOerr(0, err_val, 'aMpN', &
            ' dimensions mismatch in aMpN ')
    end if

    
    

    select type (matrix_M)
    type is (spmat)
       if ( matrix_N%nterm < matrix_M%nterm) then
          call matrix_M%info(0)
          call matrix_N%info(0)
          rc = IOerr(0, err_val, 'aMpN', &
               ' matrix_N has less non-zeros than matrix_M ')
       end if
       !
       ! use riderector for nonzeros of M into N
       !
       if ( present(Mcoeff2Ncoeff)) then
          do i=1,matrix_M%nrow
             do j=matrix_M%ia(i),matrix_M%ia(i+1)-1
                matrix_N%coeff(Mcoeff2Ncoeff(j)) = &
                     matrix_N%coeff(Mcoeff2Ncoeff(j)) +&
                     alpha*matrix_M%coeff(j)
             end do
          end do
       else
          if ( present(same_pattern)) then
             !
             ! sum nonzeros
             !
             if (same_pattern) then
                matrix_N%coeff=matrix_N%coeff+alpha*matrix_M%coeff
             end if
          else
             !
             ! slow but working procedure
             !
             do irow=1,matrix_M%nrow
                !
                ! defined destionation column indees bounds
                !
                start  = matrix_N%ia(irow)
                finish = matrix_N%ia(irow+1)-1
                nel = finish-start+1
                pos = 1

                !
                ! search column indeces position
                !
                !write(*,*) irow
                !write(*,*) matrix_M%ja(matrix_M%ia(irow):matrix_M%ia(irow+1)-1)
                !write(*,*) matrix_N%ja(matrix_N%ia(irow):matrix_N%ia(irow+1)-1)
                do j=matrix_M%ia(irow),matrix_M%ia(irow+1)-1
                   jcol = matrix_M%ja(j)
                   nel = finish-start+1 
                   pos = ifind(nel , matrix_N%ja(start:finish),jcol)
                   if (pos .eq. 0) rc = IOerr(0, err_val, 'aMpN', &
                        ' dimensions mismatch in aMpN ')
                   ind = start + pos - 1
                   matrix_N%coeff(ind) = &
                        matrix_N%coeff(ind) +&
                        alpha*matrix_M%coeff(j)
                   start = start + pos
                end do
             end do
          end if
       end if
    type is (diagmat)
       do i=1,matrix_N%nrow
          j=matrix_N%idiag(i)
          matrix_N%coeff(j)=matrix_N%coeff(j)+alpha*Matrix_M%diagonal(i)
       end do
    type is (scalmat)
       do i=1,matrix_N%nrow
          j=matrix_N%idiag(i)          
          matrix_N%coeff(j)=matrix_N%coeff(j)+alpha*Matrix_M%scalar
       end do
    end select

  end subroutine aMpN


  !>-------------------------------------------------------------
  !> Procedure to form an explicit sparse matrix from a block
  !> matrix. 
  !>
  !> (public procedure for type spmat)
  !> 
  !> usage:
  !>     call 'var'%(info,lun_err,block_matrix)
  !>
  !> where 
  
  !> \param[in   ] info         -> integer. Flag for errors 
  !> \param[in   ] lun_err      -> integer. Logical I/O unit
  !> \param[in   ] block_matrix -> type(block_linop)
  !>                               Block linear operator
  !>                               made by components of type
  !>                               spmat, diagmat, scalmat
  !<-------------------------------------------------------------
  recursive subroutine form_pair_linop(this,info,lun_err, matrix)
    use Globals
    class(spmat),      intent(inout) :: this
    integer,           intent(inout) :: info
    integer,           intent(in   ) :: lun_err
    class(pair_linop),  intent(in   ) :: matrix
    !local
    integer :: i

    if ( matrix%type .eq. 'LC') then
       call form_linear_combination(this,info,lun_err, matrix)
    else if ( matrix%type .eq. 'LP') then
       call form_composed(this,info,lun_err, matrix)
    else
       
    end if
    call this%is_like(matrix)
    
  contains
    !>-------------------------------------------------------------
    !> Procedure to form an explicit sparse matrix from a linear
    !> combination of sparse operator.
    !>
    !> M = a1*M1+a2*M2+..+an*Mn
    !>
    !> M   = sparse matrix
    !> RHS = linear_combination_linop
    !>
    !> (public procedure for type spmat)
    !> 
    !> usage:
    !>     call 'var'%(info,lun_err,linear_combination)
    !>
    !> where:
    !> \param[in   ] info       -> integer. Flag for errors 
    !> \param[in   ] lun_err    -> integer. Logical I/O unit
    !> \param[in   ]
    !>       linear_combination -> type(linear_combination_linop)
    !>                             Block linear operator
    !>                             made by components of type
    !>                             spmat, diagmat, scalmat
    !<-------------------------------------------------------------
    recursive subroutine form_linear_combination(this,&
         info,lun_err, linear_combination,is_diagonal)
      use Globals
      class(spmat),      intent(inout) :: this
      integer,           intent(inout) :: info
      integer,           intent(in   ) :: lun_err
      type(pair_linop),  intent(in   ) :: linear_combination
      logical, optional, intent(inout) :: is_diagonal
      ! local
      logical :: rc
      integer :: res
      integer :: i,j,k,m,irow,icol,istart,iend,ilinop
      integer :: ntemp,nrow,ncol,nnz,nterm,ndiagonals
      integer, allocatable :: nnzinrow(:)
      integer, allocatable :: ia_temp(:)
      integer, allocatable :: ja_temp(:)
      real(kind=double), allocatable :: coeff_temp(:)
      real(kind=double) :: alpha,leading_alpha
      type(spmat), target, allocatable :: temp_spmats(:)
      type(array_linop), allocatable :: local_list(:)
      character(len=256) :: msg
      type(diagmat) :: diagonal_matrix
      logical :: debug

      debug=.False.

      if (debug) write(*,*) '  '
      if (debug) write(*,*) 'BEGIN LIN COMB COMPOSED'
      if (debug) write(*,*) '  '
      if (debug) call linear_combination%info(6)

      info=0
      nrow = linear_combination%nrow
      ncol = linear_combination%ncol
      leading_alpha=linear_combination%alpha
      
      !
      ! count the number of sparse matrices to be formed
      ! and redirect point to temporary sparse matrices
      !
      allocate( nnzinrow(nrow),stat=res)
      if(res .ne. 0) rc = IOerr(lun_err,&
           err_alloc,  'block2sparse', &
           ' work array nnzinrow ')
      local_list=linear_combination%linop_list
      ntemp=0
      ndiagonals=0
      do ilinop=1, linear_combination%nlinops
         alpha  = linear_combination%alphas(ilinop)
         if (debug) call linear_combination%linop_list(ilinop)%linop%info(6)
         if (alpha .ne. zero) then         
            select type (mat=>linear_combination%linop_list(ilinop)%linop)
            type is (spmat)
               ! nothing to do
            type is (block_linop)
               ntemp=ntemp+1
            type is (new_linop)
               ntemp=ntemp+1
            type is (pair_linop)
               ntemp=ntemp+1
            type is (diagmat)
               ndiagonals=ndiagonals+1
            type is (scalmat)
               ndiagonals=ndiagonals+1
            class default
               rc = IOerr(lun_err,&
                    wrn_val,  'form_linear_combination', &
                    ' type not supported in linear operator number: ',ilinop)
               return
            end select
         end if
      end do

      !
      ! shorthand to build diagonal matrices
      !
      if (ndiagonals .eq. linear_combination%nlinops ) then
         if (present(is_diagonal)) is_diagonal=.True.
         !
         ! init diagonal matrix
         !
         call diagonal_matrix%init(lun_err, nrow,ncol)

         !
         ! form diagonal matrix
         !
         do ilinop=1, linear_combination%nlinops
            alpha  = linear_combination%alphas(ilinop)
            if (alpha .ne. zero) then
               select type( mat => linear_combination%linop_list(ilinop)%linop )
               type is (diagmat)
                  do irow = 1, nterm
                     diagonal_matrix%diagonal(irow) = diagonal_matrix%diagonal(irow)+alpha*mat%diagonal(irow)
                  end do
                  nnzinrow=nnzinrow+1
               type is(scalmat)
                  do irow = 1, nterm
                     diagonal_matrix%diagonal(irow) = diagonal_matrix%diagonal(irow)+alpha*mat%scalar
                  end do
               end select
            end if
         end do
         diagonal_matrix%diagonal=leading_alpha*diagonal_matrix%diagonal

         !
         ! convert to sparse matrix
         !
         call this%form_diag(info,lun_err,diagonal_matrix)

         !
         ! free memory 
         !
         call diagonal_matrix%kill(lun_err)
         return
      end if
      
      
      if (ntemp>0) then
         allocate(temp_spmats(ntemp),stat=res)
         if(res .ne. 0) rc = IOerr(lun_err,&
              err_alloc,  'block2sparse', &
              ' work temp_spmats ')
      end if
      ntemp=0
      do ilinop=1, linear_combination%nlinops
         alpha  = linear_combination%alphas(ilinop)
         if (alpha .ne. zero) then         
            select type (mat=>linear_combination%linop_list(ilinop)%linop)
            type is (block_linop)
               ntemp=ntemp+1
               call temp_spmats(ntemp)%form_block(info,lun_err,mat)
               local_list(ilinop)%linop => temp_spmats(ntemp)
            type is (new_linop)
               ntemp=ntemp+1
               call temp_spmats(ntemp)%form_new_linop(info,lun_err,mat)
               local_list(ilinop)%linop => temp_spmats(ntemp)
            type is (pair_linop)
               ntemp=ntemp+1
               call temp_spmats(ntemp)%form_pair_linop(info,lun_err,mat)
               local_list(ilinop)%linop => temp_spmats(ntemp)
            end select
         end if
      end do



      !
      ! count nonzeros in each row
      !
      ntemp    = 0
      nnzinrow = 0
      do ilinop = 1 , linear_combination%nlinops
         alpha  = linear_combination%alphas(ilinop)
         if (alpha .ne. zero) then         
            select type (mat=>local_list(ilinop)%linop)
            type is (spmat)
               do k=1,mat%nrow
                  nnzinrow(k)=&
                       nnzinrow(k) + &
                       mat%ia(k+1)-mat%ia(k)
               end do
            type is(diagmat)
               nnzinrow(1:nrow) = nnzinrow(1:nrow)+1
            type is(scalmat)
               if ( abs( mat%scalar)> small ) then 
                  nnzinrow(1:nrow) = nnzinrow(1:nrow)+1
               end if
            class default
               info=1
               write (msg,*) 'Linear operator ', ilinop , &
                    ' in block (', irow,',', icol, ') not supported'
               if(res .ne. 0) rc = IOerr(lun_err,&
                    err_val, 'form', &
                    etb(msg))
            end select
         end if
      end do

      allocate(ia_temp(nrow+1),stat=res)
      if (res.ne.0) rc=IOerr(lun_err, err_alloc, 'kill_spmat', &
           ' work array ia_temp',res)

      ia_temp(1)=1
      do i=1,nrow
         !write(6,*) i,nnzinrow(i) 
         ia_temp(i+1)=ia_temp(i)+nnzinrow(i)
      end do

      !write(6,*) sum(nnzinrow) , ia_temp(nrow+1)-1


      allocate(ja_temp(ia_temp(nrow+1)-1),stat=res)
      if (res.ne.0) rc=IOerr(lun_err, err_alloc, 'kill_spmat', &
           ' work array ja_temp',res)

      allocate(coeff_temp(ia_temp(nrow+1)-1),stat=res)
      if (res.ne.0) rc=IOerr(lun_err, err_alloc, 'kill_spmat', &
           ' work array ia_temp',res)


      ! We simply concatenate the non-zeros indeces
      ! and coefficients. We use a csr-like storage.
      ! ia_temp contains the info to get lower and
      ! upper bound limit for each row.
      ! We use nnz_row to keep track of the added
      ! non-zeros.
      ! Example:
      !
      ! ja        coeff
      ! ( 1 2 4 ) ( 1.0 2.0 4.0 )
      ! ( 3 4 )   ( 3.0 4.0) 
      ! becomes
      ! ( 1   2   4   3   4  ) 
      ! ( 1.0 2.0 4.0 3.0 4.0)
      !
      ntemp=0
      nnzinrow=0
      do ilinop = 1 , linear_combination%nlinops
         alpha  = linear_combination%alphas(ilinop)      
         if (alpha .ne. zero) then         
            select type (mat=>local_list(ilinop)%linop)
            type is (spmat)
               do irow = 1, nrow
                  istart = ia_temp(irow)+nnzinrow(irow)-1
                  !write(6,*) ilinop, irow, istart,ia_temp(irow),nnzinrow(irow)
                  i=0
                  do j = mat%ia(irow), mat%ia(irow+1)-1
                     i = i + 1 
                     ja_temp   (istart+i) = mat%ja(j)
                     coeff_temp(istart+i) = alpha*mat%coeff(j)
                  end do
                  nnzinrow(irow)=&
                       nnzinrow(irow) + &
                       mat%ia(irow+1)-mat%ia(irow)
               end do
            type is (diagmat)
               do irow = 1, nrow
                  istart = nnzinrow(irow)
                  ja_temp   (istart+irow) = irow
                  coeff_temp(istart+irow) = alpha*mat%diagonal(irow)
               end do
               nnzinrow=nnzinrow+1
            type is(scalmat)
               do irow = 1, nrow
                  istart = nnzinrow(irow)
                  ja_temp   (istart+irow) = irow
                  coeff_temp(istart+irow) = alpha*mat%scalar
               end do
               nnzinrow=nnzinrow+1
            !type is (zeromat)
               ! nothing to do in this case
            class default
               info=1
               write (msg,*) 'Linear operator ', ilinop , &
                    ' in block (', irow,',', icol, ') not supported'
               if(res .ne. 0) rc = IOerr(lun_err,&
                    err_val, 'form', &
                    etb(msg))
            end select
         end if
      end do

      nterm=0
      do irow=1,nrow
         !
         ! sort so that column index and the corresponding nnz zeros
         ! are in increasing order. We use sort_matrix_line since
         ! it work also in the case
         ! EXAMPLE
         ! ja        coeff
         ! ( 1 2 4 ) ( 1.0 2.0 4.0 )
         ! ( 3 4 )   ( 3.0 4.0) 
         !
         ! More efficient algorithm should exploit that
         ! we are merging and sorting arrays that are already sorted
         !
         !write(6,*) irow, &
         !     ia_temp(irow), ia_temp(irow+1)-1,&
         !     nnzinrow(irow),ia_temp(irow)+nnzinrow(irow)-1
         istart= ia_temp(irow)
         iend  = ia_temp(irow+1) - 1
         nnz   = iend - istart + 1
         !write(6,*) ja_temp(istart:iend)
         !write(6,*) coeff_temp(istart:iend)

         call sort_matrix_line(nnz,&
              ja_temp(istart:iend),&
              coeff_temp(istart:iend))
         !write(6,*) ja_temp(istart:iend)
         !write(6,*) coeff_temp(istart:iend)

         !
         ! we merge the duplicate columns indeces
         ! and we sum the corresponding non-zeros
         !
         m = ia_temp(irow)
         icol = ja_temp(istart)
         ! count actual non-zeros
         k=1
         ! starting from the second element in the list
         do j = istart+1,iend
            if ( ja_temp(j) .eq. icol ) then
               !  same index thus 
               !  sum the corresponding non-zeros
               coeff_temp(m) = coeff_temp(m) + coeff_temp(j)
            else
               !
               ! new index. 
               ! 
               icol = ja_temp(j)
               m=m+1
               k=k+1
               ja_temp(m)    = ja_temp(j)
               coeff_temp(m) = coeff_temp(j)
            end if
         end do
         nnzinrow(irow)=k
         !write(6,*) ja_temp(istart:istart+k-1)
         !write(6,*) coeff_temp(istart:istart+k-1)
         nterm=nterm+k
      end do

      !
      ! initialize the matrix
      !
      call this%init(lun_err,nrow,ncol,nterm,'csr')
      call this%is_like(linear_combination)
      
      !
      ! assign ia, ja,coeff
      !
      this%ia(1)=1
      do irow=1,nrow
         this%ia(irow+1)=this%ia(irow)+nnzinrow(irow)
         k=this%ia(irow)-1
         do j=ia_temp(irow),ia_temp(irow)+nnzinrow(irow)-1
            k=k+1
            this%ja(k)    =    ja_temp(j)
            this%coeff(k) = coeff_temp(j)
         end do
      end do
      this%coeff=leading_alpha*this%coeff
      

      
      !
      ! free memory
      !
      if (ntemp>0) then
         do i=1,ntemp
            call temp_spmats(i)%kill(lun_err)
         end do
         deallocate(temp_spmats,stat=res)
         if(res .ne. 0) rc = IOerr(lun_err,&
              err_dealloc,  'form_linear_combination', &
              ' work temp_spmats ')
      end if

      deallocate(nnzinrow,ia_temp,ja_temp,coeff_temp,stat=res)
      if(res .ne. 0) rc = IOerr(lun_err,&
           err_dealloc,  'form_linear_combination', &
           ' work nnzinrow,ia_temp,ja_temp,coeff_temp ')

      if (debug) call this%info(6)
      if (debug) write(*,*) '  '
      if (debug) write(*,*) 'END LIN COMB COMPOSED'
      if (debug) write(*,*) '  '

    end subroutine form_linear_combination


    
    subroutine sort_matrix_line(nnz,ja,coeff)
      use Globals
      implicit none
      integer,           intent(in   ) :: nnz
      integer,           intent(inout) :: ja(nnz)
      real(kind=double), intent(inout) :: coeff(nnz)
      !local 
      integer :: i,j, indx,isgn, itemp
      real(kind=double) :: rtemp

      !  Initialize.
      i = 0
      indx = 0
      isgn = 0
      j = 0
      do 
         call global_heapsort(nnz, indx, i,j,isgn)
         if (indx .gt. 0 ) then
            ! SWAP ELEMENT 

            ! swap column indeces
            itemp = ja(i)
            ja(i) = ja(j)
            ja(j) = itemp

            ! swap real nnzero coeff
            rtemp    = coeff(i)
            coeff(i) = coeff(j)
            coeff(j) = rtemp
         else if ( indx .lt. 0) then
            ! COMPARE
            isgn = 0
            if ( ja(i) .lt.  ja(j) ) isgn = -1
            if ( ja(i) .gt.  ja(j) ) isgn = 1
         else if ( indx .eq. 0 ) then
            exit
         end if
      end do
    end subroutine sort_matrix_line

    !
    ! This subroutine should be faster that using heapsort
    ! It has not been tested
    !
!!$    subroutine merge_sorted(nlinops,max_nnz,nnnz,nnz,max_nnz_out,merged)
!!$      implicit none
!!$      integer, intent(in   ) :: nlinops
!!$      integer, intent(in   ) :: nnnz(nlinops)
!!$      integer, intent(in   ) :: nnz(max_nnz,nlinops)
!!$      integer, intent(in   ) :: max_nnz
!!$      integer, intent(inout) :: max_nnz_out
!!$      integer, intent(inout) :: merged(max_nnz_out)
!!$
!!$      begin(:)=1
!!$      mininline=1
!!$
!!$      imin=max(nnz(1,1),nnz(1,2))
!!$      completed=.False.
!!$      do while (completed) 
!!$         do i=1,nlinops
!!$            if ( begin(i) .le. nnz(i) ) then               
!!$               if ( nnz(begin(i),i) .lt. imin ) then
!!$                  imin=nnz(begin(i),i)
!!$                  begin(i)=begin(i)+1
!!$               else if ( nnz(begin(i),1) .eq. imin ) then
!!$                  begin(i)=begin(i)+1
!!$               end if
!!$            end if
!!$         end do
!!$         m=m+1
!!$         merged(m)=imin
!!$
!!$         !
!!$         ! compelted if we visit all the nnz 
!!$         !
!!$         completed = .True.
!!$         do i=1,nlinops
!!$            completed = compelted .and. ( begin(i) .eq. nnnz(i) )
!!$         end do
!!$      end do
!!$    end subroutine merge_sorted
    
    !>-------------------------------------------------------------
    !> Procedure to form an explicit sparse matrix from a block
    !> matrix. 
    !>
    !> (public procedure for type spmat)
    !> 
    !> usage:
    !>     call 'var'%(info,lun_err,block_matrix)
    !>
    !> where 
    
    !> \param[in   ] info         -> integer. Flag for errors 
    !> \param[in   ] lun_err      -> integer. Logical I/O unit
    !> \param[in   ] block_matrix -> type(block_linop)
    !>                               Block linear operator
    !>                               made by components of type
    !>                               spmat, diagmat, scalmat,
    !<-------------------------------------------------------------
    recursive subroutine form_composed(this,info,lun_err, composed_matrix)
    use Globals
    use SimpleMatrix
    type(spmat),        intent(inout) :: this
    integer,             intent(inout) :: info
    integer,             intent(in   ) :: lun_err
    type(pair_linop), intent(in   ) :: composed_matrix

    !local
    logical :: rc
    integer :: res
    integer :: i,j,k
    integer :: nrow_loc, ncol_loc,shift,nterm,dim,nlinops
    integer :: istart,iend,ilinop,ntemp,last_sparse
    integer :: irow, icol,jcol, pos, start, finish,nel,ind
    type(spmat)  :: temp_spmat,temp_spmat2, temp_right
    type(diagmat) :: diagonal_right
    real(kind=double) :: alpha,leading_alpha
    integer, allocatable :: nnzinrow(:)
    integer, allocatable :: row_list(:)
    character(len=256) :: msg    
    logical :: debug

    debug=.False.

    if (debug) write(*,*) '  '
    if (debug) write(*,*) 'BEGIN FORM COMPOSED'
    if (debug) write(*,*) '  '
    if (debug) call composed_matrix%info(6)

    
    info=0

    leading_alpha=composed_matrix%alpha   
    nlinops=composed_matrix%nlinops
    do ilinop=nlinops,1,-1
       select type (mat=>composed_matrix%linop_list(ilinop)%linop)
       type is (scalmat)
          if (abs(mat%scalar)<1e-15) then
             call this%init(lun_err, composed_matrix%nrow,composed_matrix%ncol,0,'csr',&
                  is_symmetric=composed_matrix%nrow .eq. composed_matrix%ncol)
             this%ia(:)=1
          return
       end if
    end select
    end do

    !
    ! find the right-most operator that will be sparse
    ! cand copy it into temp_spmat
    !
    last_sparse=0
    do ilinop=nlinops,1,-1
       if (debug) call composed_matrix%linop_list(ilinop)%linop%info(6)
       select type (mat=>composed_matrix%linop_list(ilinop)%linop)
       type is (block_linop)
          last_sparse=ilinop
          call temp_spmat%form_block(info,lun_err,mat)
          !call mat%info(6)
          exit
       type is (new_linop)
          last_sparse=ilinop
          call temp_spmat%form_new_linop(info,lun_err,mat)
          exit
       type is (pair_linop)
          last_sparse=ilinop
          call temp_spmat%form_pair_linop(info,lun_err,mat)
          exit
       type is (spmat)
          last_sparse=ilinop
          temp_spmat=mat
          exit
       class default
          info=1
          write (msg,*) 'Linear operator ', ilinop , &
               ' in composed linear operator not supported'
          if(res .ne. 0) rc = IOerr(lun_err,&
               err_val, 'form', &
               etb(msg))
       end select
    end do


    if (debug) write(*,*)  'last_sparse, nlinops',last_sparse, nlinops
    !
    ! form a diagonal matrix with the product of the
    ! operator afetr the the last sparse linear operator
    !
    if ( last_sparse .ne. nlinops) then
       !
       ! set dimension of diagonal matrix
       ! wiht thif (debug)if (debug)if (debug)if (debug)if (debug)e minimal dimension active
       !
       
       dim=composed_matrix%linop_list(nlinops)%linop%nrow
       do ilinop=nlinops,last_sparse+1,-1         
          dim=min(dim, composed_matrix%linop_list(ilinop)%linop%nrow)
          dim=min(dim, composed_matrix%linop_list(ilinop)%linop%ncol)
       end do
       call diagonal_right%init(lun_err,&
            dim,&
            dim)
       if (debug) write(*,*) 'RIGHT MOST DIAGONAL MATRIX'
       if (debug) call diagonal_right%info(6)

       !
       ! build the right most diagonal matrix
       !
       diagonal_right%diagonal(1:dim)=one
       do ilinop=nlinops,last_sparse+1,-1
          select type (mat=>composed_matrix%linop_list(ilinop)%linop)
          type is(diagmat)
             diagonal_right%diagonal(1:dim) =  diagonal_right%diagonal(1:dim) * mat%diagonal(1:dim) 
          type is(scalmat)
             call dscal(dim,mat%scalar,diagonal_right%diagonal(1:dim),1)
          end select
       end do

       if ( last_sparse .eq. 0  ) then
          !
          ! convert 2 sparse matrix
          !
          call dscal(dim,leading_alpha,diagonal_right%diagonal(1:dim),1)
          if (debug) call this%info(6)
          if (debug) write(*,*) minval(this%coeff),maxval(this%coeff)

          !
          ! convert into sparse matrix
          !
          call this%form_diag(info,lun_err,diagonal_right)
          
          !
          ! free memory
          !
          call diagonal_right%kill(lun_err)
          return
       end if

       

       !
       ! select only active rows and column 
       !
       if ( ( dim < composed_matrix%linop_list(last_sparse+1)%linop%nrow ) .or. &
            ( dim < composed_matrix%linop_list(nlinops)%linop%ncol       ) ) then
          allocate (row_list(dim)) 
          do i=1,dim
             row_list(i)=i
          end do
          if ( dim < composed_matrix%linop_list(last_sparse+1)%linop%nrow ) then
             temp_spmat2=temp_spmat
             call temp_spmat2%select_permute_columns(dim,row_list, temp_spmat, lun_err)
          end if
          if  ( dim < composed_matrix%linop_list(nlinops)%linop%ncol ) then
             temp_spmat2=temp_spmat
             call temp_spmat2%select_permute_rows(dim,row_list ,temp_spmat,lun_err)
          end if
       end if


       !
       ! compute A * D
       !
       call temp_spmat%MxD(lun_err,diagonal_right%diagonal(1:dim))

       
       !
       ! add zero below and at right
       ! A becomes A 0
       !           0 0
       ! fitting the dimensions
       if ( dim < composed_matrix%linop_list(last_sparse+1)%linop%nrow ) then
          temp_spmat%ncol=composed_matrix%linop_list(nlinops)%linop%ncol
       end if
       if ( dim < composed_matrix%linop_list(nlinops)%linop%ncol ) then
          call add_tailing_zeros(temp_spmat,&
               lun_err,composed_matrix%linop_list(last_sparse)%linop%nrow-dim)
       end if

       !
       ! free memory
       !
       call diagonal_right%kill(lun_err)
       if ( allocated(row_list) ) deallocate (row_list)
       
    end if
    

    
    do ilinop=last_sparse-1,1,-1
       !
       ! copy last result into temp_right
       !
       temp_right=temp_spmat

       !
       ! selcet type
       !
       select type ( mat => composed_matrix%linop_list(ilinop)%linop )
       type is (block_linop)
          call temp_spmat2%form_block(info,lun_err,mat)
          call temp_spmat%mult_MDN( lun_err, &
               temp_spmat2,& ! lhs
               temp_right, & ! rhs
               100, 100*temp_spmat2%nrow)
          call temp_spmat2%kill(lun_err)
       type is (new_linop)
          !
          ! create a temporary matrix to form current operator
          !
          call temp_spmat2%form_new_linop(info,lun_err,mat)
          !
          ! muliply LHS RHS
          ! where RHS=(what we multiply so far by
          ! the definition temp_right=temp_spmat)
          !
          call temp_spmat%mult_MDN( lun_err, temp_spmat, &
               temp_right, 100, 100*temp_spmat2%nrow)
          !
          ! free memory
          !
          call temp_spmat2%kill(lun_err)
       type is (spmat)
          !
          ! just multiply
          !
          call temp_spmat%mult_MDN(lun_err, mat, temp_right, 100, 100*mat%nrow)
       type is(diagmat)
          call temp_spmat%DxM(lun_err,mat%diagonal(1:temp_right%nrow))
          if (mat%nrow .ne. mat%ncol) call add_tailing_zeros(temp_spmat, lun_err,mat%nrow-this%nrow)
       type is(scalmat)
          temp_spmat%coeff = temp_spmat%coeff*mat%scalar
          if (mat%nrow .ne. mat%ncol) call add_tailing_zeros(temp_spmat, lun_err,mat%nrow-this%nrow)
       end select
    end do

    !
    ! multiply leading alpha
    !
    temp_spmat%coeff = leading_alpha*temp_spmat%coeff

    !
    ! defined matrix
    !
    this = temp_spmat
    call this%is_like(composed_matrix)

    !
    ! free memory
    !
    call temp_spmat%kill(lun_err)
    
    if (debug) call this%info(6)
    if (debug) write(*,*) '  '
    if (debug) write(*,*) 'END FORM COMPOSED'
    if (debug) write(*,*) '  '

  end subroutine form_composed
    subroutine add_tailing_zeros(this,lun_err,nextra_lines)
      use Globals
      type(spmat), intent(inout) :: this
      integer, intent(in   ) :: lun_err
      integer, intent(in   ) :: nextra_lines
      !local
      type(spmat) :: temp

      temp=this
      call this%init(lun_err, &
           temp%nrow+nextra_lines,temp%ncol,temp%nterm,'csr',&
           is_symmetric=temp%is_symmetric,&
           triangular=temp%triangular)
      this%ia(1:temp%nrow+1)=temp%ia(1:temp%nrow+1)
      this%ia(temp%nrow+2:this%nrow+1)=temp%ia(temp%nrow+1)
      this%ja   =temp%ja
      this%coeff=temp%coeff
      call temp%kill(lun_err)
      
    end subroutine add_tailing_zeros
  
    
    
  end subroutine form_pair_linop


  recursive subroutine form_new_linop(this,info,lun_err, matrix)
    use Globals
    class(spmat),             intent(inout) :: this
    integer,                  intent(inout) :: info
    integer,                  intent(in   ) :: lun_err
    class(new_linop), target, intent(in   ) :: matrix
    !local
    logical :: rc
    integer :: res
    integer :: i,j,nlinops,index
    type(spmat), target, allocatable :: temp_spmats(:)
    type(new_linop)  :: copy

    !
    ! short-hand for faster procedure if the is only one
    ! pair_linop in new_linop 
    !
    if ( matrix%npairs .eq. 1) then
       call this%form_pair_linop(info,lun_err, matrix%pair_list(matrix%npairs))
       call this%is_like(matrix)
    else  
       !
       ! allocate space for temporary sparse matrices 
       !
       allocate (temp_spmats(matrix%npairs-1),stat=res)
       if(res .ne. 0) rc = IOerr(lun_err,&
            err_alloc,  'form_new_linop', &
            'work array temp_spmats')
       !
       ! create a copy in which will redirects pointers
       ! to pair_linop type to formed sparse matrices
       !
       copy = matrix

       !
       ! Form pairs but the last.
       ! If there is one one pair in new_linop these parts is ignored
       !
       do i=1,copy%npairs-1
          !
          ! Form sparse matrices. All matrices are "concrete" type
          ! diag_mat , spmat, full mat etc.
          ! The first pair should be always made by this type of matrices. 
          !
          call temp_spmats(i)%form_pair_linop(info,lun_err, matrix%pair_list(i) )

          !
          ! redirect pointers stored in the next pair_linop to formed matrices
          !
          do j=1,copy%pair_list(i+1)%nlinops
             select type( mat => copy%pair_list(i+1)%linop_list(j)%linop )
             type is (pair_linop)
                index = mat%index
                copy%pair_list(i+1)%linop_list(j)%linop => temp_spmats(index)
             end select
          end do
       end do

       !
       ! form last pair
       !
       call this%form_pair_linop(info,lun_err, copy%pair_list(copy%npairs))
       call this%is_like(matrix)
       !
       ! free memory
       !
       call copy%kill(lun_err)
       deallocate (temp_spmats,stat=res)
       if(res .ne. 0) rc = IOerr(lun_err,&
            err_alloc,  'form_new_linop', &
            'work array temp_spmats')
    end if
    
  end subroutine form_new_linop


  
  !>-------------------------------------------------------------
  !> Procedure to form an explicit sparse matrix from a block
  !> matrix. 
  !>
  !> (public procedure for type spmat)
  !> 
  !> usage:
  !>     call 'var'%(info,lun_err,block_matrix)
  !>
  !> where 
 
  !> \param[in   ] info         -> integer. Flag for errors 
  !> \param[in   ] lun_err      -> integer. Logical I/O unit
  !> \param[in   ] block_matrix -> type(block_linop)
  !>                               Block linear operator
  !>                               made by components of type
  !>                               spmat, diagmat, scalmat
  !<-------------------------------------------------------------
  recursive subroutine form_block(this,info,lun_err, block_matrix)
    use Globals
    class(spmat),      intent(inout) :: this
    integer,           intent(inout) :: info
    integer,           intent(in   ) :: lun_err
    type(block_linop), intent(in   ) :: block_matrix

    !local
    logical :: rc
    integer :: res
    integer :: i,j,k
    integer :: nrow_loc, ncol_loc,shift,nterm
    integer :: istart,iend,ilinop,ntemp
    integer :: irow, icol,jcol, pos, start, finish,nel,ind
    type(spmat), target, allocatable :: temp_spmats(:)
    real(kind=double) :: alpha
    integer, allocatable :: nnzinrow(:)
    integer, allocatable :: transpose_structure(:,:)
    type(array_linop), allocatable :: local_list(:)
    character(len=256) :: msg

    info=0
    
    allocate( nnzinrow(block_matrix%nrow),stat=res)
    if(res .ne. 0) rc = IOerr(lun_err,&
         err_alloc,  'block2sparse', &
         ' work array nnzinrow ')

    local_list=block_matrix%linop_list
    ntemp=0
    do i=1, block_matrix%nnzblock
       ilinop = block_matrix%block_structure(1,i)
       alpha  = block_matrix%alphas(i)

       if (alpha .ne. zero) then         
          select type (mat=>block_matrix%linop_list(ilinop)%linop)
          type is (block_linop)
             ntemp=ntemp+1
          type is (new_linop)
             ntemp=ntemp+1
          end select
       end if
    end do
    if (ntemp>0) then
       
       allocate(temp_spmats(ntemp),stat=res)
       if(res .ne. 0) rc = IOerr(lun_err,&
            err_alloc,  'block2sparse', &
            ' work temp_spmats ')
       ntemp=0
       do i=1, block_matrix%nnzblock
          ilinop = block_matrix%block_structure(1,i)
          alpha  = block_matrix%alphas(i)

          if (alpha .ne. zero) then         
             select type (mat=>block_matrix%linop_list(ilinop)%linop)
             type is (block_linop)
                ntemp=ntemp+1
                call temp_spmats(ntemp)%form_block(info,lun_err,mat)
                local_list(ilinop)%linop => temp_spmats(ntemp)
             type is (new_linop)
                ntemp=ntemp+1
                call temp_spmats(ntemp)%form_new_linop(info,lun_err,mat)
                local_list(ilinop)%linop => temp_spmats(ntemp)
             end select
          end if
       end do
    end if
    !
    ! count nonzeros in each row
    !
    ntemp    = 0
    nnzinrow = 0
    do i=1, block_matrix%nnzblock
       ilinop = block_matrix%block_structure(1,i)
       irow   = block_matrix%block_structure(2,i)
       icol   = block_matrix%block_structure(3,i)
       alpha  = block_matrix%alphas(i)
       
       if (alpha .ne. zero) then         
          nrow_loc = block_matrix%nrow_vec(irow)
          ncol_loc = block_matrix%ncol_vec(icol)


          istart = sum(block_matrix%nrow_vec(1:irow-1))+1
          iend   = istart + nrow_loc - 1

          select type (mat=>local_list(ilinop)%linop)
          type is (spmat)
             do k=1,mat%nrow
                nnzinrow(istart+k-1)=&
                     nnzinrow(istart+k-1) + &
                     mat%ia(k+1)-mat%ia(k)
             end do
          type is(diagmat)
             nnzinrow(istart:iend)= nnzinrow(istart:iend)+1
          type is(scalmat)
             if ( abs(mat%scalar)>small ) then
                nnzinrow(istart:iend)= nnzinrow(istart:iend)+1
             end if
          class default
             info=1
             write (msg,*) 'Linear operator ', ilinop , &
                  ' in block (', irow,',', icol, ') not supported'
             if(res .ne. 0) rc = IOerr(lun_err,&
                  err_val, 'form', &
                  etb(msg))
          end select
       end if
    end do

    !
    ! allocate memory and form ia
    !
    nterm=sum(nnzinrow)
    call this%init(lun_err, &
         block_matrix%nrow, block_matrix%ncol, nterm,&
         'csr')
    call this%is_like(block_matrix)
    this%nterm = nterm
    this%ia(1)=1
    do i=1,this%nrow
       this%ia(i+1)=this%ia(i)+nnzinrow(i)
    end do

    !
    ! use nnzinrow as scratch array to store offset 
    !
    ntemp=0
    nnzinrow=this%ia(1:block_matrix%nrow)
    do i=1,block_matrix%nnzblock
       
       ! we read the component by column
       ilinop = block_matrix%block_structure(1,i)
       irow   = block_matrix%block_structure(2,i)
       icol   = block_matrix%block_structure(3,i)
       alpha  = block_matrix%alphas(i)
       
       if (alpha .ne. zero) then
          nrow_loc = block_matrix%nrow_vec(irow)
          ncol_loc = block_matrix%ncol_vec(icol)

          ! column index = local cloumn index + shift
          shift=sum(block_matrix%ncol_vec(1:icol-1))

          ! portion of offeset to be considered
          istart = sum(block_matrix%nrow_vec(1:irow-1))+1
          iend   = istart + nrow_loc - 1

          ! select different type and fix in nnzinrow(istart:iend) the
          ! next starting index in ja 
          select type (mat=>local_list(ilinop)%linop)
          type is (spmat)
             do k=1,mat%nrow
                do j = 1,mat%ia(k+1)-mat%ia(k)
                   this%ja   (nnzinrow(istart+k-1)+j-1)=&
                        mat%ja(mat%ia(k)+j-1)+shift
                   this%coeff(nnzinrow(istart+k-1)+j-1)=&
                        alpha*mat%coeff(mat%ia(k)+j-1)
                end do
                ! add the number of non-zeros includedx
                nnzinrow(istart+k-1) = nnzinrow(istart+k-1) + &
                     mat%ia(k+1) - mat%ia(k)
             end do
          type is (diagmat)
             do j=1,min(mat%nrow,mat%ncol)
                this%ja   (nnzinrow(istart+j-1))= &
                     shift + j 
                this%coeff(nnzinrow(istart+j-1))= &
                     alpha*mat%diagonal(j)
             end do
             nnzinrow(istart:iend) = nnzinrow(istart:iend) + 1
          type is (scalmat)
             do j=1,min(mat%nrow,mat%ncol)
                this%ja   (nnzinrow(istart+j-1))= shift + j 
                this%coeff(nnzinrow(istart+j-1))= alpha*mat%scalar
             end do
             nnzinrow(istart:iend) = nnzinrow(istart:iend) + 1
          end select
       end if
    end do

    !
    ! set properties
    !
    call this%is_like(block_matrix)

    !
    ! free memory
    !
    deallocate(nnzinrow,stat=res)
    if(res .ne. 0) rc = IOerr(lun_err,&
         err_dealloc, 'block2sparse', &
         ' work array nnzinrow ')
    
  end subroutine form_block


  !>---------------------------------------------------------------
  !> Procedure to transform a symmetric matrix from csr into ssr
  !> format. Initialize the varaible mat_out if it is not
  !> 
  !> usage: call var%ssr2csr(lun_err,mat_out)
  !>
  !!> where:
  !> \param[in ] lun_err            -> integer. I\O unit for error message
  !> \param[out] (optional) mat_out -> type(spmat). Matrix in csr format
  !<---------------------------------------------------------------------
  subroutine form_diag(this, info, lun_err,diagonal_matrix)
    use Globals
    implicit none
    class(spmat),  intent(inout) :: this
    integer,       intent(inout) :: info
    integer,       intent(in   ) :: lun_err
    type(diagmat), intent(in   ) :: diagonal_matrix
    ! local
    integer :: nterm,i
    
    info=0
    !
    ! initialize the matrix
    !
    nterm=min(diagonal_matrix%ncol,diagonal_matrix%nrow)
    call this%init(lun_err,&
         diagonal_matrix%nrow,&
         diagonal_matrix%ncol,&
         nterm,'csr')
    call this%is_like(diagonal_matrix)
    

    !
    ! assign ia, ja,coeff
    !
    this%ia(1)=1
    do i=1,this%nrow
       this%ia(i+1)=this%ia(i)+1
    end do
    do i=1,nterm
       this%ja(i)=i
    end do
    do i=1,nterm
       this%coeff(i)=diagonal_matrix%diagonal(i)
    end do
    
  end subroutine form_diag
  
  
end module SparseMatrix

