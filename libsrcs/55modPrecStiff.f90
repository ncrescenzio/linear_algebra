module StiffPrecondioner
  use Globals
  use LinearOperator
  use Matrix
  use SparseMatrix
  use StdSparsePrec
  use DenseMatrix
  use StiffnessMatrix
  use DensePreconditioner
  use Scratch
  implicit none
  private
  !>-------------------------------------------------------------
  !> Structure variable containg all the varaibles
  !> and the procedure to build and apply  preconditioner
  !> for matrix given a sparse matrix.
  !> Scratch arrays for the construction and application of the
  !> preconditioner are included in this type.
  !> The preconditioner are
  !>  0) P=ID
  !>  1) P=Diag(A)
  !>------------------------------------------------------------- 
  type, extends(abs_linop), public :: stiffprec
     !> Flag prec. has been initialized
     logical :: is_built = .false.
     !> Controls of prec. construction and application
     type(input_prec) :: ctrl
     !> Number of sparse matrices preconditioner
     integer :: nsp_prec=0
     !> Sparse triangular preconditioner
     !> Dimension(nsp_prec)
     type(stdprec), allocatable :: sparse_precs(:)
     !> Number of diagonal matrices
     integer :: ndiag=0
     !> Diagonal matrices
     !> Dimension(nequ,ndiag)
     real(kind=double), allocatable :: diagonals(:,:)
     !> Number of sparse matrices preconditioner
     integer :: ndenseprec=0
     !> Preconditioners for dense matrices
     !> Dimension(ndenseprec)
     type(denseprec), allocatable :: dense_precs(:)
   contains
     !> Static constructor 
     !> (procedure public for type stiffprec)
     procedure, public, pass :: init => init_stiffprec
     !> Static constructor 
     !> (procedure public for type stiffprec)
     procedure, public, pass :: assembly => assembly_stiffprec
     !> static destructor
     !> (procedure public for type stiffprec)
     procedure, public, pass :: kill => kill_stiffprec
     !> Info procedure  sparse matrix
     !> (public procedure for type stiffprec)
     procedure, public, pass :: info => info_stiffprec
     !> Given vector (vec) compute (pvec) = Prec . vec  
     !> (public procedure for type stiffprec)
     procedure, public, pass :: matrix_times_vector => Mxv_stiffprec
  end type stiffprec
contains
  !>-------------------------------------------------------------
  !> Info Procedure.
  !> (procedure public for type input_prec)
  !>
  !> usage:
  !>     call 'var'%info(lun_out)
  !>
  !> where:
  !> \param[in] lun_out -> integer. I/O unit for output mesg.
  !<-------------------------------------------------------------
  function should_I_assembly(this,matrix2prec) result(assembly)
    use Globals
    implicit none
    class(input_prec), intent(in) :: this
    type(stiffmat), target, intent(in) :: matrix2prec
    logical :: assembly
    !local
    class(abs_matrix), pointer :: matrix
    

    assembly = .False.

    select case  (this%prec_type)
    case ('identity')
       assembly = .False.
    case ('diag')
       assembly = .False.
    case ('svd')
       assembly = .True.
    case ('C')
       assembly = .True.
    end select
       
  end function should_I_assembly
  

  !>-------------------------------------------------------------
  !> Static constructor.
  !> (procedure public for type stiffprec)
  !> Instantiate and initiliaze variable of type stiffprec
  !>
  !> usage:
  !>     call 'var'%init(lun_err,ctrl, matrix)
  !> where:
  !> \param[in] lun_err -> integer. Error logical unit
  !> \param[in] ctrl    -> type(input_prec) Controls of prec. construction
  !> \param[in] matrix  -> type(spmat). Matrix that will be used to build
  !>                                       the preconditier
  !<-------------------------------------------------------------
  subroutine init_stiffprec(this, lun_err, matrix, ctrl,info) 
    use Globals
    implicit none
    !var
    class(stiffprec), intent(inout) :: this
    integer,          intent(in   ) :: lun_err
    type(stiffmat),   intent(in   ) :: matrix
    type(input_prec), intent(in   ) :: ctrl
    integer,          intent(inout) :: info
    ! local vars
    integer :: res
    logical :: rc
    integer :: iprec,nrow, nterm, n_fillin, niaux, nraux,ntermp,nequ
    real(kind=double) :: tol_fillin

    nrow  = matrix%nrow

    ! clean if necessary
    if (this%is_initialized) call this%kill(lun_err)

    nequ       = matrix%nrow
    this%nrow  = matrix%nrow
    this%ncol  = matrix%ncol

    n_fillin   = ctrl%n_fillin
    tol_fillin = ctrl%tol_fillin
    ! Make a copy
    this%ctrl  = ctrl
    

    if ( ctrl%prec_type .eq. 'identity') then
       this%is_symmetric = .True.
    end if
    if (ctrl%prec_type .eq. 'diag') then
       this%is_symmetric = .True.
       ! allocation 
       allocate(this%diagonals(nequ,1),stat=res) 
       if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'init_stiffprec', &
            ' type prec member diagonal',res)
    end if
    
    if ( ctrl%prec_type .eq. 'C') then
       this%is_symmetric = .True.
       select type ( mat => matrix%matrix_A)
       type is (densemat) 
          this%ndenseprec = 1
          if (.not. allocated (this%dense_precs)) then
             allocate (this%dense_precs(1),stat=res)
             if(res .ne. 0) rc = IOerr(lun_err, err_alloc, &
                  ' init_stiffprec', &
                  ' type prec member diagonal',res)
          end if
          call this%dense_precs(1)%init(lun_err,matrix%dense,ctrl)
       end select
    end if
    if ( ctrl%prec_type .eq. 'IC') then
       this%is_symmetric = .True.
       select type ( mat => matrix%matrix_A) 
!!$       type is (adjmat)
!!$          if (.not. allocated (this%sparse_precs)) then
!!$             this%nsp_prec = 1
!!$             allocate (this%sparse_precs(1),stat=res)
!!$             if(res .ne. 0) rc = IOerr(lun_err, err_alloc, &
!!$                  ' init_stiffprec', &
!!$                  ' type prec member sparse_precs',res)
!!$          end if
       type is (spmat)
          !
          ! allocate one sparse matrix for 
          !
          if (.not. allocated (this%sparse_precs)) then
             this%nsp_prec = 1
             allocate (this%sparse_precs(1),stat=res)
             if(res .ne. 0) rc = IOerr(lun_err, err_alloc, &
                  ' init_stiffprec', &
                  ' type prec member sparse_precs',res)
          end if          
       end select
    end if
    
    if ( ctrl%prec_type .eq.'svd') then 
       this%is_symmetric = .True.
       select type (mat => matrix%matrix_A)
       type is (densemat)
          this%ndenseprec = 1
          if (.not. allocated (this%dense_precs)) then
             allocate (this%dense_precs(1),stat=res)
             if(res .ne. 0) rc = IOerr(lun_err, err_alloc, &
                  ' init_stiffprec', &
                  ' type prec member dense_prec',res)
          end if
          call this%dense_precs(1)%init(lun_err,matrix%dense,ctrl)
       end select
    end if

    this%is_initialized = .true.
    if (info .eq. 0) this%is_built = .true.

  end subroutine init_stiffprec

  !>-------------------------------------------------------------
  !> Static constructor.
  !> (procedure public for type stiffprec)
  !> Instantiate and initiliaze variable of type stiffprec
  !>
  !> usage:
  !>     call 'var'%init(lun_err,ctrl, matrix)
  !> where:
  !> \param[in] lun_err -> integer. Error logical unit
  !> \param[in] ctrl    -> type(input_prec) Controls of prec. construction
  !> \param[in] matrix  -> type(spmat). Matrix that will be used to build
  !>                                       the preconditier
  !<-------------------------------------------------------------
  subroutine assembly_stiffprec(this,&
       lun_err, matrix, ctrl,info) 
    use Globals
    use Timing
    
    implicit none
    !var
    class(stiffprec),      intent(inout) :: this
    integer,               intent(in   ) :: lun_err
    type(stiffmat),        intent(inout) :: matrix
    type(input_prec), intent(in   ) :: ctrl
    integer,               intent(inout) :: info
    ! local vars
    integer :: res
    logical :: rc
    integer :: iprec,nrow, nterm, n_fillin, niaux, nraux,ntermp,nequ
    real(kind=double) :: tol_fillin
    type(Tim) :: cholesky,svd
    type(spmat) :: matT

    
    call cholesky%init()
    call svd%init()

    nrow  = matrix%nrow
    nequ  = matrix%nrow

    if ( .not. this%ctrl%prec_type .eq. ctrl%prec_type) &
         rc = IOerr(lun_err, err_inp, &
         ' assembly_stiffprec' //&
         ' ctrls not compatible'// &
         ' init     : '//etb(this%ctrl%prec_type)//&
         ' assembly : '//etb(ctrl%prec_type))

    select case( this%ctrl%prec_type) 
    case  ('identity')
       this%is_symmetric = .True.
       this%is_built     = .True. 
       info = 0
    case  ('diag')
       this%is_symmetric = .True.
       !
       ! construction
       !
       call matrix%get_diagonal(this%diagonals(:,1))
       this%diagonals(:,1) = one / &
            this%diagonals(:,1)
       
       this%is_built     = .True.
       info = 0
    case ('C') 
       this%is_symmetric = .True.
       select type (mat => matrix%matrix_A)
       type is (densemat)         
          if ( .not. matrix%is_assembled ) then 
             rc = IOerr(lun_err, wrn_out, &
                  ' assembly_stiffprec' ,&
                  ' stiff matrix not assembled explicitely')
             info = -1
          end if
          call cholesky%set('start')
          call this%dense_precs(1)%assembly(lun_err,matrix%dense,info)
          call cholesky%set('stop')
          if ( info .eq. 0) then 
             this%is_built = .True.
          else
             rc = IOerr(lun_err, wrn_out, &
                  ' assembly_stiffprec' ,&
                  ' assembly of dense cholesky'//&
                  ' factorization failed with info =', info)
          end if
          call cholesky%info(6, 'CHOLESKY ASSEMBLY ' )
          call cholesky%kill()
       end select
          
    case ('IC')     
        select type(mat => matrix%matrix_A)
        type is (spmat) 
           if ( .not. matrix%is_assembled ) then 
              matT=mat
              call matT%transpose(lun_err)
              call matrix%sparse%mult_MDN(lun_err,&
                   mat,matT,&
                   nrow,100*nrow,&
                   matrix%tdens*matrix%inv_weight)
              call matT%kill(lun_err)
           end if
           !
           ! build preconditioner
           !
           call this%sparse_precs(1)%init(lun_err, info,&
                ctrl,nequ, matrix%sparse)
           if ( info .eq. 0) then 
              this%is_built = .True.
           else
              rc = IOerr(lun_err, wrn_out, &
                   ' assembly_stiffprec' ,&
                   ' assembly of incomplete cholesky'//&
                   ' factorization failed with info =', info)
           end if
!!$       type is (adjmat)
!!$          if ( .not. matrix%is_assembled ) then 
!!$             rc = IOerr(lun_err, wrn_out, &
!!$                  ' assembly_stiffprec' ,&
!!$                  ' stiff matrix not assembled explicitely')
!!$             info = -1
!!$          end if
!!$          !
!!$          ! build preconditioner
!!$          !
!!$          call this%sparse_precs(1)%init(lun_err, info,&
!!$               ctrl,nequ, matrix%sparse)
!!$
       end select
    case ('svd')
       this%is_symmetric = .True.
       select type (mat => matrix%matrix_A)
       type is (densemat)
          if ( .not. matrix%is_assembled ) then 
             rc = IOerr(lun_err, wrn_out, &
                  ' assembly_stiffprec' ,&
                  ' stiff matrix not assembled explicitely')
             info = -1
          end if
          
          call svd%set('start')
          call this%dense_precs(1)%assembly(lun_err,matrix%dense,info)
          call svd%set('stop')
          if ( info .eq. 0) then 
             this%is_built = .True.
          else
             rc = IOerr(lun_err, wrn_out, &
                  ' assembly_stiffprec' ,&
                  ' assembly of dense cholesky'//&
                  ' factorization failed with info =', info)
          end if
          call svd%info(6, ' svd ' )
          call svd%kill()
       end select
       
       
    case default
       write(lun_err,*) 'wrong flag', ctrl%prec_type
    end select



  end subroutine assembly_stiffprec


  !>-------------------------------------------------------------
  !> Static desconstructor.
  !> (procedure public for type ctrl_precond)
  !> Free memory used  by stiffprec members
  !>
  !> usage:
  !>     call 'var'%kill(lun_err)
  !<-------------------------------------------------------------
  subroutine kill_stiffprec(this, lun_err) 
    use Globals
    implicit none
    !var
    class(stiffprec), intent(inout) :: this
    integer,        intent(in   ) :: lun_err
    ! local
    logical :: rc
    integer :: res,i

    if ( allocated(this%sparse_precs) ) then
       do i=1,this%nsp_prec
          call this%sparse_precs(i)%kill(lun_err)
       end do
       deallocate (this%sparse_precs,stat=res) 
       if(res .ne. 0) rc = IOerr(lun_err, err_dealloc, 'init_mat', &
            '  type stiffprec member sparse matrices',res)
       this%nsp_prec = 0       
    end if
    if ( allocated( this%diagonals ) ) then
       deallocate (this%diagonals,stat=res) 
       if(res .ne. 0) rc = IOerr(lun_err, err_dealloc, 'init_mat', &
            '  type stiffprec member diagonals',res)
       this%ndiag = 0 
    end if

    call this%ctrl%kill()

    this%is_built       = .false.

    !
    ! reset to default
    !
    call this%to_default()
    

  end subroutine kill_stiffprec


  !>-------------------------------------------------------------
  !> Info Procedure.
  !> (procedure public for type stiffprec)
  !>
  !> usage:
  !>     call 'var'%info(lun_out)
  !<-------------------------------------------------------------
  subroutine info_stiffprec(this, lun) 
    use Globals
    implicit none
    !var
    class(stiffprec), intent(in) :: this
    integer,        intent(in) :: lun

    if ( this%is_initialized) then
       call this%ctrl%info(lun)      
    else
       write(lun,*) 'Preconditioner not initialized '
    end if
  end subroutine info_stiffprec

  !>----------------------------------------------------------
  !> Precedure to apply preconditioner
  !> (procedure public for type stiffprec)
  !> Prodecure computing the vector w=result of        
  !>              w=P^{-1}v
  !> for a given vector v (arrays' dimension must match)
  !>
  !> usage: call 'var'apply%(vec_in,vec_out,[info])
  !> 
  !> \param[in   ] vec_in  -> real, dimension(this%nequ)
  !>                            Vector v where stiffprec will be applied
  !> \param[inout] vec_out  -> real, dimension(this%nequ)
  !>                            Vector w=P^{-1}v 
  !> \param[inout] (optional) info -> integer. Flag for errors 
  !>---------------------------------------------------------------------
  subroutine Mxv_stiffprec(this, vec_in, vec_out,info,lun_err) 
    use Globals
    implicit none
    !var
    class(stiffprec),  intent(inout) :: this
    real(kind=double), intent(in   ) :: vec_in(this%nrow)     
    real(kind=double), intent(inout) :: vec_out(this%ncol)
    integer,           intent(inout) :: info
    integer,           intent(in   ) :: lun_err
           
    ! local
    integer i,j,k,m,mm,n1,info_loc,nequ
    real(kind=double) :: a,dnrm2

    nequ = this%nrow
    
    select case (this%ctrl%prec_type)
    case ('identity')
       call dcopy(nequ,vec_in,1,vec_out,1)
       info_loc=0
    case ('diag')
       info_loc=0
       vec_out = vec_in * this%diagonals(:,1) 
    case ('C')
       info_loc=0
       if (allocated(this%dense_precs)) then
          call this%dense_precs(1)%Mxv(vec_in,vec_out,info_loc,lun_err)
       end if
       if (allocated(this%sparse_precs)) then
          call this%sparse_precs(1)%Mxv(vec_in,vec_out,info_loc,lun_err)
       end if
    case ('IC')
       info_loc=0
       if (allocated(this%dense_precs)) then
          call this%dense_precs(1)%Mxv(vec_in,vec_out,info_loc,lun_err)
       end if
       if (allocated(this%sparse_precs)) then
          call this%sparse_precs(1)%Mxv(vec_in,vec_out,info_loc,lun_err)
       end if
    case ('svd')
       info_loc=0
       call this%dense_precs(1)%Mxv(vec_in,vec_out,info_loc,lun_err)
    end select
       
    info = info_loc

  end subroutine Mxv_stiffprec

  
end module StiffPrecondioner

