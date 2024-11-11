module StdSparsePrec
  use Globals
  use LinearOperator
  use Matrix
  use SparseMatrix
  use Scratch
  use LinearSolver
  implicit none
  private
  !> Structure variable containing the controls
  !> for the construction the type stdprec
  !> in this module
  type, public :: input_prec
     !> Integer identifying the Preconditinoer
     !>     'identity'  P=Id (nothing is done)
     !>     'diag'      P=diag(A)^{-1}
     !>     'IC'        P=(U^{T}U)^{-1} with A~(U^T)U
     !>     'ILU'       P=(LU)^{-1} with A~LU
     !>     'C'         P=(M^{T}M)^{-1} with A=(U^T)U
     !>     'LU'        P=(LU)^{-1} with A=LU
     character(len=20) :: prec_type='identity'
     !> Number of maximal extra-non-zero tem for
     !> row in the construction of IC_{fillin}
     integer :: n_fillin=0
     !> Tolerance for the dual drop procedure
     !> in the construction of IC_{fillin} (iprec=4)
     real(kind=double) :: tol_fillin=0.0d0
     !> Control for factorization breakdown
     !> 0 = stop factorization 
     !> 1 = try to recover factorization 
     !> 2 = try recovering in regime of "normal" number
     integer :: factorization_job=2
   contains
     !> static constructor
     !> (procedure public for type input_prec)
     procedure, public, pass :: init => init_input_prec
     !> static destructor
     !> (procedure public for type input_prec)
     procedure, public, pass :: kill => kill_input_prec
     !> Info procedure.
     !> (public procedure for type input_prec)
     procedure, public, pass :: info => info_input_prec
  end type input_prec

  !>-------------------------------------------------------------
  !> Structure variable containg all the varaibles
  !> and the procedure to build and apply standard preconditioner
  !> given a sparse matrix.
  !> Scratch arrays for the construction and application of the
  !> preconditioner are included in this type.
  !> The preconditioner are
  !>  1) P=Diag(A)
  !>  3) P=IC_0(A)
  !>  4) P=IC_{fillin}(A)
  !>--------------------------------------------------------
  type, extends(abs_linop), public :: stdprec
  !type,  public :: stdprec
     !> Integer identifying the Preconditinoer
     !>     'identity'     P=Id (nothing is done)
     !>     'inv_diag'     P=diag(A)^{-1}
     !>     'inv_sqrtdiag' P=diag(A)^{-1/2}
     !>     'inv_L'        P=L^{-1} with A~LU
     !>     'inv_U'        P=L^{-1} with A~LU
     !>     'inv_LU'       P=L^{-1} with A~LU
     character(len=20) :: label='identity'
     !> Flag prec. has been initialized
     logical :: is_built = .false.
     !> Controls of prec. construction and application
     type(input_prec) :: ctrl
     !>------------------------------------------------------------
     !> Number of sparse matrices used 
     integer :: nspmats=0
     !> Sparse matrices stored
     type(spmat), allocatable :: spmats(:)
     !>-------------------------------------------------------------
     !> Number of diagonal matrices used 
     integer :: ndiagmats=0
     !> Diagonal matrix
     real(kind=double), allocatable :: diagmats(:,:)
     !> Auxilary variable for the application of the precondtioner
     type(scrt) :: aux_apply
     !> Auxilary variable for the construction of the precondtioner
     type(scrt) :: aux_build
   contains
     !> Static constructor 
     !> (procedure public for type stdprec)
     procedure, public, pass :: init => init_stdprec
     !> Static constructor 
     !> (procedure public for type stdprec)
     procedure, public, pass :: assembly => assembly_stdprec
     !> static destructor
     !> (procedure public for type stdprec)
     procedure, public, pass :: kill => kill_stdprec
     !> Info procedure  sparse matrix
     !> (public procedure for type stdprec)
     procedure, public, pass :: info => info_stdprec
     !> Compute pvec = P^{-1} vec 
     !> (public procedure for type stdprec)
     procedure, public,  pass :: matrix_times_vector => Mxv_stdprec
     !> Procedure to build the pair of preconditioners
     !> P1, P2 from one original preconditioner P
     !> Three optional are avaliable
     !> 'left' => precs(1) = P, Id
     !> 'right' => precs(1) = Id, P
     !> 'split' => precs(1) = sqrt(P)_1, sqrt(P)_2
     !> Where sqrt(P)_i depend on the type of preconditioner P
     !> (procedure public for type stdprec)
     procedure, public, pass :: one2two => init_stdprecs 
  end type stdprec
  
 


contains
  !>-------------------------------------------------------------
  !> Static constructor.
  !> (procedure public for type ctrl_precond)
  !> Instantiate and initilize variable of type ctrl_precond
  !>
  !> usage:
  !>     call 'var'%init(prec_type,[
  !>                    [n_fillin=..],&
  !>                    [tol_fillin =..])
  !>
  !> where:
  !> \param[in] prec_type             -> string. Id. for standard preconditier
  !>                                     'identity' => P = ID
  !>                                     'diag'     => P is Diag(A)
  !>                                     'IC'       => P is IC(A)
  !>                                     'ILU'      => P is ILU(A) 
  !> \param[in] (optional) n_fillin    -> integer. (Used only for 'LU' prec )
  !>                                      Number of extra non-zero element
  !>                                      for row of Incomplete Cheoleski
  !>                                      and LU decomposition 
  !> \param[in] (optional) tol_fillin -> real. (Used only for 'IC'
  !>                                            and 'ILU' prec )
  !>                                        Drop tolerance for Incomplete
  !>                                        Cheoleski and LU decomposition 
  !<-------------------------------------------------------------
  subroutine init_input_prec(this,&
       err_unit,&
       prec_type, n_fillin, tol_fillin)
    use Globals

    implicit none
    class(input_prec),           intent(inout) :: this
    integer,                     intent(in   ) :: err_unit
    character(len=*),            intent(in   ) :: prec_type
    integer,           optional, intent(in   ) :: n_fillin
    real(kind=double), optional, intent(in   ) :: tol_fillin
    ! local
    logical :: rc
    !this%approach = 'STDPREC'
    if ( ( prec_type .ne. 'identity' ) .and. &
         ( prec_type .ne. 'diag'     ) .and. &
         ( prec_type .ne. 'IC'       ) .and. &
         ( prec_type .ne. 'ILU'      ) .and. & 
         ( prec_type .ne. 'C'        ) .and. &
         ( prec_type .ne. 'SVD'      ) ) then
       rc = IOerr(err_unit, wrn_inp, 'init_input_prec', &
            'prec_type = '//etb(prec_type)//' not supperted')
    end if
    this%prec_type   = etb(prec_type)
    
    !
    ! optional arguments for incomplete factorization based precs.
    !
    if (present(n_fillin))   this%n_fillin   = n_fillin
    if (present(tol_fillin)) this%tol_fillin = tol_fillin
    
  end subroutine init_input_prec


  !>-------------------------------------------------------------
  !> Static desconstructor.
  !> (procedure public for type ctrl_precond)
  !> Reset to default values type input_prec
  !>
  !> usage:
  !>     call 'var'%kill()
  !<-------------------------------------------------------------
  subroutine kill_input_prec(this)
    use Globals
    implicit none
    class(input_prec), intent(inout) :: this

    this%prec_type   = 'identity'
    this%n_fillin    = 0
    this%tol_fillin  = 0.0d0

  end subroutine kill_input_prec


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
  subroutine info_input_prec(this,lun)
    use Globals
    implicit none
    class(input_prec), intent(in) :: this
    integer,          intent(in) :: lun

    if ( (this%prec_type .ne. 'IC') .and.(this%prec_type .ne. 'ILU')) then
       write(lun,'(a,a)') &
            'prec_type = ',etb(this%prec_type)
    else
       write(lun,'(a,a,a,a,I3,a,1pe12.5,a)') &
            'prec_type = ', etb(this%prec_type),&
            ' ( ', &
            'nfill =', this%n_fillin,&
            ', drop. tol.=',this%tol_fillin,&
            ')'
    end if

  end subroutine info_input_prec
  
  !>-------------------------------------------------------------
  !> Static constructor.
  !> (procedure public for type stdprec)
  !> Instantiate and initiliaze variable of type stdprec
  !>
  !> usage:
  !>     call 'var'%init(lun_err,ctrl, matrix)
  !> where:
  !> \param[in] lun_ersr -> integer. Error logical unit
  !> \param[in] ctrl    -> type(input_prec) Controls of prec. construct.
  !> \param[in] matrix  -> type(spmat). Matrix that will be used to build
  !>                                       the preconditier
  !<-------------------------------------------------------------  
  subroutine init_stdprec(this,lun_err,info,ctrl,nequ, matrix) 
    use Globals
    use SparseMatrix
    use LinearSolver, only : input_solver
    implicit none
    !var
    class(stdprec),        intent(inout) :: this
    integer,               intent(in   ) :: lun_err
    integer,               intent(inout) :: info
    type(input_prec),      intent(in   ) :: ctrl
    integer,               intent(in   ) :: nequ
    class(spmat), optional,intent(in   ) :: matrix

    
    ! local vars
    integer :: res
    logical :: rc
    integer :: i
    integer :: info_build
    integer :: nrow, nterm, n_fillin, niaux, nraux,ntermp
    real(kind=double) :: tol_fillin
    type(input_solver) :: ctrl_solver

    
    ! 
    ! initialized defualt case prec=Id
    !
    ! clean if necessary
    if (this%is_initialized) call this%kill(lun_err)
    
    !
    ! copy controls
    !
    this%ctrl = ctrl
    
    if ( .not. ( present(matrix) ) ) then
       ! 
       ! default case prec = identity
       ! 
       this%is_symmetric   = .true.
       this%nrow = nequ
       this%ncol = nequ
       info = 0
       this%is_initialized = .true.
       this%is_built       = .true.
       
    else
       if ( this%nrow .ne. this%ncol ) then
          rc = IOerr(lun_err, wrn_val, 'init_stdprec', &
               ' matrix is not squared',res)
          return
       end if

       !
       ! prec=diag(A)^{-1} 
       !       or   
       ! prec=(LU)^{-1} A~LU
       !
       nrow  = matrix%nrow
       nterm = matrix%nterm
       this%is_symmetric   = .true.
       this%nrow = nequ
       this%ncol = nequ
       n_fillin   = ctrl%n_fillin
       tol_fillin = ctrl%tol_fillin

       select case( ctrl%prec_type ) 
       case  ('identity')
          ! 
          ! default case prec = identity
          ! 
          info      = 0
          this%is_symmetric   = .true.
          this%nrow = matrix%nrow
          this%ncol = matrix%ncol
          this%is_initialized = .true.
          this%is_built       = .true.
       case  ('diag')
          !
          ! set appliation
          !
          this%label = 'inv_diag'          
          info      = 0
          this%is_symmetric   = .true.
          this%nrow = matrix%nrow
          this%ncol = matrix%ncol
          !
          ! allocate memory
          ! sparse matrix not rquired
          this%nspmats   = 0
          ! inverse of the diagonal 
          this%ndiagmats = 1
          allocate(this%diagmats(nequ,this%ndiagmats),stat=res) 
          if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'init_stdprec', &
               ' type prec member diagonal',res)

          ! build
          call matrix%get_diagonal(this%diagmats(:,1))
          this%diagmats(:,1) = one / this%diagmats(:,1)

          this%is_initialized = .true.
          this%is_built       = .true.
          
          info = 0
       case ('IC')
          !
          ! set application
          !
          this%label = 'inv_IC'
          info      = 0
          this%is_symmetric   = .true.
          this%nrow = matrix%nrow
          this%ncol = matrix%ncol

          if ( matrix%is_symmetric) then
             ! 
             ! for and optimazed computation and application 
             ! of the preconditioner we compute the sparse
             ! matrix diag(U)^{-1} U and D(U)^2.
             ! Thus the inversion will required less computaion
             !
             ! 1 - Allocation for factor diag(U)^{-1} U 
             !     and D(U)^2
             ! 2 - Build diag(U)^{-1} U factor and inverse of D(U)^2  
             ! 3 - scratch vector y for solving U^t U x = b 
             !     U^t y = b
             !     U   x = y

             ! 1.1 - allocate memory
             this%nspmats   = 1
             allocate(this%spmats(this%nspmats),stat = res)
             if(res .ne. 0) rc = IOerr(lun_err, err_alloc,&
                  'init_stdprec', &
                  ' type prec member spmats',res) 
             !
             this%ndiagmats = 1
             allocate(this%diagmats(nequ,this%ndiagmats),stat = res)
             if(res .ne. 0) rc = IOerr(lun_err, err_alloc,&
                  'init_stdprec', &
                  ' type prec member spmats',res)   
             this%is_initialized = .true.


             ! 2.1 : build U factor
             ! this%diagmats(1:nequ,1) =
             ! this%spmats(1)=
             info_build = 0
             call matrix%incomplete_cholesky(lun_err,&
                  n_fillin, tol_fillin, ctrl%factorization_job,&
                  info_build,&
                  this%spmats(1),&                   ! diag(U)^{-1} U
                  diagonal=this%diagmats(1:nequ,1))  ! diag(U)

             if ( info_build .ne. 0 ) then
                info = info_build
                return
             else
                this%is_built       = .true.
                info = 0 
             end if

             ! compute invert (D(U)^2)
             this%diagmats(1:nequ,1) = one / this%diagmats(:,1)**2 

             ! 3 - scratch vector y
             call this%aux_apply%init(lun_err,0,nrow)


          else
             rc = IOerr(lun_err, err_inp,&
                  'init_stdprec', &
                  ' prec_type not consistent with non-symmetric matrix'&
                  //etb(ctrl%prec_type))
          end if
       case('ILU')          
          !
          ! set application
          !
          this%label = 'inv_ILU'
          this%is_symmetric   = .false.
          info      = 0
          this%nrow = matrix%nrow
          this%ncol = matrix%ncol

          ! we compute 
          ! 1 - U and L factors 
          ! 2 - inverse of the diagonal of U and L factor
          ! 3 - scratch vector y for solving L U x = b 
          !     L y = b
          !     U x = y

          ! 1.1 - U and L factor
          this%nspmats   = 2
          allocate(this%spmats(this%nspmats),stat = res)
          if(res .ne. 0) rc = IOerr(lun_err, err_alloc,&
               'init_stdprec', &
               ' type prec member spmats',res)
          this%is_initialized = .true.


          call matrix%incomplete_lu(lun_err,&
               n_fillin, tol_fillin,ctrl%factorization_job, &
               info_build,&
               this%spmats(1),& !lower
               this%spmats(2))  !upper
          

          if ( info_build .ne. 0 ) then
             info = info_build
             rc =  IOerr(lun_err, err_out,&
                  'init_stdprec ', &
                  'incomplete_lu should return info_build=0, info_build = ',&
                  info_build)
          else
             this%is_built       = .true.  
             info = 0 
          end if


          ! 2.1 - inverse of diag( U )
          this%ndiagmats = 1
          allocate(this%diagmats(this%nrow,this%ndiagmats),stat = res)
          if(res .ne. 0) rc = IOerr(lun_err, err_alloc,&
               'init_stdprec', &
               ' type prec member diagmats',res)

          ! 2.2 - compute inverse of the diagonal 
          call this%spmats(2)%get_diagonal(&
               this%diagmats(:,1))
          ! inverse of diag(U)
          do i=1,this%nrow
             if ( this%diagmats(i,1) .ne. zero) then
                this%diagmats(i,1) = one / this%diagmats(i,1) 
             else
                info = -1
                rc =  IOerr(lun_err, err_out,&
                     'init_stdprec ', &
                     ' factor U in incomplete_lu has'//&
                     ' null diagonal entry at line ' ,&
                     i)
                return
             end if
          end do
          ! 3 - scratch vector y
          call this%aux_apply%init(lun_err,0,nequ)
          

       case default
          info = -1
          rc =  IOerr(lun_err, err_out,&
               ' init_stdprec ', &
               ' prec type not supported passed = '//etb(ctrl%prec_type))
          return

       end select

    end if
  end subroutine init_stdprec

  !>-------------------------------------------------------------
  !> Static constructor.
  !> (procedure public for type stdprec)
  !> Instantiate and initiliaze variable of type stdprec
  !>
  !> usage:
  !>     call 'var'%init(lun_err,ctrl, matrix)
  !> where:
  !> \param[in] lun_ersr -> integer. Error logical unit
  !> \param[in] ctrl    -> type(input_prec) Controls of prec. construct.
  !> \param[in] matrix  -> type(spmat). Matrix that will be used to build
  !>                                       the preconditier
  !<-------------------------------------------------------------  
  subroutine assembly_stdprec(this,lun_err,info,ctrl,matrix) 
    use Globals
    use SparseMatrix
    implicit none
    !var
    class(stdprec),        intent(inout) :: this
    integer,               intent(in   ) :: lun_err
    integer,               intent(inout) :: info
    type(input_prec),      intent(in   ) :: ctrl
    class(spmat),          intent(in   ) :: matrix

    
    ! local vars
    integer :: res
    logical :: rc
    integer :: i
    integer :: info_build
    integer :: nequ,nrow, nterm, n_fillin, niaux, nraux,ntermp
    real(kind=double) :: tol_fillin
    type(input_solver) :: ctrl_solver

    
    
    !
    ! copy controls
    !
    this%ctrl = ctrl

    
    !
    ! prec=diag(A)^{-1} 
    !       or   
    ! prec=(LU)^{-1} A~LU
    !
    nrow  = matrix%nrow
    nterm = matrix%nterm
    nequ  = matrix%ncol
    n_fillin   = ctrl%n_fillin
    tol_fillin = ctrl%tol_fillin
    
    select case( ctrl%prec_type ) 
    case  ('identity')
       ! 
       ! default case prec = identity
       ! 
       this%label          = 'identity'   
       this%is_symmetric   = .true.
       this%is_built       = .true.
       
       info = 0
    case  ('diag')
       !
       ! set appliation
       !
       this%label = 'inv_diag'    
       this%is_built       = .true.
       
       ! build 
       call matrix%get_diagonal(this%diagmats(:,1))
       this%diagmats(:,1) = one / this%diagmats(:,1)
          
       info = 0
    case ('IC')
       !
       ! set application
       !
       this%label          = 'inv_IC'
       this%is_symmetric   = .true.
       if ( .not. matrix%is_symmetric)  then
          rc = IOerr(lun_err, wrn_inp,&
               'init_stdprec', &
               ' prec_type not consistent with non-symmetric matrix'&
               //etb(ctrl%prec_type))
          info = 1
       else
          !
          ! Build diag(U)^{-1} U factor and inverse of D(U)^2 
          ! Thus the inversion will required less computaion
  
          ! 1 : build U factor
          info_build = 0
          call matrix%incomplete_cholesky(lun_err,&
               n_fillin, tol_fillin, ctrl%factorization_job,&
               info_build,&
               this%spmats(1),&                   ! diag(U)^{-1} U
               diagonal=this%diagmats(1:nequ,1))  ! diag(U)

          if ( info_build .ne. 0 ) then
             info = info_build
             return
          else
             this%is_built       = .true.
             info = 0 
          end if

          ! compute invert (D(U)^2)
          this%diagmats(1:nequ,1) = one / this%diagmats(:,1)**2 

       end if
    case('ILU')          
       !
       ! set application
       !
       this%label          = 'inv_ILU'
       this%is_symmetric   = .false.

       ! we compute 
       ! 2 - inverse of the diagonal of U and L factor
       ! 3 - scratch vector y for solving L U x = b 
       !     L y = b
       !     U x = y
       call matrix%incomplete_lu(lun_err,&
            n_fillin, tol_fillin,  ctrl%factorization_job, &
            info_build,&
            this%spmats(1),& !lower
            this%spmats(2))  !upper

       if ( info_build .ne. 0 ) then
          info = info_build
          rc =  IOerr(lun_err, err_out,&
               'init_stdprec ', &
               'incomplete_lu should return info_build=0, info_build = ',&
               info_build)
       else
          this%is_built       = .true.  
          info = 0 
       end if

       ! 2.2 - compute inverse of the diagonal 
       call this%spmats(2)%get_diagonal(&
            this%diagmats(:,1))

       ! inverse of diag(U)
       do i=1,this%nrow
          if ( this%diagmats(i,1) .ne. zero) then
             this%diagmats(i,1) = one / this%diagmats(i,1) 
          else
             info = -1
             rc =  IOerr(lun_err, err_out,&
                  'init_stdprec ', &
                  ' factor U in incomplete_lu has'//&
                  ' null diagonal entry at line ' ,&
                  i)
             return
          end if
       end do


    end select

  end subroutine assembly_stdprec


  
  !>-------------------------------------------------------------
  !> Static desconstructor.
  !> (procedure public for type ctrl_precond)
  !> Free memory used  by stdprec members
  !>
  !> usage:
  !>     call 'var'%kill(lun_err)
  !<-------------------------------------------------------------
  subroutine kill_stdprec(this, lun_err) 
    use Globals
    implicit none
    !var
    class(stdprec), intent(inout) :: this
    integer,        intent(in   ) :: lun_err
    ! local
    logical :: rc
    integer :: res,i

    if ( allocated(this%spmats) ) then
       do i = 1,this%nspmats
          if ( this%spmats(i)%is_initialized ) then
             call this%spmats(i)%kill(lun_err)
          end if
       end do
       deallocate (this%spmats,stat=res) 
       if(res .ne. 0) rc = IOerr(lun_err, err_dealloc, 'init_mat', &
            '  type stdprec member spmats',res)
    end if
    this%nspmats=0

    if ( allocated( this%diagmats ) ) then
       deallocate (this%diagmats,stat=res) 
       if(res .ne. 0) rc = IOerr(lun_err, err_dealloc, 'init_mat', &
            '  type stdprec member diagmats',res)
    end if
    this%ndiagmats=0


    if ( this%aux_build%is_initialized ) call this%aux_build%kill(lun_err)
    if ( this%aux_apply%is_initialized ) call this%aux_apply%kill(lun_err)

    call this%ctrl%kill()


    call this%to_default()

  end subroutine kill_stdprec


  !>-------------------------------------------------------------
  !> Info Procedure.
  !> (procedure public for type stdprec)
  !>
  !> usage:
  !>     call 'var'%info(lun_out)
  !<-------------------------------------------------------------
  subroutine info_stdprec(this, lun) 
    use Globals
    implicit none
    !var
    class(stdprec), intent(in) :: this
    integer,        intent(in) :: lun
    ! local vars
    integer :: i

    if ( this%is_initialized) then
       call this%ctrl%info(lun)
       write(lun,'(a,a,I1,a,I1)') etb(this%label), &
            '| nspamts=', this%nspmats,&
            '| ndiagmts=',this%ndiagmats
       do i=1,this%nspmats
          call this%spmats(i)%info(lun)
       end do
       
    else
       write(lun,*) 'Preconditioner not initialized '
    end if
  end subroutine info_stdprec

  !>----------------------------------------------------------
  !> Procedure to apply preconditioner
  !> (procedure public for type stdprec)
  !> Prodecure computing the vector w=result of        
  !>              w=P^{-1}v
  !> for a given vector v (arrays' dimension must match)
  !>
  !> usage: call 'var'apply%(vec_in,vec_out,[info])
  !> 
  !> \param[in   ] vec_in  -> real, dimension(this%nequ)
  !>                            Vector v where stdprec will be applied
  !> \param[inout] vec_out  -> real, dimension(this%nequ)
  !>                            Vector w=P^{-1}v 
  !> \param[inout] (optional) info -> integer. Flag for errors 
  !> \param[in   ] (optional) info -> integer. Logic unit for error msg.
  !>---------------------------------------------------------------------
  subroutine Mxv_stdprec(this, vec_in, vec_out,info,lun_err) 
    use Globals
    use Timing
    implicit none
    !var
    class(stdprec),    intent(inout) :: this
    real(kind=double), intent(in   ) :: vec_in(this%ncol)     
    real(kind=double), intent(inout) :: vec_out(this%nrow)
    integer,           intent(inout) :: info
    integer,           intent(in   ) :: lun_err
    ! local
    integer :: i,j,k,m,mm,n1,nequ,lun
    integer :: error_unit=0
    real(kind=double) :: a,dnrm2
    character(len=256) :: flag
    !> Time in scalar product
    type(Tim) :: trisol1,trisol2
    integer :: info_prec


    nequ  = this%nrow

    select case (this%label)
    case ('identity')
       !
       ! Identity (default precnditioner)
       !
       call dcopy(nequ,vec_in,1,vec_out,1)

       info_prec=0
    case( 'inv_diag') 
       ! 
       ! apply D(A)^{-1}
       !
       vec_out = vec_in * this%diagmats(1:nequ,1)

       info_prec=0
    case( 'inv_sqrtdiag') 
       ! 
       ! apply D(A)^{-1/2}
       !
       vec_out = vec_in * this%diagmats(1:nequ,1)
       
       info_prec=0

    case( 'inv_IC')
       !
       ! vec_out =  (U^{T}U)^{-1} vec_in
       !           ==
       ! solve U^{T}U (vec_out)  = vec_in
       !          
       ! only upper triangular factor is stored

       ! we stored D=(diag(U)^2)^{-1} and tilde(U)=diag(U)^{-1} U

       !       
       ! solve title{U} y = rhs
       !
       call this%spmats(1)%solve_triangular_unitdiag(&
            error_unit,&
            vec_in,&
            this%aux_apply%raux(1:nequ),&
            transpose='T')

      
       ! apply  diag(U)^2)^{-1} to y
       this%aux_apply%raux(1:nequ) = &
            this%aux_apply%raux(1:nequ) * this%diagmats(1:nequ,1)

       !       
       ! solve U^{T} x = y
       !
       ! inverse_diagonal = inverse of diag(U)
       call this%spmats(1)%solve_triangular(&
             error_unit,&
            this%aux_apply%raux(1:nequ),&
            vec_out)
       info_prec=0
       
    case( 'inv_ILU' ) 
       !
       ! 
       ! vec_out = ( L U )^{-1} vec_in
       !           ==
       ! solve L U (vec_out)  = vec_in
       !          
       !
       ! compute y = L^{-1} vec_in  
       !          ==
       ! solve L y = vec_in
       !
       ! ( L has only one on the diagonals)
       call this%spmats(1)%solve_triangular_unitdiag(&
            error_unit,&
            vec_in,&
            this%aux_apply%raux(1:nequ))
       !
       
       !
       ! compute vec_out = U^{-1} y = U^{-1} L^{-1} vec_in  
       !          ==
       ! solve U vec_out = y = L^{-1} vec_in 
       !
       ! (use the inverse of the diagonal of U in this%diagmats)
       call this%spmats(2)%solve_triangular(&
            error_unit,&
            this%aux_apply%raux(1:nequ),&
            vec_out)!,&
            !inverse_diagonal=this%diagmats(1:nequ,1))  
       
       info_prec=0
       
    case ('inv_L')
       !
       ! 
       ! vec_out = ( L )^{-1} vec_in
       !           ==
       ! solve L (vec_out)  = vec_in
       !          
       ! ( L has only one on the  diagonals)
       call this%spmats(1)%solve_triangular_unitdiag(&
            error_unit,&
            vec_in,&
            vec_out)
      ! call this%spmats(1)%info(6)
       info_prec=0
    case ('inv_U')
       !
       ! compute vec_out = U^{-1} vec_in  
       !          ==
       ! solve U vec_out  = vec_in
       !
       ! solve diag(U)^{-1} U y = vec_in
       !call this%spmats(1)%info(6)
       call this%spmats(1)%solve_triangular(&
            error_unit,&
            vec_in,&
            vec_out,&
            inverse_diagonal=this%diagmats(1:nequ,1))
       info_prec=0
       

    end select

    info = info_prec
    
  end subroutine Mxv_stdprec

  subroutine init_stdprecs(this,lun_err,ctrl_split,these)
    use Globals
    implicit none
    class(stdprec),   intent(in   ) :: this
    integer,          intent(in   ) :: lun_err
    character(len=*), intent(in   ) :: ctrl_split
    type(stdprec),    intent(inout) :: these(2)
    ! local
    logical :: rc
    integer :: res
    integer :: nequ
    integer :: info
    type(input_prec) :: ctrl_loc
    
    nequ = this%nrow
    
    select case (this%label)
    case('identity') 
       select type(this)
       type is (stdprec)
          these(1) = this
          these(2) = this
       end select
    case('inv_diag') 
       select type(this)
       type is (stdprec)
          these(1) = this
       end select
       these(1)%diagmats=sqrt(these(1)%diagmats)
       these(1)%label = 'inv_sqrtdiag'
       these(2) = these(1)
    case('inv_IC')
       rc = IOerr(lun_err, err_inp,&
            'init_stdprecs', &
            ' requested precs non valid in symmetric mode: IC') 
    case('inv_ILU')
       if ( this%nspmats .eq. 1) &
            rc = IOerr(lun_err, err_inp,&
            'init_stdprecs', &
            ' requested precs non valid in symmetric mode: ILU')  
       
       select case ( ctrl_split )
       case ('left') 
          ! P1=(LU)^{-1}
          select type(this)
          type is (stdprec)
             these(1)       = this       
          end select
          these(1)%is_symmetric = .false.
          

          ! P2= Id
          call ctrl_loc%init(lun_err,'identity')
          call these(2)%init(lun_err, info, ctrl_loc,nequ)

       case ('right')
          ! P1 = Id
          call ctrl_loc%init(lun_err,'identity') 
          call these(1)%init(lun_err, info, ctrl_loc,nequ)

          ! P2 = (LU)^{-1}
          select type(this)
          type is (stdprec)
             these(2)       = this
          end select
          these(2)%is_symmetric = .false.

       case('split')
          !
          ! both factors are computed 
          !

          !
          ! allocation and assignment for these(1)
          !
          these(1)%label         = 'inv_L'
          these(1)%nrow          = this%nrow
          these(1)%ncol          = this%ncol
          these(1)%is_symmetric  = .false.
          !
          these(1)%nspmats   = 1
          allocate(these(1)%spmats(these(1)%nspmats),stat = res)
          if(res .ne. 0) rc = IOerr(lun_err, err_alloc,&
               'init_stdprec', &
               ' type prec member spmats',res) 
          !
          ! inverse of the diagonal not required
          ! aux not required

          !
          ! copy of L 
          !
          these(1)%spmats(1) = this%spmats(1)

          !
          ! allocation and assignment for these(2)
          !
          these(2)%label = 'inv_U'
          these(2)%nrow          = this%nrow
          these(2)%ncol          = this%ncol
          these(2)%is_symmetric  = .false.
          !
          these(2)%nspmats   = 1
          allocate(these(2)%spmats(these(2)%nspmats),stat = res)
          if(res .ne. 0) rc = IOerr(lun_err, err_alloc,&
               'init_stdprec', &
               ' type prec member spmats',res) 
          !
          these(2)%ndiagmats = 1
          allocate(these(2)%diagmats(nequ,these(2)%ndiagmats),stat = res)
          if(res .ne. 0) rc = IOerr(lun_err, err_alloc,&
               'init_stdprec', &
               ' type prec member spmats',res)  

          ! aux not required

          !
          ! copy of U and its diagonal
          !
          these(2)%spmats(1) = this%spmats(2) ! U factor
          these(2)%diagmats  = this%diagmats  ! diag(U)^{-1}

       end select
    end select

    these(1)%is_initialized = .true.
    these(2)%is_initialized = .true.
    these(1)%ctrl = this%ctrl
    these(2)%ctrl = this%ctrl
    these(1)%nrow          = this%nrow
    these(1)%ncol          = this%ncol
    these(1)%is_symmetric  = .false.
    these(2)%nrow          = this%nrow
    these(2)%ncol          = this%ncol
    these(2)%is_symmetric  = .false.

  end subroutine init_stdprecs

  




end module StdSparsePrec

