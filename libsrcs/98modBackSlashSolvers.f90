module BackSlashSolver
  use Globals
  use LinearSolver
  use LinearOperator
  use Matrix
  use SimpleMatrix
  use SparseMatrix
  use StdSparsePrec
  use dagmg_wrap
  implicit none
  !> The goal is being able to solve a linear system
  !> Ax=b writing x=b.\.A., like the Matlab backslash.
  !> Thus, this module provides shortcuts
  !> to the differrent linear solver, involving
  !> different type of matrices.
  private
  type, public, extends(abs_linop) :: inv
     !> Label for inverse action
     character(len=50) :: approach='UNDEFINED'
     !> Pointer to matrix
     class(abs_linop), pointer :: matrixA
     !> controls of outer solver
     type(input_solver) :: ctrl_solver
     !> controls of inner solver (precondioner)
     type(input_prec) :: ctrl_prec
     !> Info od application
     type(output_solver) :: info_solver
     !> Comulative iterations:
     integer :: cumulative_iterations=0
     !> Comulative applications:
     integer :: cumulative_applications=0
     !> Inverse via Iterative solver
     type(inverse) :: iterative_solver
     !> Preconditoners for Iterative solver
     type(stdprec) :: prec_left
     type(stdprec) :: prec_right
     !> Multigrid solver
     type(agmg_inv) :: multigrid_solver
     !> Diagonal inverse
     type(diagmat) :: invdiag
   contains
     !> static constructor 
     !> (procedure public for type inv)
     procedure, public, pass :: init => init_inv
     !> static destructor
     !> (procedure public for type inv)
     procedure, public, pass :: kill => kill_inv
     !> static constructor 
     !> (procedure public for type inv)
     procedure, public, pass :: info => info_inv
     !> Procedure to compute 
     !>         y = M * x 
     !> (public procedure for type inv)
     procedure, public, pass :: matrix_times_vector => Mxv_inv
  end type inv

!!$  type, public, extends(abs_linop) :: multigrid
!!$     !> Label for initialization
!!$     logical :: is_initialized = .False.
!!$     !> integer
!!$     integer :: nlevel
!!$     !> Matrices
!!$     type(array_linop), allocatable :: A_matrices (:)
!!$     !> Smoothers
!!$     type(array_linop), allocatable :: two_level_multigrids (:)
!!$     !> Solvers
!!$     type(array_linop), allocatable :: smoothers (:)
!!$     !> Resctricitons
!!$     type(array_linop), allocatable :: restrictions (:)
!!$     !> Interpoaltion
!!$     type(array_linop), allocatable :: restrictions (:)
!!$     !> Local Scratch varaible
!!$     type(scrt) :: aux
!!$   contains
!!$     !> static constructor 
!!$     !> (procedure public for type inv)
!!$     procedure, public, pass :: init => init_multigrid
!!$     !> static destructor
!!$     !> (procedure public for type inv)
!!$     procedure, public, pass :: kill => kill_multigrid
!!$     !> Procedure to compute 
!!$     !>         y = M * x 
!!$     !> (public procedure for type inv)
!!$     procedure, public, pass :: Mxv=> Mxv_multigrid
!!$  end type multigrid
!!$
!!$   type, public, extends(abs_linop) :: twolevel_multigrid
!!$     !> Label for initialization
!!$     logical :: is_initialized = .False.
!!$     !> Matrix finer level
!!$     class(abs_linop), pointer :: matrixA_fineh
!!$     !> Matrix coarser level
!!$     class(abs_linop), pointer :: matrixA_coarse2h
!!$     !> Matrix coarser level
!!$     class(abs_linop), pointer :: restriction_fineh_to_coarse2h
!!$     !> Matrix coarser level
!!$     class(abs_linop), pointer :: interpolation_coarse2h_to_fineh
!!$     !> Matrix coarser level
!!$     class(abs_linop), pointer :: smootherh
!!$     !> Matrix coarser level
!!$     class(abs_linop), pointer :: solver_coarseh
!!$     !> Local Scratch varaible
!!$     type(scrt) :: aux
!!$   contains
!!$     !> static constructor 
!!$     !> (procedure public for type inv)
!!$     procedure, public, pass :: init => init_twolevel_multigrid
!!$     !> static destructor
!!$     !> (procedure public for type inv)
!!$     procedure, public, pass :: kill => kill_twolevel_multigrid
!!$     !> Procedure to compute 
!!$     !>         y = M * x 
!!$     !> (public procedure for type inv)
!!$     procedure, public, pass :: Mxv=> Mxv_twolevel_multigrid
!!$  end type twolevel_multigrid
!!$  
  
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
  !> \param[in] lun_err               -> integer. Error logical unit
  !> \param[in] nrow                  -> integer. Number of rows 
  !> \param[in] nrow                  -> integer. Number of columns
  !> \param[in] nterm                 -> integer. Number of non-zero term
  !>                                     stored in the matrix
  !> \param[in] (optional) is_sym     -> Logical. T/F flag for symmetric matrix
  !> \param[in] (optional) triangular -> Character. 
  !>                                       'N' = not triangular ( the default)
  !>                                       'U' = uppertriangular
  !>                                       'L' = lower_triangular' 
  !<-------------------------------------------------------------
  subroutine init_inv(this, &
       matrix,ctrl_solver,ctrl_prec,&
       ortogonalization_matrix,lun)
    use Globals
    implicit none
    !var
    class(inv),                            intent(inout) :: this
    type(spmat),                   target, intent(in   ) :: matrix
    class(input_solver), optional, target, intent(in   ) :: ctrl_solver
    type(input_prec),    optional, target, intent(in   ) :: ctrl_prec
    class(abs_linop),    optional, target, intent(in   ) :: ortogonalization_matrix
    integer, optional, intent(in   ) :: lun
    !local
    logical :: rc
    integer :: lun_err
    integer :: info_prec
    type(input_prec) :: ctrl_loc
    
    !
    if (present(lun)) then
       lun_err = lun
    else
       lun_err = 0
    end if

    !
    ! free memory 
    !
    if (this%is_initialized) call this%kill(lun_err)

    !
    ! set properteis of linear operator
    !
    this%nrow=matrix%nrow
    this%ncol=matrix%ncol
    this%is_symmetric=matrix%is_symmetric
    this%triangular = 'N'
    
    this%matrixA=> matrix

    
    
    
    !
    ! set controls 
    !
    if ( present(ctrl_solver) ) then
       this%approach=ctrl_solver%approach ! ITERATIVE or AGMG
       !
       ! copy locally the outer solver
       !
       this%ctrl_solver=ctrl_solver
       if ( present (ctrl_prec ) )then
          !
          ! copy
          !
          this%ctrl_prec=ctrl_prec
       else
          !
          ! set default prec. with ITERATIVE + incomplete factorization
          !
          if ( ctrl_solver%approach .eq. 'ITERATIVE') then
             if ( matrix%is_symmetric ) then
                this%ctrl_prec%prec_type='IC'
             else
                this%ctrl_prec%prec_type='ILU'
             end if
             ! double the number of non zeros in the preconditonier
             this%ctrl_prec%n_fillin=2*(matrix%nterm/matrix%nrow)
             this%ctrl_prec%tol_fillin=1.0d-4
          end if
       end if
    else
       if ( present ( ctrl_prec ) )then
          !
          ! default :: solver is multigrid
          !
          this%approach='STDPREC'
          this%ctrl_prec=ctrl_prec
       else
          this%ctrl_solver%approach='AGMG'
          this%approach='AGMG'
       end if
    end if

    !
    ! default solver is Multigrid
    !
    select case(this%approach)
    case ('AGMG')
       call this%multigrid_solver%init(lun_err,matrix,this%ctrl_solver)
    case ('ITERATIVE')
       call this%iterative_solver%init(lun_err,matrix%nrow)
       call this%prec_left%init(lun_err, info_prec,&
            this%ctrl_prec, matrix%nrow, matrix)
       if (info_prec .ne. 0) then
          rc = IOerr(lun_err, wrn_inp, 'init_inv', &
               ' error build preconditioner, info=',info_prec)
       end if
       call this%iterative_solver%set(matrix,&
            ctrl_solver=this%ctrl_solver,&
            prec_left=this%prec_left,&
            ortogonalization_matrix=ortogonalization_matrix)
    case ('STDPREC')
       call this%prec_left%init(lun_err, info_prec,&
            this%ctrl_prec, matrix%nrow, matrix)
       if (info_prec .ne. 0) then
          rc = IOerr(lun_err, wrn_inp, 'init_inv', &
               ' error build preconditioner, info=',info_prec)
       end if

    end select
    this%cumulative_iterations = 0
    this%cumulative_applications = 0
    this%is_initialized = .True.
    

  end subroutine init_inv

  subroutine info_inv(this,lun)
    implicit none
    class(inv),  intent(in   ) :: this
    integer,           intent(in   ) :: lun
    character(len=256) :: str1,str2
    !local
    integer :: i,j
    character(len=3) :: str

    write(lun,*) 'APPROXIMATE INVERSE FOR MATRIX :'
    call this%matrixA%info(lun)
    write(lun,*) 'APPROACH:',etb(this%approach)
    if ( this%approach .eq. 'STDPREC') then
       call this%ctrl_prec%info(lun)
    else
       call this%ctrl_solver%info(lun)
    end if
    
  end subroutine info_inv
  

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
  !> \param[in   ] info -> integer. Info number
  !>                                  in case of error   
  !<-------------------------------------------------------------
  subroutine Mxv_inv(this,vec_in,vec_out, info,lun_err)
    use Globals
    implicit none
    class(inv),      intent(inout) :: this
    real(kind=double), intent(in   ) :: vec_in(this%ncol)
    real(kind=double), intent(inout) :: vec_out(this%nrow)
    integer,           intent(inout) :: info
    integer,           intent(in   ) :: lun_err

    select case (this%approach)
    case ('AGMG')
       call this%multigrid_solver%Mxv(vec_in,vec_out, info,lun_err)
       this%info_solver=this%multigrid_solver%info_solver
       this%cumulative_iterations = this%cumulative_iterations +&
            this%info_solver%iter
       this%cumulative_applications = this%cumulative_applications +1
    case ('ITERATIVE')
       !
       ! reset controls in case ued change it
       !
       this%iterative_solver%ctrl_solver=this%ctrl_solver
       !
       ! apply iterative method to solve A vec_out = vec_in
       !
       call  this%iterative_solver%Mxv(vec_in,vec_out, info,lun_err)

       !
       ! store information
       !
       this%info_solver=this%iterative_solver%info_solver
       this%cumulative_iterations = this%cumulative_iterations +&
            this%info_solver%iter
       this%cumulative_applications = this%cumulative_applications +1
    case ('STDPREC')
       call  this%prec_left%Mxv(vec_in,vec_out, info,lun_err)
    end select

  end subroutine Mxv_Inv

  subroutine kill_inv(this, lun_err)
    implicit none
    ! vars
    class(inv),intent(inout) :: this
    integer,     intent(in   ) :: lun_err
    ! local vars
    integer :: res
    logical :: rc

    !
    ! Free memory
    !
    
    this%matrixA=> null()
    call this%iterative_solver%kill(lun_err)
    call this%prec_left%kill(lun_err)
    call this%multigrid_solver%kill(lun_err)
    call this%invdiag%kill(lun_err)

    !
    ! reset to default
    !
    call this%to_default()
    

  end subroutine kill_inv

!!$  !>-------------------------------------------------------------
!!$  !> Static constructor.
!!$  !> (procedure public for type spmat)
!!$  !> Instantiate (allocate if necessary)
!!$  !> and initilize (by also reading from input file)
!!$  !> variable of type spmat
!!$  !>
!!$  !> usage:
!!$  !>     call 'var'%init(lun_err, nrow, nterm )
!!$  !>
!!$  !> where:
!!$  !> \param[in] lun_err               -> integer. Error logical unit
!!$  !> \param[in] nrow                  -> integer. Number of rows 
!!$  !> \param[in] nrow                  -> integer. Number of columns
!!$  !> \param[in] nterm                 -> integer. Number of non-zero term
!!$  !>                                     stored in the matrix
!!$  !> \param[in] (optional) is_sym     -> Logical. T/F flag for symmetric matrix
!!$  !> \param[in] (optional) triangular -> Character. 
!!$  !>                                       'N' = not triangular ( the default)
!!$  !>                                       'U' = uppertriangular
!!$  !>                                       'L' = lower_triangular' 
!!$  !<-------------------------------------------------------------
!!$  subroutine init_two_level_multigrid(this, &
!!$       info,lun_err,matrixA_fine,&
!!$       restriction_fine_to_coarse2,interpolation_coarse_to_fine,&
!!$       smoother,solver_coarse,aux)
!!$    class(two_level_multigrid), intent(inout) :: this
!!$    integer,                    intent(inout) :: info
!!$    integer,                    intent(in   ) :: lun_err
!!$    class(abs_linop), target,   intent(in   ) :: matrixA_fine
!!$    class(abs_linop), target,   intent(in   ) :: restriction_fine_to_coarse2
!!$    class(abs_linop), target,   intent(in   ) :: interpolation_coarse_to_fine
!!$    class(abs_linop), target,   intent(in   ) :: smoother
!!$    class(abs_linop), target,   intent(in   ) :: solver_coarse
!!$    type(scrt), target,optional,intent(inout) :: aux  
!!$    ! local
!!$    integer :: nh,n2h
!!$    type(scrt), pointer :: aux_loc
!!$
!!$    !
!!$    ! check dimension
!!$    !
!!$    nh=matrixA_fine%nrow
!!$    n2h=restriction_fineh%ncol
!!$
!!$    !
!!$    ! assign pointers
!!$    !
!!$    this%matrixA_fineh0=>matrixA_fineh
!!$    this%restriction_fineh_to_coarse2h=>restriction_fineh_to_coarse2h
!!$    this%interpolation_coarse2h_to_fineh=>interpolation_coarse2h_to_fineh
!!$    this%smootherh=>smootherh
!!$    this%solver_coarseh=>solver_coarseh
!!$
!!$
!!$    !
!!$    ! prepare work space
!!$    !
!!$    if (present(aux)) then
!!$       if ( present(aux) ) then
!!$          if( .not. aux%check(0,9*nequ) ) &
!!$               rc = IOerr(lun_err, wrn_inp, 'minres_solver', &
!!$               ' aux array too small',res)
!!$          aux_loc => aux_minres
!!$       else
!!$          call this%aux%init(lun_err,0,9*nequ)
!!$          aux_loc => this%aux
!!$       end if
!!$    end if
!!$    iend=0
!!$    call aux_loc%range(nh,ibegin,iend)
!!$    sol_h =>aux_loc%raux(ibegin:iend)
!!$    call aux_loc%range(nh,ibegin,iend)
!!$    scr_fine =>aux_loc%raux(ibegin:iend)
!!$    call aux_loc%range(nh,ibegin,iend)
!!$    err_fine =>aux_loc%raux(ibegin:iend)
!!$    call aux_loc%range(n2h,ibegin,iend)
!!$    scr_coarse =>aux_loc%raux(ibegin:iend)
!!$
!!$
!!$
!!$  end subroutine init_two_level_multigrid
!!$
!!$  subroutine kill_two_level_multigrid(this, lun_err)
!!$    class(two_level_multigrid), intent(inout) :: this
!!$    integer,                    intent(in   ) :: lun_err
!!$   
!!$    !
!!$    ! deassociate pointers
!!$    !
!!$    this%matrixA_fine => null()
!!$    this%restriction_fine_to_coarse => null()
!!$    this%interpolation_coarse_to_fine => null()
!!$    this%smoother => null()
!!$    this%solver_coarse => null()
!!$
!!$    !
!!$    ! free memory
!!$    !
!!$    if ( this%aux%is_initialized) call this%aux%kill(lun_err)
!!$    
!!$  end subroutine kill_two_level_multigrid
!!$
!!$  !>-------------------------------------------------------------
!!$  !> Procedure to compute Matrix vector product
!!$  !>         vec_out = (M) times (vec_in)
!!$  !> (public procedure for type spmat)
!!$  !> 
!!$  !> usage:
!!$  !>     call 'var'%Mxv(vec_in,vec_out,[info])
!!$  !>
!!$  !> where 
!!$  !> \param[in   ] vec_in          -> real. dimension('var'%ncol)
!!$  !>                                  vector to be multiplied
!!$  !> \param[inout] vec_out         -> real. dimension('var'%nrow)
!!$  !>                                  vector (M) times (vec_in) 
!!$  !> \param[in   ] (optional) info -> integer. Info number
!!$  !>                                  in case of error   
!!$  !<-------------------------------------------------------------
!!$  recursive subroutine Mxv_two_level_multigrid(this,vec_in,vec_out, info,lun_err)
!!$    use Globals
!!$    implicit none
!!$    class(inv),      intent(inout) :: this
!!$    real(kind=double), intent(in   ) :: vec_in(this%ncol)
!!$    real(kind=double), intent(inout) :: vec_out(this%nrow)
!!$    integer, optional, intent(inout) :: info
!!$    integer, optional, intent(in   ) :: lun_err
!!$
!!$    !
!!$    ! compute u_h
!!$    !
!!$    call this%smoother%mxv(vec_in,this%sol_h,info,lun_err)
!!$
!!$    !
!!$    ! compute residuum r_h=b-A_h*u_n
!!$    !
!!$    call this%matrixA_fine%Mxv(this%sol_h,this%scr_fine,info,lun_err)
!!$    this%scr_fine = vec_in - this%scr_fine
!!$
!!$    !
!!$    ! project to coarser grid r_2h = R r_h
!!$    !
!!$    call this%restriction_fineh_to_coarse2h%Mxv(this%scr_fine,&
!!$         this%scr_coarse,info,lun_err)
!!$
!!$    !
!!$    ! solve at the coarser level
!!$    !
!!$    ! A_2h err_2h = r_2h
!!$    !
!!$    call this%solver_coarse%Mxv(this%scr_coarse,this%scr_fine(1:n2h),info,lun_err)
!!$
!!$    !
!!$    ! project to finer grid  err_h = I err_2h
!!$    !
!!$    call this%interpolation_coarse_to_fine%Mxv(this%scr_fine(1:n2h),&
!!$         this%err_fine,info,lun_err)
!!$
!!$    !
!!$    ! update u_h
!!$    !
!!$    this%sol_h=this%sol_h+this%err_fine
!!$    call this%smoother%mxv(this%sol_h,vec_out,info,lun_err)
!!$
!!$  end subroutine Mxv_two_level_multigrid
!!$
!!$
!!$  subroutine init_multigrid(this, &
!!$       info,lun_err,matrixA_fine,&
!!$       restriction_fine_to_coarse2,interpolation_coarse_to_fine,&
!!$       smoother,solver_coarse,aux)
!!$    class(two_level_multigrid), intent(inout) :: this
!!$    integer,                    intent(inout) :: info
!!$    integer,                    intent(in   ) :: lun_err
!!$    class(abs_linop), target,   intent(in   ) :: matrixA_fine
!!$    class(abs_linop), target,   intent(in   ) :: restriction_fine_to_coarse2
!!$
!!$
!!$  end subroutine init_multigrid
!!$
!!$  subroutine init_V_cycle(this, &
!!$       info,lun_err,&
!!$       nlevel,A_matrices,&
!!$       restrictions,interpolations,&
!!$       smoothers,solver_coarse,aux)
!!$    class(two_level_multigrid), intent(inout) :: this
!!$    integer,                    intent(inout) :: info
!!$    integer,                    intent(in   ) :: lun_err
!!$    class(abs_linop), target,   intent(in   ) :: A_matrices(nlevel)
!!$    class(abs_linop), target,   intent(in   ) :: restrictions(nlevel-1)
!!$    class(abs_linop), target,   intent(in   ) :: interpolations(nlevel-1)
!!$    class(abs_linop), target,   intent(in   ) :: smoothers(nlevel-1)
!!$    class(abs_linop), target,   intent(in   ) :: coarser_solver
!!$    
!!$    allocate( two_level_multigrid(nlevel) ) 
!!$    
!!$    do i=1,nlevel
!!$       call this%two_level_multigrid(i)%init(info,lun_err,&
!!$            A_matrices(i),&
!!$            restrictions(ilevel),&
!!$            interpolations(ilevel),&
!!$            smoothers(ilevel),&
!!$            this%two_level_multigrid(i-1),&
!!$            this%aux)
!!$    end do
!!$
!!$  end subroutine init_two_level_multigrid
  
end module BackSlashSolver
