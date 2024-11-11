module dagmg_wrap
  use Globals
  use LinearOperator
  use SparseMatrix
  use Scratch
  use LinearSolver
  implicit none
  private
  !> Derived type for the definition of the 
  !> inverse of sparse matrix via
  !> algebraic multigrid solver
  !> AGgregation-based Multigrid solver
  type, extends(abs_linop),public :: agmg_inv
     !> Controls for AGMG linear solver
     type(input_solver) :: ctrl_solver
     !> Information on linear solver application
     type(output_solver) :: info_solver
     !> Pointer to the original matrix
     class(abs_linop) , pointer :: original_matrix=>null()
     !> Copy of the sparse matrix ( agmg will modified its structure)
     !> after initialization
     type(spmat)  :: mat
     !> Scratch array for temp array 
     type(scrt) :: aux
     !> Array for storing several info_solver
     integer :: napplication=0
     integer :: max_application=0
     type(output_solver), allocatable :: history(:)
   contains
     !> Procedure to set the inverse application 
     !> Setting tolerance, maxium number of iterations
     !> and the matrix to be inverted
     !> (public procedure for type spmat)
     procedure, public, pass :: init => init_agmg_inv
     !> Static descructor
     procedure, public, pass :: kill => kill_agmg_inv
     !> Procedure to compute 
     !>         y = M * x 
     !> (public procedure for type spmat)
     procedure, public, pass :: matrix_times_vector => Mxv_agmg_inv
    
  end type agmg_inv
contains
  !>--------------------------------------------------------------
  !> Statics cosntructor. Given a sparse matrix A defines
  !> an approximate linear operator acting as A^{-1}
  !> using the AGMG multigrid solver
  !> (public procedure for class abs_linop)
  !> 
  !> usage:
  !>     call 'var'%kill(lun_err,matrix,[ctrl_solver[])
  !>
  !> where 
  !> \param[in    ]       lun_err   -> integer. Info number in case of error
  !> \param[in    ]       matrix    -> type(spmat) Matrix A to be "inverted"
  !> \param[in,opt] ctrl_solver     -> type(input_solver) AGMG solver controls
  !> \param[in,opt] max_application -> Integer. Number of apllications stored
  !<--------------------------------------------------------------
  subroutine init_agmg_inv(this,lun_err,matrix,ctrl_solver,max_application)
    implicit none
    class(agmg_inv),              intent(inout) :: this
    integer,                      intent(in   ) :: lun_err
    type(spmat),target,           intent(in   ) :: matrix
    type(input_solver), optional, intent(in   ) :: ctrl_solver
    integer,            optional, intent(in   ) :: max_application
    !local
    integer :: iprint
    
    !
    ! set properties linear operator
    !
    this%nrow         = matrix%nrow
    this%ncol         = matrix%ncol
    this%is_symmetric = matrix%is_symmetric
    this%is_initialized = .True.

    
    !
    ! create local copy of sparse matrix
    !
    this%mat = matrix
    this%original_matrix => matrix

    !
    ! set work arrays
    !    
    call this%aux%init(lun_err,0,2*matrix%ncol)
    
    !
    ! set control agmg via 
    ! 
    if ( present(ctrl_solver) ) then
       this%ctrl_solver = ctrl_solver   
    end if

    
    
    if (this%ctrl_solver%iprt .eq. 0) iprint = 0
    if (this%ctrl_solver%iprt .lt. 0) iprint = -1
    if (this%ctrl_solver%iprt == 1) iprint = 0 ! print only after 
    if (this%ctrl_solver%iprt > 1) iprint = this%ctrl_solver%lun_out

    call dagmg(&
         this%mat%ncol, & ! n
         this%mat%coeff,   & ! n
         this%mat%ja,   & ! n
         this%mat%ia,   & ! n
         this%aux%raux(1:this%nrow),        & ! f
         this%aux%raux(1+this%nrow:2*this%nrow),       & ! x
         1,             & !ijob setup 
         -1,  & ! iprint
         1,& !nrestart or CG if nrestart=1
         0,& ! iter
         this%ctrl_solver%tol_sol) ! tol

    if (present(max_application)) then
       this%max_application=max_application
       allocate(this%history(max_application))
    end if

    



  end subroutine init_agmg_inv

  !>--------------------------------------------------------------
  !> Statics destructor
  !> (public procedure for class abs_linop)
  !> 
  !> usage:
  !>     call 'var'%kill(lun_err)
  !>
  !> where 
  !> \param[in   ] lun_err -> integer. Info number in case of error 
  !<--------------------------------------------------------------
  subroutine kill_agmg_inv(this,lun_err)
    implicit none
    class(agmg_inv),intent(inout) :: this
    integer,        intent(in   ) :: lun_err
    

    if ( this%is_initialized ) then
       call dagmg(&
            this%mat%ncol, & ! n
            this%mat%coeff,   & ! n
            this%mat%ja,   & ! n
            this%mat%ia,   & ! n
            this%mat%coeff,        & ! f
            this%mat%coeff,       & ! x
            -1,             & !ijob free memory
            -1,  & ! iprint
            1,& !nrestart or CG if nrestart=1
            0,& ! iter
            this%ctrl_solver%tol_sol) ! tol

       !
       ! free memory
       !
       call this%mat%kill(lun_err)
       this%original_matrix => null()
       call this%aux%kill(lun_err)
       if (allocated(this%history)) then
          deallocate(this%history)
       end if
       this%napplication=0

       !
       ! reset to default
       !
       call this%to_default()
    end if
  end subroutine kill_agmg_inv

  !>-------------------------------------------------------------
  !> Abstract procedure defining the interface for a general
  !> matrix - vector multiplication
  !>         vec_out = (M) times (vec_in)
  !>
  !> For this tpye M~=A^{-1} where A is given in init procedure.
  !> The operator is only approximately linear if
  !> rough tolerance are given.
  !> (public procedure for class abs_linop)
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
  subroutine Mxv_agmg_inv(this,vec_in,vec_out, info,lun_err)
    use Globals
    implicit none
    class(agmg_inv),   intent(inout) :: this
    real(kind=double), intent(in   ) :: vec_in(this%ncol)
    real(kind=double), intent(inout) :: vec_out(this%nrow)
    integer,           intent(inout) :: info
    integer,           intent(in   ) :: lun_err
    ! local
    integer :: iter,nrest,ijob
    real(kind=double), allocatable :: work(:)
    real(kind=double) :: dnrm2,rhs_norm
    integer :: iprint
    
    !
    ! count applications
    !
    this%napplication=this%napplication+1
    
    call this%info_solver%init()
    call this%info_solver%tot%set('start')

    iter = this%ctrl_solver%imax
    if (this%ctrl_solver%iprt .eq. 0) iprint = 0
    if (this%ctrl_solver%iprt .lt. 0) iprint = -1
    if (this%ctrl_solver%iprt == 1) iprint = 0 ! print only after 
    if (this%ctrl_solver%iprt > 1) iprint = this%ctrl_solver%lun_out
    

    !
    ! compute rhs norm
    !
    rhs_norm=dnrm2(this%ncol, vec_in,1)

    !
    ! compute initial residual
    !
    if (this%ctrl_solver%compute_initial_residue .eq. 1) then
       call this%original_matrix%Mxv(vec_out, this%aux%raux(1:this%ncol),info,lun_err)
       this%aux%raux(1:this%ncol)=vec_in-this%aux%raux(1:this%ncol)
       this%info_solver%resini=dnrm2(this%ncol, this%aux%raux(1:this%ncol),1)
       if (this%ctrl_solver%iexit .eq. 1) this%info_solver%resini=this%info_solver%resini/rhs_norm
    end if

    !
    ! copy vec_in (agmg use it in inout mode)
    !
    this%aux%raux(1:this%ncol)=vec_in
    if (this%mat%is_symmetric) then
       nrest=1
    else
       nrest=this%ctrl_solver%nrestart
    end if
    ijob=2
    if (this%ctrl_solver%iexit .eq. 0) then 
       this%ctrl_solver%tol_sol=this%ctrl_solver%tol_sol/rhs_norm
    end if

    call dagmg(&
         this%mat%ncol, & ! n
         this%mat%coeff,   & ! n
         this%mat%ja,   & ! n
         this%mat%ia,   & ! n
         this%aux%raux(1:this%ncol),        & ! f
         vec_out,       & ! x
         ijob,             & !ijob solve after setup
         iprint,  & ! iprint
         nrest,& !nrestart or CG if nrestart=1
         iter,& ! iter
         this%ctrl_solver%tol_sol) ! tol
    this%info_solver%scheme_used='AGMG'

    if ( iter < 0) then
       if ( this%ctrl_solver%iprt<0) then
          this%info_solver%ierr = 0
       else
          this%info_solver%ierr = 1
       end if
       this%info_solver%iter = -iter
    else
       this%info_solver%ierr = 0
       this%info_solver%iter = iter
    end if

    
    !
    ! compute real residue
    !
    if ( this%ctrl_solver%compute_real_residue .eq. 1) then
       call this%original_matrix%Mxv(vec_out,this%aux%raux(1:this%ncol),info,lun_err)
       this%aux%raux(1:this%ncol)=this%aux%raux(1:this%ncol)-vec_in
       this%info_solver%resreal=dnrm2(this%ncol, this%aux%raux(1:this%ncol),1)
       if (this%ctrl_solver%iexit .eq. 1) this%info_solver%resreal=this%info_solver%resreal/rhs_norm
    end if

    if (allocated(this%history)) then
       this%history(mod(this%napplication,this%max_application))=this%info_solver
    end if

    info=this%info_solver%ierr
    call this%info_solver%tot%set('stop')
    if (this%ctrl_solver%iprt == 1) call this%info_solver%info(this%ctrl_solver%lun_out)

    ! set info to zero. Used when agmg is used as preconditioner
    if (this%ctrl_solver%ignore_errors) info=0
    
  end subroutine Mxv_agmg_inv

end module Dagmg_wrap
