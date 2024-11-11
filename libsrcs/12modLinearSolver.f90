module LinearSolver
    use Globals
    use Timing
    use LinearOperator
    use Matrix
    use Scratch
    !use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, stdout=>output_unit, stderr=>error_unit
    !>-----------------------------------------------------
    !> Library for the solution of Linear System 
    !>        A x = b
    !> with method Kryol.
    !>
    !> Method currently supported 
    !> PCG     : Preconditioned Conjugate Gradient 
    !> BICGSTAB: Bi-Conjugate Gradient Stabilized
    !> MINRES  : Minimal residual
    !> GMRES   : Generalized minimal residua
    !>
    !> When preconditioners are used we solve
    !>
    !>    Pl^{-1} A PR^{-1} PR x = PL^{-1} b
    !>
    !> where PL and PR, are left and right preconditioner
    !> respectively.
    !>
    !> The main procedure is called "linear_solver"
    !> and it will handle :
    !>  - any matrix A described any type extending
    !>    the abstract type "abs_linop", 
    !>    for which the matrix-vector multiplication must be defined;
    !>  - any left and right preconditioner extending the abstract type
    !>    "abs_linop"
    !> TODO: 
    !> 1-  check gmres  is wrong
    !> 2 - Define a derived type made only by procedures
    !>     for vector operation
    !>     ddot, dnorm, euclidian norm, daxpy, dxapy, sum, 
    !>     dscal, normres
    !>     (in particular normres to removed ndir, nnodir)
    !>  We defined a defualt type with the classical opration
    !>  used if not vector-vector operation type is passed.
    implicit none
    private
    public :: input_solver, output_solver,linear_solver,uzawa_arrow_hurwics
    !>-----------------------------------------------------
    !> Derived Type containing all control parameters
    !> for the solution of linear system
    !>   A x = b
    !> All values are set to a default value
    !>------------------------------------------------------
    type, public :: input_solver
       !> Solver approach. Supported
       !>  ITERATE (includes all solver that requires only matrix vector operation)
       !>  AGMG    (used for sparse matrix)
       !> TODO: include 
       !>  DIRECT
       character(len=20) :: approach='ITERATIVE'
       !> Solver scheme. Now supported
       !>  PCG
       !>  BICGSTAB
       !>  GMRES
       !>  MINRES
       !> TODO: include :
       !>  LSQR
       character(len=20) :: scheme='BICGSTAB'
       !> I/O err msg. unit
       integer :: lun_err=0!stderr
       !> I/O out msg. unit
       integer :: lun_out=6!stdout
       !> Integer identifying the Preconditinoer
       !>     iexit=0 => exit on absolute residual               
       !>     iexit=1 => exit on |r_k|/|b|     
       integer :: iexit=1
       !> Number of maximal iterations for iterative solver
       integer :: imax=1000
       !> Flag for printing option for itereative solver
       integer :: iprt=0
       !> Flag for starting solution
       !>  isol=0 => SOL      is the initial solution for pcg
       !>  isol=1 => PREC*RHS is the initial solution for pcg
       integer :: isol=0
       !> Tolerance for the scheme
       real(kind=double) :: tol_sol=1.0d-13
       !> Ortogonalized w.r.t. to given vectors 
       !> (usally the kernel of the matrix)
       integer :: iort=0
       !> Number of krylov directions saved in GMRES algorithm
       !> Used only for GMRES
       integer :: nrestart=20
       !> Flag to active to print steps for debugging
       !> debug=0 => no print 
       !> debug=1 => print steps
       integer :: debug=0
       !> Flag to compute real residual 
       !> (absolute or relative according to iexit)
       !> compute_real_residua = 0 => not do extra computation
       !> compute_real_residua = 1 => compute residua 
       integer :: compute_initial_residue=1
       !> Flag to compute real residual 
       !> (absolute or relative according to iexit)
       !> compute_real_residua = 0 => not do extra computation
       !> compute_real_residua = 1 => compute residua 
       integer :: compute_real_residue=1
       !> Flexible PCG controls. Number of restart
       !> See "Flexible Conjugate Gradients " YVAN NOTAY
       integer :: fpcg_max_restart =1
       !> Variable used to print info where
       !> linear solver defines a linear operator.
       !> Used by inverse type
       logical :: print_info_solver = .False.
       !> Ignore errrors. setup info = 0
       logical :: ignore_errors = .False.
     contains
       !> Static constructor
       !> (procedure public for type input_solver)
       procedure, public, pass :: init => init_input_solver
       !> static destructor
       !> (procedure public for type input_solver)
       procedure, public, pass :: kill => kill_input_solver
       !> Info procedure.
       !> (public procedure for type input_solver)
       procedure, public, pass :: info => info_input_solver
       !> Info procedure.
       !> (public procedure for type input_solver)
       procedure, public, pass :: info2str => info_input_solver_string
       !> Read from file procedure
       !> (procedure public for type input_solver)
       procedure, public, pass :: read => read_input_solver
    end type input_solver

   


    !>-----------------------------------------------------
    !> Derived Type containing all information parameters
    !> for iterative method used in linear_solver, like
    !> number of iterations required, residuals, timing, etc
    !>------------------------------------------------------
    type, public :: output_solver
       !> Integer flag for error
       !>   ierr=0  => solution achievied                
       !>   ierr=1  => iterations exceeded max number(for iterative solver)
       !>   ierr=2  => BICGSTAB rhopk1<small
       !>   ierr=3  => Precondtioner application error
       !>   ierr=-1 => BICGSTAB : NaN found
       integer :: ierr=-999
       !> Integer flag for preconditioenr application
       integer :: info_prec = -999
       !> Number of iterations required
       integer :: iter=-999
       !> Norm of initial residum vector \|A x_0 -b\|
       real(kind=double) :: resini=-huge
       !> Norm of final real residum \|A x -b\|
       real(kind=double) :: resreal=-huge 
       !> Norm of final residum vector
       real(kind=double) :: resnorm=-huge
       !> Euclidain norm of rhs term 
       !> (Dirichlet nodes are not taken into account)
       ! TODO this shold be removed and computed outside
       real(kind=double) :: bnorm=-huge
       !> Total Time solver
       type(Tim) :: tot
       !> Overhead time
       type(Tim) :: ovh
       !> Time in matrix vector operations
       type(Tim) :: matvec
       !> Time in vector vector operations
       type(Tim) :: vecvec
       !> Time in scalar product
       type(Tim) :: scalprod
       !> Time in preconditioner application
       type(Tim) :: prec
       !> Time in ortgonalization
       type(Tim) :: ort
       !> Integer flag for solver used
       character(len=20) :: scheme_used='NotDefined'
       !> Estimate of conditioning number ( MINRES only)
       real(kind=double) :: matrix_cond=huge
     contains
       !> static constructor
       !> (procedure public for type output_solver)
       procedure, public, pass :: init => init_output_solver
       !> static destructor
       !> (procedure public for type output_solver)
       procedure, public, pass :: kill => kill_output_solver
       !> Info procedure.
       !> (public procedure for type output_solver)
       procedure, public, pass :: info => info_output_solver
       !> Info procedure.
       !> (public procedure for type output_solver)
       procedure, public, pass :: time => infotiming_output_solver
       !> Info procedure.
       !> (public procedure for type output_solver)
       procedure, public, pass :: info2str => info_output_solver_2string
       !> Info procedure.
       !> (public procedure for type output_solver)
       procedure, public, pass :: time2str => str_timing_output_solver
    end type output_solver


    

   

    !>-----------------------------------------------------------
    !> Vector-Vector operation abstract type
    !> TODO : 
    !> scalar product
    !> eucliadian norm
    !> ....
    !>----------------------------------------------------------
    type, abstract :: abs_vec_operation
       !>  Vector operation abstract interface
     contains
       !>-----------------------------------------------------------
       !> Procedure for y=ax*y
       procedure(alpha_x_plus_y), deferred :: axpy
    end type abs_vec_operation
    abstract interface
       !>-------------------------------------------------------------
       !>
       !<-------------------------------------------------------------
       subroutine alpha_x_plus_y(this,nvec,alpha,xvec,xinc,yvec,yinc)
         use Globals
         import abs_vec_operation
         implicit none
         class(abs_vec_operation), intent(in   ) :: this
         integer,           intent(in   ) :: nvec
         real(kind=double), intent(in   ) :: alpha
         real(kind=double), intent(in   ) :: xvec(nvec)
         integer,           intent(in   ) :: xinc
         real(kind=double), intent(inout) :: yvec(nvec)
         integer,           intent(in   ) :: yinc
       end subroutine alpha_x_plus_y
    end interface

    

    !>-------------------------------------------------------------
    !> Structure variable for storing
    !> coefficints alpha beta presnorm 
    !> and vectors pres
    !> TODO: remove this structure and use a larger scratch varible to
    !> store additional informattions
    !<--------------------------------------------------------------
    type, public :: pcgcoeff
       !> Logical Check is type has been initialized
       logical :: is_initialized=.false.
       !> Array length of eigvector
       integer :: length_pres
       !> Number of preconditionted vector res
       integer :: npres
       !> Number of preconditionted vector res actually stored in WMAT
       integer :: npres_stored=0
       !> Dimension(npres)
       !> Real array with Ddiag
       real(kind=double), allocatable :: alpha(:)
       !> Dimension(npres)
       !> Real array of eigenvectors
       real(kind=double), allocatable :: beta(:)
       !> Dimension(npres)
       !> Real array of eigenvectors
       real(kind=double), allocatable :: pres_dot_res(:)
       !> Dimension(length_pres,npres)
       !> Real array vectro P^{-1} res
       real(kind=double), allocatable :: pres(:,:)
     contains
       !> static constructor
       !> (public for type eig)
       procedure, public, pass :: init => init_pcgcoeff
       !> static destructor
       !> (public for type eig)
       procedure, public, pass :: kill => kill_pcgcoeff
       !> info procedure
       !> (public for type eig)
       procedure, public, pass :: info => info_pcgcoeff
       !> Fill the type with the information of the PCG procedure
       procedure, public, pass :: fill => fill_pcgcoeff
    end type pcgcoeff

    !>-----------------------------------------------------
    !> Derived Type containing all control parameters
    !> for iterative method of linear system
    !>   A x = b
    !> with A abstract type abs_linop, continaing the
    !> matrix-vector.
    !> All values are set to default
    !>------------------------------------------------------
    type, extends(abs_linop), public :: inverse
       !> Copy and pointer to all variable
       !> called by linear_solver
       type(input_solver)  :: ctrl_solver
       type(output_solver) :: info_solver
       class(abs_linop), pointer :: matrix
       class(abs_linop), pointer :: prec_left
       class(abs_linop), pointer :: prec_right
       class(abs_linop), pointer :: ortogonalization_matrix
       type(scrt) :: aux
       integer :: total_iter = 0
       !type(
     contains
       !> Static constructor
       !> (procedure public for type inverse)
       procedure, public, pass :: init => init_inverse
       !> static destructor
       !> (procedure public for type inverse)
       procedure, public, pass :: kill => kill_inverse
       !> Info procedure
       !> (procedure public for type inverse)
       procedure, public, pass :: info => info_inverse
       !> Set-up procedure
       !> (public procedure for type inverse)
       procedure, public, pass :: set => set_inverse
       !> Application procedure.
       !> (public procedure for type input_solver)
       procedure, public, pass :: matrix_times_vector => apply_inverse      
    end type inverse

    

  contains
    !>-------------------------------------------------------------
    !> static constructor.
    !> (procedure public for type input_solver)
    !> Instantiate variable of type input_solver
    !> All parameters are optional and, if not passed, they
    !> are set to default values 
    !>
    !> usage:
    !>     call 'var'%init()
    !>
    !> where:
    !> \param[in] (optional) err_unit   -> integer. I/O error for not supported
    !>                                     parameters
    !> \param[in] (optional) scheme     -> string :: .Id. solver  used
    !> \param[in] (optional) lun_err    -> integer. I/O error  unit  
    !> \param[in] (optional) lun_out    -> integer. I/O output unit for 
    !> \param[in] (optional) imax       -> integer. Max nmb of iterations 
    !> \param[in] (optional) iprt       -> integer. Flag for printing info
    !>                         iprt = -1 ( no print, neither errors)
    !>                         iprt = 0 (no print)
    !>                         iprt = 1 (print only end)
    !>                         iprt = 2 (print at each iteration)
    !> \param[in] (optional) isol     -> integer. Flag for initial solution
    !>                         isol = 0 (start from sol)
    !>                         isol = 1 (start from P^-1 rhs)
    !> \param[in] (optional) tol_sol -> real. Tolerance to ve achieved
    !> TODO : make scheme not key sensitive
    !<-------------------------------------------------------------
    subroutine init_input_solver(this,&
         err_unit,&
         approach,&
         scheme,&
         lun_err,lun_out,&
         iexit,&
         imax,&
         iprt,&
         isol,&
         tol_sol,&
         iort,&
         nrestart,&
         fpcg_max_restart)
      use Globals
      implicit none
      class(input_solver),         intent(inout) :: this 
      integer,                     intent(in   ) :: err_unit
      character(len=*), optional,  intent(in   ) :: approach      
      character(len=*), optional,  intent(in   ) :: scheme
      integer, optional,           intent(in   ) :: lun_err 
      integer, optional,           intent(in   ) :: lun_out
      integer, optional,           intent(in   ) :: iexit
      integer, optional,           intent(in   ) :: imax
      integer, optional,           intent(in   ) :: iprt
      integer, optional,           intent(in   ) :: isol
      real(kind=double), optional, intent(in   ) :: tol_sol
      integer, optional,           intent(in   ) :: iort
      integer, optional,           intent(in   ) :: nrestart
      integer, optional,           intent(in   ) :: fpcg_max_restart
      ! local 
      logical :: rc
      character(len=256) :: work_str  
      
      if( present(approach ) ) then
         if ( ( etb(approach) .eq. 'ITERATIVE'     ) .or. &
              ( etb(approach) .eq. 'AGMG'  )       ) &
              then
            this%approach   = etb(approach)
         else
            rc = IOerr(err_unit, err_inp , 'init_input_solver', &
                 'approach = '//etb(approach)// 'not supported' )
         end if
      end if

      if( present(scheme  ) ) then
         if ( ( scheme .eq. 'PCG'     ) .or. &
              ( scheme .eq. 'FPCG'  ) .or. &
              ( scheme .eq. 'MINRES'  ) .or. &
              ( scheme .eq. 'BICGSTAB') .or. &
              ( scheme .eq. 'FGMRES') .or. &
              ( scheme .eq. 'GMRES') ) then
            this%scheme   = etb(scheme)
         else
            rc = IOerr(err_unit, err_inp , 'init_input_solver', &
                 'scheme = '//etb(scheme)// 'not supported' )
         end if
      end if            
      if( present(tol_sol ) ) this%tol_sol  = tol_sol
      if( present(lun_err ) ) this%lun_err  = lun_err
      if( present(lun_out ) ) this%lun_out  = lun_out
      if( present(iexit   ) ) this%iexit    = iexit
      if( present(imax    ) ) this%imax     = imax 
      if( present(iprt    ) ) this%iprt     = iprt
      if( present(isol    ) ) this%isol     = isol
      if( present(tol_sol ) ) this%tol_sol  = tol_sol
      if( present(iort    ) ) this%iort     = iort 
      if( present(nrestart) ) this%nrestart = nrestart
      if( present(fpcg_max_restart)) this%fpcg_max_restart = fpcg_max_restart 
      

    end subroutine init_input_solver

    !>-------------------------------------------------------------
    !> Read linear solver from file.
    !> (procedure public for type input_solver)
    !> Reset to default values all ctrls
    !>
    !> usage:
    !>     call 'var'%kill()
    !>
    !<-------------------------------------------------------------
    subroutine read_input_solver(this,lun_err,file_ctrl)
      use Globals
      implicit none
      class(input_solver), intent(inout) :: this 
      integer, intent(in) :: lun_err
      type(file), intent(in) :: file_ctrl
      !local
      integer:: lun,info
      character(len=256) :: fname,input
      
      lun = file_ctrl%lun
      fname = file_ctrl%fn
      
      ! linear solver method controls
      call find_nocomment(lun,lun_err,input,fname,'approach',info)
      if (info.eq.0) read(input,*) this%approach

      call find_nocomment(lun,lun_err,input,fname,'scheme',info)
      if (info.eq.0) read(input,*) this%scheme
      
      call find_nocomment(lun,lun_err,input,fname,'file_err',info)
      if (info.eq.0) read(input,*) this%lun_err

      call find_nocomment(lun,lun_err,input,fname,'file_out',info)
      if (info.eq.0) read(input,*) this%lun_out
 
      call find_nocomment(lun,lun_err,input,fname,'iexit',info)
      if (info.eq.0) read(input,*) this%iexit

      call find_nocomment(lun,lun_err,input,fname,'imax',info)
      if (info.eq.0) read(input,*) this%imax
 
      call find_nocomment(lun,lun_err,input,fname,'iprt',info)
      if (info.eq.0) read(input,*) this%iprt

      call find_nocomment(lun,lun_err,input,fname,'isol',info)
      if (info.eq.0) read(input,*) this%isol

      call find_nocomment(lun,lun_err,input,fname,'tol',info)
      if (info.eq.0) read(input,*) this%tol_sol
 
      call find_nocomment(lun,lun_err,input,fname,'iort',info)
      if (info.eq.0) read(input,*) this%iort

      
    contains
      subroutine find_nocomment(lun,lun_err,input,fname,var_name,info)
        use Globals
        integer,            intent(in   ) :: lun
        integer,            intent(in   ) :: lun_err
        character(len=256), intent(in   ) :: fname
        character(len=256), intent(inout) :: input
        character(len=*),   intent(in   ) :: var_name
        integer, optional,intent(inout) :: info
        !local
        logical :: rc
        integer :: res

        character(len=256) clean
        character(len=1) first
        logical :: found

        if ( present(info)) info=0

        found = .false.
        do while( .not. found ) 
           read(lun,'(a)',iostat = res) input
           if(res .ne. 0) then
              rc = IOerr(lun_err, wrn_inp , 'Ctrl_read', &
                   fname//'input text for member'//etb(var_name),res)
              if ( present(info)) info=-1
              return
           end if

           clean = etb(input)
           read(clean,*) first
           if ( ( first .eq. '!') .or. &
                ( first .eq. '%') .or. &
                ( first .eq. '#') ) then
              found=.false.
           else
              found=.true.
           end if
        end do
      end subroutine find_nocomment
    end subroutine read_input_solver

    

    !>-------------------------------------------------------------
    !> Static constructor.
    !> (procedure public for type input_solver)
    !> Reset to default values all ctrls
    !>
    !> usage:
    !>     call 'var'%kill()
    !>
    !<-------------------------------------------------------------
    subroutine kill_input_solver(this)
      use Globals
      implicit none
      class(input_solver), intent(inout) :: this 

      this%scheme     = 'BICGSTAB'
      this%lun_err    = 6
      this%lun_out    = 6
      this%iexit      = 1
      this%imax       = 1000 
      this%iprt       = 0
      this%isol       = 0 
      this%tol_sol    = 1.0d-13
      this%iort       = 0 
      this%nrestart   = 0

    end subroutine kill_input_solver

    !>-------------------------------------------------------------
    !> Info procedure
    !> (procedure public for type star)
    !>
    !> usage:
    !>     call 'var'%info()
    !>
    !> \param[in] lun -> integer. I/O logical unit
    !<-------------------------------------------------------------
    subroutine info_input_solver(this,lun)
      use Globals
      implicit none
      class(input_solver), intent(in) :: this 
      integer,             intent(in) :: lun
      !local
      character(len=256) :: solver_used 
      
      write(lun,'(a,1x,a,a,I2,a,I2,a,I1,a,I6,a,I1,a,1pe7.1,a,I2,a)') &
              etb(this%approach),etb(this%scheme),'(err,out=',&
              this%lun_err,',',this%lun_out,&
              ' / iexit=',this%iexit,&
              ' / imax=', this%imax,&
              ' / isol=', this%isol,&
              ' / tol=', this%tol_sol,&
              ' / iort=', this%iort,')'
    
    end subroutine info_input_solver

     !>-------------------------------------------------------------
    !> Info procedure
    !> (procedure public for type star)
    !>
    !> usage:
    !>     call 'var'%info()
    !>
    !> \param[in] lun -> integer. I/O logical unit
    !<-------------------------------------------------------------
    subroutine info_input_solver_string(this,string)
      use Globals
      implicit none
      class(input_solver), intent(in   ) :: this 
      character(len=256),  intent(inout) :: string

      !local
      character(len=3) :: str1
      character(len=30) ::str2
      string=' ' 

      if ( this%iexit .eq. 0 ) then
         str1='abs'
      else
         str1='rel'
      end if

      
      write(string,'(a,1x,a,a,a,a,1pe7.1,a,I4)') &
              etb(this%approach),etb(this%scheme),' ',etb(str1),' tol=' , this%tol_sol,&
              ' - imax=', this%imax

      str2=''
      if ( this%iort .ne. 0 ) write(str2,'(a,I2)') ' ort. freq.=',this%iort

      string=etb(etb(string)//etb(str2))
    
    end subroutine info_input_solver_string
    
    
    !>-------------------------------------------------------------
    !> Static constructor.
    !> (procedure public for type output_solver)
    !> Instantiate variable set to zero all varaible output_pcg
    !> as in defautl
    !>
    !> usage:
    !>     call 'var'%init()
    !<-------------------------------------------------------------
    subroutine init_output_solver(this)
      use Globals
      implicit none
      class(output_solver) ,intent(inout) :: this

      this%ierr            = 0
      this%iter            = 0
      this%resini          = zero
      this%resreal         = zero
      this%resnorm         = zero
      this%bnorm           = zero
     
      call this%tot%init()
      call this%ovh%init()
      call this%matvec%init()
      call this%vecvec%init()      
      call this%scalprod%init()
      call this%prec%init()
      call this%ort%init()
      
    end subroutine init_output_solver

    !>-------------------------------------------------------------
    !> Static deconstructor.
    !> (procedure public for type output_pcg)
    !> Reset to default all varaible output_solver
    !>
    !> usage:
    !>     call 'var'%init()
    !<-------------------------------------------------------------   
    subroutine kill_output_solver(this)
      use Globals
      implicit none
      class(output_solver) ,intent(inout) :: this

      this%ierr            = -999
      this%iter            = -999
      this%resini          = -huge
      this%resreal         = -huge
      this%resnorm         = -huge
      this%bnorm           = -huge
      
      call this%tot%kill()
      call this%ovh%kill()
      call this%matvec%kill()
      call this%vecvec%kill()      
      call this%scalprod%kill()
      call this%prec%kill()
      call this%ort%kill()


    end subroutine kill_output_solver

    !>-------------------------------------------------------------
    !> Static constructor.
    !> (procedure public for type output_pcg)
    !> Instantiate variable of type output_pcg
    !>
    !> usage:
    !>     call 'var'%init()
    !>
    !<-------------------------------------------------------------
    subroutine info_output_solver(this,lun)
      use Globals
      implicit none
      class(output_solver), intent(in) :: this  
      integer,              intent(in) :: lun
      !local
      character(len=256) :: msg=' '
      call this%info2str(msg)
      write(lun,*) etb(msg)
    end subroutine info_output_solver

    !>-------------------------------------------------------------
    !> Static constructor.
    !> (procedure public for type output_pcg)
    !> Instantiate variable of type output_pcg
    !>
    !> usage:
    !>     call 'var'%init()
    !>
    !<-------------------------------------------------------------
    subroutine infotiming_output_solver(this,lun)
      use Globals
      implicit none
      class(output_solver), intent(in) :: this  
      integer,              intent(in) :: lun
      
      call this%tot%info      (lun,'Total          =')
      call this%ovh%info      (lun,'Overhead       =')
      call this%matvec%info   (lun,'Matrix-vector  =')
      call this%vecvec%info   (lun,'Vector-vector  =')
      call this%scalprod%info (lun,'Scalar product =')
      call this%prec%info     (lun,'Prec. applic.  =')
      call this%ort%info      (lun,'Ortogon.       =')

      
    end subroutine infotiming_output_solver


    

    !>-------------------------------------------------------------
    !> Static constructor.
    !> (procedure public for type output_pcg)
    !> Instantiate variable of type output_pcg
    !>
    !> usage:
    !>     call 'var'%init()
    !>
    !<-------------------------------------------------------------
    function str_timing_output_solver(this) result(str)
      use Globals
      implicit none
      class(output_solver), intent(in) :: this  
      character(len=256)  :: str
      !
      character(len=30) :: loc
      write(loc,'(f5.2,a)') this%tot%CUMwct,'s:'
      str = etb(loc)
      
      
      
      write(loc,'(a,f4.1,a)') 'PREC=',100*this%prec%CUMwct/this%tot%CUMwct,'%,'
      str=etb(etb(str)//etb(loc))
      
      write(loc,'(a,f4.1,a)') 'MXV=',100*this%matvec%CUMwct/this%tot%CUMwct,'%,'
      str=etb(etb(str)//etb(loc))

      write(loc,'(a,f4.1,a)') 'VEC=',1.0d2*this%vecvec%CUMwct/this%tot%CUMwct,'%,'
      str=etb(etb(str)//etb(loc))

      write(loc,'(a,f4.1,a)') 'SCAL=',1.0d2*this%scalprod%CUMwct/this%tot%CUMwct,'%,'
      str=etb(etb(str)//etb(loc))
      
      write(loc,'(a,f4.1,a)') 'ORT=',1.0d2*this%ort%CUMwct/this%tot%CUMwct,'%,'
      str=etb(etb(str)//etb(loc))

      write(loc,'(a,f4.1,a)') 'OVH=',1.0d2*this%ovh%CUMwct/this%tot%CUMwct,'%'
      str=etb(etb(str)//etb(loc))
      
      
    end function str_timing_output_solver

    !>------------------------------------------------------------------
    !> Info procedure 
    !> (procedure public for type output_pcg)
    !> Generate a string summarizing output data
    !>
    !> usage:
    !>     call 'var'%info2str()
    !>
    !> \param[out]            string  -> character(len=70). String with output data
    !> \param[in ] (optional) iformat -> character(len=*). Integer variables format
    !> \param[in ] (optional) rformat -> character(len=*). Real variables format
    !<-------------------------------------------------------------------
    subroutine info_output_solver_2string(this,string,iformat,rformat)
      use Globals
      implicit none
      class(output_solver),       intent(in   ) :: this 
      character(len=256),         intent(inout) :: string
      character(len=*), optional, intent(in   ) :: iformat
      character(len=*), optional, intent(in   ) :: rformat
      !local
      character(len=256) :: outformat
      character(len=256) :: timeinfo=' '
      if (present(iformat) .and. present(rformat) ) then
         outformat=etb(&
              '('//etb(iformat)//',1x,'//&
              etb(iformat)//',1x,'//&
              etb(rformat)//',1x,'//&
              etb(rformat)//',1x,'//&
              etb(rformat)//',1x,'//&
              'a,a,a,1pe9.2,a)')
      else
         outformat='(I2,1x,I4,1x,3(1pe7.1,1x),a,8a,a)'
      end if
      if ( abs(this%tot%CUMwct-this%tot%wct) <1e-10) then
         write(timeinfo,'(a,1pe7.1,a)') 'WCT:',this%tot%CUMwct,' s'
      else
         write(timeinfo,'(a,1pe7.1,1x,a,1pe7.1,a)') 'CUMWCT',this%tot%CUMwct,'WCT:',this%tot%wct,' s'
      end if
      string=' '
      write(string,etb(outformat)) &
           this%ierr,&
           this%iter,&
           this%resini,&
           this%resnorm,&
           this%resreal,&
           ' <<',etb(this%scheme_used),'>>',&
           etb(timeinfo)
      
    end subroutine info_output_solver_2string


    !>-------------------------------------------------------------
    !> Linear solver procedure.
    !> Solves linear sysstem in the form A x = b
    !> where A is any type of matrix for which the procedure
    !> Mxv        => matrix vector multiplication
    !> diag_scale => scaling by a diagonal matrix
    !>
    !> usage:
    !>     call linear_solver(matrix, rhs,  sol, &
    !>     info_solver, &
    !>     ctrl_solver,&
    !>     prec_left,prec_right,&
    !>     aux,&
    !>     abpres)
    !> 
    !> usage:
    !>     call 'var'%Mxv(vec_in,vec_out,[info])
    !>
    !> where 
    !> \param[inout] matrix            -> Abstact matrix, with 
    !>                                    matrix-vector multiplication and 
    !>                                    matrix(transpose)-vector multipl.
    !> \param[in   ] rhs               -> real, dimension(matrix%nrow)
    !>                                    RHS of the system A x = b
    !> \param[inout] sol               -> real, dimension(matrix%ncol)
    !>                                    solution of the system A x = b 
    !> \param[inout] info_solver       -> type(output_solver)
    !>                                    Info parameters of linear solver 
    !> \param[in   ] ctrl_solver       -> type(input_solver)
    !>                                    Controls parameters of linear solver 
    !> \param[inout] (option.) prec_left-> Abstrac precoditioner
    !>                                    with apply operation. 
    !>                                    Prec. used in PCG procedure.
    !>                                    Left prec. used in other procedures.
    !> \param[inout] (option.) prec_right-> Abstrac precoditioner
    !>                                    with apply operation. 
    !>                                    Not used in PCG procedure.
    !>                                    right prec. used in other procedures.
    !> \param[inout] (optional) abpres -> type((pcgcoeff)
    !>                                    It stores alpha, beta and 
    !>                                    preconditioned residum in PCG
    !
    !
    recursive subroutine linear_solver(&
         matrix, rhs,  sol, &
         info_solver, &
         ctrl_solver,&
         prec_left,prec_right,&
         aux,&
         abpres,&
         ortogonalization_matrix)
      use Scratch
      use SimpleMatrix
      
      class(abs_linop),                intent(inout) :: matrix
      real(kind=double),               intent(in   ) :: rhs(matrix%ncol)
      real(kind=double),               intent(inout) :: sol(matrix%nrow)
      type(output_solver),             intent(inout) :: info_solver
      type(input_solver),              intent(in   ) :: ctrl_solver
      class(abs_linop),target,optional, intent(inout) :: prec_left
      class(abs_linop),target,optional, intent(inout) :: prec_right
      type(scrt),            optional, intent(inout) :: aux
      type(pcgcoeff),        optional, intent(inout) :: abpres
      class(abs_linop),    optional,  intent(inout) :: ortogonalization_matrix
          
      ! local
      logical :: rc
      integer :: lun_err, lun_out,nequ
      character(len=20) ::  scheme
      ! local pointers to handle preconditioner assignment
      class(abs_linop), pointer :: prec_left_loc=>null()
      class(abs_linop), pointer :: prec_right_loc=>null()
      ! local identity preconditioner
      type(scalmat),    target :: identity
      logical :: exit_test
      real(kind=double) :: dnrm2
      
      
      !
      ! set logical unit for err and output
      !
      lun_err = ctrl_solver%lun_err
      lun_out = ctrl_solver%lun_out

      !
      !
      !
      nequ = matrix%nrow
      

      !
      ! clean info solver
      !
      call info_solver%init()

      !
      ! compute rhs norm and set inverse_residum_weight
      !
      info_solver%bnorm = dnrm2(nequ,rhs,1)
      

      !
      ! handles case of zero RHS
      !
      call zerorhs(nequ,info_solver%ierr,&
           info_solver%bnorm,info_solver%resnorm,sol,rhs,exit_test)
      if (exit_test) return

      
      !
      ! check ortogonalization matrix
      !
      if (  present( ortogonalization_matrix) ) then 
         if ( ( ( ortogonalization_matrix%ncol .ne. nequ) .or. &
              ( ortogonalization_matrix%nrow .ne. nequ) ) ) then
            rc = IOerr(lun_err, wrn_val, ' linear_solver', &
                 'ortogonalization matrix does not fit matrix dimension')
            info_solver%ierr = 1
            return
         end if
      end if

      !
      ! set identity as preconditioners if no user-defined is passed
      !
      if ( (.not. present(prec_left)) .or. (.not. present(prec_right) ) ) then
         call identity%eye(matrix%nrow)
      end if
      if ( present(prec_left)) then
         prec_left_loc => prec_left
      else
         prec_left_loc => identity
      end if

      if ( present(prec_right)) then
         prec_right_loc => prec_right
      else
         prec_right_loc => identity
      end if

      !
      ! optional application of preconditoner to the rhs
      ! use x0=PL^{-1} rhs
      !
      if ( ctrl_solver%isol .eq. 1 ) then
         info_solver%info_prec = 0
         call info_solver%prec%set('start')
         call prec_left_loc%matrix_times_vector(rhs,sol,info_solver%info_prec,ctrl_solver%lun_err)
         call info_solver%prec%set('stop')
         if ( info_solver%info_prec .ne. 0 ) then
            rc = IOerr(ctrl_solver%lun_err, wrn_val, ' linear_solver', &
                 ' error computating sol = prec_left  rhs ( isol =1),'//&
                 ' info prec = ',info_solver%info_prec)
            info_solver%ierr = 3
            return
         end if
      end if
      


      !
      ! set scheme
      !
      scheme = ctrl_solver%scheme
      !if ( matrix%is_symmetric .and. prec_left_loc%is_symmetric ) scheme = 'PCG' 
      if ( ctrl_solver%debug .gt. 0 ) write( lun_out, *) ' Solve Linear solved via ', scheme 

      
      
      select case (etb(scheme))  
      case ('PCG')
         !
         ! PCG
         !
         call pcg_solver(matrix, rhs, sol,&
              info_solver, ctrl_solver,&
              prec=prec_left_loc,&
              aux_pcg=aux,&
              abpres_pcg=abpres,&
              ortogonalization_matrix= ortogonalization_matrix)
      case ('FPCG')
         !
         ! FPCG
         !
         call flexible_pcg_solver(matrix, rhs, sol,&
              info_solver, ctrl_solver,&
              prec=prec_left_loc,&
              aux_pcg=aux,&
              abpres_pcg=abpres,&
              ortogonalization_matrix= ortogonalization_matrix)


      case ('MINRES')
         !
         ! MINRES
         !
         call minres_solver(matrix, rhs, sol,&
              info_solver, ctrl_solver,&
              prec=prec_left_loc,&
              aux_minres=aux,&
              ortogonalization_matrix= ortogonalization_matrix)


      case ('BICGSTAB')
         call bicgstab_solver(matrix, rhs, sol,&
              info_solver, ctrl_solver,&
              prec_left_loc,& ! prec left
              prec_right_loc,& ! prec right
              aux_bicgstab=aux, &
              ortogonalization_matrix= ortogonalization_matrix)

      case ('GMRES')
         call gmres_solver(matrix, rhs, sol,&
              info_solver, ctrl_solver,&
              prec_left_loc,& ! prec left
              prec_right_loc,& ! prec right
              aux_gmres=aux)!, &
         !ortogonalization_matrix= ortogonalization_matrix)

      case ('FGMRES')
         call fgmres_solver(matrix, rhs, sol,&
              info_solver, ctrl_solver,&
              prec_left_loc,& ! prec left
              prec_right_loc,& ! prec right
              aux_fgmres=aux)

         
      case ('JACOBI_SEIDEL')
         call jacobi_seidel_solver(&
              matrix, rhs, sol,&
              info_solver, ctrl_solver,&
              prec_left_loc,& ! prec left
              aux=aux )
      case default
         rc = IOerr(lun_err, err_val, 'linear_solver', &
              ' scheme '//etb(scheme)//' not supported ')
       
      end select

      
      !
      ! store in info scheme used 
      !
      info_solver%scheme_used=etb(scheme)

      ! optional printing of convergence rate
      if (ctrl_solver%iprt.ge.1) call info_solver%info(lun_out)



      !
      ! clean
      !
      prec_left_loc=>null()
      prec_right_loc=>null()
      call identity%kill()

    end subroutine linear_solver


    

    !>-------------------------------------------------------------
    !> Static constructor.
    !> (procedure public for type lan)
    !> Instantiate variable of type lan
    !>
    !> usage:
    !>
    !> where:
    !> \param[in] lun_err      -> integer. Error mess. logic unit
    !> \param[in] length_pres -> integer. Eigenvector dimension
    !> \param[in] npres        -> integer. Number of prec. res. vectors
    !<-------------------------------------------------------------
    subroutine init_pcgcoeff(this, lun_err, length_pres, npres)
      use Globals
      implicit none
      !vars
      class(pcgcoeff), intent(inout) :: this
      integer,         intent(in   ) :: lun_err
      integer,         intent(in   ) :: length_pres
      integer,         intent(in   ) :: npres
      !local 
      integer :: res
      logical :: rc

      this%is_initialized = .true.    
      this%length_pres  = length_pres
      this%npres        = npres
      this%npres_stored = 0


    allocate (&
         this%alpha(npres),&
         this%beta(npres),&
         this%pres_dot_res(npres),&
         this%pres(length_pres,npres),&
         stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'init_pcgcoeff', &
         'type pcgcoeff member Ddiag, Tdiag, WMAT',res)

    this%alpha         = zero
    this%beta          = zero
    this%pres_dot_res  = zero
    this%pres          = zero

  end subroutine init_pcgcoeff

  !>-------------------------------------------------------------
  !> Static destructor.
  !> (procedure public for type pcgcoeff)
  !> Kill the variable of type pcgcoeff
  !>
  !> usage:
  !>     call 'var'%init(lun_err)
  !>
  !> where:
  !> \param[in] lun_err -> integer. Error mess. logic unit
  !<-------------------------------------------------------------
  subroutine kill_pcgcoeff(this, lun_err )
    use Globals
    implicit none
    !vars
    class(pcgcoeff), intent(inout) :: this
    integer,    intent(in   ) :: lun_err
    !local 
    integer :: res
    logical ::rc

    deallocate (&
         this%alpha,&
         this%beta,&
         this%pres_dot_res,&
         this%pres,&
         stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'kill_eig', &
         'type pcgcoeff member alpaha beta pres_dot_res pres',res)
    this%is_initialized = .false.    
  end subroutine kill_pcgcoeff

  !>-------------------------------------------------------------
  !> Info procedure.
  !> (procedure public for type eig)
  !> Kill the variable of type eig
  !>
  !> usage:
  !>     call 'var'%init(id_prt,lun)
  !>
  !> where:
  !> \param[in] id_prt -> integer. Flag dor printing eigval.
  !> \param[in] lun    -> integer. Logic unit
  !<-------------------------------------------------------------
  subroutine info_pcgcoeff(this, lun)
    use Globals
    implicit none
    !vars
    class(pcgcoeff), intent(in) :: this
    integer,     intent(in) :: lun

    write(lun,*) 'Eigen- length= ', this%length_pres
    write(lun,*) 'Npres        = ', this%npres
    write(lun,*) 'Npres_stored = ', this%npres_stored
   
  end subroutine info_pcgcoeff

  !>-----------------------------------------------------
  !> Procedure to fill arrays 
  !> Ddiag,TDiag and WMAT 
  !> used by for pcgcoeff decomposition
  !> with the alpha, beta, pres produced by PCG
  !>----------------------------------------------------
  subroutine fill_pcgcoeff(this,&
       iter,alpha, beta, pres_dot_res,pres)
    use Globals
    implicit none
    class(pcgcoeff),   intent(inout) :: this
    integer,           intent(in   ) :: iter
    real(kind=double), intent(in   ) :: alpha
    real(kind=double), intent(in   ) :: beta
    real(kind=double), intent(in   ) :: pres_dot_res
    real(kind=double), intent(in   ) :: pres(this%length_pres)

    if ( iter .le. this%npres) then
       this%alpha(iter)        = alpha
       this%beta(iter)         = beta
       this%pres_dot_res(iter) = pres_dot_res
       this%pres(:,iter)       = pres(:)
       this%npres_stored       = iter
    end if
    
  end subroutine fill_pcgcoeff


  recursive subroutine pcg_solver(matrix, rhs, sol,&
       info, ctrl, &
       prec,&
       aux_pcg,&
       abpres_pcg,&
       residual_norm,&
       ortogonalization_matrix)
    use Globals
    use Matrix
    use Scratch
    use Norms
    use Timing
    class(abs_linop),           intent(inout) :: matrix
    real(kind=double),           intent(in   ) :: rhs(matrix%nrow)
    real(kind=double),           intent(inout) :: sol(matrix%ncol)
    type(output_solver),         intent(inout) :: info
    type(input_solver),          intent(in   ) :: ctrl
    class(abs_linop),   optional, intent(inout) :: prec
    type(scrt), optional,target, intent(inout) :: aux_pcg
    type(pcgcoeff),    optional, intent(inout) :: abpres_pcg
    class(abs_norm), optional,target, intent(inout) :: residual_norm
    class(abs_linop), optional,  intent(inout) :: ortogonalization_matrix

    !
    ! local vars
    !
    ! logical
    logical :: rc
    logical :: exit_test
    ! string
    character(len=256) :: msg
    ! integer
    integer :: i
    integer :: lun_err, lun_out 
    integer :: iprt,iexit,imax,isol
    integer :: nequ,dim_ker,ndir=0
    integer :: iort
    integer :: info_prec
    ! real
    real(kind=double) ,pointer :: axp(:), pres(:), ainv(:), resid(:), pk(:),scr(:)
    real(kind=double) :: alpha, beta, presnorm
    real(kind=double) :: ptap
    real(kind=double) :: dnrm2,ddot,normres
    real(kind=double) :: tol,rort
    real(kind=double) :: inverse_residum_weight
    ! derivced type
    type(scrt), target  :: aux_work
    type(scrt), pointer :: aux_loc
    class(abs_norm), pointer :: local_norm
    type(euc_norm), target   :: euclidian_norm
    
    lun_err = ctrl%lun_err
    lun_out = ctrl%lun_out
    
    !
    ! symmtretry check
    !
    if (.not. matrix%is_symmetric) then
       rc= IOerr(lun_err, wrn_inp, 'pcg_solver', &
            'non symmetric matrix ')
       info%ierr = 999
       return
    end if

    

    if  (.not. prec%is_symmetric)  then
       msg='non symmetric preconditioner '
       rc= IOerr(lun_err, wrn_inp, 'pcg_solver', &
            msg)
       info%ierr = 999
       return
    end if


    call info%tot%set('start')
    call info%ovh%set('start')
    !
    ! set dimensions
    !
    nequ    = matrix%ncol

    !
    ! local copy of controls
    !
 
    iprt    = ctrl%iprt
    imax    = ctrl%imax
    isol    = ctrl%isol
    iexit   = ctrl%iexit
    tol     = ctrl%tol_sol
    iort    = ctrl%iort

    !
    ! set work arrays
    !    
    if ( present(aux_pcg) ) then
       aux_loc => aux_pcg
    else
       call aux_work%init(lun_err,0,6*nequ)
       aux_loc => aux_work  
    end if

    aux_loc%raux = zero
    axp          => aux_loc%raux(       1 :   nequ)
    pres         => aux_loc%raux(  nequ+1 : 2*nequ)
    ainv         => aux_loc%raux(2*nequ+1 : 3*nequ)
    resid        => aux_loc%raux(3*nequ+1 : 4*nequ)
    pk           => aux_loc%raux(4*nequ+1 : 5*nequ)
    scr          => aux_loc%raux(5*nequ+1 : 6*nequ)

    info%ierr = 0
    exit_test = .false.
    info%iter = 0

    !
    ! set metric to test the residuum and the norm of rhs
    !
    if( present ( residual_norm ) ) then
       local_norm => residual_norm
    else
       local_norm => euclidian_norm
    end if
    

    !
    ! In case of ortogonalization w.r.t to the matrix kernel 
    ! control if rhs is ortogonal to the kernel 
    ! up to tolerance required 
    !
    if (iort .ne. 0) then 
       if ( present(ortogonalization_matrix) ) then
          call ortogonalization_matrix%matrix_times_vector(rhs,scr,info%ierr,lun_err)

          if ( info%ierr .ne. 0) then
             call info%ovh%set('stop')
             call info%tot%set('stop ')
             return
          end if
          
          rort=dnrm2(nequ, scr,1)
          if  ( rort >  tol ) then
             write(msg,'(1e12.5)') rort
             rc = IOerr(lun_err, wrn_val, ' pcg_solver', &
                  'RHS not ortogonal to ortogonalization matrix'&
                  //etb(msg))
             info%ierr = 1
             call info%ovh%set('stop')
             call info%tot%set('stop ')
             return
          end if
       end if
    end if
    
    !
    ! print head table for convergence profile
    !
    if (iprt .ge. 2) write(lun_out,'(a)') '      iter        resnorm'

    !
    ! compute rhs norm and set inverse_residum_weight
    !
    info%bnorm = dnrm2(nequ,rhs,1)


    !
    ! handles case of zero RHS
    !
    call zerorhs(nequ,info%ierr,&
         info%bnorm,info%resnorm,sol,rhs,exit_test)

    !
    ! set inverse_residum_weight
    ! inverse_residum_weight = one       => resnorm = |Ax-b|
    ! inverse_residum_weight = one/bnrom => resnorm = |Ax-b|/|b|
    !
    if (ctrl%iexit.eq.0) then
       inverse_residum_weight = one
    else
       inverse_residum_weight = one / info%bnorm
    end if
    call info%ovh%set('stop')
    

    !
    ! calculate initial residual (res = b-M*x)
    !      
    call info%matvec%set('start')
    call matrix%matrix_times_vector(sol,resid,info%ierr,lun_err)
    call info%matvec%set('stop')
    if (info%ierr .ne. 0 ) return 
    
    

    call info%vecvec%set('start')
    resid = rhs - resid
    call info%vecvec%set('stop')
        
    !
    ! compute initial norm of residum
    !
    call info%scalprod%set('start')
    info%resini = dnrm2(nequ,resid,1)*inverse_residum_weight
    call info%scalprod%set('stop')

    if ( isnan(info%resini) ) then
       rc = IOerr(lun_err, wrn_val, ' pcg_solver', &
               'resini = NaN')
       info%ierr =-1
       call info%ovh%set('stop')
       call info%tot%set('stop ')
       return
    end if

    !
    ! exit if the initial residuum is already below the required tolerance 
    !
    if ( info%resini .lt. tol) then
       !
       info%resreal = info%resini
       info%resnorm = zero
       !
       ! free memory
       !
       call info%ovh%set('start')
       if ( aux_work%is_initialized ) call aux_work%kill(lun_err)  
       aux_loc => null() 
       call info%ovh%set('stop')
       call info%tot%set('stop ')
       return
    end if
    
    
    !
    ! cycle
    ! 
    do while (.not. exit_test)
       info%iter = info%iter + 1   
       !
       ! compute  pres = PREC  (r_k+1)
       !
       call info%prec%set('start')
       info_prec = 0
       call prec%matrix_times_vector(resid,pres,info_prec,lun_err)
       call info%prec%set('stop')
       if ( info_prec .ne. 0 ) then
          rc = IOerr(lun_err, wrn_val, ' pcg_solver', &
               ' error computating sol = P^{-1} r_{k+1},'//&
               ' info prec = ',info_prec )
          info%ierr = 3
          info%info_prec = info_prec
          return
       end if
       
       
       
       !
       ! optional (compute (pres)^T res for abpres variables)
       !
       if (present (abpres_pcg) ) then
          call info%ovh%set('start')
          presnorm = ddot(nequ,pres,1,resid,1)
          call info%ovh%set('start')
       end if

       !
       !  calculates \beta_k
       !
       call info%scalprod%set('start')
       if (info%iter.eq.1) then
          beta = zero
       else
          beta = -ddot(nequ,pres,1,axp,1)/ptap
       end if
       call info%scalprod%set('stop')

       !
       !  calculates p_k+1:=pres_k+beta_k*p_k
       !
       call info%vecvec%set('start')
       call dxpay(nequ,pres,1,beta,pk,1)
       call info%vecvec%set('stop')
              
          
       !
       !  calculates axp_k+1:= matrix * p_k+1
       !      
       call info%matvec%set('start')
       call matrix%matrix_times_vector(pk,axp,info%ierr,lun_err)
       call info%matvec%set('stop')

       
       !
       !  calculates \alpha_k
       !
       call info%scalprod%set('start')
       ptap  = ddot(nequ,pk,1,axp,1)
       alpha = ddot(nequ,pk,1,resid,1)/ptap
       call info%scalprod%set('stop')

       
       !
       ! optional storage of pcg variables (alpha, beta, etc.)
       ! contained in type pcg_coeff
       !
       call info%ovh%set('start')
       if (present(abpres_pcg)) then
          call abpres_pcg%fill(&
               info%iter,alpha, beta, presnorm,pres)
       end if
       call info%ovh%set('stop')

       !
       !  calculates x_k+1 and r_k+1
       !
       call info%vecvec%set('start')
       call daxpy(nequ,alpha,pk,1,sol,1)
       call daxpy(nequ,-alpha,axp,1,resid,1)
       call info%vecvec%set('stop')

       
       !
       ! optional ortogonalization 
       !
       if ( ( ctrl%iort  .gt. 0       ) .and. &
            ( mod(info%iter,ctrl%iort) .eq. 0 ) ) then
          if ( present(ortogonalization_matrix) ) then
             call info%ort%set('start')
             ! compute Qx=AA^T
             call ortogonalization_matrix%matrix_times_vector(&
                  sol,scr,info%ierr,lun_err)
             sol = sol - scr
             ! compute Qx=AA^T
             call ortogonalization_matrix%matrix_times_vector(&
                  resid,scr,info%ierr,lun_err)
             resid = resid - scr
             call info%ort%set('stop')
          end if
       end if

       !
       !  compute residum
       !
       call info%scalprod%set('start')       
       info%resnorm = dnrm2(nequ,resid,1)*inverse_residum_weight
       call info%scalprod%set('stop')

       if ( isnan(info%resnorm) ) then
          rc = IOerr(lun_err, wrn_val, ' pcg_solver', &
               'resnorm = NaN')
          info%ierr =-1 
          return
       end if
       
       if ( ( info%resnorm > 1.0d6 * info%resini)   ) then
          info%ierr =-1 
          return
       end if

       ! optional printing of convergence rate
       if (iprt.ge.2) write(lun_out,'(i10,e15.7)') &
            info%iter,info%resnorm


       exit_test = (info%iter.gt.imax .or. info%resnorm .le. tol)
       if (info%iter.ge.imax) info%ierr = 1
    end do

    !
    ! optional ortogonalization 
    !
    if ( ctrl%iort .gt. 0  ) then
       if ( present(ortogonalization_matrix) ) then
          call info%ort%set('start')
          ! compute Qx=AA^T
          call ortogonalization_matrix%matrix_times_vector(&
               sol,scr,info%ierr,lun_err)
          sol = sol - scr
       end if
    end if

    !
    ! compute final residum
    !
    call info%matvec%set('start')
    call matrix%matrix_times_vector(sol,resid,info%ierr,lun_err)
    call info%matvec%set('stop')

    call info%vecvec%set('start')
    resid = rhs - resid
    call info%vecvec%set('stop')
    
    !
    call info%scalprod%set('start')
    info%resreal = dnrm2(nequ,resid,1)*inverse_residum_weight
    call info%scalprod%set('stop')
    if ( isnan(info%resreal) ) then
       info%ierr =-1
       return
    end if


    !
    ! free memory if required
    !
    call info%ovh%set('start')
    if ( aux_work%is_initialized ) call aux_work%kill(lun_err)  
    aux_loc => null() 
    call info%ovh%set('stop')
    call info%tot%set('stop ')


  end subroutine pcg_solver

  recursive subroutine flexible_pcg_solver(matrix, rhs, sol,&
       info, ctrl, &
       prec,&
       aux_pcg,&
       abpres_pcg,&
       residual_norm,&
       ortogonalization_matrix)
    use Globals
    use Matrix
    use Scratch
    use Norms
    use Timing
    class(abs_linop),           intent(inout) :: matrix
    real(kind=double),           intent(in   ) :: rhs(matrix%nrow)
    real(kind=double),           intent(inout) :: sol(matrix%ncol)
    type(output_solver),         intent(inout) :: info
    type(input_solver),          intent(in   ) :: ctrl
    class(abs_linop),   optional, intent(inout) :: prec
    type(scrt), optional,target, intent(inout) :: aux_pcg
    type(pcgcoeff),    optional, intent(inout) :: abpres_pcg
    class(abs_norm), optional,target, intent(inout) :: residual_norm
    class(abs_linop), optional,  intent(inout) :: ortogonalization_matrix

    !
    ! local vars
    !
    ! logical
    logical :: rc
    logical :: exit_test
    ! string
    character(len=256) :: msg
    ! integer
    integer :: i,j,nres
    integer :: lun_err, lun_out 
    integer :: iprt,iexit,imax,isol
    integer :: nequ,dim_ker,ndir=0
    integer :: iort
    integer :: info_prec
    integer :: start,finish,slot,last_slot
    integer :: mj
    ! real
    real(kind=double) ,pointer :: axp(:), pres(:), ainv(:), resid(:), pk(:),scr(:)
    real(kind=double) :: alpha, beta, presnorm
    real(kind=double) :: ptap
    real(kind=double) :: dnrm2,ddot,normres
    real(kind=double) :: tol,rort
    real(kind=double) :: inverse_residum_weight
    ! derivced type
    type(scrt), target  :: aux_work
    type(scrt), pointer :: aux_loc
    class(abs_norm), pointer :: local_norm
    type(euc_norm), target   :: euclidian_norm
    
    lun_err = ctrl%lun_err
    lun_out = ctrl%lun_out
    
    !
    ! symmtretry check
    !
    if (.not. matrix%is_symmetric) then
       rc= IOerr(lun_err, wrn_inp, 'pcg_solver', &
            'non symmetric matrix ')
       info%ierr = 999
       return
    end if

    

    if  (.not. prec%is_symmetric)  then
       msg='non symmetric preconditioner '
       rc= IOerr(lun_err, wrn_inp, 'pcg_solver', &
            msg)
       info%ierr = 999
       return
    end if


    call info%tot%set('start')
    call info%ovh%set('start')
    !
    ! set dimensions
    !
    nequ    = matrix%ncol

    !
    ! local copy of controls
    !
 
    iprt    = ctrl%iprt
    imax    = ctrl%imax
    isol    = ctrl%isol
    iexit   = ctrl%iexit
    tol     = ctrl%tol_sol
    iort    = ctrl%iort

    !
    ! set work arrays
    !    
    if ( present(aux_pcg) ) then
       aux_loc => aux_pcg
    else
       call aux_work%init(lun_err,0,&
            6*nequ+2*nequ*(ctrl%fpcg_max_restart))
       ! 6*nequ for classical pcg       
       ! nequ*(ctrl%fpcg_max_restart-1) directions p_i
       ! nequ*(ctrl%fpcg_max_restart-1) Conjugate direction Ap_i 
       aux_loc => aux_work  
    end if

    aux_loc%raux = zero
    axp          => aux_loc%raux(       1 :   nequ)
    pres         => aux_loc%raux(  nequ+1 : 2*nequ)
    ainv         => aux_loc%raux(2*nequ+1 : 3*nequ)
    resid        => aux_loc%raux(3*nequ+1 : 4*nequ)
    pk           => aux_loc%raux(4*nequ+1 : 5*nequ)
    scr          => aux_loc%raux(5*nequ+1 : 6*nequ)

    info%ierr = 0
    exit_test = .false.
    info%iter = 0

    !
    ! set metric to test the residuum and the norm of rhs
    !
    if( present ( residual_norm ) ) then
       local_norm => residual_norm
    else
       local_norm => euclidian_norm
    end if
    

    !
    ! In case of ortogonalization w.r.t to the matrix kernel 
    ! control if rhs is ortogonal to the kernel 
    ! up to tolerance required 
    !
    if (iort .ne. 0) then 
       if ( present(ortogonalization_matrix) ) then
          call ortogonalization_matrix%matrix_times_vector(rhs,scr,info%ierr,lun_err)
          rort=dnrm2(nequ, scr,1)
          if  ( rort >  tol ) then
             write(msg,'(1e12.5)') rort
             rc = IOerr(lun_err, wrn_val, ' pcg_solver', &
                  'RHS not ortogonal to ortogonalization matrix'&
                  //etb(msg))
             info%ierr = 1
             return
          end if
       end if
    end if
    
    !
    ! print head table for convergence profile
    !
    if (iprt .ge. 1) write(lun_out,'(a)') '      iter        resnorm'

    !
    ! compute rhs norm and set inverse_residum_weight
    !
    info%bnorm = dnrm2(nequ,rhs,1)


    !
    ! handles case of zero RHS
    !
    call zerorhs(nequ,info%ierr,&
         info%bnorm,info%resnorm,sol,rhs,exit_test)

    !
    ! set inverse_residum_weight
    ! inverse_residum_weight = one       => resnorm = |Ax-b|
    ! inverse_residum_weight = one/bnrom => resnorm = |Ax-b|/|b|
    !
    if (ctrl%iexit.eq.0) then
       inverse_residum_weight = one
    else
       inverse_residum_weight = one / info%bnorm
    end if
    call info%ovh%set('stop')
    

    !
    ! calculate initial residual (res = b-M*x)
    !      
    call info%matvec%set('start')
    call matrix%matrix_times_vector(sol,resid,info%ierr,lun_err)
    call info%matvec%set('stop')

    call info%vecvec%set('start')
    resid = rhs - resid
    call info%vecvec%set('stop')
        
    !
    ! compute initial norm of residum
    !
    call info%scalprod%set('start')
    info%resini = dnrm2(nequ,resid,1)*inverse_residum_weight
    call info%scalprod%set('stop')

    if ( isnan(info%resini) ) then
       rc = IOerr(lun_err, wrn_val, ' pcg_solver', &
               'resini = NaN')
       info%ierr =-1  
       return
    end if

    !
    ! exit if the initial residuum is already below the required tolerance 
    !
    if ( info%resini .lt. tol) then
       !
       info%resreal = info%resini
       info%resnorm = zero
       !
       ! free memory
       !
       call info%ovh%set('start')
       if ( aux_work%is_initialized ) call aux_work%kill(lun_err)  
       aux_loc => null() 
       call info%ovh%set('stop')
       call info%tot%set('stop ')
       return
    end if
    
    
    !
    ! cycle
    ! 
    do while (.not. exit_test)
       info%iter = info%iter + 1   
       !
       ! compute  pres = PREC  (r_k+1)
       !
       call info%prec%set('start')
       info_prec = 0
       call prec%matrix_times_vector(resid,pres,info_prec,lun_err)
       call info%prec%set('stop')
       if ( info_prec .ne. 0 ) then
          rc = IOerr(lun_err, wrn_val, ' pcg_solver', &
               ' error computating sol = P^{-1} r_{k+1},'//&
               ' info prec = ',info_prec )
          info%ierr = 3
          info%info_prec = info_prec
          return
       end if
       
       
       
       !
       ! optional (compute (pres)^T res for abpres variables)
       !
       if (present (abpres_pcg) ) then
          call info%ovh%set('start')
          presnorm = ddot(nequ,pres,1,resid,1)
          call info%ovh%set('start')
       end if

       !
       !  calculates p_k+1:=pres_k+\sum_{j=m_j}^{k}beta_j*p_j
       !  with 
       !  mi=k+1-max(1,mod(k-1,fcg_max_restart+1)
       !  Algorithm 2.1 in "Flexible Conjugate Gradient" Notay
       !  (i in paper = k - 1 in this algorithm)
       ! 
       call info%vecvec%set('start')
       !
       ! pk=pres
       !
       call dcopy(nequ,pres,1,pk,1)

       !
       ! add beta_j*p_j
       !
       i=info%iter-1
       if (info%iter > 1) then
          nres = max (1, mod(info%iter-1, ctrl%fpcg_max_restart+1))
          do j= 1,nres
             !
             !  calculates \beta_k
             !
             ! compute beta
             ! find the interval where Ap_{i} is stored (see after)
             slot = mod(ctrl%fpcg_max_restart+last_slot-j,ctrl%fpcg_max_restart)+1
             start =6*nequ+(slot-1)*2*nequ+nequ+1              
             finish=6*nequ+(slot-1)*2*nequ+nequ+nequ
             call info%scalprod%set('start')
             beta = -ddot(nequ,pres,1,aux_loc%raux(start:finish),1)/ptap
             call info%scalprod%set('stop')
             ! find the interval where p_{i} is stored (see after)
             
             start =6*nequ+(slot-1)*2*nequ+1
             finish=6*nequ+(slot-1)*2*nequ+nequ
             call info%vecvec%set('start')
             call daxpy(nequ,beta, aux_loc%raux(start:finish),1,pk,1)
             start =6*nequ+(slot-1)*2*nequ+nequ+1              
             finish=6*nequ+(slot-1)*2*nequ+nequ+nequ
             
             call info%vecvec%set('stop')
          end do
          
          
       end if
       !
       ! store pk 
       !
       slot=mod(info%iter-1, ctrl%fpcg_max_restart)+1
       
       start =6*nequ+(slot-1)*2*nequ+1
       finish=6*nequ+(slot-1)*2*nequ+nequ
       aux_loc%raux(start:finish) = pk
       
       
       !
       !  calculates axp_k+1:= matrix * p_k+1
       !      
       call info%matvec%set('start')
       call matrix%matrix_times_vector(pk,axp,info%ierr,lun_err)
       call info%matvec%set('stop')

       !
       ! store Apk^{k}
       !
       slot=mod(info%iter-1, ctrl%fpcg_max_restart)+1
       last_slot=slot
       start =6*nequ+(slot-1)*2*nequ+nequ+1              
       finish=6*nequ+(slot-1)*2*nequ+nequ+nequ
       aux_loc%raux(start:finish) = axp

       !write(*,*) 'store',slot,info%iter,aux_loc%raux(start)



       
       !
       !  calculates \alpha_k
       !
       call info%scalprod%set('start')
       ptap  = ddot(nequ,pk,1,axp,1)
       alpha = ddot(nequ,pk,1,resid,1)/ptap
       call info%scalprod%set('stop')

       !
       ! optional storage of pcg variables (alpha, beta, etc.)
       ! contained in type pcg_coeff
       !
       call info%ovh%set('start')
       if (present(abpres_pcg)) then
          call abpres_pcg%fill(&
               info%iter,alpha, beta, presnorm,pres)
       end if
       call info%ovh%set('stop')

       !
       !  calculates x_k+1 and r_k+1
       !
       call info%vecvec%set('start')
       call daxpy(nequ,alpha,pk,1,sol,1)
       call daxpy(nequ,-alpha,axp,1,resid,1)
       call info%vecvec%set('stop')

       !
       ! optional ortogonalization 
       !
       if ( ( ctrl%iort  .gt. 0       ) .and. &
            ( mod(info%iter,ctrl%iort) .eq. 0 ) ) then
          if ( present(ortogonalization_matrix) ) then
             call info%ort%set('start')
             ! compute Qx=AA^T
             call ortogonalization_matrix%matrix_times_vector(&
                  sol,scr,info%ierr,lun_err)
             sol = sol - scr
             ! compute Qx=AA^T
             call ortogonalization_matrix%matrix_times_vector(&
                  resid,scr,info%ierr,lun_err)
             resid = resid - scr
             call info%ort%set('stop')
          end if
       end if

       !
       !  compute residum
       !
       call info%scalprod%set('start')
       
       info%resnorm = dnrm2(nequ,resid,1)*inverse_residum_weight
       call info%scalprod%set('stop')

       if ( isnan(info%resnorm) ) then
          rc = IOerr(lun_err, wrn_val, ' pcg_solver', &
               'resnorm = NaN')
          info%ierr =-1 
          return
       end if
       
       if ( ( info%resnorm > 1.0d6 * info%resini)   ) then
          info%ierr =-1 
          return
       end if

       ! optional printing of convergence rate
       if (iprt.ge.1) write(lun_out,'(i10,e15.7)') &
            info%iter,info%resnorm


       exit_test = (info%iter.gt.imax .or. info%resnorm .le. tol)
       if (info%iter.ge.imax) info%ierr = 1
    end do


    !
    ! compute final residum
    !
    call info%matvec%set('start')
    call matrix%matrix_times_vector(sol,resid,info%ierr,lun_err)
    call info%matvec%set('stop')

    call info%vecvec%set('start')
    resid = rhs - resid
    call info%vecvec%set('stop')
    
    !
    call info%scalprod%set('start')
    info%resreal = dnrm2(nequ,resid,1)*inverse_residum_weight
    call info%scalprod%set('stop')
    if ( isnan(info%resreal) ) then
       info%ierr =-1
       return
    end if


    !
    ! free memory if required
    !
    call info%ovh%set('start')
    if ( aux_work%is_initialized ) call aux_work%kill(lun_err)  
    aux_loc => null() 
    call info%ovh%set('stop')
    call info%tot%set('stop ')


  end subroutine flexible_pcg_solver


  
  !>-----------------------------------------------------------------
  !> Subroutine that handles case of zero RHS
  !>----------------------------------------------------------------
  subroutine zerorhs(nequ,info, bnorm,resnorm,sol,rhs,exit_test)
    use Globals
    implicit none
    integer,           intent(in   ) :: nequ
    integer,           intent(inout) :: info
    real(kind=double), intent(inout) :: bnorm
    real(kind=double), intent(inout) :: resnorm
    real(kind=double), intent(inout) :: sol(nequ)
    real(kind=double), intent(in   ) :: rhs(nequ)


    logical, intent(inout) ::   exit_test
    !local 

    real(kind=double) :: dnrm2

    exit_test=.FALSE.
    if (bnorm .lt. small) then
       sol = zero
       resnorm = zero
       info = 0
       exit_test=.true.
    end if
  end subroutine zerorhs


  subroutine gmres_solver(&
       matrix, rhs, sol,&
       info_solver, ctrl, &
       prec_left,&
       prec_right,&
       aux_gmres)
    use Globals
    use LinearOperator
    use Scratch
    class(abs_linop),            intent(inout) :: matrix
    real(kind=double),           intent(in   ) :: rhs(matrix%nrow)
    real(kind=double),           intent(inout) :: sol(matrix%ncol)
    class(output_solver),        intent(inout) :: info_solver
    class(input_solver),         intent(in   ) :: ctrl
    class(abs_linop),            intent(inout) :: prec_left
    class(abs_linop),            intent(inout) :: prec_right
    type(scrt), optional,target, intent(inout) :: aux_gmres

    ! local
    ! errors
    logical :: rc
    logical :: exit_test
    character(len=256) :: msg    
    integer :: INFOBLAS,res
   
    ! integer controls
    integer :: lun_err, lun_out, iprt,iexit,imax,isol, iort, nrestart
    ! indeces
    integer :: i,k,ikryl,iter
    integer :: igmres1,igmres2,igmres3,igmres4,igmres5,igmres6,igmres7,end
    integer :: info_prec
    ! lenght
    integer :: niaux, nraux, nequ,dim_ker

    ! real controls
    real(kind=double) :: tol,rort
    ! functions 
    real(kind=double) :: dnrm2,ddot,normres   
    ! work arrays
!!$    real(kind=double) ,pointer :: scr(:), resid(:), wvec(:), ykryl(:), zvec(:)
!!$    real(kind=double) ,pointer :: hess(:,:), vkryl(:,:)
    real(kind=double) ,allocatable :: scratch(:), pvec(:), vykryl(:)
    real(kind=double) ,allocatable :: cos1(:), sin1(:), rhs_trired(:)
    real(kind=double), allocatable :: resid(:), wvec(:), ykryl(:), zvec(:)
    real(kind=double) ,allocatable :: hess(:,:), vmatkryl(:,:)
    real(kind=double) :: bnorm, beta, resnorm,scale
    real(kind=double) :: inverse_residum_weight
    type(scrt), target  :: aux_work
    type(scrt), pointer :: aux_loc
    integer :: ndir=0

    !
    ! start algorithm 
    !
    call info_solver%tot%set('start')
    call info_solver%ovh%set('start')
    
    !
    ! set dimensions
    !
    nequ    = matrix%ncol

    !
    ! local copy of controls (EF remove?)
    !
    lun_err  = ctrl%lun_err
    lun_out  = ctrl%lun_out
    iprt     = ctrl%iprt
    imax     = ctrl%imax
    isol     = ctrl%isol
    iexit    = ctrl%iexit
    tol      = ctrl%tol_sol
    iort     = ctrl%iort
    
    nrestart = ctrl%nrestart   


    !
    ! setup for work arrays
    !
    igmres1 = 1
    igmres2 = igmres1 + nequ + 3*nrestart + 1
    igmres3 = igmres2 + nequ
    igmres4 = igmres3 + nequ*(nrestart+1)
    igmres5 = igmres4 + (nrestart+1)*nrestart
    igmres6 = igmres5 + nrestart
    igmres7 = igmres6 + nequ
    end     = igmres7 + nequ

!!$    if ( .not. present(aux_gmres) ) then
!!$       niaux = 0
!!$       nraux = nequ * ( nrestart + 5 ) + nrestart**2 + 5*nrestart + 1
!!$       call aux_work%init(lun_err, niaux,nraux)
!!$       aux_loc => aux_work   
!!$    else
!!$       aux_loc => aux_gmres
!!$    end if
!!$    
!!$    scr                             => aux_loc%raux(igmres1:igmres2-1)
!!$    resid                           => aux_loc%raux(igmres2:igmres3-1)
!!$    vmatkryl(1:nequ,     1:nrestart+1) => aux_loc%raux(igmres3:igmres4-1)
!!$    hess(1:nrestart+1,1:nrestart  ) => aux_loc%raux(igmres4:igmres5-1)
!!$    ykryl                           => aux_loc%raux(igmres5:igmres6-1)
!!$    wvec                            => aux_loc%raux(igmres6:igmres7-1)
!!$    zvec                            => aux_loc%raux(igmres7:end-1    )


!    write(*,*) 'info gmres_solv ',nrestart, nequ
    allocate(&
         rhs_trired(nrestart+1),& ! scr(3*nrestart+nequ+1)
         cos1(nrestart),&         ! 
         sin1(nrestart),&         !
         scratch(nequ),&          !
         pvec(nequ),&
         vykryl(nequ),&
         resid(nequ),&
         vmatkryl(nequ,nrestart+1),&
         hess(nrestart+1,nrestart),&
         ykryl(nrestart),&
         wvec(nequ),&
         zvec(nequ),&
         stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'gmres_solver', &
         'works arrays',res)
        
    !
    ! normal exit flag
    !
    info_solver%ierr = 0
    !
    ! set up initial constants and start loop
    !
    exit_test = .false.
    if (iprt.gt.0) write(lun_out,'(a)') '      iter        resnorm'
    info_solver%bnorm = dnrm2(nequ,rhs,1)
    
    
    !
    ! handles case of zero RHS and properly initializes exit_test
    !
    if (info_solver%bnorm.eq.zero) then
       call zerorhs(nequ,info_solver%ierr,info_solver%bnorm,info_solver%resnorm,sol,rhs,exit_test)
       if (exit_test) return
    end if


    !
    ! define the scaling factor of the residum
    ! iexit = 0 => 1.0
    ! iexit = 1 => 1.0 / || rhs ||
    !
    if (ctrl%iexit.eq.0) then
       inverse_residum_weight = one
    else
       inverse_residum_weight = one/info_solver%bnorm
    end if

    iter = 0
    !
    ! loop for the restarted version
    !
    do while (.not. exit_test .and. iter.lt.imax)
       !write(6,*) '-------------------------',iter

       !
       !  calculate initial residual b-m*x
       !
       call info_solver%matvec%set('start')
       call matrix%matrix_times_vector(sol,resid,info_solver%ierr,lun_err)
       call info_solver%matvec%set('stop')
       !
       call info_solver%vecvec%set('start')
       wvec=rhs
       resid = wvec - resid
       call info_solver%vecvec%set('stop')
       if(iter.eq.0) then
          !
          ! compute initial residum
          !
          call info_solver%scalprod%set('start')
          info_solver%resini = &
               dnrm2(nequ,resid,1) * inverse_residum_weight
          call info_solver%scalprod%set('stop')

          if(info_solver%resini.le.ctrl%tol_sol) then
             info_solver%resreal = info_solver%resini
             info_solver%iter = 0
             return
          end if          
       end if

       ! old librykrysol
       ! precondition if necessary
       !
       !  call precnsy(iprec,.true.,nequ,ntermp,iap,jap,idiag,
       !  1                prec,scr(2+3*nrestart),res,vkryl(1,1))
       !

       info_prec = 0
       call prec_left%matrix_times_vector( resid, vmatkryl(1:nequ,1),info_prec,lun_err)
       if ( info_prec .ne. 0 ) then
          rc = IOerr(lun_err, wrn_val, ' gmres_solver', &
               ' error computating sol = p^{-1} rhs ( isol =1),'//&
               ' info prec = ',info_prec)
          info_solver%ierr = 3
          info_solver%info_prec = info_prec
          return
       end if


!!$       write(100,*) resid
!!$       write(100,*) 
!!$       write(100,*) vmatkryl(:,1)
!!$       write(100,*) 

       beta = dnrm2(nequ,vmatkryl(:,1),1)

       call drscl(nequ,beta,vmatkryl(:,1),1)

!!$       write(100,*) beta
!!$       write(100,*)
!!$       write(100,*) vmatkryl(:,1)
!!$       write(100,*)
       !
       ! loop on the Krylov subspace
       !
       ikryl = 0

       do while ( .not. exit_test .and. ikryl.lt.nrestart)
          ikryl = ikryl+1
                
          !
          ! apply P_right v
          !
          call prec_right%matrix_times_vector(vmatkryl(:,ikryl),pvec,&
               info_solver%info_prec,ctrl%lun_err)
          if (info_solver%info_prec .ne. 0) then
             info_solver%ierr = 3
             return
          end if
          
          
          !
          ! compute zvec =A P_{2}  v
          !
          call matrix%matrix_times_vector(pvec,zvec,info_solver%ierr,lun_err)

          !
          ! compute wvec = P_{1} A P_{2} v
          !
          call prec_left%matrix_times_vector(zvec,wvec,&
               info_solver%info_prec,ctrl%lun_err)
          if (info_solver%info_prec .ne. 0) then
             info_solver%ierr = 3
             return
          end if
          
!!$          write(100,*) wvec
!!$          write(100,*) 
!!$          write(100,*) vmatkryl(:,ikryl)
!!$          write(100,*)
          !
          !  modified Gram-Schmidt orthogonalization
          !
          do k = 1,ikryl
             hess(k,ikryl) = ddot(nequ,wvec,1,vmatkryl(1,k),1)
             call daxpy(nequ,-hess(k,ikryl),vmatkryl(1,k),1,wvec,1)
          end do

!!$          write(100,*) iter,ikryl,nrestart
!!$          write(100,*) 
!!$          write(100,*) hess
!!$          write(100,*) 
          
          hess(ikryl+1,ikryl) = dnrm2(nequ,wvec,1)
          call dcopy(nequ,wvec,1,vmatkryl(1:nequ,ikryl+1),1)
          call drscl(nequ,hess(ikryl+1,ikryl),vmatkryl(1:nequ,ikryl+1),1)
          
          !
          !  triangularization of ||\beta e_1 - \bar{H_m} y||_2
          !
          call local_trired(ikryl,&
               nrestart,&
               beta,&
               hess,&
               rhs_trired,&  ! rhs
               cos1,&        ! c1
               sin1)         ! s1

!!$          write(6,*) iter,ikryl,nrestart
!!$          write(6,*) 
!!$          write(6,*) hess
!!$          write(6,*)
!!$          write(6,*) 'scr', info_solver%bnorm, ikryl, iter
!!$          write(6,*) scr(1:nrestart+1)
!!$          write(6,*) 


          !
          ! compute weighted residual norm
          !          
          info_solver%resnorm = abs(rhs_trired(ikryl+1)) * inverse_residum_weight

          if(iprt.ne.0)&
               write(lun_out,'(i4,1pe15.5)') iter+ikryl,info_solver%resnorm

          exit_test=( info_solver%resnorm .lt. ctrl%tol_sol )
       end do
       !
       ! solve the linear system
       !       
       call dcopy(ikryl,rhs_trired,1,ykryl,1)
       call dlatrs( 'U', 'N', 'N', 'N', ikryl, hess, nrestart+1, ykryl,&
            scale, scratch, infoblas ) 
       if (infoblas .ne. 0) then
          rc = IOerr(lun_err, wrn_inp, 'gmres_solver', &
               ' in sub. lsqsol: sub. dlatrs info=',infoblas)
       end if

       !
       !  calculate the approximate solution x_m = x_0 + (K^{-1}) V_m y_m
       !
       vykryl = zero
       
       call dgemv('N',nequ,ikryl,one,vmatkryl,nequ,&
            ykryl,1,one,& ! ! EF put zero instead of one
            vykryl,1)
       !
       ! restore the original variable (if necessary)
       ! here zvec is a work array containing U^{-1} V_m y_m
       !
       ! in libkrysol
       ! call precnsy(iprec,.false.,nequ,ntermp,iap,jap,idiag,
       !  prec,wvec,scr(2+3*mmax),zvec)

       !
       ! apply PR v
       !
       call prec_right%matrix_times_vector(vykryl,zvec,&
            info_solver%info_prec,ctrl%lun_err)
       if (info_solver%info_prec .ne. 0) then
          info_solver%ierr = 3
          return
       end if
       
       call daxpy(nequ,one,zvec,1,sol,1)

       iter = iter + ikryl

    end do
    info_solver%iter = iter
    
    if(info_solver%iter.ge.imax) info_solver%ierr=1

    !
    ! compute final residum
    !
    call info_solver%matvec%set('start')
    call matrix%matrix_times_vector(sol,resid,info_solver%ierr,lun_err)
    call info_solver%matvec%set('stop')

    call info_solver%vecvec%set('start')
    wvec = rhs
    resid = wvec - resid
    call info_solver%matvec%set('stop')
    !
    call info_solver%scalprod%set('start')
    info_solver%resreal=dnrm2(nequ,resid,1)*inverse_residum_weight
    call info_solver%scalprod%set('stop')

!!$    if ( aux_work%is_initialized ) then
!!$       aux_loc =>null()
!!$       !call aux_work%info(0)
!!$       call aux_work%kill(lun_err)
!!$    end if

    deallocate(&
         rhs_trired,&
         cos1,&
         sin1,&
         scratch,&
         pvec,&
         vykryl,&
         resid,&
         vmatkryl,&
         hess,&
         ykryl,&
         wvec,&
         zvec,&
         stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_dealloc, 'gmres_solver', &
         'works arrays',res)
    
  contains
    
    !************************* TRIRED ***********************************
    !
    subroutine local_trired(ndim,nrestart,beta,hess,rhs,c1,s1)
      !
      !  reduces matrix \bar{H_m}  to triangular form by 
      !  successive application of Givens rotations and updates the RHS vector
      !
      use Globals
      implicit none
      integer,           intent(in   ) :: ndim,nrestart
      real(kind=double), intent(in   ) :: beta
      real(kind=double), intent(inout) :: hess(nrestart+1,nrestart)
      real(kind=double), intent(inout) :: rhs(nrestart+1)
      real(kind=double), intent(inout) :: c1(nrestart)
      real(kind=double), intent(inout) :: s1(nrestart)
      
      ! local
      integer ::  i,j
      real(kind=double):: hyp,save_hess,save_rhs

!!$      write(*,*) 'in_trired ', ndim, nrestart, beta
!!$      write(300,'(a,/,(5e13.4))') 'c1 ', c1
!!$      write(300,'(a,/,(5e13.4))') 's1 ', s1
!!$      write(300,'(a,/,(5e13.4))') 'hess ', hess
      !
      !   update the new elements of \bar{H_m} with the old Givens rotation
      !
      if (ndim.eq.1) then
         rhs(1)=beta
      else
         do i=1,ndim-1
            save_hess      =  c1(i)*hess(i,ndim) + s1(i)*hess(i+1,ndim)
            hess(i+1,ndim) = -s1(i)*hess(i,ndim) + c1(i)*hess(i+1,ndim)
            hess(i,ndim)   = save_hess
         end do
      end if

      !
      !  construct the Givens rotation
      !
      hyp  = sqrt(hess(ndim,ndim)**2 + hess(ndim+1,ndim)**2)
      s1(ndim) = hess(ndim+1,ndim)/hyp
      c1(ndim) = hess(ndim,ndim)/hyp
      !        |1           |
      !        | 1          |
      !        |  ....      |
      !  W_i = |   c_i  s_i |   W_i is an (i+1)x(i+1) matrix (i = ndim)
      !        |  -s_i  c_i |
      !
      !  applies the latest Givens rotation
      !
      !  to \bar{H_m}...
      !
      hess(ndim,ndim)   = &
           c1(ndim)*hess(ndim,ndim) + &
           s1(ndim)*hess(ndim+1,ndim)
      hess(ndim+1,ndim) = zero
      !
      !  and to rhs
      !
      save_rhs    =  c1(ndim)*rhs(ndim)
      rhs(ndim+1) = -s1(ndim)*rhs(ndim)
      rhs(ndim)   =  save_rhs

!!$      write(200,*)ndim
!!$      do j = 1,ndim
!!$         write(200,*)j,(hess(j,i),i=1,ndim),rhs(j)
!!$      end do

    end subroutine local_trired
    
  end subroutine gmres_solver
  
   recursive subroutine bicgstab_solver(&
       matrix, rhs, sol,&
       info, ctrl, &
       prec_left,&
       prec_right,&
       aux_bicgstab,&
       ortogonalization_matrix&
       )
    use Globals
    use LinearOperator
    use Scratch
    use Timing
    class(abs_linop),           intent(inout) :: matrix
    real(kind=double),           intent(in   ) :: rhs(matrix%nrow)
    real(kind=double),           intent(inout) :: sol(matrix%ncol)
    class(output_solver),        intent(inout) :: info
    class(input_solver),         intent(in   ) :: ctrl
    class(abs_linop),             intent(inout) :: prec_left
    class(abs_linop),             intent(inout) :: prec_right
    type(scrt), optional,target, intent(inout) :: aux_bicgstab
    class(abs_linop), optional,  intent(inout) :: ortogonalization_matrix

    
    ! local
    ! errors
    logical :: rc
    logical :: exit_test,test
    character(len=256) :: msg    
    integer :: INFOBLAS,res
   
    ! integer controls
    integer :: lun_err, lun_out, iprt,iexit,imax,isol, iort, nrestart
    integer :: info_prec
    ! indeces
    integer :: iter, i
    ! lenght
    integer :: nequ,ndir=0

    ! real controls
    real(kind=double) :: tol
    ! functions 
    real(kind=double) :: dnrm2,ddot,normres   
    ! work arrays
    real(kind=double), pointer :: resid(:), scr(:), res0(:)
    real(kind=double), pointer :: pk(:), pkprime(:)
    real(kind=double), pointer :: apk(:), skprime(:), ask(:)
    real(kind=double) :: inverse_residum_weight
    real(kind=double) :: bnorm,resnorm,snorm
    real(kind=double) :: rhok, rhokp1, alpha, beta, omega 
    real(kind=double) :: inv_res_weight,rort
    type(scrt), target  :: aux_work
    type(scrt), pointer :: aux_loc
    

    ! allocatable arrays
!!$    real(kind=double), allocatable :: resid(:), scr(:), res0(:)
!!$    real(kind=double), allocatable :: pk(:), pkprime(:)
!!$    real(kind=double), allocatable :: apk(:), skprime(:), ask(:)

    !
    ! start algorithm 
    !
    call info%tot%set('start')
    call info%ovh%set('start')
    
    !
    ! set dimensions
    !
    nequ    = matrix%ncol

    !
    ! local copy of controls (EF remove?)
    !
    lun_err  = ctrl%lun_err
    lun_out  = ctrl%lun_out
    iprt     = ctrl%iprt
    imax     = ctrl%imax
    isol     = ctrl%isol
    iexit    = ctrl%iexit
    tol      = ctrl%tol_sol
    iort     = ctrl%iort
    nrestart = ctrl%nrestart 

    !
    ! set work arrays
    !
    if ( present(aux_bicgstab) ) then
       aux_loc => aux_bicgstab
       test=aux_bicgstab%check(0,8*nequ)
       if( .not. aux_loc%check(0,8*nequ) ) &
            rc = IOerr(lun_err, err_inp, 'bicgstab_solver', &
            ' aux array too small',res)
    else
       call aux_work%init(lun_err,0,9*nequ)
       aux_loc => aux_work  
    end if
    resid   => aux_loc%raux(       1 :   nequ)
    scr     => aux_loc%raux(  nequ+1 : 2*nequ)
    res0    => aux_loc%raux(2*nequ+1 : 3*nequ)
    pk      => aux_loc%raux(3*nequ+1 : 4*nequ)
    pkprime => aux_loc%raux(4*nequ+1 : 5*nequ)
    apk     => aux_loc%raux(5*nequ+1 : 6*nequ)
    ask     => aux_loc%raux(6*nequ+1 : 7*nequ)
    skprime => aux_loc%raux(7*nequ+1 : 8*nequ)

!!$ THESE LINES ARE LEFT FOR DEBUG PURPOSE
!!$    allocate(&
!!$         resid(nequ),&
!!$         scr(nequ),&
!!$         res0(nequ),&
!!$         pk(nequ),&
!!$         pkprime(nequ),&
!!$         apk(nequ),&
!!$         ask(nequ),&
!!$         skprime(nequ),&
!!$         stat=res)
!!$    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'bicgstab_solver', &
!!$         'works arrays',res)

    !
    ! In case of ortogonalization w.r.t to the matrix kernel 
    ! control if rhs is ortogonal to the kernel 
    ! up to tolerance required 
    !
    if ( (iort .ne. 0           ) .and. &
         present( ortogonalization_matrix ) ) then 
       if (matrix%is_symmetric   ) then
          call ortogonalization_matrix%matrix_times_vector(rhs,scr,info%ierr,lun_err)
          rort=dnrm2(nequ, scr,1)
          if  ( rort >  tol ) then
             write(msg,'(1e12.5)') rort
             rc = IOerr(lun_err, wrn_val, ' pcg_solver', &
                  'RHS not ortogonal to ortogonalization matrix'&
                  //etb(msg))
             info%ierr = 1
             return
          end if
       end if
    end if
   
    !
    ! normal exit flag
    !
    info%ierr = 0
    
    !
    ! set up initial constants and start loop
    !
    exit_test = .false.
    if (iprt.gt.0) write(lun_out,'(a)') 'iter        resnorm'
    
    !
    ! eval rhs norm
    !
    bnorm = dnrm2(nequ,rhs,1)
    info%bnorm = bnorm

    !
    ! handles case of zero RHS
    !
    call zerorhs(nequ,info%ierr,&
         info%bnorm,info%resnorm,sol,rhs,exit_test)

    !
    ! handles case of zero RHS and properly initializes exit_test
    !
    if (bnorm.eq.zero) then
       call zerorhs(nequ,info%ierr,bnorm,info%resnorm,sol,rhs,exit_test)
       if (exit_test) return
    end if


    !
    ! set inverse_residum_weight
    ! inverse_residum_weight = one       => resnorm = |Ax-b|
    ! inverse_residum_weight = one/bnrom => resnorm = |Ax-b|/|b|
    !
    if (ctrl%iexit.eq.0) then
       inverse_residum_weight = one
    else
       inverse_residum_weight = one / info%bnorm
    end if
    

    !
    ! preliminaries ended
    !
    call info%ovh%set('stop')


    ! 
    !  start main loop   
    !
    iter = 0

    !
    !  calculate initial residual b-M*x
    !
    if (ctrl%debug .gt. 0 ) write(lun_out,*) 'Matrix res'
    call info%matvec%set('start')
    call matrix%matrix_times_vector(sol,resid,info%ierr,lun_err)
    call info%matvec%set('stop')
    !
    call info%vecvec%set('start')
    scr=rhs
    resid = scr - resid
    call info%vecvec%set('stop')

    !
    ! compute initial residual
    !
    call info%scalprod%set('start')
    info%resini = &
         dnrm2(nequ,resid,1) * inverse_residum_weight
    call info%scalprod%set('stop')
    
    if ( isnan(info%resini) ) then
       info%ierr =-1 
       return
    end if


    if(info%resini.le.ctrl%tol_sol) then
       info%resreal = info%resini
       info%iter = 0
       return
    end if


 

    ! 
    !  precondition if neccessary
    !  (here res0 is used as scratch vector)
    !
    !  call precnsy(iprec,.true.,nequ,ntermp,iap,jap,idiag,&
    !               prec,scr,res0,res)
    call dcopy(nequ,resid,1,res0,1)
    info_prec = 0
    if (ctrl%debug .gt. 0 ) write(lun_out,*) 'Precleft^{-1} resid'
    call prec_left%matrix_times_vector(res0,resid,info_prec,lun_err)
    if ( info_prec .ne. 0 ) then
       rc = IOerr(lun_err, wrn_val, ' bicgstab_solver ', &
            ' error computing resid = pleft^{-1} res0,'//&
            ' info prec = ',info_prec)
       info%ierr = 3
       info%info_prec = info_prec
       return
    end if
    call dcopy(nequ,resid,1,pk,1)
    call dcopy(nequ,resid,1,res0,1)
    rhok = ddot(nequ,resid,1,res0,1)  ! EF = dnrm2(nequ,resid,1) ?

    !
    !  loop on the Krylov subspace
    !
    do while (.not. exit_test .and. iter.le.imax)
       iter = iter + 1
       !
       !  apply the preconditioned matrix: p' := U^-1 p
       !          (e.g. if split prec.:  apk <-- L^-1 A U^-1)
       !
       !  call aprecdotx(iprec,nequ,nterm,ntermp,ia,ja,iap,jap,idiag,
       !  1      sysmat,prec,scr,pk,pkprime,apk)

       !
       ! apply v' = P_right v
       !
       call info%prec%set('start')
       if (ctrl%debug .gt. 0 ) write(lun_out,*) '    Precright^{-1} pk'
       info_prec = 0
       call prec_right%matrix_times_vector(pk,pkprime,info_prec,lun_err)
       if ( info_prec .ne. 0 ) then
          rc = IOerr(lun_err, wrn_val, ' bicgstab_solver', &
               ' error computing pk = prec_right^{-1} pkprim,'//&
               ' info prec = ',info_prec)
          info%ierr = 3
          info%info_prec = info_prec
          return
       end if
       call info%prec%set('stop')  
!!$       if ( ( matrix%dim_kernel    .gt. 0       ) .and. &
!!$            ( ctrl%iort  .gt. 0       ) .and. &
!!$            ( mod(info%iter,ctrl%iort) .eq. 0 ) ) then
!!$          call info%ort%set('start')
!!$          call ortogonalize(nequ,matrix%dim_kernel,matrix%kernel,pkprime)
!!$          call info%ort%set('stop')
!!$       end if


       !
       ! compute scr =A P_{2}  v
       !
       if (ctrl%debug .gt. 0 ) write(lun_out,*) '              Matrix PrecRight^{-1} pk'
       call info%matvec%set('start')
       call matrix%matrix_times_vector(pkprime,scr,info%ierr,lun_err)
       call info%matvec%set('stop')
!!$       if ( ( matrix%dim_kernel    .gt. 0       ) .and. &
!!$            ( ctrl%iort  .gt. 0       ) .and. &
!!$            ( mod(info%iter,ctrl%iort) .eq. 0 ) ) then
!!$          call info%ort%set('start')
!!$          call ortogonalize(nequ,matrix%dim_kernel,matrix%kernel,scr)
!!$          call info%ort%set('stop')
!!$       end if


       !
       ! compute apk = P_{1} A P_{2} v
       !
       if (ctrl%debug .gt. 0 ) write(lun_out,*) 'PrecLeft^{-1} Matrix PrecRight^{-1} pk'
       call info%prec%set('start')
       info_prec = 0
       call prec_left%matrix_times_vector(scr,apk,info_prec,lun_err)
       if ( info_prec .ne. 0 ) then
          rc = IOerr(lun_err, wrn_val, ' bicgstab_solver', &
               ' error computing apk = prec_left^{-1} scr, info prec = ',info_prec)
          info%ierr = 3
          info%info_prec = info_prec
          return
       end if
       call info%prec%set('stop')
!!$       if ( ( matrix%dim_kernel    .gt. 0       ) .and. &
!!$            ( ctrl%iort  .gt. 0       ) .and. &
!!$            ( mod(info%iter,ctrl%iort) .eq. 0 ) ) then
!!$          call info%ort%set('start')
!!$          call ortogonalize(nequ,matrix%dim_kernel,matrix%kernel,skprime)
!!$          call info%ort%set('stop')
!!$       end if


       !
       !  alpha := r^T r_0/r_0^T A p'
       !
       call info%scalprod%set('start')
       alpha = rhok/ddot(nequ,apk,1,res0,1)
       if (ctrl%debug .eq. 1) write(lun_out,*) 'alpha =', alpha
       call info%scalprod%set('stop')
       
       !
       !  s := r - alpha A p
       !
       call info%vecvec%set('start')
       call daxpy(nequ,-alpha,apk,1,resid,1)
       call info%vecvec%set('stop')

       call info%scalprod%set('start')
       snorm = dnrm2(nequ,resid,1)
       if (ctrl%debug .eq. 1) write(lun_out,*) 'snorm =', snorm
       resnorm = snorm * inverse_residum_weight
       call info%scalprod%set('stop')

       if ( isnan(info%resnorm) ) then
          info%ierr = -1  
          return
       end if

       if(resnorm.lt.tol) then
          !
          !  convergence achieved
          !
          !  x := x + alpha p'
          !
          call info%vecvec%set('start')
          call daxpy(nequ,alpha,pkprime,1,sol,1)
          call info%vecvec%set('stop')
          exit_test=.true.
       else
          !
          !  s' := K^-1 s
          !
          !  call aprecdotx(iprec,nequ,nterm,ntermp,ia,ja,iap,jap,idiag,
          !  1                     sysmat,prec,scr,res,skprime,ask)
          
          !
          ! apply P_right v
          !
          info_prec = 0
          call info%prec%set('start')
          call prec_right%matrix_times_vector(resid,skprime,info_prec,lun_err)
          call info%prec%set('stop')
          if ( info_prec .ne. 0 ) then
             rc = IOerr(lun_err, wrn_val, ' bicgstab_solver', &
                  ' error computing skprime = prec_right^{-1} resid,'// &
                  'info prec = ',info_prec)
             info%ierr = 3
             info%info_prec = info_prec
             return
          end if

                 
          !
          ! compute zvec =A P_{2}  v
          !
          call info%matvec%set('start')
          call matrix%matrix_times_vector(skprime,scr,info%ierr,lun_err)
          call info%matvec%set('stop')
          !
          ! compute wvec = P_{1} A P_{2} v
          !
          info_prec = 0
          call info%prec%set('start')
          call prec_left%matrix_times_vector(scr,ask,info%ierr,lun_err)
          call info%prec%set('stop')
          
          if ( info_prec .ne. 0 ) then
             rc = IOerr(lun_err, wrn_val, ' bicgstab_solver', &
                  ' error computing ask = prec_left^{-1} scr,'// &
                  'info prec = ',info_prec)
             info%ierr = 3
             info%info_prec = info_prec
             return
          end if
          

          !
          !  omega := s^T A s' / (A s')^T A s'
          !
          call info%scalprod%set('start')
          omega = ddot(nequ,resid,1,ask,1)/ddot(nequ,ask,1,ask,1)
          if (ctrl%debug .eq. 1) write(lun_out,*) 'omega =', omega
          call info%scalprod%set('stop')
          
          !
          !  x := x + alpha p' + omega s'
          !
!!$          if ( present(ortogonalization_matrix) .and. &
!!$               ( ctrl%iort  .gt. 0            ) .and. &
!!$               ( mod(info%iter,ctrl%iort) .eq. 0 ) ) then
!!$             call info%ort%set('start')
!!$             ! compute Qx=AA^Tx
!!$             call ortogonalization_matrix%matrix_times_vector(&
!!$                  pkprime,scr)
!!$             pkprime = pkprime - scr
!!$             
!!$             ! compute Qx=AA^Tx
!!$             call ortogonalization_matrix%matrix_times_vector(&
!!$                  skprime,scr)
!!$             skprime = skprime - scr
!!$             call info%ort%set('stop')
!!$             call info%ort%set('stop')
!!$          end if  
          call info%vecvec%set('start')
          call daxpy(nequ,alpha,pkprime,1,sol,1)
          call daxpy(nequ,omega,skprime,1,sol,1)
          call info%vecvec%set('stop')
!!$          if ( present(ortogonalization_matrix) .and. &
!!$               ( ctrl%iort  .gt. 0            ) .and. &
!!$               ( mod(info%iter,ctrl%iort) .eq. 0 ) ) then
!!$             call info%ort%set('start')
!!$             ! compute Qx=AA^Tx
!!$             call ortogonalization_matrix%matrix_times_vector(&
!!$                  sol,scr)
!!$             sol = sol - scr
!!$          end if
          
          !
          !  r := s - omega A s'
          !
          call info%vecvec%set('start')
          call daxpy(nequ,-omega,ask,1,resid,1)
          call info%vecvec%set('stop')

          call info%scalprod%set('start')
          snorm = dnrm2(nequ,resid,1)
          resnorm = snorm * inverse_residum_weight 
          if (ctrl%debug .eq. 1) write(lun_out,*) 'snorm =', snorm
          call info%scalprod%set('stop')
          if ( isnan(resnorm) ) then
             info%ierr =-1  
             return
          end if
          
          if ( present(ortogonalization_matrix) .and. &
               ( ctrl%iort  .gt. 0            ) .and. &
               ( mod(info%iter,ctrl%iort) .eq. 0 ) ) then
             call info%ort%set('start')
             ! compute Qx=AA^Tx
             call ortogonalization_matrix%matrix_times_vector(&
                  sol,scr,info%ierr,lun_err)
             sol = sol - scr
             call info%ort%set('stop')
          end if  
               

          exit_test=(resnorm.lt.tol)
          if ( resnorm > 1.0d5*info%resini ) then
             info%ierr = -1
             rc = IOerr(lun_err, wrn_val, ' bicgstab_solver', &
                  ' error resnorm > 1.e5 initial residum ',info%ierr)
             return
          end if

          if (.not. exit_test) then
             !
             !  rho := r^T r_0
             !
             call info%scalprod%set('start')
             rhokp1 = ddot(nequ,resid,1,res0,1)
             call info%scalprod%set('stop')

             if (abs(rhokp1).le.verysmall) then
                exit_test = .true.
                info%ierr = 2
             end if

             beta = (rhokp1/rhok)*(alpha/omega)
             rhok = rhokp1
             !
             !  p := r + beta (p - omega A p')
             !
             call info%vecvec%set('start')
             call daxpy(nequ,-omega,apk,1,pk,1)
             call dxpay(nequ,resid,1,beta,pk,1)
             call info%vecvec%set('stop')
          end if
       end if
       !
       ! optional print convergence profile
       !
       if(iprt.ne.0) write(lun_out,'(i4,1pe22.14)') iter,resnorm
    end do

    info%resnorm = resnorm
    info%iter = iter
    
    if(info%iter.ge.imax) info%ierr=1
    
    !
    ! compute final residum
    !
    call info%matvec%set('start')
    call matrix%matrix_times_vector(sol,resid,info%ierr,lun_err)
    call info%matvec%set('stop')

    call info%vecvec%set('start')
    scr = rhs
    resid = scr - resid
    call info%vecvec%set('stop')
    !
    call info%scalprod%set('start')
    info%resreal=dnrm2(nequ,resid,1)*inverse_residum_weight
    call info%scalprod%set('stop')

    !
    ! free memory if required
    !
    call info%ovh%set('start')
    if ( aux_work%is_initialized ) call aux_work%kill(lun_err)  
    aux_loc => null() 
    call info%ovh%set('stop')
    call info%tot%set('stop')




!!$ THESE LINES ARE LEFT FOR DEBUG PURPOSE
!!$    deallocate(&
!!$         resid,&
!!$         scr,&
!!$         res0,&
!!$         pk,&
!!$         pkprime,&
!!$         apk,&
!!$         ask,&
!!$         skprime,&
!!$         stat=res)
!!$    if(res .ne. 0) rc = IOerr(lun_err, err_dealloc, 'bicgstab_solver', &
!!$         'works arrays',res)

  
  end subroutine bicgstab_solver
    
  function loc_normres(nequ,ndir,noddir,vec,scr) result(norm)
    use Globals
    implicit none
    integer,           intent(in   ) :: nequ
    integer,           intent(in   ) :: ndir
    integer,           intent(in   ) :: noddir(ndir)
    real(kind=double), intent(in   ) :: vec(nequ)
    real(kind=double), intent(inout) :: scr(nequ)
    real(kind=double) :: norm
    !local 
    integer :: i
    real(kind=double) :: dnrm2


    if (ndir.eq.0) then
       norm=dnrm2(nequ,vec,1)
    else
       call dcopy(nequ,vec,1,scr,1)
       do i=1,ndir 
          scr(noddir(i))=zero
       end do
       norm=dnrm2(nequ,scr,1)
    end if
  end function loc_normres

  subroutine dxpay(n,dx,incx,da,dy,incy)
    use Globals
    !
    !     constant times a vector plus a vector.
    !     y:=x+a*y
    !     uses unrolled loops for increments equal to one.
    !     jack dongarra, linpack, 3/11/78.
    !     modified 12/3/93, array(1) declarations changed to array(*)
    !
    
    implicit none
    integer :: i,incx,incy,ix,iy,m,mp1,n
    real(kind=double) ::  dx(n),dy(n),da
    
    if ( ( incx .ne. 1) .or. (incy .ne. 1) )  then
       write(*,*) 'x or y increment not equal 1'
       stop
    end if
    dy=dx+da*dy

       
!!$    !
!!$    if(n.le.0)return
!!$    if (da .eq. 0.0d0) then
!!$       call dcopy(n,dx,1,dy,1)
!!$       return
!!$    end if
!!$    if(incx.eq.1.and.incy.eq.1) go to 20
!!$    !
!!$    ! code for unequal increments or equal increments
!!$    ! not equal to 1
!!$    !
!!$    ix = 1
!!$    iy = 1
!!$    if(incx.lt.0)ix = (-n+1)*incx + 1
!!$    if(incy.lt.0)iy = (-n+1)*incy + 1
!!$    do i = 1,n
!!$       dy(iy) = dx(iy) + da*dy(ix)
!!$       ix = ix + incx
!!$       iy = iy + incy
!!$    end do
!!$    return
!!$    !
!!$    ! code for both increments equal to 1
!!$    !
!!$    !
!!$    ! clean-up loop
!!$    !
!!$20  m = mod(n,4)
!!$    if( m .eq. 0 ) go to 40
!!$    do i = 1,m
!!$       dy(i) = dx(i) + da*dy(i)
!!$    end do
!!$    if( n .lt. 4 ) return
!!$40  mp1 = m + 1
!!$    do 50 i = mp1,n,4
!!$       dy(i) = dx(i) + da*dy(i)
!!$       dy(i + 1) = dx(i + 1) + da*dy(i + 1)
!!$       dy(i + 2) = dx(i + 2) + da*dy(i + 2)
!!$       dy(i + 3) = dx(i + 3) + da*dy(i + 3)
!!$50     continue
!!$       return
!!$    end do
    
  end subroutine dxpay

  subroutine jacobi_seidel_solver(matrix, rhs, sol,&
       info, ctrl, &
       prec,&
       aux,&
       residual_norm)
    use Globals
    use LinearOperator
    use Scratch
    use Norms
    class(abs_linop),           intent(inout) :: matrix
    real(kind=double),           intent(in   ) :: rhs(matrix%nrow)
    real(kind=double),           intent(inout) :: sol(matrix%ncol)
    class(output_solver),        intent(inout) :: info
    class(input_solver),         intent(in   ) :: ctrl
    class(abs_linop),   optional, intent(inout) :: prec
    type(scrt), optional,target, intent(inout) :: aux
    class(abs_norm), optional,target, intent(inout) :: residual_norm
    !local
    logical :: rc
    logical :: exit_test
    integer :: nequ,ndir=0
    integer :: info_prec
    class(abs_norm), pointer :: local_norm
    type(euc_norm), target   :: euclidian_norm
    type(scrt), target  :: aux_work
    type(scrt), pointer :: aux_loc
    real(kind=double) ,pointer :: pres(:), resid(:)
    
    real(kind=double) :: inverse_residum_weight
    real(kind=double) :: dnrm2


    !*********************************************************
    ! Setup parts 
    !*********************************************************
    !
    ! Set metric to test the residuum and the norm of rhs
    !
    if( present ( residual_norm ) ) then
       local_norm => residual_norm
    else
       local_norm => euclidian_norm
    end if

    !
    ! Set work arrays
    !    
    if ( present(aux) ) then
       aux_loc => aux
    else
       call aux_work%init(ctrl%lun_err,0,2*nequ)
       aux_loc => aux_work  
    end if

    aux_loc%raux = zero
    resid        => aux_loc%raux(       1 :   nequ)
    pres         => aux_loc%raux(  nequ+1 : 2*nequ)

    !
    ! Set to zero info and counter
    !
    info%ierr = 0
    exit_test = .false.
    info%iter = 0

   
    !
    ! Set up initial constants and start loop
    !
    exit_test = .false.
    if (ctrl%iprt.gt.0) write(ctrl%lun_out,'(a)') '      iter        resnorm'
    
    !****************************************************
    ! Start computation 
    !***************************************************

    !
    ! eval rhs norm
    !
    info%bnorm = dnrm2(nequ,rhs,1)

    !
    ! handles case of zero RHS
    !
    call zerorhs(nequ,info%ierr,&
         info%bnorm,info%resnorm,sol,rhs,exit_test)
           
    !
    ! handles case of zero RHS and properly initializes exit_test
    !
    if (info%bnorm.eq.zero) then
       call zerorhs(nequ,info%ierr,info%bnorm,&
            info%resnorm,sol,rhs,exit_test)
       if (exit_test) return
    end if

    !
    ! set inverse_residum_weight
    ! inverse_residum_weight = one       => resnorm = |Ax-b|
    ! inverse_residum_weight = one/bnrom => resnorm = |Ax-b|/|b|
    !
    if (ctrl%iexit.eq.0) then
       inverse_residum_weight = one
    else
       inverse_residum_weight = one / info%bnorm
    end if
    
    !
    ! calculate initial residual (res = b-M*x)
    !      
    call info%matvec%set('start')
    call matrix%matrix_times_vector(sol,resid,info%ierr,ctrl%lun_err)
    call info%matvec%set('stop')

    call info%vecvec%set('start')
    resid = rhs - resid
    call info%vecvec%set('stop')
        
    !
    ! compute initial norm of residum
    !
    call info%scalprod%set('start')
    info%resini = dnrm2(nequ,resid,1)*inverse_residum_weight
    call info%scalprod%set('stop')

    !
    ! exit if the initial residuum is already below the required tolerance 
    !
    if ( info%resini .lt. ctrl%tol_sol) then
       !
       info%resreal = info%resini
       info%resnorm = zero
       !
       ! free memory
       !
       call info%ovh%set('start')
       if ( aux_work%is_initialized ) call aux_work%kill(ctrl%lun_err)  
       aux_loc => null() 
       call info%ovh%set('stop')
       call info%tot%set('stop ')
       return
    end if
    
    !
    ! cycle
    ! 
    do while (.not. exit_test)
       info%iter = info%iter + 1   
       !
       ! compute  pres = PREC  (r_k+1)
       !
       info_prec = 0
       call info%prec%set('start')       
       call prec%matrix_times_vector(resid,pres,info_prec,ctrl%lun_err)
       call info%prec%set('stop')
       if ( info_prec .ne. 0 ) then
          rc = IOerr(ctrl%lun_err, wrn_val, ' jacobi_seidel_solver', &
               ' error computating P^{-1} r_{k+1},'//&
               ' info prec = ',info_prec )
          info%ierr = 3
          info%info_prec = info_prec
          return
       end if

       
       
       !
       !  calculates x_k+1:=x_k + pres_k
       !
       call info%vecvec%set('start')
       call dxpay(nequ,pres,1,one,sol,1)
       call info%vecvec%set('stop')
          
       !
       ! calculate initial residual (res = b-M*x)
       !      
       call info%matvec%set('start')
       call matrix%matrix_times_vector(sol,resid,info%ierr,ctrl%lun_err)
       call info%matvec%set('stop')

       call info%vecvec%set('start')
       resid = rhs - resid
       call info%vecvec%set('stop')

       !
       !  compute residum
       !       
       call info%scalprod%set('start')
       info%resnorm = dnrm2(nequ,resid,1)*inverse_residum_weight
       call info%scalprod%set('stop')

       ! optional printing of convergence rate
       if (ctrl%iprt.ge.1) write(ctrl%lun_out,'(i10,e15.7)') &
            info%iter,info%resnorm

       exit_test = (info%iter.gt.ctrl%imax .or. info%resnorm .le. ctrl%tol_sol)
       if (info%iter.ge. ctrl%imax) info%ierr = 1
    end do

    !
    ! fill info
    !
    info%resreal = info%resnorm


    !
    ! free memory if required
    !
    call info%ovh%set('start')
    if ( aux_work%is_initialized ) call aux_work%kill(ctrl%lun_err)  
    aux_loc => null() 
    call info%ovh%set('stop')
    call info%tot%set('stop ')

  end subroutine jacobi_seidel_solver
  

!!$  !>-------------------------------------------------------------
!!$  !> Linear solver procedure.
!!$  !> Solves linear sysstem in the form A x = b
!!$  !> where A is any type of matrix for which the procedure
!!$  !> Mxv        => matrix vector multiplication
!!$  !> diag_scale => scaling by a diagonal matrix
!!$  !>
!!$  !> usage:
!!$  !>     call linear_solver(matrix, rhs,  sol, &
!!$  !>     info_solver, &
!!$  !>     ctrl_solver,&
!!$  !>     prec_left,prec_right,&
!!$  !>     aux,&
!!$  !>     abpres)
!!$  !> 
!!$  !> usage:
!!$  !>     call 'var'%matrix_times_vector(vec_in,vec_out,[info])
!!$  !>
!!$  !> where 
!!$  !> \param[inout] matrix            -> Abstact matrix, with 
!!$  !>                                    matrix-vector multiplication and 
!!$  !>                                    matrix(transpose)-vector multipl.
!!$  !> \param[in   ] rhs               -> real, dimension(matrix%nrow)
!!$  !>                                    RHS of the system A x = b
!!$  !> \param[inout] sol               -> real, dimension(matrix%ncol)
!!$  !>                                    solution of the system A x = b 
!!$  !> \param[inout] info_solver       -> type(output_solver)
!!$  !>                                    Info parameters of linear solver 
!!$  !> \param[in   ] ctrl_solver       -> type(input_solver)
!!$  !>                                    Controls parameters of linear solver 
!!$  !> \param[inout] (optin.) prec_left-> Abstrac precoditioner
!!$  !>                                    with apply operation. 
!!$  !>                                    Prec. used in PCG procedure.
!!$  !>                                    Left prec. used in other procedures.
!!$  !> \param[inout] (option.) prec_right-> Abstrac precoditioner
!!$  !>                                    with apply operation. 
!!$  !>                                    Not used in PCG procedure.
!!$  !>                                    right prec. used in other procedures.
!!$  !> \param[inout] (option.) abpres  -> type((pcgcoeff)
!!$  !>                                    It stores alpha, beta and 
!!$  !>                                    preconditioned residum in PCG
!!$  !
!!$  !<-------------------------------------------------------------
!!$  subroutine multigrid_linear_solver(&
!!$       nlevel,&
!!$       matrices, rhss,  sol, &
!!$       info_solvers, &
!!$       ctrl_solvers,&
!!$       precs_left,precs_right,&
!!$       aux,&
!!$       abpres)
!!$    use Scratch
!!$    integer,                         intent(in   ) :: nlevel
!!$    class(abs_linop),                intent(inout) :: matrix(nlevel)
!!$    real(kind=double),               intent(in   ) :: rhs(matrices%(level)%ncol)
!!$    real(kind=double),               intent(inout) :: sol(matrices%(level)%ncol)
!!$    type(array_mat),                 intent(inout) :: pre_smoother(nlevel-1)
!!$    type(array_mat),                 intent(inout) :: post_smoother(nlevel-1)
!!$    type(array_mat),                 intent(inout) :: restriction_operators(nlevel-1)
!!$    type(array_mat),                 intent(inout) :: prolongator_operators(nlevel-1)
!!$    class(abs_linop),                intent(inout) :: coaster_matrix_solver
!!$    type(output_solver),             intent(inout) :: info_solver
!!$    type(input_solver),              intent(in   ) :: ctrl_solver
!!$    class(array_prec),target,optional, intent(inout) :: precs_left
!!$    class(array_prec),target,optional, intent(inout) :: precs_right
!!$    type(scrt),            optional, intent(inout) :: aux
!!$    type(pcgcoeff),        optional, intent(inout) :: abpres
!!$    ! local
!!$    integer :: lun_err, lun_out
!!$    character(len=20) ::  scheme
!!$
!!$    !
!!$    ! set logical unit for err and output
!!$    !
!!$    lun_err = ctrl_solver%lun_err
!!$    lun_out = ctrl_solver%lun_out
!!$
!!$    !
!!$    ! V cycle
!!$    !
!!$    call matrices(1)%mat%matrix_times_vector(sol,res)
!!$    res=rhs-res
!!$    
!!$    !
!!$    ! cycle to the coaser matrix
!!$    !
!!$    do i=1, level-1
!!$       call linear_solver(&
!!$            matrices%i%mat, rhss(i),  sol, &
!!$            info_solver, &
!!$            ctrl_solver,&
!!$            precs_left,precs_right,&
!!$            aux)
!!$       call projectors(i)%matrix_times_vector(sol, proj_sol)
!!$       
!!$       sol=proj_sol
!!$       
!!$    end do
!!$    
!!$    !
!!$    ! exact solution at the coastest level   
!!$    !
!!$    call linear_solver(&
!!$            matrices(nlevel)%mat, rhss(nlevel),  sol, &
!!$            infos_solver(nlevel, &
!!$            ctrl_solver,&
!!$            precs_left,precs_right,&
!!$            aux)
!!$
!!$    !
!!$    ! cycle to the finest  matrix
!!$    !
!!$    do i=level, 1,-1 
!!$
!!$       call linear_solver(&
!!$            matrices%i%mat, rhss(i),  sol, &
!!$            info_solver, &
!!$            ctrl_solver,&
!!$            precs_left,precs_right,&
!!$            aux)
!!$
!!$       call interpolator(i)%matrix_times_vector(sol, proj_sol)
!!$       
!!$       sol=proj_sol
!!$       
!!$    end do
!!$
!!$
!!$
!!$  end subroutine multigrid_linear_solver

!>------------------------------------------------------------------
!>
!>     MINRES algorithm developed from code written by Michael A. Saunders
!>     and other. Below we report the credits and the description of
!>     the original source code.
!>
!>------------------------------------------------------------------
!>
!>     MINRES is an implementation of the algorithm described in
!>     the following reference:
!>
!>     C. C. Paige and M. A. Saunders (1975),
!>     Solution of sparse indefinite systems of linear equations,
!>     SIAM J. Numer. Anal. 12(4), pp. 617-629.
!>------------------------------------------------------------------
!>
!>
!>     MINRES development:
!>            1972: First version, similar to original SYMMLQ.
!>                  Later lost @#%*!
!>        Oct 1995: Tried to reconstruct MINRES from
!>                  1995 version of SYMMLQ.
!>     30 May 1999: Need to make it more like LSQR.
!>                  In middle of major overhaul.
!>     19 Jul 2003: Next attempt to reconstruct MINRES.
!>                  Seems to need two vectors more than SYMMLQ.  (w1, w2)
!>                  Lanczos is now at the top of the loop,
!>                  so the operator Aprod is called in just one place
!>                  (not counting the initial check for symmetry).
!>     22 Jul 2003: Success at last.  Preconditioning also works.
!>                  minres.f added to http://www.stanford.edu/group/SOL/.
!>
!>     FUTURE WORK: A stopping rule is needed for singular systems.
!>                  We need to estimate ||Ar|| as in LSQR.  This will be
!>                  joint work with Sou Cheng Choi, SCCM, Stanford.
!>                  Note that ||Ar|| small => r is a null vector for A.
!>
!>
!>     Michael A. Saunders           na.msaunders@na-net.ornl.gov
!>     Department of MS&E            saunders@stanford.edu
!<     Stanford University
!>     Stanford, CA 94305-4026       (650) 723-1875
!>------------------------------------------------------------------





  recursive subroutine minres_solver(matrix, rhs, sol,&
       info, ctrl, &
       prec,&
       aux_minres,&
       residual_norm,&
       ortogonalization_matrix)
    use Globals
    use LinearOperator
    use Scratch
    use Norms
    class(abs_linop),           intent(inout) :: matrix
    real(kind=double),           intent(in   ) :: rhs(matrix%nrow)
    real(kind=double),           intent(inout) :: sol(matrix%ncol)
    class(output_solver),        intent(inout) :: info
    class(input_solver),         intent(in   ) :: ctrl
    class(abs_linop),   optional, intent(inout) :: prec
    type(scrt), optional,target, intent(inout) :: aux_minres
    class(abs_norm), optional,target, intent(inout) :: residual_norm
    class(abs_linop), optional,  intent(inout) :: ortogonalization_matrix

    !
    ! local vars
    !
    
    ! logical
    logical :: rc
    logical :: exit_test
    
    ! strings
    character(len=256) :: message


    ! integer 
    integer :: res
    integer :: nequ
    integer :: dim_ker
    integer :: info_prec
    integer :: lun_err
    integer :: lun_out
    integer :: iort
    integer :: ndir ! useless varible

    ! real
    real(kind=double) :: inverse_residum_weight
    real(kind=double) :: rort
    
    ! derived type
    type(scrt), target  :: aux_work
    type(scrt), pointer :: aux_loc
    class(abs_norm), pointer :: local_norm
    type(euc_norm), target   :: euclidian_norm
    
    ! local copy of in/out arguemnts 
    ! argument in orignal code
    integer           :: n, nout, itnlim, istop, itn
    logical           :: checkA, precon
    real(kind=double) :: shift, rtol, Anorm, Acond, rnorm, ynorm
    real(kind=double), pointer :: x(:),y(:),r1(:), r2(:)
     real(kind=double), pointer :: v(:), w(:),w1(:), w2(:),my_work(:)
    ! local copy of local arguments of original code
    real(kind=double)  alfa  , beta  , beta1 , cs 
    real(kind=double)  dbar  , delta , denom , diag  
    real(kind=double)  eps   , epsa  , epsln , epsr  , epsx
    real(kind=double)  gamma , gbar  , gmax  , gmin 
    real(kind=double)  oldb  , oldeps, qrnorm, phi   , phibar
    real(kind=double)  rhs1  , rhs2  , s     , sn    , t     
    real(kind=double)  tnorm2, ynorm2, z
    integer            i,ibegin,iend
    logical            debug, prnt

    real(kind=double) :: ten = 1.0d1
    real(kind=double) :: ddot
    real(kind=double) :: dnrm2

    character*16       enter, exit
    character*52       msg(-1:8)
    character(len=256)  :: err_msg

    !
    ! local copy of controls
    !
    lun_err = ctrl%lun_err
    lun_out = ctrl%lun_out
    
    !
    ! symmtretry check
    !
    if (.not. matrix%is_symmetric) then
       rc= IOerr(lun_err, wrn_val, 'minres_solver', &
            ' non symmetric matrix ')
       info%ierr = 999
       return
    end if

    if  (.not. prec%is_symmetric)  then
       rc= IOerr(lun_err, wrn_val, 'minres_solver', &
            ' non symmetric preconditioner ')
       info%ierr = 999
       return
    end if


    call info%tot%set('start')
    call info%ovh%set('start')
    !
    ! set dimensions
    !
    nequ    = matrix%ncol

    
    iort    = ctrl%iort

    !
    ! set work arrays
    !    
    if ( present(aux_minres) ) then
       if( .not. aux_minres%check(0,9*nequ) ) &
            rc = IOerr(lun_err, wrn_inp, 'minres_solver', &
            ' aux array too small',res)
       aux_loc => aux_minres
    else
       call aux_work%init(lun_err,0,9*nequ)
       aux_loc => aux_work  
    end if
    iend=0
    call aux_loc%range(nequ,ibegin,iend)
    x =>aux_loc%raux(ibegin:iend)
    call aux_loc%range(nequ,ibegin,iend)
    y =>aux_loc%raux(ibegin:iend)
    call aux_loc%range(nequ,ibegin,iend)
    r1 =>aux_loc%raux(ibegin:iend)
    call aux_loc%range(nequ,ibegin,iend)
    r2 =>aux_loc%raux(ibegin:iend)
    call aux_loc%range(nequ,ibegin,iend)
    v =>aux_loc%raux(ibegin:iend)
    call aux_loc%range(nequ,ibegin,iend)
    w =>aux_loc%raux(ibegin:iend)
    call aux_loc%range(nequ,ibegin,iend)
    w1 =>aux_loc%raux(ibegin:iend)
    call aux_loc%range(nequ,ibegin,iend)
    w2 =>aux_loc%raux(ibegin:iend)
    call aux_loc%range(nequ,ibegin,iend)
    my_work =>aux_loc%raux(ibegin:iend)
    
    !
    ! set to zero 
    !

    info%ierr = 0
    info%iter = 0

    !
    ! set metric to test the residuum and the norm of rhs
    !
    if( present ( residual_norm ) ) then
       local_norm => residual_norm
    else
       local_norm => euclidian_norm
    end if
    

    !
    ! In case of ortogonalization w.r.t to the matrix kernel 
    ! control if rhs is ortogonal to the kernel 
    ! up to tolerance required 
    !
    if ( (iort .ne. 0           ) .and. &
         present( ortogonalization_matrix ) ) then 
       if (matrix%is_symmetric   ) then
          call ortogonalization_matrix%matrix_times_vector(rhs,my_work,info%ierr,lun_err)
          rort=dnrm2(nequ, my_work,1)
          if  ( rort >  ctrl%tol_sol ) then
             write(err_msg,'(1e12.5)') rort
             rc = IOerr(lun_err, wrn_val, ' pcg_solver', &
                  'RHS not ortogonal to ortogonalization matrix'&
                  //etb(err_msg))
             info%ierr = 1
             return
          end if
       end if
    end if
    
    !
    ! compute rhs norm and set inverse_residum_weight
    !
    info%bnorm = dnrm2(nequ,rhs,1)

    !
    ! handles case of zero RHS
    !
    call zerorhs(nequ,info%ierr,&
         info%bnorm,info%resnorm,sol,rhs,exit_test)

    !
    ! set inverse_residum_weight
    ! inverse_residum_weight = one       => resnorm = |Ax-b|
    ! inverse_residum_weight = one/bnrom => resnorm = |Ax-b|/|b|
    !
    if (ctrl%iexit.eq.0) then
       inverse_residum_weight = one
    else
       inverse_residum_weight = one / info%bnorm
    end if

    !
    call info%ovh%set('stop')

    !
    ! optional application of preconditoner to the rhs
    !
    info_prec = 0
    if ( ctrl%isol .eq. 1 ) then
       info_prec = 0
       call info%prec%set('start')
       call prec%matrix_times_vector(rhs,sol,info_prec,lun_err)
       call info%prec%set('stop')
       if ( info_prec .ne. 0 ) then
          rc = IOerr(lun_err, wrn_val, ' pcg_solver', &
               ' error computating sol = p^{-1} rhs ( isol =1),'//&
               ' info prec = ',info_prec)
          info%ierr = 3
          info%info_prec = info_prec
          return
       end if

    end if

    
    !
    ! calculate initial residual (res = b-M*x)
    !      
    call info%matvec%set('start')
    call matrix%matrix_times_vector(sol,r1,info%ierr,lun_err)
    call info%matvec%set('stop')

    call info%vecvec%set('start')
    r1 = rhs - r1
    call info%vecvec%set('stop')
        
    !
    ! compute initial norm of residum
    !
    call info%scalprod%set('start')
    info%resini = dnrm2(nequ,r1,1)*inverse_residum_weight
    call info%scalprod%set('stop')

    if ( isnan(info%resini) ) then
       rc = IOerr(lun_err, wrn_val, ' minres_solver', &
            'resini = NaN')
       info%ierr =-1 
       return
    end if


    !
    ! exit if the initial residuum is already below the required tolerance 
    !
    if ( info%resini .lt. ctrl%tol_sol) then
       !
       info%resreal = info%resini
       info%resnorm = zero
       !
       ! free memory
       !
       call info%ovh%set('start')
       if ( aux_work%is_initialized ) call aux_work%kill(lun_err)  
       aux_loc => null() 
       call info%ovh%set('stop')
       call info%tot%set('stop ')
       return
    end if




!!$      subroutine MINRES( n, b, r1, r2, v, w, w1, w2, x, y,
!!$     $                   Aprod, Msolve, checkA, precon, shift,
!!$     $                   nout , itnlim, rtol,
!!$     $                   istop, itn, Anorm, Acond, rnorm, ynorm )
    !
    
    n=nequ
    x=sol


    !inputs
    checkA = (ctrl%debug>1)
    precon = .True.
    shift  = zero
    
    nout = 0
    itnlim = ctrl%imax
    rtol   = ctrl%tol_sol

!     
!    see outputs  after main loop
!
!    info%info = istop
!    info%iter = itn
!    info%resnorm = rnorm
!    ! Anorm, Acond, ynorm useful but not passed
!    info%matrix_cond = Acond
    
! ------------------------------------------------------------------
!
!     MINRES  is designed to solve the system of linear equations
!
!                Ax = b
!
!     or the least-squares problem
!
!         min || Ax - b ||_2,
!
!     where A is an n by n symmetric matrix and b is a given vector.
!     The matrix A may be indefinite and/or singular.
!
!     1. If A is known to be positive definite, the Conjugate Gradient
!        Method might be preferred, since it requires the same number
!        of iterations as MINRES but less work per iteration.
!
!     2. If A is indefinite but Ax = b is known to have a solution
!        (e.g. if A is nonsingular), SYMMLQ might be preferred,
!        since it requires the same number of iterations as MINRES
!        but slightly less work per iteration.
!
!     The matrix A is intended to be large and sparse.  It is accessed
!     by means of a subroutine call of the form
!     SYMMLQ development:
!
!                call Aprod ( n, x, y )
!
!     which must return the product y = Ax for any given vector x.
!
!
!     More generally, MINRES is designed to solve the system
!
!                (A - shift*I) x = b
!     or
!         min || (A - shift*I) x - b ||_2,
!
!     where  shift  is a specified scalar value.  Again, the matrix
!     (A - shift*I) may be indefinite and/or singular.
!     The work per iteration is very slightly less if  shift = 0.
!
!     Note: If  shift  is an approximate eigenvalue of  A
!     and  b  is an approximate eigenvector,  x  might prove to be
!     a better approximate eigenvector, as in the methods of
!     inverse iteration and/or Rayleigh-quotient iteration.
!     However, we're not yet sure on that -- it may be better
!     to use SYMMLQ.
!
!     A further option is that of preconditioning, which may reduce
!     the number of iterations required.  If M = C C' is a positive
!     definite matrix that is known to approximate  (A - shift*I)
!     in some sense, and if systems of the form  My = x  can be
!     solved efficiently, the parameters precon and Msolve may be
!     used (see below).  When  precon = .true., MINRES will
!     implicitly solve the system of equations
!
!             P (A - shift*I) P' xbar  =  P b,
!
!     i.e.                  Abar xbar  =  bbar
!     where                         P  =  C**(-1),
!                                Abar  =  P (A - shift*I) P',
!                                bbar  =  P b,
!
!     and return the solution       x  =  P' xbar.
!     The associated residual is rbar  =  bbar - Abar xbar
!                                      =  P (b - (A - shift*I)x)
!                                      =  P r.
!
!     In the discussion below, eps refers to the machine precision.
!     eps is computed by MINRES.  A typical value is eps = 2.22d-16
!     for IEEE double-precision arithmetic.
!
!     Parameters
!     ----------
!
!     n       input      The dimension of the matrix A.
!
!     b(n)    input      The rhs vector b.
!
!     r1(n)   workspace
!     r2(n)   workspace
!     v(n)    workspace
!     w(n)    workspace
!     w1(n)   workspace
!     w2(n)   workspace
!
!     x(n)    output     Returns the computed solution  x.
!
!     y(n)    workspace
!
!     Aprod   external   A subroutine defining the matrix A.
!                        For a given vector x, the statement
!
!                              call Aprod ( n, x, y )
!
!                        must return the product y = Ax
!                        without altering the vector x.
!
!     Msolve  external   An optional subroutine defining a
!                        preconditioning matrix M, which should
!                        approximate (A - shift*I) in some sense.
!                        M must be positive definite.
!                        For a given vector x, the statement
!
!                              call Msolve( n, x, y )
!
!                        must solve the linear system My = x
!                        without altering the vector x.
!
!                        In general, M should be chosen so that Abar has
!                        clustered eigenvalues.  For example,
!                        if A is positive definite, Abar would ideally
!                        be close to a multiple of I.
!                        If A or A - shift*I is indefinite, Abar might
!                        be close to a multiple of diag( I  -I ).
!
!                        NOTE.  The program calling MINRES must declare
!                        Aprod and Msolve to be external.
!
!     checkA  input      If checkA = .true., an extra call of Aprod will
!                        be used to check if A is symmetric.  Also,
!                        if precon = .true., an extra call of Msolve
!                        will be used to check if M is symmetric.
!
!     precon  input      If precon = .true., preconditioning will
!                        be invoked.  Otherwise, subroutine Msolve
!                        will not be referenced; in this case the
!                        actual parameter corresponding to Msolve may
!                        be the same as that corresponding to Aprod.
!
!     shift   input      Should be zero if the system Ax = b is to be
!                        solved.  Otherwise, it could be an
!                        approximation to an eigenvalue of A, such as
!                        the Rayleigh quotient b'Ab / (b'b)
!                        corresponding to the vector b.
!                        If b is sufficiently like an eigenvector
!                        corresponding to an eigenvalue near shift,
!                        then the computed x may have very large
!                        components.  When normalized, x may be
!                        closer to an eigenvector than b.
!
!     nout    input      A file number.
!                        If nout .gt. 0, a summary of the iterations
!                        will be printed on unit nout.
!
!     itnlim  input      An upper limit on the number of iterations.
!
!     rtol    input      A user-specified tolerance.  MINRES terminates
!                        if it appears that norm(rbar) is smaller than
!                              rtol * norm(Abar) * norm(xbar),
!                        where rbar is the transformed residual vector,
!                              rbar = bbar - Abar xbar.
!
!                        If shift = 0 and precon = .false., MINRES
!                        terminates if norm(b - A*x) is smaller than
!                              rtol * norm(A) * norm(x).
!
!     istop   output     An integer giving the reason for termination...
!
!              -1        beta2 = 0 in the Lanczos iteration; i.e. the
!                        second Lanczos vector is zero.  This means the
!                        rhs is very special.
!                        If there is no preconditioner, b is an
!                        eigenvector of A.
!                        Otherwise (if precon is true), let My = b.
!                        If shift is zero, y is a solution of the
!                        generalized eigenvalue problem Ay = lambda My,
!                        with lambda = alpha1 from the Lanczos vectors.
!
!                        In general, (A - shift*I)x = b
!                        has the solution         x = (1/alpha1) y
!                        where My = b.
!                        
!               0        b = 0, so the exact solution is x = 0.
!                        No iterations were performed.
!
!               1        Norm(rbar) appears to be less than
!                        the value  rtol * norm(Abar) * norm(xbar).
!                        The solution in  x  should be acceptable.
!
!               2        Norm(rbar) appears to be less than
!                        the value  eps * norm(Abar) * norm(xbar).
!                        This means that the residual is as small as
!                        seems reasonable on this machine.
!
!               3        Norm(Abar) * norm(xbar) exceeds norm(b)/eps,
!                        which should indicate that x has essentially
!                        converged to an eigenvector of A
!                        corresponding to the eigenvalue shift.
!
!               4        Acond (see below) has exceeded 0.1/eps, so
!                        the matrix Abar must be very ill-conditioned.
!                        x may not contain an acceptable solution.
!
!               5        The iteration limit was reached before any of
!                        the previous criteria were satisfied.
!
!               6        The matrix defined by Aprod does not appear
!                        to be symmetric.
!                        For certain vectors y = Av and r = Ay, the
!                        products y'y and r'v differ significantly.
!
!               7        The matrix defined by Msolve does not appear
!                        to be symmetric.
!                        For vectors satisfying My = v and Mr = y, the
!                        products y'y and r'v differ significantly.
!
!               8        An inner product of the form  x' M**(-1) x
!                        was not positive, so the preconditioning matrix
!                        M does not appear to be positive definite.
!
!                        If istop .ge. 5, the final x may not be an
!                        acceptable solution.
!
!     itn     output     The number of iterations performed.
!
!     Anorm   output     An estimate of the norm of the matrix operator
!                        Abar = P (A - shift*I) P',   where P = C**(-1).
!
!     Acond   output     An estimate of the condition of Abar above.
!                        This will usually be a substantial
!                        under-estimate of the true condition.
!
!     rnorm   output     An estimate of the norm of the final
!                        transformed residual vector,
!                           P (b  -  (A - shift*I) x).
!
!     ynorm   output     An estimate of the norm of xbar.
!                        This is sqrt( x'Mx ).  If precon is false,
!                        ynorm is an estimate of norm(x).
!
!
!     Subroutines and functions
!
!     USER       Aprod , Msolve
!     BLAS1      daxpy , dcopy , ddot  , dnrm2  } These are all in
!     Utilities  daxpy2, dload2, dscal2         } the file minresblas.f




!!$      data               enter /' Enter MINRES.  '/,
!!$     $                   exit  /' Exit  MINRES.  '/

!!$      data               msg
!!$      / 'beta2 = 0.  If M = I, b and x are eigenvectors of A',&
!!$        'beta1 = 0.  The exact solution is  x = 0',&
!!$        'Requested accuracy achieved, as determined by rtol',&
!!$        'Reasonable accuracy achieved, given eps',&
!!$        'x has converged to an eigenvector',&
!!$        'Acond has exceeded 0.1/eps',&
!!$        'The iteration limit was reached',&
!!$        'Aprod  does not define a symmetric matrix',&
!!$        'Msolve does not define a symmetric matrix',&
!!$        'Msolve does not define a pos-def preconditioner' /
!     ------------------------------------------------------------------
    debug = .false.
    if ( ctrl%debug .eq. 1) debug = .True.

!     ------------------------------------------------------------------
!     Compute eps, the machine precision.  The call to daxpy is
!     intended to fool compilers that use extra-length registers.
!     31 May 1999: Hardwire eps so the debugger can step thru easily.
!     ------------------------------------------------------------------
      eps    = 2.22d-16    ! Set eps = zero here if you want it computed.
      if (eps .le. zero) then
         eps    = two**(-12)
   10       eps    = eps / two
            x(1)   = eps
            y(1)   = one
            call daxpy ( 1, one, x, 1, y, 1 )
            if (y(1) .gt. one) go to 10
         eps    = eps * two
      end if

!     ------------------------------------------------------------------
!     Print heading and initialize.
!     ------------------------------------------------------------------
      if (nout .gt. 0) then
         write(nout, 1000) enter, n, checkA, precon,&
                          itnlim, rtol, shift
      end if
      istop  = 0
      itn    = 0
      Anorm  = zero
      Acond  = zero
      rnorm  = zero
      ynorm  = zero
      
      ! original code
      !
      ! call dload2( n, zero, x )
      x=zero

!     ------------------------------------------------------------------
!     Set up y and v for the first Lanczos vector v1.
!     y  =  beta1 P' v1,  where  P = C**(-1).
!     v is really P' v1.
!     ------------------------------------------------------------------
      call dcopy ( n, rhs, 1, y , 1 )         ! y  = b
      call dcopy ( n, rhs, 1, r1, 1 )         ! r1 = b
      
      ! original code
      !
      ! if ( precon ) call Msolve( n, b, y )

      !
      ! apply preconditioner
      !
      if ( precon) then
         info_prec = 0
         call info%prec%set('start')
         call prec%matrix_times_vector(rhs,y,info_prec,ctrl%lun_err)
         call info%prec%set('stop')
         if ( info_prec .ne. 0 ) then
            rc = IOerr(ctrl%lun_err, wrn_val, ' pcg_solver', &
                 ' error computating y = p^{-1} rhs,'//&
                 ' info prec = ',info_prec)
            info%ierr = 3
            info%info_prec = info_prec
            return
         end if
         
      end if
      
      beta1  = ddot  ( n, rhs, 1, y, 1 )
      if (beta1 .lt. zero) then    ! M must be indefinite.
         istop = 8
         go to 900
      end if
      
      ! original code
      !
      !if (beta1 .eq. zero) then    ! b = 0 exactly.  Stop with x = 0.
      !   istop = 0
      !   go to 900
      !end if
      !
      ! Nothing to be done, already handled before 

      beta1  = sqrt( beta1 )       ! Normalize y to get v1 later.
      
!     ------------------------------------------------------------------
!     See if Msolve is symmetric.
!     ------------------------------------------------------------------
      if (checkA  .and.  precon) then
         ! original code
         !
         !call Msolve( n, y, r2 )
         info_prec = 0
         call prec%matrix_times_vector(y,r2,info_prec,ctrl%lun_err)
         if ( info_prec .ne. 0 ) then
            rc = IOerr(ctrl%lun_err, wrn_val, ' pcg_solver', &
                 ' error computating sol = p^{-1} rhs ( isol =1),'//&
                 ' info prec = ',info_prec)
            info%ierr = 3
            info%info_prec = info_prec
            return
         end if
         s      = ddot  ( n, y, 1, y, 1 )
         t      = ddot  ( n,r1, 1,r2, 1 )
         z      = abs( s - t )
         epsa   = (s + eps) * eps**0.33333d+0
         if (z .gt. epsa) then
            istop = 7
            rc = IOerr(ctrl%lun_err, wrn_val, ' minres_solver', &
                 ' preconditioner passed is not symmetric ')
            info%ierr = -1
            return
         end if
      end if

!     ------------------------------------------------------------------
!     See if Aprod  is symmetric.
!     ------------------------------------------------------------------
      if (checkA) then
         ! original code
         !
         !call Aprod ( n, y, w )
         !call Aprod ( n, w, r2 )
         !
         call info%matvec%set('start')
         call matrix%matrix_times_vector(y,w,info%ierr,lun_err)
         call matrix%matrix_times_vector(w,r2,info%ierr,lun_err)
         call info%matvec%set('stop')
         
         s      = ddot  ( n, w, 1, w, 1 )
         t      = ddot  ( n, y, 1,r2, 1 )
         z      = abs( s - t )
         epsa   = (s + eps) * eps**0.33333d+0
         if (z .gt. epsa) then
            istop = 6
            rc = IOerr(ctrl%lun_err, wrn_val, ' minres_solver', &
                 ' Matrix is not symmetric ')
            info%ierr = -1
            return
            ! go to 900
         end if
      end if

!     ------------------------------------------------------------------
!     Initialize other quantities.
!     ------------------------------------------------------------------
      oldb   = zero
      beta   = beta1
      dbar   = zero
      epsln  = zero
      qrnorm = beta1
      phibar = beta1
      rhs1   = beta1
      rhs2   = zero
      tnorm2 = zero
      ynorm2 = zero
      cs     = - one
      sn     = zero
      ! original code
      !
      !call dload2( n, zero, w  )        ! w  = 0
      !call dload2( n, zero, w2 )        ! w2 = 0
      w=0
      w2=0
      call dcopy ( n, r1, 1, r2, 1 )    ! r2 = r1

      if (debug) then
         write(*,*) ' '
         write(*,*) 'beta ', beta
         write(*,*) ' '
      end if

!     ------------------------------------------------------------------
!     Main iteration loop.
!     ------------------------------------------------------------------
      ! oginal code 
      !100   itn    = itn   +  1               ! k = itn = 1 first time through
      !if (istop .ne. 0) go to 900
      !
      if (ctrl%iprt.gt.0) write(ctrl%lun_out,'(a)') '      iter        resnorm'
      do while (istop .eq. 0) 
         itn    = itn   +  1   
         !-----------------------------------------------------------------
         ! Obtain quantities for the next Lanczos vector vk+1, k = 1, 2,...
         ! The general iteration is similar to the case k = 1 with v0 = 0:
         !
         !   p1      = Operator * v1  -  beta1 * v0,
         !   alpha1  = v1'p1,
         !   q2      = p2  -  alpha1 * v1,
         !   beta2^2 = q2'q2,
         !   v2      = (1/beta2) q2.
         !
         ! Again, y = betak P vk,  where  P = C**(-1).
         ! .... more description needed.
         !-----------------------------------------------------------------
         s      = one / beta            ! Normalize previous vector (in y).
         call dscal2( n, s, y, v )      ! v = vk if P = I
         
         ! original code
         !
         !call Aprod ( n, v, y )
         call info%matvec%set('start')
         call matrix%matrix_times_vector(v,y,info%ierr,lun_err)
         call info%matvec%set('stop')

         call daxpy ( n, (- shift), v, 1, y, 1 )
         if (itn .ge. 2) then
            call daxpy ( n, (- beta/oldb), r1, 1, y, 1 )
         end if

         alfa   = ddot  ( n, v, 1, y, 1 )     ! alphak

         call daxpy ( n, (- alfa/beta), r2, 1, y, 1 )
         call dcopy ( n, r2, 1, r1, 1 )
         call dcopy ( n,  y, 1, r2, 1 )

         ! original code
         !
         ! if ( precon ) call Msolve( n, r2, y )
         if ( precon) then
            call prec%matrix_times_vector(r2,y,info_prec,ctrl%lun_err)
            if ( info_prec .ne. 0 ) then
               rc = IOerr(ctrl%lun_err, wrn_val, ' pcg_solver', &
                    ' error computating sol = p^{-1} rhs ( isol =1),'//&
                    ' info prec = ',info_prec)
               info%ierr = 3
               info%info_prec = info_prec
               return
            end if
         end if

         oldb   = beta                        ! oldb = betak
         beta   = ddot  ( n, r2, 1, y, 1 )    ! beta = betak+1^2
         if (beta .lt. zero) then
            istop = 6
            go to 900
         end if

         beta   = sqrt( beta )                ! beta = betak+1
         tnorm2 = tnorm2 + alfa**2 + oldb**2 + beta**2

         if (itn .eq. 1) then                 ! Initialize a few things.
            if (beta/beta1 .le. ten*eps) then ! beta2 = 0 or ~ 0.
               istop = -1                     ! Terminate later.
            end if
            !tnorm2 = alfa**2
            gmax   = abs( alfa )              ! alpha1
            gmin   = gmax                     ! alpha1
         end if

         ! Apply previous rotation Qk-1 to get
         !   [deltak epslnk+1] = [cs  sn][dbark    0   ]
         !   [gbar k dbar k+1]   [sn -cs][alfak betak+1].

         oldeps = epsln
         delta  = cs * dbar  +  sn * alfa ! delta1 = 0         deltak
         gbar   = sn * dbar  -  cs * alfa ! gbar 1 = alfa1     gbar k
         epsln  =               sn * beta ! epsln2 = 0         epslnk+1
         dbar   =            -  cs * beta ! dbar 2 = beta2     dbar k+1

         ! Compute the next plane rotation Qk

         gamma  = sqrt( gbar**2 + beta**2 )   ! gammak
         cs     = gbar / gamma                ! ck
         sn     = beta / gamma                ! sk
         phi    = cs * phibar                 ! phik
         phibar = sn * phibar                 ! phibark+1

         if (debug) then
            write(*,*) ' '
            write(*,*) 'alfa ', alfa
            write(*,*) 'beta ', beta
            write(*,*) 'gamma', gamma
            write(*,*) 'delta', delta
            write(*,*) 'gbar ', gbar
            write(*,*) 'epsln', epsln
            write(*,*) 'dbar ', dbar
            write(*,*) 'phi  ', phi
            write(*,*) 'phiba', phibar
            write(*,*) ' '
         end if

         ! Update  x.

         denom = one/gamma

         do i = 1, n
            w1(i) = w2(i)
            w2(i) = w(i)
            w(i)  = ( v(i) - oldeps*w1(i) - delta*w2(i) ) * denom
            x(i)  =   x(i) +   phi * w(i)
         end do

         ! optional ortogonalization w.r.t. the kernel     
         if ( present(ortogonalization_matrix) .and. &
               ( ctrl%iort  .gt. 0            ) .and. &
               ( mod(info%iter,ctrl%iort) .eq. 0 ) ) then
             call info%ort%set('start')
             ! compute (I-QQ^T)x
             call ortogonalization_matrix%matrix_times_vector(&
                  w1,my_work,info%ierr,lun_err)
             w1 = w1 - my_work
             
             call ortogonalization_matrix%matrix_times_vector(&
                  w2,my_work,info%ierr,lun_err)
             w2 = w2 - my_work

             call ortogonalization_matrix%matrix_times_vector(&
                  w,my_work,info%ierr,lun_err)
             w = w - my_work

             call ortogonalization_matrix%matrix_times_vector(&
                  x,my_work,info%ierr,lun_err)
             x = x - my_work
             
             call info%ort%set('stop')
          end if 
         


         ! Go round again.

         gmax   = max( gmax, gamma )
         gmin   = min( gmin, gamma )
         z      = rhs1 / gamma
         ynorm2 = z**2  +  ynorm2
         rhs1   = rhs2  -  delta * z
         rhs2   =       -  epsln * z

         ! Estimate various norms and test for convergence.

         Anorm  = sqrt( tnorm2 )
         ynorm  = sqrt( ynorm2 )
         epsa   = Anorm * eps
         epsx   = Anorm * ynorm * eps
         epsr   = Anorm * ynorm * rtol
         diag   = gbar
         if (diag .eq. zero) diag = epsa

         qrnorm = phibar
         rnorm  = qrnorm

         if ( isnan(rnorm) ) then
            rc = IOerr(lun_err, wrn_val, ' minres_solver', &
                 'resnorm = NaN')
            info%ierr =-1 
            return
         end if

        

         ! Estimate  cond(A).
         ! In this version we look at the diagonals of  R  in the
         ! factorization of the lower Hessenberg matrix,  Q * H = R,
         ! where H is the tridiagonal matrix from Lanczos with one
         ! extra row, beta(k+1) e_k^T.

         Acond  = gmax / gmin

         ! See if any of the stopping criteria are satisfied.
         ! In rare cases, istop is already -1 from above (Abar = const*I).

         if (istop .eq. 0) then
            if (itn    .ge. itnlim    ) istop = 5
            if (Acond  .ge. 0.1d+0/eps) istop = 4
            if (epsx   .ge. beta1     ) istop = 3
            if (qrnorm .le. epsx      ) istop = 2
            if (qrnorm .le. epsr      ) istop = 1
            if (qrnorm *inverse_residum_weight  .le. rtol    ) istop = 1
         end if



         ! See if it is time to print something.

         if (nout .gt. 0) then
            prnt   = .false.
            if (n      .le. 40         ) prnt = .true.
            if (itn    .le. 10         ) prnt = .true.
            if (itn    .ge. itnlim - 10) prnt = .true.
            if (mod(itn,10)  .eq.     0) prnt = .true.
            if (qrnorm .le.  ten * epsx) prnt = .true.
            if (qrnorm .le.  ten * epsr) prnt = .true.
            if (Acond  .ge. 1.0d-2/eps ) prnt = .true.
            if (istop  .ne.  0         ) prnt = .true.

            if ( prnt ) then
               if (    itn     .eq. 1) write(nout, 1200)
               write(nout, 1300) itn, x(1), qrnorm, Anorm, Acond
               if (mod(itn,10) .eq. 0) write(nout, 1500)
            end if
         end if

         ! optional printing of convergence rate
         if (ctrl%iprt.gt.0) then
            write(ctrl%lun_out,*)&!'(i10,e15.7)') &
                 itn,qrnorm,rnorm,epsx,epsr

            call matrix%matrix_times_vector(x,my_work,info%ierr,lun_err)
            my_work = rhs - my_work
            write(ctrl%lun_out,*)  itn, dnrm2(nequ,my_work,1)*inverse_residum_weight,dnrm2(nequ,my_work,1)
         end if
         ! orignal code go to 100
         !
      end do
      !------------------------------------------------------------------
      ! End of main iteration loop.
      !------------------------------------------------------------------

      !
      ! convert outputs from original code into outputs fr this code
      !

      !
      ! copy results
      !
      sol = x
     
      ! errors
      if  (istop .le. 2) then 
         info%ierr = 0
      else
         info%ierr = -1
      end if
      info%iter = itn
      info%resnorm = rnorm
      if ( isnan(info%resnorm) ) then
         rc = IOerr(lun_err, wrn_val, ' minres_solver', &
              'resnorm = NaN')
         info%ierr =-1 
         return
      end if


      ! Anorm, Acond, ynorm useful but not passed
      info%matrix_cond = Acond

      !
      ! compute final residum
      !
      call info%matvec%set('start')
      call matrix%matrix_times_vector(sol,r1,info%ierr,lun_err)
      call info%matvec%set('stop')

      call info%vecvec%set('start')
      r1 = rhs - r1
      call info%vecvec%set('stop')

      call info%scalprod%set('start')
      info%resreal = dnrm2(nequ,r1,1)*inverse_residum_weight
      call info%scalprod%set('stop')


      !
      ! free memory if required
      !
      call info%ovh%set('start')

      x => null()
      y => null()
      r1 => null()
      r2 => null()
      v => null()
      w => null()
      w1 => null()
      w2 => null()
      if ( aux_work%is_initialized ) call aux_work%kill(lun_err)  
      aux_loc => null() 
      call info%ovh%set('stop')
      call info%tot%set('stop ')
      

      ! Display final status.

  900 if (nout  .gt. 0) then
         write(nout, 2000) exit, istop, itn,&
              exit, Anorm, Acond,&
              exit, rnorm, ynorm
         write(nout, 3000) exit, msg(istop)
      end if

      return


 1000 format(// 1p,    a, 5x, 'Solution of symmetric   Ax = b'&
           / ' n      =', i7, 5x, 'checkA =', l4, 12x,&
           'precon =', l4 &
           / ' itnlim =', i7, 5x, 'rtol   =', e11.2, 5x,&
           'shift  =', e23.14)
1200  format(// 5x, 'itn', 8x, 'x(1)', 10x,&
           'norm(r)', 3x, 'norm(A)', 3X, 'cond(A)')
1300  format(1p, i8, e19.10, 3e10.2)
1500  format(1x)
2000  format(/ 1p, a, 5x, 'istop =', i3,   14x, 'itn   =', i8 &
           /     a, 5x, 'Anorm =', e12.4, 5x, 'Acond =', e12.4 &
           /     a, 5x, 'rnorm =', e12.4, 5x, 'ynorm =', e12.4)
3000  format(      a, 5x, a )
      
    contains
      subroutine dscal2( n, a, x, y )
        implicit           none
        integer            n
        double precision   a, x(n), y(n)
        !-----------------------------------------------------------------
        ! dscal2 sets y = a*x.
        !------------------------------------------------------------------
        integer :: i
        
        do i = 1, n
           y(i) = a*x(i)
        end do
        
    end subroutine dscal2


  end  subroutine minres_solver

  

    subroutine init_inverse(this,lun_err,nequ)
    implicit none
    class(inverse),    intent(inout) :: this
    integer,           intent(in   ) :: lun_err
    integer, intent(in   ) :: nequ
    ! local 
    logical :: rc
    
    this%nrow         = nequ
    this%ncol         = nequ
    call this%aux%init(6,0, 10*nequ)
    call this%info_solver%init()

    


  end subroutine init_inverse

  subroutine kill_inverse(this,lun_err)
    implicit none    
    class(inverse),    intent(inout) :: this
    integer, intent(in  ) ::  lun_err

    this%matrix => null()
    this%prec_left  => null()
    this%prec_right  => null()
    this%ortogonalization_matrix =>null()
    call this%aux%kill(lun_err)

    call this%info_solver%kill()

  end subroutine kill_inverse

  subroutine info_inverse(this,lun)
    implicit none    
    class(inverse),    intent(in) :: this
    integer, intent(in  ) ::  lun

    write(lun, *) 'Total internal iterations=', this%total_iter
  end subroutine info_inverse

  subroutine set_inverse(this,matrix,ctrl_solver,prec_left,prec_right,&
       ortogonalization_matrix)
    implicit none
    class(inverse),   intent(inout) :: this
    class(abs_linop), target, intent(in  ) :: matrix
    type(input_solver),  intent(in  ) :: ctrl_solver
    class(abs_linop),   optional, target, intent(in  ) :: prec_left
    class(abs_linop),   optional, target, intent(in  ) :: prec_right
    class(abs_linop),   optional, target, intent(in  ) :: ortogonalization_matrix
    ! local
    logical :: rc
    integer :: lun_err
    
    lun_err = ctrl_solver%lun_err
    
    this%is_symmetric = matrix%is_symmetric
    if ( matrix%nrow .ne. matrix%ncol ) then
       rc = IOerr(lun_err, wrn_val, 'init_solver', &
            ' non-squared matrix passed')
       return
    end if
    

    this%matrix => matrix
    this%ctrl_solver = ctrl_solver
    if (present(prec_left)) then
       this%prec_left  => prec_left
    else
       this%prec_left  => null()
    end if

    if (present(prec_right)) then
       this%prec_right => prec_right
    else
       this%prec_right => null()
    end if
    
    if( present(ortogonalization_matrix)) then
       this%ortogonalization_matrix => ortogonalization_matrix
       this%ctrl_solver%iort = 1
    else
       this%ortogonalization_matrix =>null()
    end if
       
    this%total_iter = 0
    this%is_symmetric = this%is_symmetric

  end subroutine set_inverse


!!$  function inv(matrix,lun_err,ctrl_solver,prec_left) result(inv)
!!$    implicit none
!!$    class(abs_matrix), intent(in   ) :: matrix
!!$    integer, optional, intent(in   ) :: lun_err
!!$    type(input_solver), optional, intent(in  ) :: ctrl_solver
!!$    class(abs_linop),  intent(in   ) :: prec_left
!!$    type(inverse) :: inv
!!$    ! local 
!!$    logical :: rc
!!$    integer :: lun_err_loc
!!$    type(input_solcer) :: ctrl_solver_loc
!!$    type(input_prec) :: ctrl_prec_loc
!!$
!!$    lun_err = 0
!!$    if (present(lun_err)) lun_err_loc=lun_err
!!$    
!!$    !
!!$    ! allocate memory
!!$    !
!!$    call inv%init(lun_err,matrix)
!!$
!!$    select type(mat=>matrix)
!!$       !
!!$       ! sparse matrix via iterative methods
!!$       !
!!$    type is (spmat) 
!!$       !
!!$       ! select interative scheme
!!$       !
!!$       if (present(ctrl_solver)) then
!!$          ctrl_solver_loc = ctrl_solver
!!$       else
!!$          if ( matrix%is_symmetric ) then
!!$             ctrl_solver_loc%scheme = 'MINRES'
!!$          else
!!$             ctrl_solver_loc%scheme = 'BICGSTAB'
!!$          end if
!!$       end if
!!$          
!!$
!!$       if (present(prec_left)) then
!!$          prec_left_loc => prec_left
!!$       else
!!$          if ( mat%is_symmetric ) then
!!$             ctrl_prec_loc%scheme = 'IC'
!!$          else
!!$             ctrl_prec_loc%scheme = 'ILU'
!!$          end if
!!$          ctrl_prec_loc%n_fillin = 10
!!$          ctrl_prec_loc%tol_fillin = 1.0d-3
!!$          call inv%spprec%init(lun_err_loc,ctrl_prec_loc,mat)
!!$       end if
!!$    end select
!!$
!!$       
!!$
!!$   
!!$       
!!$
!!$    call inv%set(matrix,ctrl_solver,prec_left=prec_left_loc)
!!$
!!$
!!$    
!!$
!!$
!!$  end function inv

  subroutine apply_inverse(this,vec_in,vec_out,info,lun_err)
    implicit none
    class(inverse),    intent(inout) :: this
    real(kind=double), intent(in   ) :: vec_in(this%ncol)
    real(kind=double), intent(inout) :: vec_out(this%nrow)
    integer,           intent(inout) :: info
    integer,           intent(in   ) :: lun_err
    ! local
    logical :: rc
    integer :: res
    real(kind=double) :: dnrm2
    character(256) ::msg,msg1
       
    !if (vec_out=zero
    info=0
    call  linear_solver(&
         this%matrix, vec_in,  vec_out, &
         this%info_solver, &
         this%ctrl_solver,&
         prec_left=this%prec_left,&
         aux=this%aux,&
         ortogonalization_matrix=this%ortogonalization_matrix)
    this%total_iter =  this%total_iter + this%info_solver%iter
    if ( this%ctrl_solver%print_info_solver) &
         call this%info_solver%info(this%ctrl_solver%lun_out)


    if (this%info_solver%ierr .ne. 0) then
       info = this%info_solver%ierr
       rc = IOerr(lun_err, wrn_val, ' apply_inverse', &
            ' failure in internal linear solver procedure',info)
       call this%info_solver%info(lun_err)
       return
    end if
    

  end subroutine apply_inverse

    

  !(A BT  ) (x) = f 
  !(C -D  ) (y) = g
  subroutine uzawa_arrow_hurwics(&
       matrix_A, &
       matrix_BT, &
       matrix_C, &
       matrix_D,&
       vec_f, vec_g,&
       vec_x, vec_y,&
       alpha_relax, w_relax,&
       ctrl_solver, info_solver,&
       prec_A, prec_S,aux)
    use SimpleMatrix
    implicit none
    class(abs_linop),  intent(inout) :: matrix_A
    class(abs_linop),  intent(inout) :: matrix_BT
    class(abs_linop),  intent(inout) :: matrix_C
    class(abs_linop),  intent(inout) :: matrix_D
    real(kind=double), intent(in   ) :: vec_f(matrix_A%nrow)
    real(kind=double), intent(in   ) :: vec_g(matrix_D%nrow)
    real(kind=double), intent(inout) :: vec_x(matrix_A%nrow)
    real(kind=double), intent(inout) :: vec_y(matrix_D%nrow)
    type(input_solver),  intent(in  ) :: ctrl_solver
    type(output_solver),  intent(inout) :: info_solver
    real(kind=double), intent(in) :: alpha_relax
    real(kind=double), intent(in) :: w_relax
    type(scrt), optional,target, intent(inout) :: aux
    class(abs_linop), target, optional, intent(inout) :: prec_A
    class(abs_linop), target, optional, intent(inout) :: prec_S
    !local
    logical :: rc
    logical :: steady
    integer :: iter,ibegin
    integer :: na,ns,nequ,iend
    type(scrt), target  :: aux_work
    type(scrt), pointer :: aux_loc
    type(scalmat),    target :: identity_na
    type(scalmat),    target :: identity_ns
    class(abs_linop),  pointer :: prec_A_loc
    class(abs_linop),  pointer :: prec_S_loc
    real(kind=double), allocatable :: scr0(:), scr1(:), scr2(:) 
    real(kind=double) :: dnrm2, incx, incy,r1,r2

    na=matrix_A%nrow
    ns=matrix_C%nrow

    write(*,*) na,ns

    !
    ! handle prec_A
    !
    if ( present( prec_A ) ) then
       write(*,*) 'prec_A present' 
       prec_A_loc => prec_A
    else
       call identity_na%eye(na)
       prec_A_loc => identity_na
    end if

    !
    ! handle prec_S
    !
    if ( present( prec_S ) ) then
       prec_S_loc => prec_S
    else
       write(*,*) 'prec_S not present'
       call identity_ns%eye(ns)
       prec_S_loc => identity_ns
    end if

    !
    ! handle scratch variable
    !
    nequ=max(na,ns)
    if ( present(aux) ) then
       if( .not. aux%check(0,3*nequ) ) &
            rc = IOerr(ctrl_solver%lun_err, wrn_inp, 'minres_solver', &
            ' aux array too small')
       aux_loc => aux
    else
       call aux_work%init(ctrl_solver%lun_err,0,3*nequ)
       aux_loc => aux_work  
    end if
    
    
!!$    iend=0
!!$    call aux_loc%range(nequ,ibegin,iend)
!!$    scr0 =>aux_loc%raux(ibegin:iend)
!!$    call aux_loc%range(nequ,ibegin,iend)
!!$    scr1 =>aux_loc%raux(ibegin:iend)
!!$    call aux_loc%range(nequ,ibegin,iend)
!!$    scr2 =>aux_loc%raux(ibegin:iend)
    allocate(scr0(nequ), scr1(nequ),scr2(nequ))
    



    iter = 0
    steady = .False.
    info_solver%ierr=0
    do while ( ( iter < ctrl_solver%imax ) .and. ( .not. steady ))   
       iter = iter + 1
       ! compute f-Ax_{k}-B^T y_{k}
       call matrix_A%matrix_times_vector(vec_x, scr1(1:na),info_solver%ierr,ctrl_solver%lun_err)
       call matrix_BT%matrix_times_vector(vec_y,scr2(1:na),info_solver%ierr,ctrl_solver%lun_err)
       scr0(1:na) = vec_f - scr1(1:na) - scr2(1:na)
       if ( iter .eq. 1 ) r1=dnrm2(na,scr0,1)
       
       ! compute x_{k+1} = x_{k+1} + alpha * (~A)^{-1} ( scr )
       call prec_A_loc%matrix_times_vector(scr0(1:na),scr1(1:na),info_solver%ierr,ctrl_solver%lun_err)
       vec_x(1:na) = vec_x(1:na) + alpha_relax * scr1(1:na)
       incx=dnrm2(na,scr1,1)

       ! compute  C x_{k+1} - D y_{k} - g
       call matrix_C%matrix_times_vector(vec_x, scr1(1:ns),info_solver%ierr,ctrl_solver%lun_err)
       call matrix_D%matrix_times_vector(vec_y, scr2(1:ns),info_solver%ierr,ctrl_solver%lun_err)
       scr0(1:ns)= scr1(1:ns)-scr2(1:ns)-vec_g(1:ns)
       if ( iter .eq. 1 ) then 
          r2=dnrm2(na,scr0,1)
          info_solver%resini = r1 + r2
       end if
          

       ! compute x_{k+1} = x_{k+1} + alpha * (~S)^{-1} ( scr )
       call prec_S_loc%matrix_times_vector(scr0,scr1,info_solver%ierr,ctrl_solver%lun_err)
       vec_y = vec_y + w_relax * scr1

       incy=dnrm2(ns,scr1,1)

       steady = ( (incx + incy) < ctrl_solver%tol_sol )
!!$       if ( info_solver%resini < 1.0d8 *( incx+incy) ) then
!!$          exit
!!$       end if
       
       if ( ctrl_solver%iprt .eq. 1 ) write(ctrl_solver%lun_out,*) iter, incx, incy
    end do

    if (.not. steady ) info_solver%ierr=1
    info_solver%iter=iter
    info_solver%resnorm = incx+incy
    info_solver%resreal= incx+incy
    deallocate(scr0, scr1,scr2)
    
    
  end subroutine uzawa_arrow_hurwics



  recursive subroutine fgmres_solver(&
       matrix, rhs, sol,&
       info_solver, ctrl, &
       prec_left,&
       prec_right,&
       aux_fgmres)
    use Globals
    use LinearOperator
    use Scratch
    class(abs_linop),           intent(inout) :: matrix
    real(kind=double),          intent(in   ) :: rhs(matrix%nrow)
    real(kind=double),          intent(inout) :: sol(matrix%ncol)
    class(output_solver),       intent(inout) :: info_solver
    class(input_solver),        intent(in   ) :: ctrl
    class(abs_linop),           intent(inout) :: prec_left
    class(abs_linop),           intent(inout) :: prec_right
    type(scrt), optional,target,intent(inout) :: aux_fgmres

    ! varaible used in reverse comunication mi15id
    integer ::         ICNTL( 8 )
    real(kind=double) ::  CNTL( 4 )
    integer :: ISAVE(17)
    real(kind=double) :: RSAVE(9)
    logical :: LSAVE(4)
    real(kind=double) :: RESID,inverse_residum_weight
    integer :: IACT, N, M, LDW, LDZ, LOCY, LOCZ, LDH
    real(kind=double), pointer :: W( :, :), H( :, : ), Z(:, :)
    integer :: INFO( 4 )
    integer :: info_loc
    integer :: K
    real(kind=double) :: dnrm2
    
    ! interface 
    type(scrt), target  :: aux_work
    type(scrt), pointer :: aux_loc
    integer :: wbegin,wend
    integer :: hbegin,hend
    integer :: zbegin,zend
    integer :: double_memory_requirement
    logical :: rc,test
    integer :: res

    call info_solver%tot%set('start')
    call info_solver%ovh%set('start')

    !
    ! set variable dimension
    !
    N   = matrix%nrow
    M   = ctrl%nrestart
    LDW = N
    LDZ = N
    LDH = M + 1

    !
    ! Initialization of controls
    !
    CALL MI15ID( ICNTL, CNTL, ISAVE, RSAVE, LSAVE )
    ICNTL(1) = ctrl%lun_err ! error unit
    ICNTL(2) = ctrl%lun_err ! warning unit
    ICNTL(3) = 3            ! Precondition from both sides
    ICNTL(4) = 0            ! use relative residual
    if ( ctrl%iexit ==0) then
       write(*,*) 'FGMRES: absolute residua non supported'
       stop
    end if
    ICNTL(5) = 1  ! use sol passed as initial solution 
    ICNTL(6) = ctrl%imax ! maximum number of total iteration

    CNTL(1) = ctrl%tol_sol
    CNTL(2) = zero
    
    !
    ! set work arrays
    !    
    double_memory_requirement = LDW*(M+7)+LDH*(M+2)+LDZ*M
    wbegin=1
    wend  =LDW*(M+7)
    hbegin=wend+1
    hend  =wend+LDH*(M+2)
    zbegin=hend+1
    zend  =hend+LDZ*M
    if ( .not. present(aux_fgmres) )then
       call aux_work%init(ctrl%lun_err,0,double_memory_requirement)
       aux_loc => aux_work
    else
       if ( .not. aux_fgmres%is_initialized) then
          call aux_work%init(ctrl%lun_err,0,double_memory_requirement)
          aux_loc => aux_work
       else
          aux_loc => aux_fgmres
          test=aux_fgmres%check(0,double_memory_requirement)
          if( .not. aux_loc%check(0,double_memory_requirement) ) then
             write(ctrl%lun_err,*) 'required=', &
                  double_memory_requirement, 'passed=',aux_loc%nraux
             rc = IOerr(ctrl%lun_err, err_inp, 'fgmres_solver', &
                  ' aux array too small',res)
          end if
       end if
    end if

    W(1:LDW,1:M+7) => aux_loc%raux(wbegin:wend) 
    H(1:LDH,1:M+2) => aux_loc%raux(hbegin:hend)
    Z(1:LDZ,1:M  ) => aux_loc%raux(zbegin:zend)

    !
    ! set inverse_residum_weight
    ! inverse_residum_weight = one       => resnorm = |Ax-b|
    ! inverse_residum_weight = one/bnrom => resnorm = |Ax-b|/|b|
    !
    if (ctrl%iexit.eq.0) then
       inverse_residum_weight = one
    else
       inverse_residum_weight = one / info_solver%bnorm
    end if

    if (ctrl%compute_initial_residue == 1) then
       !
       ! calculate initial residual (res = b-M*x)
       !      
       call info_solver%matvec%set('start')
       call matrix%matrix_times_vector(sol,W(1:N,1),info_solver%ierr,ctrl%lun_err)
       call info_solver%matvec%set('stop')

       call info_solver%vecvec%set('start')
       W(1:N,1) = rhs - W(1:N,1)
       call info_solver%vecvec%set('stop')

       !
       ! compute initial norm of residum
       !
       call info_solver%scalprod%set('start')
       info_solver%resini = dnrm2(N,W(1:N,1),1)*inverse_residum_weight
       call info_solver%scalprod%set('stop')

       !
       ! exit if the initial residuum is already below the required tolerance 
       !
       if ( info_solver%resini .lt. ctrl%tol_sol) then
          !
          info_solver%resreal = info_solver%resini
          info_solver%resnorm = zero
          !
          ! free memory
          !
          call info_solver%ovh%set('start')
          if ( aux_work%is_initialized ) call aux_work%kill(ctrl%lun_err)  
          aux_loc => null() 
          call info_solver%ovh%set('stop')
          call info_solver%tot%set('stop ')
          return
       end if
    end if
    
    ! Set right hand side, b
    w(1:N,1)=rhs

    ! Set initial data with current sol
    w(1:N,2)=sol
    call info_solver%ovh%set('stop')
    
    
    !
    !  Solve the system using the Flex-GMRES(m) method
    !
    IACT = 0
    !write(*,*) ctrl%tol_sol,'IACT',IACT,'info 1',INFO(1),'info(2)',info(2),'cntl(:)', cntl(:),'ICNTL(4)',ICNTL(4)
    DO 80 K = 1, 2*N
       CALL MI15AD( IACT, N, M, W, LDW, Z, LDZ,  LOCY, LOCZ, H, LDH,&
                      RESID, ICNTL, CNTL, INFO, ISAVE, RSAVE, LSAVE )
       !write(*,*) 'IACT',IACT,'info 1',INFO(1),'info(2)',info(2),'cntl(1)', cntl(1),'ICNTL(4)',ICNTL(4)
       IF ( IACT .LT. 0 ) THEN
          !  ERROR
          WRITE( 6, '(A,I2,A)' ) ' Error return: INFO( 1  ) = ',&
               INFO( 1 ),' on exit '
          if (INFO(1) < 0  ) then
             info_solver%ierr = INFO(1) 
          end if
          info_solver%iter= INFO(2)
          EXIT
          !
       ELSEIF ( IACT .EQ. 1 ) THEN
          !  Solution found
          sol = W(1:N,2)
          info_solver%iter = INFO( 2 )
          if (ctrl%compute_real_residue .eq. 1) then
             call matrix%matrix_times_vector(sol,W(1:N,1),info_solver%ierr,ctrl%lun_err)
             W(1:N,1)=W(1:N,1)-rhs
             info_solver%resreal = dnrm2(N,W(1:N,1),1)*inverse_residum_weight
          end if

          EXIT

       ELSEIF ( IACT .EQ. 2 ) THEN
          !
          !  Perform the matrix-vector product
          !
          call info_solver%matvec%set('start')
          call matrix%matrix_times_vector(W(1:N,LOCZ),W(1:N,LOCY),info_solver%ierr,ctrl%lun_err)
          call info_solver%matvec%set('stop')          

          !
          if (ctrl%iprt .ge. 2) then
             write(ctrl%lun_out,*) INFO(2), RSAVE(5)
          end if

          
       ELSEIF ( IACT .EQ. 3 ) THEN
          !
          !  Perform the left preconditioning operation
          !
          call info_solver%prec%set('start')
          call prec_left%matrix_times_vector(W(1:N,LOCZ),W(1:N,LOCY),info_solver%ierr,ctrl%lun_err)
          call info_solver%prec%set('stop')

          
       ELSEIF ( IACT .EQ. 4 ) THEN
          !
          !  Perform the right preconditioning operation
          !
          call info_solver%prec%set('start')
          call prec_right%matrix_times_vector(W(1:N,LOCZ),W(1:N,LOCY),info_solver%ierr,ctrl%lun_err)
          call info_solver%prec%set('stop')
       END IF
80  END DO

    if (K==2*N ) then
       info_solver%ierr = 1
       info_solver%iter = K 
    end if

    

    !
    ! free memory
    !
    call info_solver%ovh%set('start')
    W => null()
    H => null()
    Z => null()
    if ( aux_work%is_initialized ) call aux_work%kill(ctrl%lun_err)  
    aux_loc => null() 
    call info_solver%ovh%set('stop')
    call info_solver%tot%set('stop')

  end subroutine fgmres_solver
    
    

end module LinearSolver

