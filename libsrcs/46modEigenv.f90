!> Module for work with eigvalue
module Eigenv
  use Globals
  use LinearOperator
  use Timing
  use DenseMatrix
  implicit none
  public ::  arpack_wrap,study_eigen
  private
  !> Structure variable for storing
  !> eigenvalue and eigenvectors
  type, public :: eigen
     !> Logical flag for initilaization
     logical :: is_initialized=.false.
     !> Array legnth of eigvectors
     integer :: length_eigen=0 
     !> number of eigenvalues
     integer :: nev=0
     !> number of first non-zero eigenvalue
     integer :: first_nonzero=0
     !> 0/1 flag if eigen. have been compute
     integer :: computed=0
     !> Dimension(nev)
     !> Real array of eigenvalues
     real(kind=double), allocatable :: eigval(:)
     !> Full matrix with dimension ( len_eigen,nev)
     !> containg eigvectors
     type(densemat) :: eigenvectors
     !> Dimension(nev)
     !> Real array of containg the residuum
     !> |A v_i-\lambda_i v_i|/|v_i| 
     real(kind=double), allocatable :: residua(:)
   contains
     !> static constructor
     !> (public for type eig)
     procedure, public, pass :: init => init_eig
     !> static destructor
     !> (public for type eig)
     procedure, public, pass :: kill => kill_eig
     !> info procedure
     !> (public for type eig)
     procedure, public, pass :: info => info_eig
     !> DACG procedure computating eigenvectors 
     !> and eigenvalues for sparse matrix 
     !> (public for type eig)
     procedure, public, pass :: dacg_algorithm
     !> LANCZOS procedure computating eigenvectors 
     !> and eigenvalues for sparse matrix prec^{-1}  matrix
     !> given coefficients alpaha beta pres_dot_res 
     !> and vectors pres
     !> (public for type eig) 
     procedure, public, pass :: lanczos_eigprec
  end type eigen
  !> Structure varible containg the input/controls
  !> variables for DACG procedure
  type, public :: input_dacg
     !> I/O err msg. unit    (Defualt = 6)
     integer :: lun_err=6
     !> I/O output msg. unit (Defualt = 6)
     integer :: lun_out=6
     !> Flag for residuum computation
     integer :: iexit=0
     !> Number of maximal iterations (Defualt = 200)
     integer :: imax=200 
     !> Number of maximal iterations (Defualt = 0)
     integer :: iprt=0 
     !> Flag for starting solution (Defualt = 0)
     !>  isol=0 => SOL      is the initial solution for dacg
     !>  isol=1 => PREC*RHS is the initial solution for dacg
     integer :: isol=0
     !> Tolerance for DACG scheme (Defualt = 1.0d-6)
     real(kind=double) :: tol_dacg=1.0d-6
   contains
     !> static constructor
     !> (procedure public for type input_dacg)
     procedure, public, pass :: init => init_input_dacg
     !> static destructor
     !> (procedure public for type ctrl_precond)
     procedure, public, pass :: kill => kill_input_dacg
     !> Info procedure.
     !> (public procedure for type ctrl_precond)
     procedure, public, pass :: info => info_input_dacg
  end type input_dacg
  !> Structure varible containg the output/info
  !> variables for DACG procedure
  type, public :: output_dacg
     !> Integer flag for error
     !>   ierr=0 => Convergence achievied                
     !>   ierr=1 => iterations exceeded max fixed 
     integer :: ierr=0
     !> Number of dacg-iterations
     integer :: totit=0
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
     !> Number of eigvalue calculted
     integer :: nev
     !> Dimension(nev)
     !> Number of iterations for each eigen value
     !> required dor DACG algorithm
     integer, allocatable :: iterations(:)
   contains
     !> static constructor
     !> (procedure public for type input_dacg)
     procedure, public, pass :: init => init_output_dacg
     !> static destructor
     !> (procedure public for type ctrl_precond)
     procedure, public, pass :: kill => kill_output_dacg
     !> Info procedure.
     !> (public procedure for type ctrl_precond)
     procedure, public, pass :: reset => reset_output_dacg
     !> Info procedure.
     !> (public procedure for type ctrl_precond)
     procedure, public, pass :: info => info_output_dacg
     !> info procedure
     !> (public for type eig)
     procedure, public, pass :: CPUtimes => timing_output_dacg
  end type output_dacg
contains
  !>-------------------------------------------------------------
  !> Static constructor.
  !> (procedure public for type eig)
  !> Instantiate variable of type eig
  !>
  !> usage:
  !>     call 'var'%init(lun_err,nev,length_eigen)
  !>
  !> where:
  !> \param[in] lun_err -> integer. Error mess. logic unit
  !> \param[in] nev     -> integer. Number of eigv.
  !> \param[in] length_eigen    -> integer. Eigenvector dimension
  !<-------------------------------------------------------------
  subroutine init_eig(this, lun_err, nev, length_eigen)
    use Globals
    implicit none
    !vars
    class(eigen),   intent(out) :: this
    integer,      intent(in ) :: lun_err
    integer,      intent(in ) :: nev
    integer,      intent(in ) :: length_eigen
    !local 
    integer :: res
    logical :: rc

    this%is_initialized = .true.
    this%length_eigen   = length_eigen
    this%nev            = nev

    allocate (&
         this%eigval(nev),&
         this%residua(nev),&
         stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'init_eig', &
         'type eig member eigval, residua, eigvec',res)

    call this%eigenvectors%init(lun_err,&
         length_eigen, nev,&
         is_symmetric=.False., triangular='N')  
    this%eigenvectors%coeff = zero
    

    this%eigval  = zero
    this%residua = large
    this%residua = huge
    this%computed= 0
  end subroutine init_eig

  !>-------------------------------------------------------------
  !> Static destructor.
  !> (procedure public for type eig)
  !> Kill the variable of type eig
  !>
  !> usage:
  !>     call 'var'%init(lun_err)
  !>
  !> where:
  !> \param[in] lun_err -> integer. Error mess. logic unit
  !<-------------------------------------------------------------
  subroutine kill_eig(this, lun_err )
    use Globals
    implicit none
    !vars
    class(eigen), intent(inout) :: this
    integer,    intent(in   ) :: lun_err
    !local 
    integer :: res
    logical ::rc

    deallocate (&
         this%eigval,&
         this%residua,&
         stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'kill_eig', &
         'member eigval residua ',res)
    call this%eigenvectors%kill(lun_err)
    this%computed=0
    this%is_initialized = .false.
    
  end subroutine kill_eig

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
  subroutine info_eig(this, lun)
    use Globals
    implicit none
    !vars
    class(eigen), intent(in) :: this
    integer,    intent(in) :: lun
    !local 
    integer :: i

    write(lun,*) this%nev ,         '  ! Number of eigvalues'
    write(lun,*) this%length_eigen, '  ! Number of equations'
    if ( this%computed .eq. 1) then
       write(lun,*) ' index eig   Eigenvalue residua'
       do i=1,this%nev
          write(lun,*) i, this%eigval(i), this%residua(i)
       end do
    else
       write(lun,*) 'Eigen* not computed'
    end if
   
  end subroutine info_eig

  !>-------------------------------------------------------------
  !> Static constructor.
  !> (procedure public for type input_dacg)
  !> Instantiate variable of type input_dacg
  !>
  !> usage:
  !>     call 'var'%init(iexit,imax,iprt,isol,tol_dacg)
  !>
  !> where:
  !> \param[in] lun_err  -> integer. I/O unit for error  msg 
  !> \param[in] lun_out  -> integer. I/O unit for output msg 
  !> \param[in] imax     -> integer. Max number of iterations 
  !> \param[in] iprt     -> integer. Flag for printing info
  !>                         iprt = 0 (no print)
  !>                         iprt = 1 (print only end)
  !>                         iprt = 2 (print at each iteration)
  !> \param[in] isol     -> integer. Flag for initial solution
  !>                         isol = 0 (start from sol)
  !>                         isol = 1 (start from P^-1 rhs)
  !> \param[in] tol_dacg -> real. Tolerance to ve achieved
  !<-------------------------------------------------------------
  subroutine init_input_dacg(this,&
       lun_err,lun_out,iexit,imax,iprt,isol,tol_dacg)
    use Globals
    implicit none
    class(input_dacg),           intent(inout) :: this 
    integer,           optional, intent(in   ) :: lun_err
    integer,           optional, intent(in   ) :: lun_out
    integer,           optional, intent(in   ) :: iexit
    integer,           optional, intent(in   ) :: imax
    integer,           optional, intent(in   ) :: iprt
    integer,           optional, intent(in   ) :: isol
    real(kind=double), optional, intent(in   ) :: tol_dacg

    if (present(lun_err))  this%lun_err  = lun_err
    if (present(lun_out))  this%lun_out  = lun_out
    if (present(iexit))    this%iexit    = iexit 
    if (present(imax))     this%imax     = imax 
    if (present(iprt))     this%iprt     = iprt
    if (present(isol))     this%isol     = isol 
    if (present(tol_dacg)) this%tol_dacg = tol_dacg
  end subroutine init_input_dacg

  !>-------------------------------------------------------------
  !> Static destructor.
  !> (procedure public for type star)
  !> Reset to defualt variable of type input_dacg
  !>
  !> usage:
  !>     call 'var'%kill()
  !>
  !<-------------------------------------------------------------
  subroutine kill_input_dacg(this)
    use Globals
    implicit none
    class(input_dacg),   intent(inout) :: this 

    this%lun_err  = 6
    this%lun_out  = 6
    this%iexit    = 0
    this%imax     = 200
    this%iprt     = 0
    this%isol     = 0
    this%tol_dacg = 1.0d-7

  end subroutine kill_input_dacg

  !>-------------------------------------------------------------
  !> Static destructor.
  !> (procedure public for type star)
  !> Reset to defualt variable of type input_dacg
  !>
  !> usage:
  !>     call 'var'%kill()
  !>
  !<-------------------------------------------------------------
  subroutine info_input_dacg(this,lun)
    use Globals
    implicit none
    class(input_dacg), intent(in) :: this 
    integer,           intent(in) :: lun
    

    write(lun,'(a,I2,a,I2,a,I1,a,I3,a,I1,a,1pe10.2,a)') &
         'DACG CTRL (err/out LU = ',this%lun_err,',',this%lun_out,&
         '/iexit=',this%iexit,&
         ' /imax=', this%imax,&
         ' /isol=', this%isol,&
         ' /tol=', this%tol_dacg,')'

!!$    write(lun,*) 'Error Unit     = ', this%lun_err
!!$    write(lun,*) 'Error Unit     = ', this%lun_out
!!$    write(lun,*) 'iexit          = ', this%iexit  
!!$    write(lun,*) 'Max iterations = ', this%imax  
!!$    write(lun,*) 'iprt           = ', this%iprt
!!$    write(lun,*) 'isol           = ', this%isol
!!$    write(lun,*) 'Tol. required  = ', this%tol_dacg
  end subroutine info_input_dacg


  !>-------------------------------------------------------------
  !> Static constructor.
  !> (procedure public for type output_dacg)
  !> Instantiate variable of type output_dacg
  !>
  !> usage:
  !>     call 'var'%init()
  !>
  !<-------------------------------------------------------------
  subroutine init_output_dacg(this,lun_err,nev)
    use Globals
    implicit none
    class(output_dacg), intent(inout) :: this
    integer,            intent(in   ) :: lun_err
    integer,            intent(in   ) :: nev 
    !local 
    logical :: rc
    integer :: res

    this%nev = nev
    allocate(&
         this%iterations(nev),&
         stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'init_output_dacg', &
         ' member iterations',res)

    this%iterations = 1000000

    call this%tot%init()
    call this%ovh%init()
    call this%matvec%init()
    call this%vecvec%init()      
    call this%scalprod%init()
    call this%prec%init()
    call this%ort%init()

  end subroutine init_output_dacg

  !>-------------------------------------------------------------
  !> Static constructor.
  !> (procedure public for type output_dacg)
  !> Instantiate variable of type output_dacg
  !>
  !> usage:
  !>     call 'var'%init()
  !>
  !<-------------------------------------------------------------
  subroutine reset_output_dacg(this)
    use Globals
    implicit none
    class(output_dacg), intent(inout) :: this
   
    call this%tot%init()
    call this%ovh%init()
    call this%matvec%init()
    call this%vecvec%init()      
    call this%scalprod%init()
    call this%prec%init()
    call this%ort%init()

  end subroutine reset_output_dacg

  
  !>-------------------------------------------------------------
  !> Static constructor.
  !> (procedure public for type info_dacg)
  !> Reset to default values
  !>
  !> usage:
  !>     call 'var'%kill()
  !>
  !<-------------------------------------------------------------
  subroutine kill_output_dacg(this,lun_err)
    use Globals
    implicit none
    class(output_dacg),   intent(inout) :: this 
    integer,         intent(in   ) :: lun_err
    !local
    logical :: rc
    integer :: res

    this%ierr      = 0 
    this%totit     = 0
    
    call this%tot%kill()
    call this%ovh%kill()
    call this%matvec%kill()
    call this%vecvec%kill()      
    call this%scalprod%kill()
    call this%prec%kill()
    call this%ort%kill()

    deallocate(&
         this%iterations,&
         stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_dealloc, 'kill_output_dacg', &
         ' member iterations',res)
  end subroutine kill_output_dacg


  !>-------------------------------------------------------------
  !> Static constructor.
  !> (procedure public for type info_dacg)
  !> Instantiate variable of type info_dacg
  !>
  !> usage:
  !>     call 'var'%init()
  !>
  !<-------------------------------------------------------------
  subroutine info_output_dacg(this,lun)
    use Globals
    implicit none
    class(output_dacg), intent(in) :: this 
    integer,         intent(in) :: lun
    ! local
    character(len=256) :: msg

    write(lun,'(a,I2,1x,a,I4,1x,a,f6.2,a)') &
         ' DACG: IERR: ', this%ierr, &
         ' ITER: ', this%totit, &
         ' CPU: ',this%tot%CUMwct,'s:'
    !if ( this%tot%CUMwct> zero) then
    !   write(lun,*) etb(this%CPUtimes())
    !end if
  end subroutine info_output_dacg

  !>-------------------------------------------------------------
  !> Static constructor.
  !> (procedure public for type output_dacg)
  !> Instantiate variable of type output_dacg
  !>
  !> usage:
  !>     str='var'%timing()
  !>\param[out] str :: character(len=256) Time statistic 
  !<-------------------------------------------------------------
  function timing_output_dacg(this) result(str)
    use Globals
    implicit none
    class(output_dacg), intent(in) :: this  
    character(len=256)  :: str
    !
    character(len=30) :: loc
    
    write(loc,'(f6.2,a)') this%tot%CUMwct,' s :'
    str = etb(loc)

    if ( this%tot%CUMwct > zero ) then
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
    end if


  end function timing_output_dacg


  !>-------------------------------------------------------------
  !> Procedure computing the eigvetoors and eigvalues of the
  !> for matrix M via DACG algorithm
  !> (procedure public for type eig)
  !>
  !> usage:
  !>     call 'var'%init(lun_err,lun_out,&
  !>                     matrix, prec    )
  !>
  !> where:
  !> \param[in] id_prt -> integer. Flag dor printing eigval.
  !> \param[in] lun    -> integer. Logic unit
  !<-------------------------------------------------------------
  subroutine dacg_algorithm( this, &
       ctrl,info, &
       matrix, &
       prec,&
       aux_dacg)
    use Globals
    use Matrix
    use Scratch
    implicit none

    class(eigen),       intent(inout) :: this
    type(input_dacg),   intent(in   ) :: ctrl
    type(output_dacg),  intent(inout) :: info
    class(abs_linop),   intent(inout) :: matrix
    class(abs_linop),    intent(inout) :: prec
    type(scrt), target, optional,intent(inout) :: aux_dacg
    !optional
    !local
    integer :: length_eigen,nev
    integer :: iter,indau, j 
    integer :: imax, isol, iprt
    integer :: totit,ierr
!!$    real(kind=double), allocatable :: g_old(:),xk(:),pk(:),ax(:),ap(:)
!!$    real(kind=double), allocatable :: g_new(:),scr(:),gprec(:)
!!$    real(kind=double), allocatable :: aux1(:)
    real(kind=double), pointer :: g_old(:),xk(:),pk(:),ax(:),ap(:)
    real(kind=double), pointer :: g_new(:),scr(:),gprec(:)
    real(kind=double), pointer :: aux1(:)

    type(scrt), target  :: aux_work
    type(scrt), pointer :: aux_loc
    


    real(kind=double)  pap,pbp,pax,pbx,gnormnew,cc,dd
    real(kind=double)  alpha,eig_old,eig_new,aa,bb,delta,gnormold,residum,err
    real(kind=double)  tolres,xax,xbx,beta,normax
    real(kind=double)  ddot,dnrm2,a1
    real*4   tprec,t,tdacg,totort,time,rtime
    logical :: test_exit,debug=.False.
    logical :: rc
    integer :: res,lun_err,lun_out

    lun_err = ctrl%lun_err
    lun_out = ctrl%lun_out
    imax    = ctrl%imax
    isol    = ctrl%isol
    iprt    = ctrl%iprt

    tolres = ctrl%tol_dacg
    length_eigen    = this%length_eigen
    nev     = this%nev


    ! set xero counter variables
    rtime  = 0.0e0
    tdacg  = 0.0e0
    tprec  = 0.0e0
    totort = 0.0e0
    totit  = 0
    ierr   = 0
    
    !
    ! set timing to zero
    !
    call  info%reset()
    
    !
    ! set work arrays
    !    
    if ( present(aux_dacg) ) then
       aux_loc => aux_dacg
    else
       call aux_work%init(lun_err,0,8*length_eigen+nev)
       aux_loc => aux_work  
    end if

    !> Redirection of working array
    g_old => aux_loc%raux( 1           :   length_eigen )
    xk    => aux_loc%raux( 1 +   length_eigen  : 2*length_eigen )
    pk    => aux_loc%raux( 1 + 2*length_eigen  : 3*length_eigen )
    ax    => aux_loc%raux( 1 + 3*length_eigen  : 4*length_eigen )
    g_new => aux_loc%raux( 1 + 4*length_eigen  : 5*length_eigen )
    ap    => aux_loc%raux( 1 + 5*length_eigen  : 6*length_eigen )
    scr   => aux_loc%raux( 1 + 6*length_eigen  : 7*length_eigen )
    gprec => aux_loc%raux( 1 + 7*length_eigen  : 8*length_eigen )
    aux1  => aux_loc%raux( 1 + 8*length_eigen  : 8*length_eigen +nev )

!!$    allocate(&
!!$         g_old(length_eigen),&
!!$         xk(length_eigen),&
!!$         pk(length_eigen),&
!!$         ax(length_eigen),&
!!$         g_new(length_eigen),&
!!$         ap(length_eigen),&
!!$         scr(length_eigen),&
!!$         gprec(length_eigen),&
!!$         aux1(nev),&
!!$         stat=res)
!!$    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'dacg', &
!!$         'work array g_old, xk, pk, ax, g_new, ap, scr, gprec,',res)

    ! DACG iterative solver
    if (iprt .gt. 0) then
       write(lun_out,*) ' j iter     tot      eigenvalue    ',&
            '  scarto          ||res||     CPU prog'
       write(lun_out,*) "==========================================",&
            "==============================="
    end if
    

    !  CICLO PER IL CALCOLO DI OGNI AUTOVALORE E AUTOVETTORE
    call info%tot%set('start')
    do indau = 1,nev
       !call tim(time,1)
       if (isol .eq. 0)  then
          call dcopy(length_eigen,this%eigenvectors%coeff(1,indau),1,xk,1)
       else 
          call initvec(length_eigen,indau,nev,xk)
       end if
       !if (debug) write(*,*) xk
       !
       ! B-orthogonalization of the initial 
       ! vector against the eigenvectors
       ! previously calculated
       !
       call info%ort%set('start')
       call reort(length_eigen,length_eigen,&
            nev,indau-1,this%eigenvectors%coeff,xk,totort,aux1)
       call info%ort%set('stop')

       !
       !  computation of Ax_0, Bx_0, x_0^TAx_0, x_0^TBx_0
       !
       call info%matvec%set('start')
       call matrix%matrix_times_vector(xk,ax,info%ierr,lun_err)
       call info%matvec%set('stop')
       if (info%ierr .ne. 0 ) then
          return
       end if

       call info%scalprod%set('start')
       xax = ddot(length_eigen,xk,1,ax,1)
       xbx = ddot(length_eigen,xk,1,xk,1)
       eig_new = xax/xbx
       if (debug) write(*,*) xax, xbx
       call info%scalprod%set('stop')


       iter = 0
       test_exit = .false.
       do while ( .not. test_exit )
          iter= iter + 1 
          !
          ! * g_new <-- ax - eig_new bx
          !
          !call dzaxpy(length_eigen,-eig_new,xk,1,ax,1,g_new,1)
          call info%vecvec%set('start')
          g_new = ax - eig_new * xk  
          call info%vecvec%set('stop')

          call info%vecvec%set('start')
          normax=dnrm2(length_eigen,ax,1)
          !        res = dnrm2(length_eigen,g_new,1)/dnrm2(length_eigen,ax,1)
          !  Relativo
          ! assoluto
          residum = dnrm2(length_eigen,g_new,1)
          !        write(100+indau,*) res                
          !        if(res.lt.tolres) goto 1000
          call info%vecvec%set('stop')

          if (debug) write(lun_out,*) 'res=', residum, 'normax',normax,'tol',tolres
          if (debug) write(lun_out,*) 'eig=', eig_new
          !if (debug) write(lun_out,*) 'sol=', xk

          ! exit when null eigevalue and a kernel vector is found
          if ( ( abs(eig_new) < small) .and. ( residum < tolres) ) then
             ! exit if | Ax - lambda *x |/| Ax | < tolres
             test_exit = .true.
          elseif( residum .lt. normax*tolres ) then
             test_exit = .true.
          else
             call info%vecvec%set('start')
             call dscal(length_eigen,2.d0/xbx,g_new,1)
             call info%vecvec%set('stop')

             !
             ! gprec  <--  K^(-1) * g_k
             !
             call info%prec%set('start')
             call prec%matrix_times_vector(g_new,gprec,info%ierr,lun_err)
             call info%prec%set('stop')
             if (info%ierr .ne. 0 ) then
                return
             end if

             call info%vecvec%set('start')
             call betak(iter,length_eigen,&
                  gnormnew,gnormold,beta,gprec,g_new,g_old)

             gnormold = gnormnew
             !
             ! * evaluation of the new search direction
             ! * p_k+1 = K^(-1) g_k + beta_k p_k
             !
             !call dxpay(length_eigen,gprec,1,beta,pk,1)

             pk = beta * pk + gprec
             call info%vecvec%set('stop')


             tdacg = tdacg + t
             call info%ort%set('start')
             call reort(length_eigen,length_eigen,nev,indau-1,&
                  this%eigenvectors%coeff,pk,totort,aux1)
             call info%ort%set('stop')

             !
             ! compute ap = A  pk
             !
             call info%matvec%set('start')
             call matrix%matrix_times_vector(pk,ap,info%ierr,lun_err)
             call info%matvec%set('stop')
             if (info%ierr .ne. 0 ) then
                return
             end if

             call info%scalprod%set('start')
             pax = ddot(length_eigen,pk,1,ax,1)
             pap = ddot(length_eigen,pk,1,ap,1)
             pbx = ddot(length_eigen,pk,1,xk,1)
             pbp = ddot(length_eigen,pk,1,pk,1)
             ! write(0,*) pax,pap,pbx,pbp 
             aa    = pap*pbx - pax*pbp
             bb    = xax*pbp - xbx*pap
             cc    = xbx*pax - xax*pbx

             delta = bb**2 - 4.d0*aa*cc
             dd    = sqrt(delta)
             alpha  = 0.5d0*(bb + dd)/aa
             ! write(0,*) 'beta',beta,'alfa',alpha

             eig_old = eig_new
             xbx  = xbx + 2.d0*alpha*pbx + alpha*alpha*pbp
             xax  = xax + 2.d0*alpha*pax + alpha*alpha*pap
             eig_new = xax/xbx
             err = (eig_new - eig_old)/eig_new
             !       write(6,*) err,eig_new

             !  
             ! *  x_k+1 = x_k + alpha_k p_k    
             ! 
             call info%scalprod%set('stop')

             call info%vecvec%set('start')
             call daxpy(length_eigen,alpha,pk,1,xk,1)
             call daxpy(length_eigen,alpha,ap,1,ax,1)
             call info%vecvec%set('stop')

             tdacg = tdacg + t

             !      iprt = 2
             if (iprt.eq.2) call prt1(lun_out,&
                  indau,iter,totit,&
                  eig_new,err,residum,time,rtime,info)

             ! exit if max number of iterations is reached
             if ( iter .eq. imax ) test_exit = .true.

             !        
             ! * end of inner loop
             !c
          end if
       end do
       info%iterations(indau) = iter
       this%residua(indau)    = residum
       if ( iter .ge. imax ) then
          ierr=ierr+1
          info%ierr = ierr
          info%totit = totit + iter
          this%computed=1
          call info%tot%set('stop')
          return
       end if

       !  
       ! * B-normalization of the eigenvector x
       !        
       a1 = 1.d0/sqrt(xbx)
       call info%vecvec%set('start')
       do j = 1,length_eigen
          this%eigenvectors%coeff(j,indau) = xk(j)*a1
       end do
       call info%vecvec%set('stop')

       this%eigval(indau) = eig_new
       if ( ( this%first_nonzero .eq. 0) .and. &
            (abs(this%eigval(indau)) .gt. small) )  then
          this%first_nonzero=indau
       end if
          
       if (iprt.eq.1) call prt2(lun_out,&
            indau,iter,totit,&
            eig_new,err,residum/normax,time,rtime,info)

    end do


    call info%tot%set('stop')


    info%ierr  = ierr
    info%totit = totit

    this%computed=1

!!$    deallocate(&
!!$         g_old,&
!!$         xk,&
!!$         pk,&
!!$         ax,&
!!$         g_new,&
!!$         ap,&
!!$         scr,&
!!$         gprec,&
!!$         aux1,&
!!$         stat=res)
!!$    if(res .ne. 0) rc = IOerr(lun_err, err_dealloc, 'dacg', &
!!$         'work array g_old, xk, pk, ax, g_new, ap, scr, gprec,',res)
    
    ! free memory
    aux_loc => null() 
    if ( aux_work%is_initialized ) call aux_work%kill(lun_err)  
    
  contains         
    !
    subroutine betak(iter,length_eigen,gnormnew,gnormold,&
         beta,gprec,g_new,g_old)
      use Globals
      implicit none

      integer, intent(in) :: length_eigen
      integer, intent(in) ::  iter

      real(kind=double), intent(out  ) :: g_old(length_eigen)
      real(kind=double), intent(in   ) :: gprec(length_eigen)
      real(kind=double), intent(in   ) :: g_new(length_eigen)      
      real(kind=double), intent(out  ) :: beta
      real(kind=double), intent(inout) :: gnormnew
      !local
      real(kind=double)  ddot,gnormold

      if (iter.eq.1) then
         beta = 0.d0
         call dcopy(length_eigen,g_new,1,g_old,1)
         gnormnew = ddot(length_eigen,gprec,1,g_new,1)
      else
         gnormnew = ddot(length_eigen,gprec,1,g_new,1)
         beta = gnormnew/gnormold
         !     return
         !     * beta with authomatic 'restart'   a <-- g_k+1 K^{-1} g_k
         !     a = ddot(length_eigen,gprec,1,g_old,1)
         !     call dcopy(length_eigen,g_new,1,g_old,1)
         !     beta = (gnormnew - a)/gnormold
         !     small = 1.e-30
         !     if (gnormnew.lt.small) then
         !     write(6,*) 'small',iter,small
         !     beta = 0.d0
         !     endif
      end if
    end subroutine betak
    !
    subroutine prt1(lun_out,indau,iter,totit,eigv2,err,residum,&
         time,rtime,info)
      use Globals
      implicit none
      integer,           intent(in) :: lun_out
      integer,           intent(in) :: indau,iter,totit
      real(kind=double), intent(in) :: eigv2,err,residum
      real*4,            intent(inout) :: time,rtime
      type(output_dacg), intent(inout) :: info

      call info%tot%set('stop')
      rtime=real(info%tot%usr)
      write(lun_out,1127)indau,iter,totit,eigv2,err,residum,rtime
      call info%tot%set('start')

      return 

1127  format(i3,i4,i8,2x,e16.8,2e15.5,f10.3)
      !2127  format(i3,i4,3i5,2x,e14.8,2e14.5,f10.3)
    end subroutine prt1

    !
    subroutine prt2(lun_out,indau,iter,totit,&
         eigv2,err,residum,time,rtime,info)
      use Globals
      implicit none
      integer,           intent(in) :: lun_out
      integer,           intent(in) :: indau
      integer,           intent(inout) :: iter,totit
      real(kind=double), intent(in) :: eigv2,err,residum
      real*4,            intent(inout) :: time,rtime
      type(output_dacg), intent(inout) :: info

      iter = iter - 1
      call info%tot%set('stop')

      !rtime = time  + rtime
      rtime=real(info%tot%usr)

      totit = totit + iter


      !write(6,'(i5,2e20.12,f10.3)') indau,eigv2,err,res
      write(lun_out,1138)indau,iter,totit,eigv2,err,residum,rtime
      !write(1111,1138)indau,iter,totit,eigv2,err,residum,rtime
      call info%tot%set('start')
1138  format(i3,1x,i4,i8,2x,e16.8,2e15.5,f10.3)
      !2138  format(i3,i4,3i5,2x,e14.8,2e14.5,f10.3)

    end subroutine prt2
    !
    subroutine reort(n,ndim,nev,p,autvec,x,tort,aux1)
      use Globals
      implicit none
      integer n,p,ndim,nev
      real(kind=double) autvec(ndim,nev),x(n),aux1(nev)
      real*4 tt,tort

      if (p.ge.1) then
         call dgemv('T',n,p,1.d0,autvec,ndim,x,1,0.d0,aux1,1)
         call dgemv('N',n,p,-1.d0,autvec,ndim,aux1,1,1.d0,x,1)
      end if
      tort = tort + tt

      return
    end subroutine reort
    !    
    subroutine initvec(length_eigen,indau,np,xk)
      use Globals
      implicit none
      !     
      !     *  initial vector guess
      !     
      integer  length_eigen,indau,np,j

      real(kind=double)   xk(length_eigen)

      do j = 1,length_eigen
         xk(j) = 1.d0
         if (j.eq.indau) xk(j) = 0.d0
         if (j.gt.np) xk(j) = 1.d0
      end do
    end subroutine initvec

  end subroutine dacg_algorithm


  subroutine lanczos_eigprec(this,&
       lun_err,&
       lun_out,&
       abpres,&
       matrix,&
       prec)
    use Globals
    use Matrix
    use LinearSolver
    implicit none

    ! input data structures
    class(eigen),      intent(inout) :: this
    integer,           intent(in   ) :: lun_err
    integer,           intent(in   ) :: lun_out
    type(pcgcoeff),    intent(in   ) :: abpres
    class(abs_matrix), intent(inout) :: matrix
    class(abs_linop),   intent(inout) :: prec
    
    ! local data structures for compute eigprec
    logical :: rc
    integer :: res
    integer :: i,iter,info,length_pres,nev,npres_stored,npres
    real(kind=double), allocatable :: matq(:,:)
    real(kind=double), allocatable :: aux(:)
    real(kind=double), allocatable :: AMAT(:,:)
    real(kind=double), allocatable :: Ddiag(:)
    real(kind=double), allocatable :: Tdiag(:)
    real(kind=double), allocatable :: WMAT(:,:)
    real(kind=double) :: den, normres,dnrm2

    

    length_pres  = abpres%length_pres
    npres_stored = abpres%npres_stored
    npres        = abpres%npres
    nev          = this%nev
    
    if ( npres_stored .le. nev ) return
    
    !
    ! allcocate workin arrays
    !
    allocate(&
         Ddiag(npres_stored),&
         Tdiag(npres_stored),&
         WMAT(length_pres,npres_stored),&
         matq(npres_stored,npres_stored),&
         aux(length_pres),&
         AMAT(length_pres, nev),&
         stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'lanczos_eigprec', &
         'local arrays Ddiag, Tdiag, matq, aux,AMAT',res)

    !
    ! build tri-diagonal matrix for Lanczos decomposition
    !
    do iter=1,npres_stored
       call dcopy(length_pres,abpres%pres(1,iter),1,WMAT(1,iter),1)
       call dscal(length_pres,&
            1.d0/sqrt(abpres%pres_dot_res(iter)),&
            WMAT(1,iter),1)

       if (iter .eq. 1) then
          Ddiag(iter) = 1/abpres%alpha(iter)        
       else
          Ddiag(iter)   = &
               1.d0/abpres%alpha(iter) + &
               abpres%beta(iter)/abpres%alpha(iter-1)
          Tdiag(iter-1) = -sqrt(abpres%beta(iter))/abpres%alpha(iter-1)
       end if
    end do
       

    !
    ! Compte eigenvales and eigenvectors
    !
    call dstev('V',npres_stored,Ddiag,Tdiag,matq,npres_stored,aux,info)
    call dcopy(nev,Ddiag,1,this%eigval,1)
    call dgemm('N','N',length_pres,nev,npres_stored,1.d0,&
         WMAT,length_pres,matq,npres_stored,0.d0,&
         this%eigenvectors%coeff,length_pres)
    this%first_nonzero = 1
    this%computed=1

    !
    ! compute residua
    !
    do i = 1,nev
       call matrix%matrix_times_vector(&
            this%eigenvectors%coeff(1:this%length_eigen,i) ,&
            AMAT(1:length_pres,i),info,lun_err)
       call prec%Mxv(AMAT(1:length_pres,i),aux,info,lun_err)
       den = dnrm2(length_pres,aux(1:length_pres),1)
       call daxpy(length_pres,-this%eigval(i),&
            this%eigenvectors%coeff(1,i),1,aux,1)
       normres = dnrm2(length_pres,aux(1:length_pres),1)/den
       write(lun_out,'(a,i6,f14.5,e15.3)') &
            'eig, relative eigenresidual',i,this%eigval(i),normres
       this%residua(i)=normres
    end do

    !
    ! free memory
    !

    deallocate(&
         Ddiag,&
         Tdiag,&
         WMAT,&
         matq,&
         aux,&
         AMAT,&
         stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_dealloc, 'lanczos_eigprec', &
         'local arrays Ddiad Tdiag WMAT matq, aux, AMAT',res)
 

  end subroutine lanczos_eigprec


  !>----------------------------------------------------------
  !> Procedure computing the eigvetoors and eigvalues of the
  !> for matrix M via DACG algorithm
  !> (procedure public for type eig)
  !>
  !> usage:
  !>     call 'var'%init(lun_err,lun_out,&
  !>                     matrix, prec    )
  !>
  !> where:
  !> \param[in] id_prt -> integer. Flag dor printing eigval.
  !> \param[in] lun    -> integer. Logic unit
  !<----------------------------------------------------------
  subroutine arpack_wrap( matrixA, lun_err,&
       nev,ncv,which, maxit,tol, &
       eigenvalues, eigenvectors,&
       info,matrixB,print_info)
    use Globals
    use Matrix
    use Scratch
    !include 'debug.h'
    !implicit none
    class(abs_linop), intent(inout) :: matrixA
    integer, intent(in ) ::lun_err
    integer, intent(in ) :: nev
    integer, intent(in ) :: ncv
    character(len=2), intent(in) :: which
    integer, intent(in ) :: maxit
    real(kind=double), intent(in   ) :: tol
    real(kind=double), intent(inout) :: eigenvalues(nev)
    real(kind=double), intent(inout) :: eigenvectors(matrixA%ncol,nev)
    integer, intent(inout) :: info
    class(abs_matrix), optional, intent(inout) :: matrixB
    logical, optional, intent(in ) :: print_info
    ! local
    ! work arrays
    real(kind=Double) , allocatable :: ax(:), d(:,:), resid(:)
    real(kind=Double) , allocatable :: v(:,:), workd(:)
    real(kind=Double) , allocatable :: workev(:)
    real(kind=Double) , allocatable :: workl(:),ones(:)
    ! controls array
    integer ::           iparam(11), ipntr(11)
    integer :: ndigit, logfil, msgets, msaitr, msapps, msaupd
    integer :: msaup2 , mseigt, mseupd 
    
    ! %---------------%
    ! | Local Scalars |
    ! %---------------%
    !
    character bmat*1
    integer ::ido, n, nx, lworkl,ldv,  ierr,infoAxv
    integer :: j, ishfts, maxitr, mode1, nconv,maxnev,maxn,maxncv
    real(kind=double) :: sigmar, sigmai,sigma,ddot
    logical :: first, rvec
    logical, allocatable :: select(:)


    ! %-----------------------------%
    ! | BLAS & LAPACK routines used |
    ! %-----------------------------%
    real(kind=double) :: dlapy2, dnrm2
  
    !
    ! set dimensions to the maximum
    ! since we are using allocatable arrays
    !
    maxn=matrixA%nrow
    n=maxn
    maxncv=ncv
    maxnev=nev
    ldv=maxn
    if ( matrixA%is_symmetric) then
       lworkl = ncv*(ncv+8)
    else
       lworkl  = 3*ncv**2+6*ncv
    end if
    allocate(ax(maxn), d(maxncv,2), resid(maxn), &
         v(ldv,maxncv), workd(3*maxn), &
         workev(3*maxncv), &
         workl(lworkl))
    d(:,:) = 0.d0
    workl(:) = 0.d0 
    
    allocate(select(maxncv))

    
    ndigit = -3
    logfil = 6
    msgets = 0
    msaitr = 0 
    msapps = 0
    msaupd = 1
    msaup2 = 0
    mseigt = 0
    mseupd = 0

    !
    ! set controls
    !
    if ( present(matrixB) ) then
       bmat='G'
    else
       bmat  = 'I'
    end if
    
    ishfts = 1
    maxitr = 300
    mode1 = 1
    !
    iparam(1) = ishfts
    iparam(3) = maxit
    iparam(7) = mode1

    !
    ! start loop
    !
    ido    = 0
    ierr   = 0

    allocate(ones(n))
    ones(:)=1.0d0

    call random_number(resid)
    resid=resid - (ddot(n,resid,1,ones,1)/ddot(n,ones,1,ones,1))*ones
    resid=resid/dnrm2(n,resid,1)
    
    
    !
    ! use dsaupd + dseupd
    ! ( dnaupd + dneupd for non-symmetric matrix)
    !
    if ( matrixA%is_symmetric) then
       do !while( ido .ne. 99 )
          call dsaupd ( ido, bmat, n, which, nev, tol, resid, &
               ncv, v, ldv, iparam, ipntr, workd, workl,&
               lworkl, ierr )

          
          select case (ido)
          case (1)
             !infoAxv=0
             call matrixA%matrix_times_vector(workd(ipntr(1):ipntr(1)+n-1),&
                  workd(ipntr(2):ipntr(2)+n-1),infoAxv,lun_err)
             if (infoAxv .ne.0 ) then
                ierr = -100
                exit
             end if
          case (-1)
             !infoAxv=0
             call matrixA%matrix_times_vector(workd(ipntr(1):ipntr(1)+n-1),&
                  workd(ipntr(2):ipntr(2)+n-1),infoAxv,lun_err)
             if (infoAxv .ne. 0 ) then
                ierr = -100
                exit
             end if
          case (99)
             exit
          case default
             
          end select
          
          
       end do

       info=ierr
       

       !
       ! %----------------------------------------%
       ! | Either we have convergence or there is |
       ! | an error.                              |
       ! %----------------------------------------%
       !
       if ( ierr .lt. 0 ) then
          !
          ! %--------------------------%
          ! | Error message. Check the |
          ! | documentation in DSAUPD. |
          ! %--------------------------%
          !
          print *, ' '
          print *, ' Error with _saupd, ierr = ', ierr
          print *, ' Check documentation in _saupd '
          print *, ' '
          !
       else 
          !
          ! %-------------------------------------------%
          ! | No fatal errors occurred.                 |
          ! | Post-Process using DSEUPD.                |
          ! |                                           | 
          ! | Computed eigenvalues may be extracted.    |  
          ! |                                           |
          ! | Eigenvectors may be also computed now if  |
          ! | desired.  (indicated by rvec = .true.)    | 
          ! |                                           |
          ! | The routine DSEUPD now called to do this  |
          ! | post processing (Other modes may require  |
          ! | more complicated post processing than     |
          ! | mode1.)                                   |
          ! |                                           |
          ! %-------------------------------------------%
           
          rvec = .true.

          call dseupd ( rvec, 'All', select, d, v, ldv, sigma,&
               bmat, n, which, nev, tol, resid, ncv, v, ldv,& 
               iparam, ipntr, workd, workl, lworkl, ierr )
          info=ierr
          !
          ! %----------------------------------------------%
          ! | Eigenvalues are returned in the first column |
          ! | of the two dimensional array D and the       |
          ! | corresponding eigenvectors are returned in   |
          ! | the first NCONV (=IPARAM(5)) columns of the  |
          ! | two dimensional array V if requested.        |
          ! | Otherwise, an orthogonal basis for the       |
          ! | invariant subspace corresponding to the      |
          ! | eigenvalues in D is returned in V.           |
          ! %----------------------------------------------%
          !
          if ( ierr .ne. 0) then
             !
             ! %------------------------------------%
             ! | Error condition:                   |
             ! | Check the documentation of DSEUPD. |
             ! %------------------------------------%
             !
             print *, ' '
             print *, ' Error with _seupd, info = ', ierr
             print *, ' Check the documentation of _seupd. '
             print *, ' '
             !
          else
             !
             if ( present(print_info) ) then
                if (print_info) then
                   nconv =  iparam(5)
                   do j=1, nconv
                      !     
                      ! %---------------------------%
                      ! | Compute the residual norm |
                      ! |                           |
                      ! |   ||  A*x - lambda*x ||   |
                      ! |                           |
                      ! | for the NCONV accurately  |
                      ! | computed eigenvalues and  |
                      ! | eigenvectors.  (iparam(5) |
                      ! | indicates how many are    |
                      ! | accurate to the requested |
                      ! | tolerance)                |
                      ! %---------------------------%
                      !
                      call matrixA%Mxv(v(1:n,j),ax(1:n),info,lun_err)
                      call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                      d(j,2) = dnrm2(n, ax, 1)
                      d(j,2) = d(j,2) / abs(d(j,1))
                      !     
                   end do
                end if
             end if
          end if
             
          !     
          ! %-----------------------------%
          ! | Display computed residuals. |
          ! %-----------------------------%
          !     
          !call dmout(6, nconv, 2, d, maxncv, -6, &
          !     'Ritz values and relative residuals')
       end if
       !     
       ! %-------------------------------------------%
       ! | Print additional convergence information. |
       ! %-------------------------------------------%
       !     
       if ( ierr .eq. 1) then
          print *, ' '
          print *, ' Maximum number of iterations reached.'
          print *, ' '                                          
       else if ( ierr .eq. 3) then
          print *, ' ' 
          print *, ' No shifts could be applied during implicit Arnoldi update, try increasing NCV.'
          print *, ' '
       end if
       !
       if ( present(print_info) ) then
          if (print_info) then
             print *, ' '
             print *, ' Size of the matrix is ', n
             print *, ' The number of Ritz values requested is ', nev
             print *, ' The number of Arnoldi vectors generated',&
                  ' (NCV) is ', ncv
             print *, ' What portion of the spectrum: ', which
             print *, ' The number of converged Ritz values is ', &
                  nconv 
             print *, ' The number of Implicit Arnoldi update',&
                  ' iterations taken is ', iparam(3)
             print *, ' The number of OP*x is ', iparam(9)
             print *, ' The convergence criterion is ', tol
             print *, ' '
          end if
       end if
    end if

    !
    ! copy result
    !
    eigenvalues(1:nev)=d(1:nev,1)
    do j=1,nev
       eigenvectors(1:n,j)=v(1:n,j)
    end do

    !
    ! free memory
    !
    deallocate(ax, d,resid , v, workd, workev, workl)    
    deallocate(select)
       

  end subroutine arpack_wrap

  subroutine study_eigen(use_arpack, use_lapack,nev,matrix)
    use Globals
    use LinearOperator
    use Densematrix
    use SimpleMatrix
    implicit none  
    logical,intent(in) :: use_arpack, use_lapack
    integer,intent(in) :: nev  
    class(abs_linop),intent(inout) :: matrix
    !local
    integer :: nrow,ncol,i,j
    integer :: info,stderr=6
    real(kind=double) :: dnrm2
    real(kind=double), allocatable :: scr(:)
    type(diagmat) :: eigenvalues
    type(densemat) :: eigenvectors
    type(densemat) :: matrix_full_formed

    nrow=matrix%nrow
    ncol=matrix%ncol

    if ( use_arpack ) then
       call eigenvalues%init(stderr,nev)
       call eigenvectors%init(stderr,matrix%nrow,nev) 
       allocate(scr(matrix%nrow))

       ! compute minimal eigenvalue
       !
       call arpack_wrap( &
            matrix,6,&
            nev,nev*4,'BE', 1000,1d-10,&
            eigenvalues%diagonal, eigenvectors%coeff,&
            i)
       write(*,*) 'Info ARAPCK=', i

       !
       ! print residual
       !
       write(*,*) 'ARpack Eigenvalues'
       do i=1,nev/2
          scr=zero
          call matrix%Mxv(eigenvectors%coeff(1:ncol,i),scr,info,6)
          scr=scr-eigenvalues%diagonal(i)*eigenvectors%coeff(1:ncol,i)
          write(*,*) 'Eig ', i, ' = ', eigenvalues%diagonal(i), &
               ' res= ', dnrm2(ncol,scr,1)/eigenvalues%diagonal(i)
       end do
       do i=1,nev/2
          j=ncol+i-nev/2
          scr=zero
          call matrix%Mxv(eigenvectors%coeff(1:ncol,i),scr,info,6)
          scr=scr-eigenvalues%diagonal(i)*eigenvectors%coeff(1:ncol,i)
          write(*,*) 'Eig ', j, ' = ', eigenvalues%diagonal(i), &
               ' res= ', dnrm2(ncol,scr,1)/eigenvalues%diagonal(i)
       end do
       write(*,*) 'Max/min ', &
            eigenvalues%diagonal(nev)/eigenvalues%diagonal(1)

       call eigenvalues%kill(stderr)
       call eigenvectors%kill(stderr)
       deallocate(scr)
    end if



    if (use_lapack) then
       call matrix_full_formed%form(stderr,matrix)
       call eigenvalues%init(stderr,matrix%nrow)
       call eigenvectors%init(stderr,matrix%nrow,matrix%ncol)
       allocate(scr(matrix%nrow))

       write(*,*) 'Spectral decomposition starts'
       info=0
       call matrix_full_formed%spectral_decomposition(stderr,&
            info,&
            eigenvalues%diagonal,&
            eigenvectors%coeff)

       write(*,*) 'Spectral decomposition completed. info= ',info

       write(*,*) 'Lapack Eigenvalues'
       do i=1,nev/2
          scr=zero
          call matrix%Mxv( eigenvectors%coeff(1:ncol,i),scr,info,6)
          scr=scr- &
               eigenvalues%diagonal(i)*&
               eigenvectors%coeff(1:ncol,i)
          write(*,*) 'Eig ', i, ' = ',  eigenvalues%diagonal(i), &
               ' res= ', dnrm2(ncol,scr,1)/ eigenvalues%diagonal(i)
       end do
       do i=1,nev/2
          j=ncol+i-nev/2
          scr=zero
          call matrix%Mxv( eigenvectors%coeff(1:ncol,j),scr,info,6)
          scr = scr-&
               eigenvalues%diagonal(j)* &
               eigenvectors%coeff(1:ncol,j)
          write(*,*) 'Eig ', j, ' = ',  eigenvalues%diagonal(j), &
               ' res= ', dnrm2(ncol,scr,1)/ eigenvalues%diagonal(j)
       end do
       write(*,*) 'Max/min ', &
            eigenvalues%diagonal(ncol)/eigenvalues%diagonal(1)

       call matrix_full_formed%kill(stderr)
       call eigenvalues%kill(stderr)
       call eigenvectors%kill(stderr)
       deallocate(scr)
    end if
  end subroutine study_eigen

end module Eigenv

