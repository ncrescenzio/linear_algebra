module RankOneUpdate
  use Globals
  use LinearOperator
  !use Preconditioner
  implicit none
  private
  !> Structure variable for the application of the 
  !> Sherman-Morrison formula. Given a rank-one update
  !>      P = P0 + alpha u v^T
  !> its inverse reads as
  !> P^{-1} = (P0)^{-1} - alpha *(P0)^{-1}u v^T (P0)^{-1}/(1 + alpha v^T (P0)^{-1} u)
  !>
  !> P0 = preconditioner=linear inverible operator
  !> u,v= vectors
  !> alpha=positive scalar (default =one)
  !> This type contains a pointer to the P0 preconditoner, 
  !> the vectors u, v and other auxiliary variables
  !>------------------------------------------------------------------------
  type, extends(abs_linop), public :: shmr
     !> Abract invertible linear operator
     class(abs_linop), pointer :: prec_zero => null()
     !> Vector u of rank-1 upate
     !> Dimension(prec_zero%nequ)
     real(kind=double), allocatable :: uvec(:)
     !> Vector v of rank-1 upate
     !> Dimension(prec_zero%nequ)
     real(kind=double), allocatable :: vvec(:)
     !> Vector y = prec_zero^{-1} uvec 
     !> Dimension(prec_zero%nequ)
     real(kind=double), allocatable :: prec_zero_uvec(:)
     !> Scalar 1.0/(1.0 + v^T (P_0)^{-1} uvec)  
     !> Dimension(prec_zero%nequ)
     real(kind=double)  :: den=zero
     !> Scalar alpha
     !> Dimension(prec_zero%nequ)
     real(kind=double)  :: alpha=one    
   contains
     !> static constructor
     !> (procedure public for type shmr)
     procedure, public, pass :: init => init_shmr
     !> static constructor
     !> (procedure public for type shmr)
     procedure, public, pass :: set => set_shmr
     !> static destructor
     !> (procedure public for type shmr)
     procedure, public, pass :: kill => kill_shmr
     !> Info procedure.
     !> (public procedure for type shmr)
     procedure, public, pass :: matrix_times_vector => apply_shmr
  end type shmr


  !> Structure variable for the application of the 
  !> Broyen update for preconditioning the jacobian in the
  !> Newton scheme
  !>    J_{F}(x_{k}) s_{k} = -F(x)
  !>              x_{k+1}  = x_{k} + s_{k}
  !> given a rank-one update
  !>      P_{k+1} = P_{k} + u v^T
  !> with 
  !>     u=(y_{k}-P_{k} s_k)/ (v^T s_k)
  !>     v=s_{k}/||s_k||
  !> where
  !>     y_k=F(x_{k+1})-F(x_{k+1})
  !> By the Sherman Morrison forumla 
  !>  we obtain 
  !> P_{k+}1 ^{-1} = P_{k})^{-1} 
  !>                 - ( P_{k}^{-1} y_{k} - s_{k} ) s_{k}^T P_{k}^{-1} ) /
  !>                   ( s_{k}^T P_{k}^{-1} y_{k} )
  !>
  !> We store the vectors 
  !>   up_k  = s_{k}
  !>   pres_k = (P_{k}^{-1} y_{k} - s_{k}) / ( s_{k}^T P_{k}^{-1} y_{k} )
  !>------------------------------------------------------------------------
  type, extends(abs_linop), public :: broyden
     !> Abract invertible linear operator
     class(abs_linop), pointer :: prec_k => null()
     !> Vector u of rank-1 upate
     !> Dimension(prec_zero%nequ)
     real(kind=double), allocatable :: ess_k(:)
     !> Vector (P_{k}^{-1} y_{k} - s_{k}) / ( s_{k}^T P_{k}^{-1} y_{k} )
     !> Dimension(prec_zero%nequ)
     real(kind=double), allocatable :: pres_k(:)
   contains
     !> static constructor
     !> (procedure public for type broyden)
     procedure, public, pass :: init => init_broyden
     !> Set Prec_k, ess_k, pres_k
     !> (procedure public for type broyden)
     procedure, public, pass :: set => set_broyden
     !> static destructor
     !> (procedure public for type broyden)
     procedure, public, pass :: kill => kill_broyden
     !> Info procedure.
     !> (public procedure for type broyden)
     procedure, public, pass :: matrix_times_vector => apply_broyden
  end type broyden


  type, extends(abs_linop), public :: rankk
     !> Abract invertible linear operator
     class(abs_linop), pointer :: prec_zero => null()
     !> Maximal updates
     integer :: nrankmax
     !> Maximal updates
     integer :: nrank
     !> Rank one updates 
     type(shmr), allocatable :: updates(:)
     !> wotk arrays
     real(kind=double), allocatable :: work(:)
   contains
     !> static constructor
     !> (procedure public for type rankk)
     procedure, public, pass :: init => init_rankk
     !> static constructor
     !> (procedure public for type rankk)
     procedure, public, pass :: set => set_rankk
     !> static destructor
     !> (procedure public for type rankk)
     procedure, public, pass :: kill => kill_rankk
     !> Info procedure.
     !> (public procedure for type rankk)
     procedure, public, pass :: matrix_times_vector => apply_rankk
  end type rankk
contains
  subroutine init_shmr(this,lun_err,nequ)
    use Globals
    implicit none
    class(shmr),     intent(inout) :: this
    integer,         intent(in   ) :: lun_err
    integer,         intent(in   ) :: nequ
    !local
    logical :: rc
    integer :: res

    this%nrow         = nequ
    this%ncol         = nequ

    allocate(&
         this%prec_zero_uvec(this%nrow),&
         this%uvec(this%nrow),&
         this%vvec(this%nrow),&
         stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'init_shmr', &
         ' type shmr member uvec, vvec, prec_zero_uvec',res)

  end subroutine init_shmr


  subroutine set_shmr(this,lun_err,&
       prec_zero,uvec,vvec,is_symmetric,info,alpha)
    use Globals
    implicit none
    class(shmr),             intent(inout) :: this
    integer,                 intent(in   ) :: lun_err
    class(abs_linop), target, intent(in   ) :: prec_zero
    real(kind=double),       intent(in   ) :: uvec(prec_zero%nrow)
    real(kind=double),       intent(in   ) :: vvec(prec_zero%nrow)    
    logical,                 intent(in   ) :: is_symmetric
    integer,        intent(inout) :: info
    real(kind=double),optional, intent(in   ) :: alpha    
    !local
    logical :: rc
    integer :: res,nequ
    real(kind=double) :: ddot
    
    nequ         = prec_zero%nrow
    
    this%nrow         = nequ
    this%ncol         = nequ
    
    this%is_symmetric = prec_zero%is_symmetric
    
    !
    ! store data and pointer
    !
    this%prec_zero => prec_zero
    this%uvec     = uvec
    this%vvec     = vvec
    
    !
    ! set alpha if passed (deafult=one)
    !
    if ( present(alpha) ) then
       if ( alpha .le. zero ) then
          rc = IOerr(lun_err, wrn_inp, 'init_shmr', &
               'not invertible update',res)
          info = -1
       else
          this%alpha = alpha
       end if
    end if

    call this%prec_zero%Mxv(uvec,this%prec_zero_uvec,info,lun_err)

    this%den = one + this%alpha * ddot(this%nrow,this%prec_zero_uvec,1,this%vvec,1)
    
    if ( abs( this%den) > small ) then
       this%den = this%alpha / this%den
    else
       rc = IOerr(lun_err, err_inp, 'init_shmr', &
            'not invertible update',res)
       info = -1
    end if
  end subroutine set_shmr


  subroutine kill_shmr(this,lun_err)
    use Globals
    implicit none
    class(shmr),     intent(inout) :: this
    integer,         intent(in   ) :: lun_err
    !local
    logical :: rc
    integer :: res 

    this%nrow      = -1
    this%is_symmetric = .false.
    this%prec_zero => null()

    deallocate(&
         this%prec_zero_uvec,&
         this%uvec,&
         this%vvec,&
         stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_dealloc, 'kill_shmr', &
         ' type shmr member prec_zero_uvec',res)

  end subroutine kill_shmr

  recursive subroutine apply_shmr(this,vec_in,vec_out,info,lun_err)
    use Globals
    implicit none
    class(shmr),       intent(inout) :: this
    real(kind=double), intent(in   ) :: vec_in(this%ncol)
    real(kind=double), intent(inout) :: vec_out(this%nrow)
    integer,           intent(inout) :: info
    integer,           intent(in   ) :: lun_err
    ! local
    real(kind=double) :: vt_prec_zero_u
    real(kind=double) :: ddot

    call this%prec_zero%Mxv(vec_in,vec_out,info,lun_err)     
    vt_prec_zero_u     = ddot(this%nrow,vec_out,1,this%vvec,1)

    call daxpy(this%nrow, - vt_prec_zero_u * this%den, &
         this%prec_zero_uvec,1,&
         vec_out,1)

  end subroutine apply_shmr

  subroutine init_broyden(this,lun_err,nequ)
    use Globals
    implicit none
    class(broyden), intent(inout) :: this
    integer,        intent(in   ) :: lun_err
    integer,        intent(in   ) :: nequ
    !local
    logical :: rc
    integer :: res

    this%nrow         = nequ
    this%ncol         = nequ

    ! allocate space
    allocate(&
         this%ess_k(this%nrow),&
         this%pres_k(this%nrow),&
         stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'init_broyden', &
         ' type broyden member ess_k,pres_k',res)

  end subroutine init_broyden


  subroutine set_broyden(this,lun_err,prec_k,ess_k,y_k,info)
    use Globals
    implicit none
    class(broyden),          intent(inout) :: this
    integer,                 intent(in   ) :: lun_err
    class(abs_linop), target, intent(in   ) :: prec_k
    real(kind=double),       intent(in   ) :: ess_k(this%nrow)
    real(kind=double),       intent(in   ) :: y_k(this%ncol)
    integer,       optional, intent(inout) :: info
    !local
    logical :: rc
    integer :: res,nequ
    real(kind=double) :: den
    real(kind=double) :: ddot

    ! assign prec_k
    nequ      = prec_k%nrow
    
    this%nrow         = nequ
    this%ncol         = nequ
    
    this%prec_k    => prec_k


    ! copy ess_k
    this%ess_k = ess_k

    ! build pres_k
    call this%prec_k%Mxv(y_k,this%pres_k,info,lun_err)
    den = ddot(this%nrow,this%pres_k,1,ess_k,1)
    if ( abs( den) < small ) then
       rc = IOerr(lun_err, wrn_inp, 'set_broyden', &
            'not invertible update',res)
       info = 1
    end if
    this%pres_k = ( this%pres_k - this%ess_k ) / den   

  end subroutine set_broyden

  subroutine kill_broyden(this,lun_err)
    use Globals
    implicit none
    class(broyden), intent(inout) :: this
    integer,        intent(in   ) :: lun_err
    
    ! local
    logical :: rc
    integer :: res

    ! deassociate pointer
    this%prec_k   => null()

    ! allocate space
    deallocate(&
         this%ess_k,&
         this%pres_k,&
         stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_dealloc, 'init_broyden', &
         ' type broyden member prec_k_inv_uvec',res)

  end subroutine kill_broyden

  recursive subroutine apply_broyden(this,vec_in,vec_out,info,lun_err)
    use Globals
    implicit none
    class(broyden),    intent(inout) :: this
    real(kind=double), intent(in   ) :: vec_in(this%ncol)
    real(kind=double), intent(inout) :: vec_out(this%nrow)
    integer,           intent(inout) :: info
    integer,           intent(in   ) :: lun_err
    ! local
    real(kind=double) :: scalar
    real(kind=double) :: ddot

    call this%prec_k%Mxv(vec_in,vec_out,info,lun_err)     
    scalar     = ddot(this%nrow,this%ess_k,1,vec_out,1)

    vec_out = vec_out - scalar * this%pres_k

  end subroutine apply_broyden

  subroutine init_rankk(this,lun_err,nequ,nrankmax)
    use Globals
    implicit none
    class(rankk),    intent(inout) :: this
    integer,         intent(in   ) :: lun_err
    integer,         intent(in   ) :: nequ
    integer,         intent(in   ) :: nrankmax
    !local
    logical :: rc
    integer :: res
    integer :: i


    this%nrow         = nequ
    this%ncol         = nequ
    this%nrankmax = nrankmax
    
    allocate(&
         this%updates(this%nrankmax),&
         this%work(nequ),&
         stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'init_rankk', &
         ' type rankk member updates work',res)

    do i=1,nrankmax
       call this%updates(i)%init(lun_err,nequ)
    end do

  end subroutine init_rankk


  subroutine set_rankk(this,lun_err,&
       prec_zero,nrank,&
       Uvectors,Vvectors,is_symmetric,info,alphas)
    use Globals
    implicit none
    class(rankk),             intent(inout) :: this
    integer,                 intent(in   ) :: lun_err
    class(abs_linop), target, intent(in   ) :: prec_zero
    integer, intent(in) :: nrank
    real(kind=double),       intent(in   ) :: Uvectors(prec_zero%nrow,nrank)
    real(kind=double),       intent(in   ) :: Vvectors(prec_zero%nrow,nrank)    
    logical,                 intent(in   ) :: is_symmetric
    integer,        intent(inout) :: info
    real(kind=double),optional, intent(in   ) :: alphas(nrank)
    !local
    logical :: rc
    integer :: res,i
    real(kind=double) :: ddot

    this%is_symmetric = is_symmetric
    
    this%nrank = nrank

    !
    ! store data and pointer
    !
    this%prec_zero => prec_zero

    call this%updates(1)%set(lun_err,&
         this%prec_zero,&
         Uvectors(:,1),Vvectors(:,1), is_symmetric,info,alphas(1))
    if ( info .ne. 0  ) then
       rc = IOerr(lun_err, wrn_inp, 'set_rankk', &
            ' update not inverible, nupdate = ', 1)
    end if
    do i=2,nrank
       call this%updates(i)%set(lun_err,&
            this%updates(i-1),&
            Uvectors(:,i),Vvectors(:,i), is_symmetric,info,alphas(i))
       if ( info .ne. 0  ) then
          rc = IOerr(lun_err, wrn_inp, 'set_rankk', &
               ' update not inverible, nupdate = ', i)
          exit
       end if
    end do
    
  end subroutine set_rankk


  subroutine kill_rankk(this,lun_err)
    use Globals
    implicit none
    class(rankk),     intent(inout) :: this
    integer,         intent(in   ) :: lun_err
    !local
    logical :: rc
    integer :: res 
    integer :: i

    this%nrow     = -1
    this%is_symmetric = .false.
    this%prec_zero => null()

    do i=1,this%nrankmax
       call this%updates(i)%kill(lun_err)
    end do

    deallocate(&
         this%updates,&
         this%work,&
         stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_dealloc, 'kill_rankk', &
         ' type rankk member updates work',res)

  end subroutine kill_rankk

  recursive subroutine apply_rankk(this,vec_in,vec_out,info,lun_err)
    use Globals
    implicit none
    class(rankk),      intent(inout) :: this
    real(kind=double), intent(in   ) :: vec_in(this%ncol)
    real(kind=double), intent(inout) :: vec_out(this%nrow)
    integer,           intent(inout) :: info
    integer,           intent(in   ) :: lun_err
    ! local
    integer :: i

!!$    call this%prec_zero%Mxv(vec_in,vec_out,info,lun_err)
!!$    do i=1, this%nrank
!!$       this%work = vec_out
!!$       call this%updates(i)%Mxv(this%work,vec_out,info,lun_err)
!!$    end do
    call this%updates(this%nrank)%Mxv(vec_in,vec_out,info,lun_err)
    
  end subroutine apply_rankk


end module RankOneUpdate


