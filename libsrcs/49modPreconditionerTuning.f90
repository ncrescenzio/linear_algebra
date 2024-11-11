module PreconditionerTuning
  use Globals
  use LinearOperator
  use SparseMatrix
  use Eigenv
  use DenseMatrix
  implicit none
  private
  !> Structure variable containing the controls
  !> for the tuning procedure for preconditioning
  !> P^{-1}= P0^{-1} + V^T D V
  !> with 
  !> V: collection of eigenvectors
  !> D: diagonal with the inverses of the eigevalues 
  type, public :: input_tuning
     !> Flag to control building/use of spectral informations
     !>     ieig = 0 => No tuning
     !>     ieig = 1 => Dacg
     !>     ieig = 2 => Lanczos
     integer :: ieig=0  
     !> Number of eigen vectors for tuned
     integer :: nev=0
     !> Flag controlling the tuning
     !>     ituning = 1 => 
     !>     ituning = 2 => 
     !>     ituning = 3 => 
     integer :: ituning=0         
     !> Controls for Dacg Procedure
     !> Used for ieig=1
     type(input_dacg) :: ctrl_dacg
     !> Number of prec. residual stored
     integer :: npres=0
   contains
     !> static constructor
     !> (procedure public for type input_tuning)
     procedure, public, pass :: init => init_input_tuning
     !> static destructor
     !> (procedure public for type input_tuning)
     procedure, public, pass :: kill => kill_input_tuning
     !> Info procedure.
     !> (public procedure for type input_tuning)
     procedure, public, pass :: info => info_input_tuning
  end type input_tuning  
  !> Structure variable containing the preconditioners 
  !> and the procedures to build and apply them
  type, extends(abs_linop), public :: tunedprec
     !> Standard sparse preconditioner
     class(abs_linop), pointer :: prec_zero => null()
     !> Controls for tunign procedure of prec_zero
     type(input_tuning) :: ctrl_tuning
     !> Spectral information
     type(eigen),  pointer :: spectral_info => null()
     !> Working arrays
     !> Dimension given by input_tuned
     real(kind=double), allocatable :: qaq (:)
   contains
     !> static constructor
     !> (public for type eig)
     procedure, public, pass :: init => init_tunedprec
     !> static destructor
     !> (public for type eig)
     procedure, public, pass :: kill => kill_tunedprec
     !> Info procedure
     !> (public for type eig)
     procedure, public, pass :: info => info_tunedprec
     !> Procedure for the application of the preconditier
     !> (public for type eig)
     procedure, public, pass :: matrix_times_vector => apply_tunedprec
  end type tunedprec

   !> Structure variable containing the preconditioners 
  !> and the procedures to build and apply them
  type, extends(abs_linop), public :: bfgs
     !> Standard sparse preconditioner
     class(abs_linop), pointer :: prec_zero => null()
     !> Standard sparse preconditioner
     class(abs_linop), pointer :: H => null()
     !> Standard sparse preconditioner
     class(densemat), pointer :: W => null()
     !> Standard sparse preconditioner
     class(densemat), pointer :: AV => null()
     !> Working arrays
     real(kind=double), allocatable :: scr (:)
     real(kind=double), allocatable :: u_scr(:)
     real(kind=double), allocatable :: y_scr (:)
     real(kind=double), allocatable :: w_scr (:)
   contains
     !> static constructor
     !> (public for type eig)
     procedure, public, pass :: init => init_bfgs
     !> static destructor
     !> (public for type eig)
     procedure, public, pass :: kill => kill_bfgs
     !> static destructor
     !> (public for type eig)
     procedure, public, pass :: set => set_bfgs
     !> Procedure for the application of the preconditier
     !> (public for type eig)
     procedure, public, pass :: matrix_times_vector => apply_bfgs
  end type bfgs


contains
  !>-------------------------------------------------------------
  !> Static constructor.
  !> (procedure public for type input_tuning)
  !> Instantiate and initilize variable of type input_tuning
  !>
  !> usage:
  !>     call 'var'%init(ieig,nev,ituning)
  !>
  !> where:
  !> \param[in] (optional) ieig      -> integer. (default=0)
  !>                                      ieig=0 => No spectral info
  !>                                      ieig=1 => Eigev via DACG
  !>                                      ieig=2 => Eigev via LANCZOS
  !> \param[in] (optional) nev       -> integer. (Used only if ieig>0)
  !>                                      Number of EigenVectors (default=0)
  !> \param[in] (optional) ituning    -> integer. (Used only if ieig>0)
  !>                                        (default=0)
  !>                                        ituning =>  1
  !>                                        ituning =>  2
  !> \param[in] (optional) ctrl_dacg -> type(input_dacg) (used only if ieig=1)
  !>                                      Controls for DACG Procedure
  !> \param[in] (optional) npres     -> integer. (Used only if ieig=2)
  !>                                      Number of prec. residual stores  
  !<-------------------------------------------------------------
  subroutine init_input_tuning(this&
       ,ieig&
       ,nev&
       ,ituning &
       ,ctrl_dacg &
       ,npres &
       )
    use Globals
    use Eigenv
    implicit none
    class(input_tuning),         intent(inout) :: this
    integer,           optional, intent(in   ) :: ieig
    integer,           optional, intent(in   ) :: nev
    integer,           optional, intent(in   ) :: ituning 
    type(input_dacg),  optional, intent(in   ) :: ctrl_dacg !(ieig=1)
    integer,           optional, intent(in   ) :: npres     !(ieig=2) 
    
    if ( present(ieig) )      this%ieig      = ieig
    if ( present(nev) )       this%nev       = nev
    if ( present(ituning) )   this%ituning    = ituning
    if ( present(ctrl_dacg) ) this%ctrl_dacg = ctrl_dacg
    if ( present(npres) )     this%npres     = npres

  end subroutine init_input_tuning


  !>-------------------------------------------------------------
  !> Static desconstructor.
  !> (procedure public for type input_tuning)
  !> Instantiate and initilize variable of type input_tuning
  !>
  !> usage:
  !>     call 'var'%kill()
  !<-------------------------------------------------------------
  subroutine kill_input_tuning(this)
    use Globals
    implicit none
    class(input_tuning), intent(inout) :: this
    !    integer,         intent(in   ) :: lun_err

    this%ieig   = 0
    this%nev    = 0
    this%ituning = 0 
    call this%ctrl_dacg%kill()
    this%npres  = 0
  end subroutine kill_input_tuning


  !>-------------------------------------------------------------
  !> Static desconstructor.
  !> (procedure public for type input_tuning)
  !> Instantiate and initilize variable of type input_tuning
  !>
  !> usage:
  !>     call 'var'%kill(lun_out)
  !>
  !> where:
  !> \param[in] lun_out -> integer. I/O unit for output mesg.
  !<-------------------------------------------------------------
  subroutine info_input_tuning(this,lun_out)
    use Globals
    implicit none
    class(input_tuning), intent(in) :: this
    integer,           intent(in) :: lun_out 
    !local
    character(len=256) :: method
    select case (this%ieig)
    case (0) 
       write(method,'(a)') 'NO TUNING' 
       write(lun_out,'(a)') etb(method)
    case (1) 
       write(method,'(a)') 'DACG'
       write(lun_out,'(a,I1,a,a,I1,a,I3,a,a,a)')&
            'ieig  = ', this%ieig,&
            ' (tuning via matrix eigen. /',&
            ' itunued = ', this%ituning, &
            ' / N. eig. = ', this%nev, &
            ' ',etb(method),')'
       call this%ctrl_dacg%info(lun_out)
    case (2) 
       write(method,'(a)') 'LANCZOS'
       write(lun_out,'(a,I1,a,a,I1,a,I3,a,a,a)')&
            'ieig=', this%ieig,&
            ' (tuning via prec. matrix eigen. /',&
            ' itunued = ', this%ituning, &
            ' / N. eigen. = ', this%nev, &
            ' ',etb(method),')'
       write(lun_out,'(a,I3)') 'Max Npres.=', this%npres
    end select
          
       
  end subroutine info_input_tuning
  !>-------------------------------------------------------------
  !> Static constructor.
  !> (procedure public for type input_tuning)
  !> Instantiate and initilize variable of type input_tuning
  !>
  !> usage:
  !>     call 'var'%init(ieig,nev,ituning)
  !>
  !> where:
  !> \param[in] nev                  -> integer. I/O err. msg. logical unit 
  !> \param[in] ctrl_tuning           -> type(input_tuning). Controls varaibles
  !>                                     for preconditioner application
  !> \param[in] (optional) ctrl_dacg -> type(input_dacg) (used only if ieig=1)
  !>                                      Controls for DACG Procedure
  !> \param[in] (optional) npres     -> integer. (Used only if ieig=2)
  !>                                      Number of prec. residual stores  
  !<-------------------------------------------------------------
  subroutine init_tunedprec(this&
       ,lun_err&
       ,prec_zero&
       ,ctrl_tuning&
       ,spectral_info&
       )
    use Globals
    use Eigenv

    implicit none
    class(tunedprec),       intent(inout) :: this
    integer,               intent(in   ) :: lun_err
    class(abs_linop), target, intent(in   ) :: prec_zero
    type(input_tuning),    intent(in   ) :: ctrl_tuning
    type(eigen),   target, intent(in   ) :: spectral_info
     !local 
    logical :: rc
    integer :: res
    
    this%is_initialized = .true.
    !
    ! set dimension and symmetry
    !
    this%nrow        = prec_zero%nrow
    this%ncol        = prec_zero%ncol
    this%is_symmetric = prec_zero%is_symmetric
    
    !
    ! set prec zero
    !
    this%prec_zero => prec_zero
    
    ! set tuning
    this%ctrl_tuning = ctrl_tuning
    this%spectral_info    => spectral_info
    allocate ( this%qaq(spectral_info%nev),&
         stat = res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'init_tunedprec', &
         '  type tunedprec member qaq',res)

  end subroutine init_tunedprec

  !>-------------------------------------------------------------
  !> Static desconstructor.
  !> (procedure public for type input_tuning)
  !> Instantiate and initilize variable of type input_tuning
  !>
  !> usage:
  !>     call 'var'%kill()
  !<-------------------------------------------------------------
  subroutine kill_tunedprec(this,lun_err)
    use Globals
    implicit none
    class(tunedprec), intent(inout) :: this
    integer,         intent(in   ) :: lun_err
    !local 
    logical :: rc
    integer :: res

    this%prec_zero        => null()
    this%spectral_info    => null()
    call this%ctrl_tuning%kill()
    
    if (allocated(this%qaq) ) then
       deallocate ( this%qaq, stat = res)
       if(res .ne. 0) rc = IOerr(lun_err, err_dealloc, 'kill_tunedprec', &
            '  type tunedprec member qaq',res)
    end if

    !
    ! reset to default
    !
    call this%to_default()
    
    
  end subroutine kill_tunedprec
     

  !>-------------------------------------------------------------
  !> Static desconstructor.
  !> (procedure public for type input_tuning)
  !> Instantiate and initilize variable of type input_tuning
  !>
  !> usage:
  !>     call 'var'%kill(lun_out)
  !>
  !> where:
  !> \param[in] lun_out -> integer. I/O unit for output mesg.
  !<-------------------------------------------------------------
  subroutine info_tunedprec(this,lun)
    use Globals
    implicit none
    class(tunedprec), intent(in) :: this
    integer,           intent(in) :: lun
    !local
    character(len=256) :: method
    select case (this%ctrl_tuning%ieig)
    case (0) 
       write(lun,*) 'no tuning'
    case (1) 
       write(method,'(a)') 'DACG'
       write(lun,'(a,I1,a,a,I1,a,I3,a,a,a)')&
            'ieig  = ', this%ctrl_tuning%ieig,&
            ' (tuning via matrix eigen. /',&
            ' itunued = ', this%ctrl_tuning%ituning, &
            ' / N. eig. = ', this%spectral_info%nev, &
            ' ',etb(method),')'
    case (2) 
       write(method,'(a)') 'LANCZOS'
       write(lun,'(a,I1,a,a,I1,a,I3,a,a,a)')&
            'ieig=', this%ctrl_tuning%ieig,&
            ' (tuning via prec. matrix eigen. /',&
            ' itunued = ', this%ctrl_tuning%ituning, &
            ' / N. eigen. = ', this%spectral_info%nev, &
            ' ',etb(method),')'
    end select
          
       
  end subroutine info_tunedprec

  !>----------------------------------------------------------
  !> Precedure to apply preconditioner
  !> (procedure public for type tunedprec)
  !> Prodecure computing the vector w=result of        
  !>              w=P^{-1}v
  !> for a given vector v (arrays' dimension must match)
  !>
  !> usage: call 'var'Mxv%(vec_in,vec_out,[info])
  !> 
  !> \param[in   ] vec_in  -> real, dimension(this%nequ)
  !>                            Vector v where tunedprec will be applied
  !> \param[inout] vec_out  -> real, dimension(this%nequ)
  !>                            Vector w=P^{-1}v 
  !> \param[inout] (optional) info -> integer. Flag for errors  
  !>---------------------------------------------------------------------
  subroutine apply_tunedprec(this, vec_in, vec_out,info,lun_err) 
    use Globals
    implicit none
    !var
    class(tunedprec),  intent(inout) :: this
    real(kind=double), intent(in   ) :: vec_in(this%nrow)     
    real(kind=double), intent(inout) :: vec_out(this%ncol)
    integer,           intent(inout) :: info
    integer,           intent(in   ) :: lun_err

    ! local 
    integer :: nequ

    nequ = this%nrow
    select case (this%ctrl_tuning%ieig)
    case (0)
       !
       ! vec_out = Prec vec_in
       !
       call this%prec_zero%Mxv(vec_in, vec_out,info,lun_err)
    case (1)
       select case (this%ctrl_tuning%ituning)
       case (1)
          !
          ! vec_out = Prec vec_in
          !
          call this%prec_zero%Mxv(vec_in, vec_out,info,lun_err)
          
          !
          ! vec_out = ( I + E diag(alpha_i) E^T ) Prec vec_in
          !
          call tuning(this%ctrl_tuning%ituning,&
               this%spectral_info,vec_in,vec_out,this%qaq)
       end select
    case (2)
       select case (this%ctrl_tuning%ituning)
       case (1)
          ! vec_out = Prec vec_in
          call this%prec_zero%Mxv(vec_in, vec_out,info,lun_err)
          
          ! vec_out = ( I + E diag(alpha_i) E^T ) Prec vec_in
          call tuning(this%ctrl_tuning%ituning,&
               this%spectral_info,vec_in,vec_out,this%qaq)
       end select
    end select
    info = 0
    
  contains 
    subroutine tuning(ituning, spectrum, vec_in,vec_out,qaq)
      use Globals
      implicit none
      !var
      integer,           intent(in   ) :: ituning
      type(eigen),       intent(in   ) :: spectrum
      real(kind=double), intent(in   ) :: vec_in(:)     
      real(kind=double), intent(out  ) :: vec_out(:)
      real(kind=double), intent(inout) :: qaq(spectrum%nev)
      ! local 
      integer :: nequ,nev,nnzev,i

      nequ = spectrum%length_eigen
      nev  = spectrum%nev
      nnzev = spectrum%nev - spectrum%first_nonzero

      select case (ituning)
      case (1)
         if (this%spectral_info%computed .eq. 1) then
            
            ! E       = (nequ x nev) matrix with eigenvectors 
            ! alpah_i = eigenvalues
            !
            ! vec_out = (I + E diag(alpha_i) E^T )(  vec_in )
            call dgemv('T',nequ,nnzev,one,&
                 spectrum%eigenvectors%coeff(1,spectrum%first_nonzero),&
                 nequ,vec_in,1,zero,this%qaq,1)
            !
            do i=spectrum%first_nonzero,nev
               qaq(i) = qaq(i)/spectrum%eigval(i)
            end do
            !
            call dgemv('N',nequ,nnzev,one,&
                 spectrum%eigenvectors%coeff(1,spectrum%first_nonzero),&
                 nequ,qaq,1,one,vec_out,1)
         end if
      end select
    end subroutine tuning


          
  end subroutine apply_tunedprec

  !>-------------------------------------------------------------
  !> Static constructor.
  !> (procedure public for type input_tuning)
  !> Instantiate and initilize variable of type input_tuning
  !>
  !> usage:
  !>     call 'var'%init(ieig,nev,ituning)
  !>
  !> where:
  !> \param[in] nev                  -> integer. I/O err. msg. logical unit 
  !> \param[in] ctrl_tuning           -> type(input_tuning). Controls varaibles
  !>                                     for preconditioner application
  !> \param[in] (optional) ctrl_dacg -> type(input_dacg) (used only if ieig=1)
  !>                                      Controls for DACG Procedure
  !> \param[in] (optional) npres     -> integer. (Used only if ieig=2)
  !>                                      Number of prec. residual stores  
  !<-------------------------------------------------------------
  subroutine init_bfgs(&
       this,&
       lun_err,&
       nequ&
       )
    use Globals
    use Eigenv

    implicit none
    class(bfgs),       intent(inout) :: this
    integer,           intent(in   ) :: lun_err
    integer,           intent(in   ) :: nequ
     !local 
    logical :: rc
    integer :: res
    
    !
    ! set dimension and symmetry
    !
    this%nrow         = nequ
    this%ncol         = nequ

    allocate ( &
         this%scr(nequ),&
         this%u_scr(nequ),&
         this%y_scr(nequ),&
         this%w_scr(nequ),&
         stat = res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'init_bfgs', &
         '  type bfgs scratch members',res)

  end subroutine init_bfgs

    !>-------------------------------------------------------------
  !> Static constructor.
  !> (procedure public for type input_tuning)
  !> Instantiate and initilize variable of type input_tuning
  !>
  !> usage:
  !>     call 'var'%init(ieig,nev,ituning)
  !>
  !> where:
  !> \param[in] lun_err       -> integer. I/O err. msg. logical unit 
  !<-------------------------------------------------------------
  subroutine kill_bfgs(this&
       ,lun_err&
       )
    use Globals
    implicit none
    class(bfgs),       intent(inout) :: this
    integer,           intent(in   ) :: lun_err
    !local 
    logical :: rc
    integer :: res

    deallocate ( &
         this%scr,&
         this%u_scr,&
         this%y_scr,&
         this%w_scr,&
         stat = res)
    if(res .ne. 0) rc = IOerr(lun_err, err_dealloc, 'init_bfgs', &
         '  type bfgs scratch members',res)

  end subroutine kill_bfgs

   !>-------------------------------------------------------------
  !> Static constructor.
  !> (procedure public for type input_tuning)
  !> Instantiate and initilize variable of type input_tuning
  !>
  !> usage:
  !>     call 'var'%init(ieig,nev,ituning)
  !>
  !> where:
  !> \param[in] nev                  -> integer. I/O err. msg. logical unit 
  !> \param[in] ctrl_tuning           -> type(input_tuning). Controls varaibles
  !>                                     for preconditioner application
  !> \param[in] (optional) ctrl_dacg -> type(input_dacg) (used only if ieig=1)
  !>                                      Controls for DACG Procedure
  !> \param[in] (optional) npres     -> integer. (Used only if ieig=2)
  !>                                      Number of prec. residual stores  
  !<-------------------------------------------------------------
  subroutine set_bfgs(this,&
       W,H,AV,prec_zero&
       )
    implicit none
    class(bfgs),           intent(inout) :: this
    type(densemat), target, intent(in   ) :: W
    class(abs_linop), target, intent(in   ) :: H
    class(densemat),   target, intent(in   ) :: AV
    class(abs_linop), target, intent(in   ) :: prec_zero
     !local 
    logical :: rc
    integer :: res
    
    this%is_symmetric  = prec_zero%is_symmetric

    this%W => W
    this%H => H
    this%AV => AV
    this%prec_zero =>  prec_zero

  end subroutine set_bfgs
  
  
  !>----------------------------------------------------------
  !> Precedure to apply preconditioner
  !> (procedure public for type bfgs)
  !> Prodecure computing the vector w=result of        
  !>              w=P^{-1}v
  !> for a given vector v (arrays' dimension must match)
  !>
  !> usage: call 'var'Mxv%(vec_in,vec_out,[info])
  !> 
  !> \param[in   ] vec_in  -> real, dimension(this%nequ)
  !>                            Vector v where bfgs will be applied
  !> \param[inout] vec_out  -> real, dimension(this%nequ)
  !>                            Vector w=P^{-1}v 
  !> \param[inout] (optional) info -> integer. Flag for errors  
  !>---------------------------------------------------------------------
  subroutine apply_bfgs(this, vec_in, vec_out,info,lun_err) 
    use Globals
    implicit none
    !var
    class(bfgs),   intent(inout) :: this
    real(kind=double), intent(in   ) :: vec_in(this%nrow)     
    real(kind=double), intent(inout) :: vec_out(this%ncol)
    integer,           intent(inout) :: info
    integer,           intent(in   ) :: lun_err

    ! local 
    integer :: nequ
    integer :: i
    integer :: info_loc

    !function z = bfgs(r,H,W,av,L)
    !u = H*(W'*r);    
    !y = r - av*u;
    !w = L'\(L\y);
    !y = H*(av'*w);
    !z = W*(u - y) + w;

    info_loc = 0
    nequ = this%nrow
    !u = H*(W'*r);
    call this%W%MTxv(vec_in,this%scr,i,lun_err)
    info_loc=info_loc+i
    call this%H%Mxv(this%scr, this%u_scr,i,lun_err)
    info_loc=info_loc+i
    
    ! y = r - av*u;
    call this%Av%Mxv(this%u_scr,this%y_scr,i,lun_err)
    info_loc=info_loc+i
    this%y_scr = vec_in - this%y_scr

    !w = L'\(L\y);
    call this%prec_zero%Mxv(this%y_scr,this%w_scr,i,lun_err)
    info_loc=info_loc+i
    
    !y = H*(av'*w);
    call this%AV%MTxv(this%w_scr, this%scr,i,lun_err)
    info_loc=info_loc+i
    call this%H%Mxv(this%scr, this%y_scr,i,lun_err)
    info_loc=info_loc+i
    
    !z = W*(u - y) + w;
    this%y_scr = this%u_scr - this%y_scr
    call this%W%Mxv(this%y_scr,vec_out,i,lun_err)
    info_loc=info_loc+i
    vec_out = vec_out + this%w_scr
        
    

    info = info_loc
          
  end subroutine apply_bfgs
  

end module PreconditionerTuning
