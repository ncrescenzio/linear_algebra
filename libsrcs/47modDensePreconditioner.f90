module DensePreconditioner
  use Globals
  use LinearOperator
  use Scratch
  use Matrix
  use StdSparsePrec
  use DenseMatrix
  use Eigenv
  use SimpleMatrix
  implicit none
  private
  !>------------------------------------------------------------------------
  !> Structure variable for the application 
  !> inverse of linear operator derived from 
  !> dense matrices.
  !> Supported option, that can be selected via input_prec
  !> type are:
  !> - identity (default)
  !> - CHOL (for dense symmetric matrices)
  !> - SVD      (for all matrices            )
  !>------------------------------------------------------------------------
  type, extends(abs_linop), public :: denseprec
     !> cheolesky
     character(len=20) :: prec_type
     !> Number of Factors 
     integer :: nfactors
     !> Full Matrix Full Matrix factorization
     type(densemat), allocatable :: factor(:)
     !> Number of diagonal matrices
     integer :: ndiags
     !> Full Matrix Full Matrix factorization
     type(diagmat), allocatable :: diagonals(:)
     !> Eigenvectors and eigenvalue of matrix 
     type(eigen) :: eig
     !> Pointer to original matrix
     class(abs_matrix), pointer :: matrix => null()
     !> Aux variables
     type(scrt) :: aux
   contains
     !> static constructor
     !> (procedure public for type denseprec)
     procedure, public, pass :: init => init_denseprec
      !> static constructor
     !> (procedure public for type denseprec)
     procedure, public, pass :: assembly => assembly_denseprec
     !> static destructor
     !> (procedure public for type denseprec)
     procedure, public, pass :: kill => kill_denseprec
     !> Procedure to compute 
     !>         y = P^{-1} * x 
     !> (public procedure for type diag_prec)
     procedure, public,  pass :: matrix_times_vector => Mxv_denseprec
  end type denseprec
  contains
    
    subroutine init_denseprec(this,lun_err,full_matrix,ctrl)
     implicit none
     class(denseprec),       intent(inout) :: this
     integer,               intent(in   ) :: lun_err
     type(densemat), target, intent(in   ) :: full_matrix
     type(input_prec),      intent(in   ) :: ctrl
     !local
     logical :: rc
     integer :: res,i,nequ
     integer :: n1,n2,n3,nb 
     integer :: ilaenv
     character(len=1) :: UPLO
     
     !
     ! check
     !
     if ( full_matrix%nrow .ne. full_matrix%ncol ) &
          rc = IOerr(lun_err, err_val, 'init_denseprec', &
          ' matrix not squared',res)


     ! set properties
     this%is_symmetric = full_matrix%is_symmetric
     this%nrow         = full_matrix%nrow
     this%ncol         = full_matrix%ncol
     nequ= this%nrow

     this%matrix => full_matrix

     select case (ctrl%prec_type) 
     case ('C')
        this%prec_type = 'C'
        this%nfactors = 1
        allocate(this%factor(1),stat=res)
        if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'init_denseprec', &
                ' type denseprec member factor',res)
        call this%factor(1)%init(lun_err, &
             nequ,nequ)
        
     case ('LU')
        this%prec_type = 'LU'
        this%nfactors = 1
        ! just one dense containg L and U
        allocate(this%factor(1),stat=res) 
        if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'init_denseprec', &
                ' type denseprec member factor',res)
        call this%factor(1)%init(lun_err, &
             nequ,nequ)

        ! init integer array for pivoting
        call this%aux%init(lun_err,nequ,0)
     case ('SVD')
        !
        ! select tpye nd action of preconditioner
        !
        this%prec_type = 'SVD'
        
        !
        ! init eigen variables and auxiliary  
        !
        call this%eig%init(lun_err, nequ, nequ)
        nb=ilaenv( 1, 'dsyev', UPLO, N1, N2, N3, nequ)
        call this%aux%init(lun_err, 0 , (nb+2)* nequ)
        
     end select
     
           

   end subroutine init_denseprec


   subroutine assembly_denseprec(this,lun_err,full_matrix,info)
     use Globals
     implicit none
     class(denseprec),       intent(inout) :: this
     integer,               intent(in   ) :: lun_err
     type(densemat),         intent(in   ) :: full_matrix
     integer,               intent(inout) :: info
          
     !local
     logical :: rc
     integer :: res,i,nequ
     character :: symmetric_storage
     
     nequ = this%nrow
     
     !
     ! check consistence or re-initi  
     !
!!$     if ( this%prec_type .ne. ctrl%prec_type) then
!!$        call this%kill(lun_err)
!!$        call this%init(lun_err,full_matrix,ctrl)
!!$     end if


     !
     ! FIXME checks consistency of data
     !
     if ( this%prec_type .eq. 'C' ) then
        this%factor(1)%coeff = full_matrix%coeff        
        CALL DPOTRF( 'U', nequ, this%factor(1)%coeff, nequ, info )
     end if

     if ( this%prec_type .eq. 'LU' ) then
        this%factor(1)%coeff = full_matrix%coeff
        call DGETRF( nequ, nequ, this%factor(1)%coeff, nequ, this%aux%iaux, &
             info )
     end if
     
     
     
     if ( this%prec_type .eq. 'SVD' ) then
        !
        ! copy matrix coeff.
        !
        this%eig%eigenvectors%coeff = full_matrix%coeff
        
        !
        ! set storage system 
        !
        symmetric_storage = full_matrix%symmetric_storage
        if (full_matrix%symmetric_storage .eq. 'F') then
           symmetric_storage='U'
        end if
        call full_matrix%spectral_decomposition(lun_err,&
             info,&
             this%eig%eigval,this%eig%eigenvectors%coeff,this%aux)
     end if


   end subroutine assembly_denseprec


   subroutine kill_denseprec(this,lun_err)
     use Globals
     implicit none
     class(denseprec),  intent(inout) :: this
     integer,          intent(in   ) :: lun_err
     !local
     logical :: rc
     integer :: res,i

     do i=1,this%nfactors
        call this%factor(i)%kill(lun_err)
     end do
     deallocate(this%factor,stat=res)
     if(res .ne. 0) rc = IOerr(lun_err, err_dealloc, 'kill_denseprec', &
          ' type denseprec member factor',res)
     
     if ( this%eig%is_initialized) call this%eig%kill(lun_err)
  
     if ( this%aux%is_initialized) call this%aux%kill(lun_err)

     
     this%is_symmetric   = .false.    

     this%matrix => null()
   end subroutine kill_denseprec


   subroutine Mxv_denseprec(this,vec_in,vec_out,info,lun_err)
     use Globals
     implicit none
     class(denseprec),  intent(inout) :: this
     real(kind=double), intent(in   ) :: vec_in(this%nrow)
     real(kind=double), intent(inout) :: vec_out(this%ncol)
     integer,           intent(inout) :: info
     integer,           intent(in   ) :: lun_err
     !local
     integer :: nequ,info_loc,i
     
     
     nequ= this%nrow
     if (this%prec_type .eq. 'C') then
        ! solve A x = b solving L L^t x = b
        vec_out = vec_in 
        CALL DPOTRS( 'U',&
             nequ, 1, this%factor(1)%coeff, nequ, vec_out, nequ, info_loc )
        
        info = info_loc
     end if

     if (this%prec_type .eq. 'LU') then
        ! solve A x = b solving L U x = b
        vec_out = vec_in 
        call DGETRS( 'N', nequ, 1, this%factor(1)%coeff, nequ,&
             this%aux%iaux, vec_out, nequ, info_loc)
        info = info_loc
     end if

     if (this%prec_type .eq. 'SVD') then
        vec_out = vec_in 
        call this%eig%eigenvectors%MTxv(vec_in,this%aux%raux(1:nequ),info,lun_err)
        do i=1,nequ
           this%aux%raux(i) = this%aux%raux(i) /this%eig%eigval(i)
        end do
        call this%eig%eigenvectors%Mxv(this%aux%raux(1:nequ),vec_out,info,lun_err)
     end if

     

   end subroutine Mxv_denseprec

 end module DensePreconditioner
