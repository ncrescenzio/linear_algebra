module SaddlePointMatrix
  use Globals
  use Matrix
  use LinearOperator
  use SimpleMatrix
  use SparseMatrix
  use LinearSolver
  use StdSparsePrec
  implicit none
  private
  !>-----------------------------------------------------------------------
  !> Structure variable for preconditioning the Saddle point Matrix
  !> M= ( A   B1^T )
  !>    ( B2 -C    )
  !> with A     := Real matrix of dimension n x n 
  !>      B1,B2 := Real matrix of dimension m x n (constraint matrix) 
  !>      C     := Real matrix of dimension m x m
  !> We follow the notation of 
  !> BMGHLJ05: Benzi, M. Golub H. and Liesen J.
  !> "Numerical solution of saddle point problems"
  !> Acta Numerica(2005)
  !>
  !> See. procedure init_saddle_point of type block_linop
  !>   
  !> BMWZ11: Benzi, M. and Wang Z.
  !> "Analysis of Augmented Lagrangian-based Preconditioners for
  !> the steady incompressible Navier-Stokes equations"
  !> SIAM Journal on Scientific Computing (2011)
  !>
  !> Most of procedure for type saddlemat are shorthands
  !> for block linear operator
  !>------------------------------------------------------------------------

  public :: augmented_saddle_point
  !>-----------------------------------------------------------------------
  !> Structure variable for application of precondiontioner for
  !> saddle point Matrix,
  !> M= ( A   B1^T )
  !>    ( B2 -C    )
  !> with A     := Real matrix of dimension n x n 
  !>      B1,B2 := Real matrix of dimension m x n (constraint matrix) 
  !>      C     := Real matrix of dimension m x m
  !> 
  !> We follow the notation and the presentation in 
  !>
  !> "Numerical solution of saddle point problems"
  !> Benzi, M. Golub H. and Liesen J.
  !> 
  !> in particular to section . First, we define the Schur complements
  !>
  !> SchurAC =    A + B1^T C^{-1} B2  )
  !> SchurCA = - (C + B2   A^{-1} B1^T)
  !>
  !> The preconditioner are bases on the following factorization of matrix M
  !>
  !> F1 :      M = ( I -B1^T C^{-1} ) ( SchurAC       0       ) ( I           0             )
  !>               ( 0  I           ) ( 0            -C       ) ( -C^{-1} B2  I             )
  !> with inverse
  !> invM1=M^{-1}= ( I         0    ) ( SchurAC^{-1}   0      ) ( I           B1^{T} C^{-1} )
  !>               ( C^{-1} B2 I    ) ( 0            -C^{-1}) ( ( 0           I             )
  !> the factorization
  !> F2=       M = ( I          0   ) ( A      0              ) ( I           A^{-1} B1^T   )
  !>               ( B2 A^{-1}  I   ) ( 0      SchurCA        ) ( 0           I             )
  !> with inverse
  !> invM2=M^{-1}= ( I -A^{-1} B1^T ) ( A^{-1} 0              ) ( I           0             )
  !>               ( 0 I            ) ( 0      SchurCA^{-1}   ) ( -B2 A^{-1}  I             )
  !>
  !> see equation (3.2 3.3)
  !> Preconditioner base on inverse of SchurAC and inverse of C 
  !> 'SAC_Full' = invM1                    
  !> 'SAC_Diag' = ( SchuAC 0  )^{-1}  
  !>              ( 0      -C )
  !> 'SAC_Full' is called constrained preconditioner
  !> 'SAC_Diag' is block diagonal  preconditioner
  !>
  !> Preconditioner base on inverse of SchurCA and inverse of A
  !> 
  !> 'SCA_Full'       = invM2
  !> 'SCA_Diag'       = ( A 0      )^{-1}
  !>                    ( 0 SchuCA )
  !> 'SCA_UpperTriang'= ( A  B1T    )^{-1}
  !>                    ( 0  SchuCA )
  !> 'SCA_LowerTriang'= ( A  0      )^{-1}
  !>                    ( B2 SchuCA )
  !> 'SCA_Full' is called constrained preconditioner
  !> 'SCA_Diag' is called block diagonal  preconditioner
  !> 'SAC_UpperTriang' is called (uppper) triangular constrained preconditioner
  !> (see BMGHLJ05)
  !> 'SAC_LowerTriang' is called (lower ) triangular constrained preconditioner
  !> (see BMWZ11)
  !>------------------------------------------------------------------------
  type, extends(abs_linop), public :: saddleprec
     !> Preconditoner type description
     !> It will take values
     !> 'SAC_Full', 'SAC_Diag' for precondtioner based on the
     !> inverse of SAC and C
     !> 'SCA_Full', 'SCA_Diag', 'SCA_UpperTriang', 'SCA_LowerTriang',
     !> for precondtioner based on the inverse of A and SCA
     character(len=20) :: prec_type='empty'
     !> Pointer to the saddle point matrix
     !> for which we want to approximate the inverse
     class(block_linop), pointer :: matrixM => null()
     !> Linear operator acting as the inverse 
     !> of the Schur complemte = A + B1^{T} C^{-1} B2.
     class(abs_linop), pointer :: invSchurAC => null()
     !> Linear operator acting as the inverse of A 
     class(abs_linop), pointer :: invC => null()
     !> Linear operator acting as  the inverse
     !> of the Schur complement = A + B1^{T} C^{-1} B2.
     class(abs_linop), pointer :: invSchurCA => null()
     !> Linear operator acting as the inverse of C
     class(abs_linop), pointer :: invA => null()
     !> Dimension(nrow)
     !> Scratch array for Matrix Vector application
     real(kind=double), allocatable :: scr_nrow(:)
     !> Dimension(ncol)
     !> Scratch array for Matrix Vector application
     real(kind=double), allocatable :: scr_ncol(:)
   contains
     !> static constructor
     !> (procedure public for type saddleprec)
     procedure, public, pass :: init => init_saddleprec
     !> static destructor
     !> (procedure public for type saddleprec)
     procedure, public, pass :: kill => kill_saddleprec
     !> Procedure to compute 
     !>         y = P^{-1} * x 
     !> (public procedure for type diag_prec)
     procedure, public,  pass :: matrix_times_vector => apply_saddleprec
  end type saddleprec
  
contains
  !>-------------------------------------------------------------
  !> Static constructor.
  !> (procedure public for type stiffprec)
  !> Instantiate and initiliaze variable of type saddleprec
  !>
  !> usage:
  !>     call 'var'%init(lun_err,matrix_B, prec_E1, prec_shur)
  !> where:
  !> \param[in] lun_err  -> integer. Error logical unit
  !> \param[in] matrix_B -> class(abs_matrix) Constraint matrix
  !> \param[in] prec_E1   -> class(prec_schur) Constraint matrix    
  !<-------------------------------------------------------------
  subroutine init_saddleprec(this,&
       lun_err,&
       matrixM,&
       prec_type,&
       invA,invSchurCA,&
       invSchurAC,invC)
    use Globals
    implicit none
    class(saddleprec), target, intent(inout) :: this
    integer,                  intent(in   ) :: lun_err
    class(block_linop), target, intent(in   ) :: matrixM
    character(len=*),         intent(in   ) :: prec_type
    class(abs_linop), target, optional, intent(in   ) :: invA
    class(abs_linop), target, optional, intent(in   ) :: invSchurCA
    class(abs_linop), target, optional, intent(in   ) :: invSchurAC
    class(abs_linop), target, optional, intent(in   ) :: invC
    !local
    logical :: rc
    integer :: res
  

    !
    ! free memory
    !
    if ( allocated(this%scr_nrow)) call this%kill(lun_err)

    this%prec_type = etb(prec_type)
    this%matrixM => matrixM

    if (.not. matrixM%is_saddle_point ) then
       rc = IOerr(lun_err, err_val,&
            ' init_saddleprec', &
            ' matrix M is not a saddle point matrix')
    end if

    
    if( &
         ( this%prec_type .eq. 'SchurACFull') .or. &
         ( this%prec_type .eq. 'SchurACDiag') .or. &
         ( this%prec_type .eq. 'SchurACTriang') &
         )then

       if ( present ( invSchurAC ) )  then
          ! check dimensions
          if ( ( invSchurAC%nrow .le. 0) .or. ( invSchurAC%nrow .le. 0)) then
             call invSchurAC%info(lun_err)
             rc = IOerr(lun_err, err_val,&
                  ' init_saddleprec', &
                  ' prec '//etb(prec_type)//&
                  ' requires invSchurSAC has zero dimension'  ,res)
          end if
          if ( ( invSchurAC%nrow .ne. this%matrixM%linop_list(1)%linop%nrow) .or. &
               ( invSchurAC%ncol .ne. this%matrixM%linop_list(1)%linop%ncol)) then
             call invSchurAC%info(lun_err)
             call this%matrixM%linop_list(1)%linop%info(lun_err)
             rc = IOerr(lun_err, err_val,&
                  ' init_saddleprec', &
                  ' prec '//etb(prec_type)//&
                  ' requires invSchurAC does not match with matrix A'  ,res)
          end if
          ! assign pointer
          this%invSchurAC => invSchurAC
       else
          rc = IOerr(lun_err, err_val,&
               ' init_saddleprec', &
               ' prec '//etb(prec_type)//&
               ' requires invSchurSAC that was not passed'  ,res)
       end if

       if ( present ( invC ) )  then
          ! check dimensions
          if ( ( invC%nrow .ne. this%matrixM%linop_list(3)%linop%nrow) .or. &
               ( invC%ncol .ne. this%matrixM%linop_list(2)%linop%ncol)) then
             call invC%info(lun_err)
             rc = IOerr(lun_err, err_val,&
                  ' init_saddleprec', &
                  ' using approach '//etb(prec_type)//&
                  ' invC does not match with matrix C')
          end if
          ! assign pointer
          this%invC => invC
       else
          rc = IOerr(lun_err, err_val,&
               ' init_saddleprec', &
               ' prec '//etb(prec_type)//&
               ' requires invC that was not passed'  ,res)

       end if

       !
       ! same dimension
       !
       this%nrow=matrixM%nrow
       this%ncol=matrixM%ncol
       
       select case (this%prec_type)
       case ('SchurACFull')
          this%is_symmetric = &
               invC%is_symmetric      .and. &
               invSchurAC%is_symmetric .and. &
               matrixM%B1equalB2 
       case ('SchurACDiag')
          ! same dimension 
          call this%is_like(matrixM)
          this%is_symmetric = &
               invC%is_symmetric      .and. &
               invSchurAC%is_symmetric
       case ('SchurACTriang') 
          this%is_symmetric = .False.
          this%triangular   = 'N'
       end select
          
       
    else if ( &
         ( this%prec_type.eq. 'SchurCAFull') .or. &
         ( this%prec_type.eq. 'SchurCADiag') .or. &
         ( this%prec_type.eq. 'SchurCAUpperTriang') .or. &
         ( this%prec_type.eq. 'SchurCALowerTriang') &
         ) then
       !
       ! check presence, otherwise return error
       !
       if ( present ( invSchurCA ) )  then
          !
          ! check size  
          !
          if ( ( invSchurCA%nrow .ne. this%matrixM%linop_list(3)%linop%nrow) .or. & ! B2 
               ( invSchurCA%ncol .ne. this%matrixM%linop_list(2)%linop%ncol)) then  ! B1T
             call invSchurCA%info(lun_err)
             call this%matrixM%info(lun_err)
             rc = IOerr(lun_err, err_val,&
                  ' init_saddleprec', &
                  ' prec '//etb(prec_type)//&
                  ' requires invSchurSAC does not match with matrix C'  ,res)
          end if
          this%invSchurCA => invSchurCA
       else
          rc = IOerr(lun_err, err_val,&
               ' init_saddleprec', &
               ' prec '//etb(prec_type)//&
               ' requires invSchurSAC that was not passed'  ,res)
       end if
       !
       ! check present of A inverse, otherwise return error
       !
       if ( present ( invA ) )  then
          !
          ! check dimension
          !
          if ( ( invA%nrow .ne. this%matrixM%linop_list(1)%linop%nrow) .or. &
               ( invA%ncol .ne. this%matrixM%linop_list(1)%linop%ncol)) then
             call invA%info(lun_err)
             call this%matrixM%linop_list(1)%linop%info(lun_err)
             rc = IOerr(lun_err, err_val,&
                  ' init_saddleprec', &
                  ' prec '//etb(prec_type)//&
                  ' requires invA does not match with matrix A'  ,res)
          end if
          this%invA => invA
       else
          rc = IOerr(lun_err, err_val,&
               ' init_saddleprec', &
               ' prec '//etb(prec_type)//&
               ' requires invC that was not passed'  ,res)

       end if

       !
       ! same dimension
       !
       this%nrow=matrixM%nrow
       this%ncol=matrixM%ncol

       select case (this%prec_type)
       case ('SchurCAfull')
          ! same symmetry, same symmetry
          call this%is_like(matrixM)
          
       case ('SchurCAdiag')
          this%is_symmetric = &
               invA%is_symmetric .and. &
               invSchurCA%is_symmetric
          this%triangular='N'
          if ( (invA%triangular .eq. 'L')  .and. &
               (invA%triangular .eq. 'L') ) then
             this%triangular='L'
          end if
          if ( (invA%triangular .eq. 'U')  .and. &
               (invA%triangular .eq. 'U') ) then
             this%triangular='U'
          end if

       case ('SchurCAtriang') 
          this%is_symmetric = .False.
          this%triangular   = 'N'
       case ('SchurCALowerTriang')
          this%is_symmetric = .False.
          this%triangular   = 'N'
       case ('SchurCAUpperTriang')
          this%is_symmetric = .False.
          this%triangular   = 'N' 
       end select
       
    else
       rc = IOerr(lun_err, err_val,&
            ' init_saddleprec', &
            ' prec '//etb(prec_type)//' not supported'  ,res)
    end if

    !
    ! allocate space 
    !
    allocate(&
         this%scr_nrow(this%nrow),&
         this%scr_ncol(this%ncol),&
         stat=res)
    if ( res .ne. 0) &
         rc = IOerr(lun_err, err_alloc, 'init_saddleprec', &
         ' type saddlemat member linops, nrow_vec, ncol_vec'//&
         ' block_structure scrt_nrow, scr_ncol')



     
  end subroutine init_saddleprec

  !>-------------------------------------------------------------
  !> Static constructor.
  !> (procedure public for type stiffprec)
  !> Instantiate and initiliaze variable of type saddleprec
  !>
  !> usage:
  !>     call 'var'%init(lun_err,matrix_B, prec_E1, prec_shur)
  !> where:
  !> \param[in] lun_err  -> integer. Error logical unit
  !> \param[in] matrix_B -> class(abs_matrix) Constraint matrix
  !> \param[in] prec_E1   -> class(prec_schur) Constraint matrix  
  
  !<-------------------------------------------------------------
  subroutine kill_saddleprec(this, lun_err)
    use Globals
    implicit none
    class(saddleprec), intent(inout) :: this
    integer,           intent(in   ) :: lun_err
    !local
    logical :: rc
    integer :: res

    ! nullfy pointers
    this%invSchurAC=> null()
    this%invSchurCA=> null()
    this%invC=> null()
    this%invA=> null()
    this%matrixM=> null()
    ! set to defalt

    call this%to_default()
    this%prec_type ='empty'

    !
    ! deallocate space 
    !
    deallocate(&
         this%scr_nrow,&
         this%scr_ncol,&
         stat=res)
    if ( res .ne. 0) &
         rc = IOerr(lun_err, err_dealloc, 'init_saddlemat', &
         ' type saddlemat member linops, nrow_vec, ncol_vec'//&
         ' block_structure scr_nrow, scr_ncol')
    
  end subroutine kill_saddleprec

 !>-------------------------------------------------------------
  !> procedure defining the interface for a general
  !> matrix-vector multiplication
  !>         vec_out = (M) times (vec_in)
  !> (public procedure for class abs_matrix)
  !> 
  !> usage:
  !>     call 'var'%matrix_times_vector(vec_in,vec_out,[info])
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
  subroutine apply_saddleprec(this,vec_in,vec_out,info,lun_err)
    use Globals

    implicit none
    class(saddleprec), intent(inout) :: this
    real(kind=double), intent(in   ) :: vec_in(this%ncol)
    real(kind=double), intent(inout) :: vec_out(this%nrow)
    integer,           intent(inout) :: info
    integer,           intent(in   ) :: lun_err
    ! local
    integer :: n,m,nm
    integer :: info_loc
    logical :: debug=.False.

    info_loc=0
    info=0
    
    m=this%matrixM%linop_list(3)%linop%nrow
    n=this%matrixM%linop_list(3)%linop%ncol
    nm=m+n

    !
    ! we use the scratch vectors in block_linear_operator
    !
    select case (this%prec_type)
    case ('SchurACFull')
       !
       ! use factorization F1
       !
       if (debug) write(6,*) 'SchurACFull'   
       
       ! compute t=C^{-1} y
       ! stored in vec_out(n+1:nm)
       call this%invC%matrix_times_vector(&
            vec_in(n+1:nm),&
            this%scr_nrow(n+1:nm),&
            info=info,lun_err=lun_err)
       if (debug) write(6,*) 'C^{-1} y, info=', info
       if (info.ne.0) return
       
       ! compute v= B1T C^{-1} y
       ! stored in this%matrixM%scr_ncol(1:n)
       call this%matrixM%linop_list(2)%linop%matrix_times_vector(&
            this%scr_nrow(n+1:nm),&
            this%scr_ncol(1:n),&
            info=info,lun_err=lun_err)
       if (debug) write(6,*) 'B1T^{-1} y, info=', info
       if (info.ne.0) return
       
       ! compute x+v=x + B1^T C^{-1} y
       !this%scr_ncol(1:n)=vec_in(1:n)+this%scr_ncol(1:n)
       call daxpy(n,one,&
            vec_in(1:n),1,&
            this%scr_ncol(1:n),1)
       
       
       ! compute vec_out(1:n) =SchurAC^{-1} ( x + B1^T C^{-1} y)
       call this%invSchurAC%matrix_times_vector(&
            this%scr_ncol(1:n),&
            vec_out(1:n),&
            info=info,lun_err=lun_err)
       if (debug) write(6,*) 'SchurAC^{-1} ( x + B1^T C^{-1} y), info=', info
       if (info.ne.0) return
          

       ! compute vec_out(1+n:nm) = C^{-1} B2 vec_out(1:n) - C^{-1} y
       ! compute B2 * vec_out(1:n) 
       call this%matrixM%linop_list(3)%linop%matrix_times_vector(&
            vec_out(1:n),&
            this%scr_ncol(1:m),&
            info=info,lun_err=lun_err)
       if (debug) write(6,*) 'C^{-1} B2 vec_out(1:n) - C^{-1} y, info=', info
       if (info.ne.0) return


       ! compute C^{-1} B2 * vec_out(1:n) 
       call this%invC%matrix_times_vector(&
            this%scr_ncol(1:m),&
            vec_out(n+1:nm),&
            info=info,lun_err=lun_err)
       if (debug) write(6,*) 'C^{-1} B2 * vec_out(1:n), info=', info
       if (info.ne.0) return
       
       !vec_out(1+n:nm)=vec_out(1+n:nm)-this%scr_nrow(n+1:nm)
       call daxpy(m,-one,&
            this%scr_nrow(n+1:nm),1,&
            vec_out(1+n:nm),1)

       
    case ('SchurACDiag')
       call this%invSchurAC%matrix_times_vector(&
            vec_in(1:n),&
            vec_out(1:n),&
            info=info,lun_err=lun_err)
       if (info.ne.0) return

       
       call this%invC%matrix_times_vector(&
            vec_in(n+1:nm),&
            vec_out(n+1:nm),&
            info=info,lun_err=lun_err)      
       if (info.ne.0) return

       
    case ('SchurCAFull')
       ! vec_in  = (f,g)
       ! vec_out = (x,y)
       ! P vec_out = vec_in
       !
       ! compute t = A^{-1} f
       ! stored in vec_out(1:n)
       call this%invA%matrix_times_vector(&
            vec_in(1:n),&
            this%scr_nrow(1:n),& ! t (do not use it we need it later)
            info=info,lun_err=lun_err)
       if (info.ne.0) return
       

       ! compute v= B2 A^{-1} f = B2 t
       call this%matrixM%linop_list(3)%linop%matrix_times_vector(&
            this%scr_nrow(1:n),&    ! t
            this%scr_nrow(1+n:nm),& ! v
            info=info,lun_err=lun_err)
       if (info.ne.0) return

       ! compute w=-B2 A^{-1}f + g = -v + g
       this%scr_ncol(n+1:nm) = &
            - this%scr_nrow(1+n:nm) + vec_in(n+1:nm) 

       ! compute y = SchurCA^{-1} ( -B2^T A^{-1}f + g)
       !           = SchurCA^{-1} w
       call this%invSchurCA%matrix_times_vector(&
            this%scr_ncol(n+1:nm),& ! w
            vec_out(n+1:nm),&       ! y
            info=info,lun_err=lun_err)
       if (info.ne.0) return

       
       ! compute w = B1T * y
       call this%matrixM%linop_list(2)%linop%matrix_times_vector(&
            vec_out(n+1:nm),&    ! y
            this%scr_ncol(1:n),& ! w (overwritten)
            info=info,lun_err=lun_err)
       if (info.ne.0) return

       
       ! compute x= A^{-1} B1T y= A^{-1} w
       call this%invC%matrix_times_vector(&
            this%scr_ncol(1:m),&
            vec_out(1:n),&
            info=info_loc,lun_err=lun_err)
       if (info.ne.0) return

       ! x=t-x=A^{-1} x - A^{-1} B1^{T} y
       vec_out(1:n) = this%scr_nrow(1:n) - vec_out(1:n) 
      
    case ('SchurCAdiag')
       ! BMGHLJ05 eq 10.1
       call this%invA%matrix_times_vector(&
            vec_in(1:n),&
            vec_out(1:n),&
            info=info,lun_err=lun_err)
       if (info.ne.0) return

       call this%invSchurCA%matrix_times_vector(&
            vec_in(n+1:nm),&
            vec_out(n+1:nm),&
            info=info,lun_err=lun_err)
       if (info.ne.0) return

       vec_out(n+1:nm)=-vec_out(n+1:nm)
     
    case ('SchurCAUpperTriang')
       ! BMGHLJ05 eq 10.11

       ! solve linear system
       ! ( A B1T) x = f
       ! ( 0 S  ) y   g
       !
       ! x=A^{-1} (f - S^{-1} g ) 
       ! y=S^{-1} g
       
       ! compute y = SchurCA^{-1} g
       call this%invSchurCA%matrix_times_vector(&
            vec_in(1+n:nm),&
            vec_out(1+n:nm),&
            info=info,lun_err=lun_err)
       if (info.ne.0) return
              
       
       ! compute v= B1 SchurCA^{-1} g
       ! stored in this%matrixM%scr(1:n)
       call this%matrixM%linop_list(2)%linop%matrix_times_vector(&
            vec_out(1+n:nm),&
            this%scr_ncol(1:n),&
            info=info,lun_err=lun_err)
       if (info.ne.0) return

       
       ! compute f+v=f-B1 SchurCA^{-1} g
       ! stroed in v
       this%scr_ncol(1+n:nm) = &
            vec_in(1:n) - &
            this%scr_ncol(1:n)

       ! compute vec_out(1:n) = SchurAC^{-1} ( g + B2 A{-1} f)
       call this%invA%matrix_times_vector(&
            this%scr_ncol(1:n),&
            vec_out(1:n),&
            info=info,lun_err=lun_err)
       if (info.ne.0) return
       
       
    case ('SchurCALowerTriang')
     
       ! equation 3.3

       ! solve linear system
       ! ( A  0) x = f
       ! ( B2 S) y   g
       !
       ! x=A^{-1} f
       ! y=S^{-1}(g-B2 A^{-1} f) 
       
       ! compute x = A^{-1} f
       call this%invA%matrix_times_vector(&
            vec_in(1:n),&
            vec_out(1:n),&
            info=info,lun_err=lun_err)
       if (info.ne.0) return

      
       ! compute v= B2 A^{-1} f
       ! stored in this%matrixM%scr(1:n)
       call this%matrixM%linop_list(3)%linop%matrix_times_vector(&
            vec_out(1:n),&
            this%scr_ncol(1+n:nm),&
            info=info,lun_err=lun_err)
       if (info.ne.0) return
       
       ! compute g+v=g-B2 A^{-1} f
       ! stroed in v
       this%scr_ncol(1+n:nm) = &
            vec_in(1+n:nm) - &
            this%scr_ncol(1+n:nm)

       ! compute vec_out(1:n) = SchurCA^{-1} ( g + B2 A{-1} f)
       call this%invSchurCA%matrix_times_vector(&
            this%scr_ncol(1+n:nm),&
            vec_out(n+1:nm),&
            info=info,lun_err=lun_err)
       if (info.ne.0) return
       
    end select

       


  end subroutine apply_saddleprec

  !>-------------------------------------------------------------
  !> Build Augmented Saddle Point component A_gamma B1T_gamma of
  !> M_gamma=( A_gamma B1T_gamma) = ( A+gammaB1TinvWCB2 B1T(I-gammainvWC) )
  !>         ( B2      -C       )   ( B2                -C                )
  !>
  !> Note that M_gamma=P_gamma * M
  !> with 
  !> P_gamma = ( I gamma B1T inv WC) 
  !>           ( 0 I               ) 
  !> 
  !> (procedure public for type stiffprec)
  !> Instantiate and initiliaze variable of type saddleprec
  !>
  !> usage:
  !>     call 'var'%init(lun_err,)
  !> where:
  !> \param[in] lun_err  -> integer. Error logical unit
  !> \param[in] matrix_B -> class(abs_matrix) Constraint matrix
  !> \param[in] prec_E1   -> class(prec_schur) Constraint matrix    
  !<-------------------------------------------------------------
  subroutine augmented_saddle_point(&
       lun_err,&
       gamma, invW,&
       matrixM,&
       matrixAgamma,matrixB1Tgamma,&
       matrixB1TinvW, augmentation_operator,eye_n,eye_m)
    use Globals
    implicit none
    integer,           intent(in   ) :: lun_err
    real(kind=double), intent(in   ) :: gamma
    class(abs_linop),  target,intent(in   ) :: invW
    class(block_linop),target,  intent(in   ) :: matrixM
    type(new_linop), target, intent(inout) :: matrixAgamma
    type(new_linop), target, intent(inout) :: matrixB1Tgamma
    type(new_linop), target,   intent(inout) :: matrixB1TinvW
    type(block_linop), intent(inout) :: augmentation_operator
    type(scalmat), target, intent(inout) :: eye_n
    type(scalmat), target, intent(inout) :: eye_m
    ! local
    logical :: rc,debug=.True.
    integer :: n,m,nm,i
    type(array_linop) :: list(3)
    real(kind=double) :: alphas(3)
    integer :: block_structure(3,3)

    if ( .not. matrixM%is_saddle_point ) stop
    
    !
    ! build Agamma=A+gamma B1T W^{-1} B2
    !
    matrixAgamma = matrixM%linop_list(1)%linop + &
         gamma*&
         matrixM%linop_list(2)%linop *&
         invW * &
         matrixM%linop_list(3)%linop    
    matrixAgamma%name=etb(etb(matrixM%linop_list(1)%linop%name)//'_gamma')

    if ( debug ) write(lun_err,*) 'Agamma'


    !
    ! build Bgamma=B1T+gamma B1T W^{-1} C=B1T
    !
    if ( matrixM%nlinops .eq. 4) then    
       matrixB1Tgamma=matrixM%linop_list(2)%linop - &
            gamma * matrixM%linop_list(2)%linop*invW*matrixM%linop_list(4)%linop
       
       matrixB1Tgamma%name=etb(etb(matrixM%linop_list(2)%linop%name)//'_gamma')
       if ( debug ) write(lun_err,*) 'done matrixB1Tgamma'
    else      
       matrixB1Tgamma=one*matrixM%linop_list(2)%linop
       if ( debug ) write(lun_err,*) 'done matrixB1T_gamma'
    end if

    matrixB1TinvW=matrixM%linop_list(2)%linop*invW
    
    call eye_n%eye(matrixB1TinvW%nrow)
    call eye_m%eye(matrixB1TinvW%ncol)


    block_structure(:,1) = (/1,1,1/)
    list(1)%linop => eye_n
    alphas(1)=one

    block_structure(:,2) = (/2,1,2/)
    list(2)%linop=> matrixB1TinvW
    alphas(2)=gamma

    block_structure(:,3) = (/3,2,2/)
    list(3)%linop=> eye_m
    alphas(3)=one


    call augmentation_operator%init(lun_err, &
         3, list,&
         2, 2, &
         3, block_structure,&
         .False.,alphas=alphas)
       

       
  end subroutine augmented_saddle_point


  !>--------------------------------------------
  !> Least Square Preconditioner from
  !> See "BLOCK PRECONDITIONERS BASED ON APPROXIMATE COMMUTATORS"
  !> (HOWARD ELMAN, VICTORIA E. HOWLE, JOHN SHADID,
  !> ROBERT SHUTTLEWORTH, AND RAY TUMINARO)
  !> for theory and notation.
  !> Implementation of the simplest version (D=G^T) and
  !> no weights. The resulting linear operatror
  !> aprroximates the inverse of  minus SchurCA = DF^{-1}G
  !<--------------------------------------------
  subroutine build_least_square_preconditioner(&
       lun_err,&
       DequalGT,&
       matrixGrad,matrixDiv,matrixF,&
       inverse_matrixDG,inverse_matrixGTG,&
       least_square_preconditioner) 
    implicit none
    integer,                intent(in   ) :: lun_err
    logical,                intent(in   ) :: DequalGT
    class(abs_linop),target,intent(in   ) :: matrixGrad
    class(abs_linop),target,intent(in   ) :: matrixDiv
    class(abs_linop),target,intent(in   ) :: matrixF
    class(abs_linop),target,intent(in   ) :: inverse_matrixDG
    class(abs_linop),target,intent(in   ) :: inverse_matrixGTG
    type(pair_linop),    intent(inout) :: least_square_preconditioner
    !local
    integer :: info
    type(array_linop) :: linops(5)

    linops(1)%linop => inverse_matrixDG
    linops(2)%linop => matrixDiv
    linops(3)%linop => matrixF
    linops(4)%linop => matrixGrad
    linops(5)%linop => inverse_matrixGTG

    call least_square_preconditioner%init(info,lun_err, &
         5, linops,'LP')

    least_square_preconditioner%is_symmetric = DequalGT
    

  end subroutine build_least_square_preconditioner

  
  
  subroutine solve_saddle_point(&
       matrix_M,&
       rhs, & ! rhs
       sol, & ! solution and initial guess
       lun_err,info,&
       prec_combination,prec_type,&
       invA,invSchurCA,invSchurAC,invC,&
       ctrl_solver, &
       info_solver,&
       aux)
    use LinearOperator
    use Scratch
    implicit none
    type(block_linop),target,           intent(inout) :: matrix_M
    real(kind=double),                  intent(in   ) :: rhs(matrix_M%nrow)
    real(kind=double),                  intent(inout) :: sol(matrix_M%ncol)
    integer, intent(in) :: lun_err
    integer,                            intent(inout) :: info
    character(len=*),                   intent(in   ) :: prec_combination
    character(len=*),                   intent(in   ) :: prec_type
    class(abs_linop), target, optional, intent(in   ) :: invA
    class(abs_linop), target, optional, intent(in   ) :: invSchurCA
    class(abs_linop), target, optional, intent(in   ) :: invSchurAC
    class(abs_linop), target, optional, intent(in   ) :: invC    
    type(input_solver), optional,       intent(inout) :: ctrl_solver
    type(output_solver), optional,       intent(inout) :: info_solver
    type(scrt), optional,target,        intent(inout) :: aux

    ! derivced type
    type(scrt), target  :: aux_work
    type(scrt), pointer :: aux_loc

    ! local
    logical :: rc,schur_AD_is_symmetric
    integer :: res
    integer :: i,j,ind
    integer :: info_prec
    integer :: ntdens, npot, nfull,n,m,nm
    type(saddleprec) :: preconditoner_M




    type(input_solver) :: ctrl_uzawa

    ! broyden updates quantities
    integer :: iter_broyden
    integer :: info_broyden
    integer :: iterm,inode

    type(file) :: fmat
    character(len=256) :: fname,directory,tail
    type(spmat),   pointer :: matrix2prec
    !
    real(kind=double) :: dnrm2,esnorm
    real(kind=double) :: max
    integer :: imax

    character(len=70) :: str
    character(len=256) :: msg,msg1,msg2



    !
    ! set work arrays
    !    
    if ( present(aux) ) then
       aux_loc => aux
    else
       call aux_work%init(lun_err,0,9*matrix_M%nrow)
       aux_loc => aux
    end if

    ! set no error flag
    !
    info = 0

    !
    ! set isol = 0 to compute only errors estimate of 
    ! | J x - b |
    select case (prec_combination)
    case ('PREC')

       call preconditoner_M%init(&
            lun_err,&
            matrix_M,&
            prec_type,&
            invA,invSchurCA,&
            invC,invSchurAC)

       call preconditoner_M%matrix_times_vector(rhs,sol,info,lun_err)

       if (present(info_solver)) then
          info_solver%ierr=info
       end if


    case ('PREC+SOLVER')
       call preconditoner_M%init(&
            lun_err,&
            matrix_M,&
            prec_type,&
            invA,invSchurCA,&
            invC,invSchurAC)

       call linear_solver(&
            matrix_M,rhs,sol,&
            info_solver, &
            ctrl_solver,&
            prec_left=preconditoner_M,&
            aux=aux)
    end select


  end subroutine solve_saddle_point


end module SaddlePointMatrix
