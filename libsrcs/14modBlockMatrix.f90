module BlockMatrix
  use Globals
  use Matrix
  implicit none
  private 
  !> Structure variable containing the variables and the procedures
  !> to handle a block matrices (e.g.)
  !> ( A   0 B  0)
  !> ( 0   C 0  D)
  !> ( B^T 0 D  0)
  !> and perform Matrix-vector and Matrix(transpose)-vector operations
  type, extends(abs_matrix), public :: blockmat
     !> Number of block matrices
     !> It is less or equal to the number of blocks
     !> In the example nmats = 4
     integer :: nmats
     !> Dimension (nmats)
     !> Array that contains the pointer to non-zero blocks
     !> In the example 
     !> mats(1)%mat => A
     !> mats(2)%mat => B
     !> mats(3)%mat => C
     !> mats(4)%mat => D
     type(array_mat), allocatable :: mats(:)
     !> Number of blocks in each row
     !> In the example nblock_row=3
     integer :: nblock_row
     !> Number of blocks in each column
     !> In the example nblock_col=4
     integer :: nblock_col
     !> Number of non-zero blocks
     !> In the example nnzblock=6
     integer :: nnzblock
     !> Dimension (3,nnzblock)
     !> It contains the list of the matrices used
     !> reading the Block Matrix line by line.
     !> Given i in 1, nnzblock
     !> - imat = block_structure(1,i)   
     !>   imat = index of the i^th-block in mats arrays
     !> - block_structure(2,i) = row    index of the i^th block 
     !> - block_structure(2,i) = column index of the i^th block 
     !> - block_structure(4,i) = 0-1 flag (0=not transpose 1=transpose)
     !> In the example
     !> block_structure(:,1)=(1,1,1)
     !> block_structure(:,2)=(2,1,3)
     !> block_structure(:,3)=(3,2,2)
     !> block_structure(:,4)=(4,2,4)
     !> block_structure(:,5)=(2,3,1)
     !> block_structure(:,6)=(4,3,3)
     integer, allocatable :: block_structure(:,:)
     !> Dimension (nnzblock)
     !> It contains the direction of the i-th component.
     !> (reading the matrix line by line)
     !> Given i in (1, nnzblock)
     !> block_structure(i) = 'N'-'T' flag
     !> 'N'= Normal Direction
     !> 'T'= Transpose direction
     !> In the example
     !> block_structure(1)='N'
     !> block_structure(2)='N'
     !> block_structure(3)='N'
     !> block_structure(4)='N
     !> block_structure(5)='T'
     !> block_structure(6)='N'
     character(len=1), allocatable :: direction(:)
     !> Dimension(nblock_row)
     !> Row dimension 
     !> In the example
     !> nrow_vec(1) = matrix_A%nrow ( must match matrix_B%nrow)
     !> nrow_vec(2) = matrix_C%nrow ( must match amtrix_D%nrow)
     !> nrow_vec(3) = matrix_B%ncol ( must match matrix_D%nrow)
     integer, allocatable :: nrow_vec(:)
     !> Dimension(nblock_row)
     !> Row dimension 
     !> In the example
     !> ncol_vec(1) = matrix_A%ncol ( must match matrix_B%nrow)
     !> ncol_vec(2) = matrix_C%ncol 
     !> ncol_vec(3) = matrix_B%ncol ( must match matrix_D%ncol)
     !> ncol_vec(4) = matrix_D%ncol 
     integer, allocatable :: ncol_vec(:)     
     !> Dimension(nrow)
     !> Scratch array for Matrix Vector application
     real(kind=double), allocatable :: scr_nrow(:)
     !> Dimension(ncol)
     !> Scratch array for Matrix Vector application
     real(kind=double), allocatable :: scr_ncol(:)
   contains
     !> static constructor
     !> (procedure public for type blockmat)
     procedure, public, pass :: init => init_blockmat
     !> static destructor
     !> (procedure public for type blockmat)
     procedure, public, pass :: kill => kill_blockmat
     !> Procedure to compute 
     !>         y = M * x 
     !> (public procedure for type blockmat)
     procedure, public, pass :: matrix_times_vector=> Mxv_block
     !> Procedure to compute 
     !>         y = M^T * x 
     !> with M^T the transposed of a matrix M
     !> (public procedure for type blockmat)
     procedure, public, pass :: matrix_transpose_times_vector => MTxv_block
  end type blockmat
contains
  !>-------------------------------------------------------------
  !> Static constructor.
  !> (procedure public for type blockmat)
  !> Instantiate variable of type blockmat
  !>
  !> usage:
  !>     call 'var'%init(lun_err, nrow, nnzblk )
  !>
  !> where:
  !> \param[in] lun_err               -> integer. Error logical unit
  !> \param[in] nmats                 -> integer. Number of rows 
  !> \param[in] nrow_block            -> integer. Number of columns
  !> \param[in] ncol_block            -> integer. Number of columns
  !> \param[in] nnzblk                -> integer. Number of non-zero term
  !>                                     stored in the matrix
  !> \param[in] block_structure       -> integer. dimension(4,nnzblock)
  !>                                      Table describing the block structure
  !> \param[in] is_symmetric          -> Logical. Flag for symmetric or not matrix
  !> \param[in] is_symmetric          -> integer. Dimension of the kernel
  !<-------------------------------------------------------------
  ! TODO automatic detection of symmetry or unsymmetric 
  subroutine init_blockmat(this, lun_err, &
       nmats, mats,&
       nblock_row, nblock_col, &
       nnzblock, block_structure, direction,&
       is_symmetric)
    use Globals
    implicit none
    !var
    class(blockmat),             intent(inout) :: this
    integer,                     intent(in   ) :: lun_err
    integer,                     intent(in   ) :: nmats
    type(array_mat),             intent(in   ) :: mats(nmats)
    integer,                     intent(in   ) :: nblock_row
    integer,                     intent(in   ) :: nblock_col
    integer,                     intent(in   ) :: nnzblock
    integer,                     intent(in   ) :: block_structure(3,nnzblock)
    character(len=1),            intent(in   ) :: direction(nnzblock)
    logical,                     intent(in   ) :: is_symmetric

    ! local vars
    logical :: rc
    integer :: res
    integer :: imat, nrow_loc, ncol_loc,irow,icol,iterm
    character(len=1) :: itrans


    this%nmats      = nmats
    this%nblock_row = nblock_row
    this%nblock_col = nblock_col
    this%nnzblock   = nnzblock

    this%is_symmetric = is_symmetric

    allocate(&
         this%mats(nmats),&
         this%nrow_vec(nblock_row),&
         this%ncol_vec(nblock_col),&
         this%block_structure(3,nnzblock),&
         this%direction(nnzblock),&
         stat=res)
    if (res .ne. 0) &
    rc = IOerr(lun_err, err_alloc, 'init_blockmat', &
         ' type blockmat member mats, nrow_vec, ncol_vec'//&
         ' block_structure direction scrt_nrow, scr_ncol')

    !
    ! copy block structure components
    !
    this%mats            = mats
    this%block_structure = block_structure
    this%direction       = direction

    !
    ! check and define block dimensions
    !
    this%nrow_vec=0
    this%ncol_vec=0
    do iterm = 1, nnzblock
       imat   = block_structure(1,iterm)
       irow   = block_structure(2,iterm)
       icol   = block_structure(3,iterm)
       itrans = direction(iterm)

       nrow_loc = mats(imat)%mat%nrow
       ncol_loc = mats(imat)%mat%ncol

       !write(*,*) imat, irow,icol, itrans, nrow_loc,ncol_loc 
       if (itrans .eq. 'N' ) then
          
          if ( this%nrow_vec(irow) .eq. 0 ) then
             ! set row dimension
             this%nrow_vec(irow) = nrow_loc
          else if (this%nrow_vec(irow) .ne. nrow_loc) then
             ! check row dimension match
             rc = IOerr(lun_err, err_inp, 'init_blockmat', &
                  '  not consistent dimensions for rows value')
          end if
          if ( this%ncol_vec(icol) .eq. 0 ) then
             ! set row dimension
             this%ncol_vec(icol) = ncol_loc
          else if ( this%ncol_vec(icol) .ne. ncol_loc) then
             ! check col dimension match
             rc = IOerr(lun_err, err_inp, 'init_blockmat', &
                  '  not consistent dimensions for cols value')
          end if
       else if (itrans .eq. 'T' ) then
          if ( this%nrow_vec(irow).eq.0 ) then
             ! set row dimension
             this%nrow_vec(irow) = ncol_loc
          else if (this%nrow_vec(irow) .ne. ncol_loc) then
             ! check row dimension match
             rc = IOerr(lun_err, err_inp, 'init_blockmat', &
                  '  not consistent dimensions for rows value')
          end if
          if ( this%ncol_vec(icol).eq.0 ) then
             ! set col dimension
             this%ncol_vec(icol) = nrow_loc
          else if (this%ncol_vec(icol) .ne. nrow_loc) then
             ! check col dimension match
             rc = IOerr(lun_err, err_inp, 'init_blockmat', &
                  '  not consistent dimensions for cols value')
          end if
       else
          ! check col dimension match
          rc = IOerr(lun_err, err_inp, 'init_blockmat', &
               ' only N or T directions allowed, passed ='//itrans)
       end if
    end do

    this%nrow = sum(this%nrow_vec)
    this%ncol = sum(this%ncol_vec)
    
    allocate(&
         this%scr_nrow(this%nrow),&
         this%scr_ncol(this%ncol),&
         stat=res)
    if ( res .ne. 0) &
         rc = IOerr(lun_err, err_alloc, 'init_blockmat', &
         ' type blockmat member mats, nrow_vec, ncol_vec'//&
         ' block_structure scrt_nrow, scr_ncol')



  end subroutine init_blockmat

  !>-------------------------------------------------------------
  !> Static destructor.
  !> (procedure public for type blockmat)
  !> Free memory and diassociate pointers
  !>
  !> usage:
  !>     call 'var'%init(lun_err)
  !>
  !> where:
  !> \param[in] lun_err               -> integer. Error logical unit
  !<-------------------------------------------------------------
  subroutine kill_blockmat(this, lun_err)
    use Globals
    implicit none
    !var
    class(blockmat),   intent(inout) :: this
    integer,           intent(in   ) :: lun_err
    ! local vars
    logical :: rc
    integer :: res

    this%nmats      = 0
    this%nblock_row = 0
    this%nblock_col = 0
    this%nnzblock   = 0
    this%nrow = 0
    this%ncol = 0


    deallocate(&
         this%mats,&
         this%nrow_vec,&
         this%ncol_vec,&
         this%block_structure,&
         this%direction,&
         stat=res)
    if (res .ne. 0) &
         rc = IOerr(lun_err, err_dealloc, 'kill_blockmat', &
         ' type blockmat member mats, nrow_vec, ncol_vec'//&
         ' block_structure scrt_nrow, scr_ncol')

    deallocate(&
         this%scr_nrow,&
         this%scr_ncol,&
         stat=res)
    if (res .ne. 0) &
    rc = IOerr(lun_err, err_dealloc, 'kill_blocksmat', &
         ' type blockmat member scr_nrow, scr_ncol')

  end subroutine kill_blockmat

  
  


  recursive subroutine Mxv_block(this,vec_in,vec_out,info,lun_err)
       use Globals
       implicit none
       class(blockmat),   intent(inout) :: this
       real(kind=double), intent(in   ) :: vec_in(this%ncol)
       real(kind=double), intent(inout) :: vec_out(this%nrow)
       integer,           intent(inout) :: info
       integer,           intent(in   ) :: lun_err

       !local
       integer :: iterm, imat, irow, icol
       character(len=1) :: itrans
       integer :: irowt, icolt,info_loc
       integer :: in_begin,in_end, out_begin, out_end

       info_loc = 0
       vec_out = zero
       this%scr_nrow = zero
       do iterm = 1, this%nnzblock
          imat   = this%block_structure(1,iterm)
          irow   = this%block_structure(2,iterm)
          icol   = this%block_structure(3,iterm)
          itrans = this%direction(iterm)
                    
          in_begin  = sum(this%ncol_vec(1:icol)) - this%ncol_vec(icol) + 1
          in_end    = sum(this%ncol_vec(1:icol))

          out_begin = sum(this%nrow_vec(1:irow)) - this%nrow_vec(irow) + 1
          out_end   = sum(this%nrow_vec(1:irow))

          if ( itrans .eq. 'N' ) then 
             call this%mats(imat)%mat%matrix_times_vector(&
                  vec_in(in_begin:in_end),&
                  this%scr_nrow(out_begin:out_end),info_loc,lun_err)
             vec_out(out_begin:out_end) = vec_out(out_begin:out_end) +&
                  this%scr_nrow(out_begin:out_end)

          else if ( itrans .eq. 'T' ) then 
             call this%mats(imat)%mat%matrix_transpose_times_vector(&
                  vec_in(in_begin:in_end),&
                  this%scr_nrow(out_begin:out_end),info_loc,lun_err)
             vec_out(out_begin:out_end) = vec_out(out_begin:out_end) +&
                  this%scr_nrow(out_begin:out_end)
          end if
          if (info_loc .ne. 0) then
             info = info_loc
             return
          end if
       end do
          
     end subroutine Mxv_block

     !>-------------------------------------------------------------
     !> Matrix transpose times vector multiplication procedure
     !>         vec_out = (M) times (vec_in)
     !> (public procedure for class abs_matrix)
     !> 
     !> usage:
     !>     call 'var'%matrix_transpose_times_vector(vec_in,vec_out,[info])
     !>
     !> where 
     !> \param[in   ] vec_in   -> real, dimension('var'%ncol)
     !>                           vector to be multiplied
     !> \param[inout] vec_out -> real, dimension('var'%nrow)
     !>                           vector (M) times (vec_in) 
     !> \param[in   ] info    -> integer. Info number
     !>                                    in case of error 
     !> \param[in   ] lun_err -> integer. Info number
     !>                                    in case of error 
     !<-------------------------------------------------------------
     subroutine MTxv_block(this,vec_in,vec_out,info,lun_err)
       use Globals
       implicit none
       class(blockmat),   intent(inout) :: this
       real(kind=double), intent(in   ) :: vec_in(this%nrow)
       real(kind=double), intent(inout) :: vec_out(this%ncol)
       integer,           intent(inout) :: info
       integer,           intent(in   ) :: lun_err
       !local
       integer :: iterm, imat, irow, icol
       character(len=1) :: itrans
       integer :: irowt, icolt
       integer :: in_begin,in_end, out_begin, out_end

       vec_out      = zero
       this%scr_ncol = zero
       do iterm = 1, this%nnzblock
          imat   = this%block_structure(1,iterm)
          irow   = this%block_structure(2,iterm)
          icol   = this%block_structure(3,iterm)
          itrans = this%direction(iterm)

          irowt = icol
          icolt = irow
          
          
          in_begin  = sum(this%ncol_vec(1:icolt)) - this%ncol_vec(icolt) + 1
          in_end    = sum(this%ncol_vec(1:icolt))

          out_begin = sum(this%nrow_vec(1:irowt)) - this%nrow_vec(irowt) + 1
          out_end   = sum(this%nrow_vec(1:irowt))
          
          if ( itrans .eq. 'N' ) then    
             call this%mats(imat)%mat%matrix_transpose_times_vector(&
                  vec_in(in_begin:in_end),&
                  this%scr_ncol(out_begin:out_end),info,lun_err)
             vec_out(out_begin:out_end) = vec_out(out_begin:out_end) +&
                  this%scr_ncol(out_begin:out_end)
             
          else if ( itrans .eq. 'T' ) then 
             call this%mats(imat)%mat%matrix_times_vector(&
                  vec_in(in_begin:in_end),&
                  this%scr_ncol(out_begin:out_end),info,lun_err)
             vec_out(out_begin:out_end) = vec_out(out_begin:out_end) +&
                  this%scr_ncol(out_begin:out_end)

          end if
       end do

       
     end subroutine MTxv_block


end module BlockMatrix
