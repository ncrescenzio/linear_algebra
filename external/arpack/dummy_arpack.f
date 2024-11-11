       subroutine dsaupd
     &   ( ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam,
     &     ipntr, workd, workl, lworkl, info )
      character  bmat*1, which*2
      integer    ido, info, ldv, lworkl, n, ncv, nev
      Double precision
     &           tol
      integer    iparam(11), ipntr(11)
      Double precision
     &     resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)

      write(*,*) 'DUMMY SUBROUTINE dsaupd in dummy_arpack.f '      
      write(*,*) 'ARPACK LIBRARIES NOT INSTALLED CORRECTLY'
      write(*,*) 'program will be stopped'
      stop

      end
      
      subroutine dnaupd
     &   ( ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam,
     &     ipntr, workd, workl, lworkl, info )

      character  bmat*1, which*2
      integer    ido, info, ldv, lworkl, n, ncv, nev
      Double precision
     &           tol
      integer    iparam(11), ipntr(14)
      Double precision
     &           resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)

      write(*,*) 'DUMMY SUBROUTINE dnaupd in dummy_arpack.f'      
      write(*,*) 'ARPACK LIBRARIES NOT INSTALLED CORRECTLY'
      write(*,*) 'program will be stopped'
      stop      
      end

    
      
      SUBROUTINE DMOUT( LOUT, M, N, A, LDA, IDIGIT, IFMT )
      CHARACTER*( * )    IFMT
      INTEGER            IDIGIT, LDA, LOUT, M, N
      DOUBLE PRECISION   A( LDA, * )

      write(*,*) 'DUMMY SUBROUTINE dseupd in DMOUT.f'      
      write(*,*) 'ARPACK LIBRARIES NOT INSTALLED CORRECTLY'
      write(*,*) 'program will be stopped'
      stop
      
      end

        subroutine dseupd (rvec  , howmny, select, d    ,
     &                   z     , ldz   , sigma , bmat ,
     &                   n     , which , nev   , tol  ,
     &                   resid , ncv   , v     , ldv  ,
     &                   iparam, ipntr , workd , workl,
     &                   lworkl, info )

      character  bmat, howmny, which*2
      logical    rvec
      integer    info, ldz, ldv, lworkl, n, ncv, nev
      Double precision sigma, tol
      integer    iparam(7), ipntr(11)
      logical    select(ncv)
      Double precision
     &           d(nev)     , resid(n)  , v(ldv,ncv),
     &     z(ldz, nev), workd(2*n), workl(lworkl)


      write(*,*) 'DUMMY SUBROUTINE dseupd in dummy_arpack.f'      
      write(*,*) 'ARPACK LIBRARIES NOT INSTALLED CORRECTLY'
      write(*,*) 'program will be stopped'
      stop
      end
