!
      SUBROUTINE MI15ID( ICNTL, CNTL, ISAVE, RSAVE, LSAVE )
      INTEGER          ICNTL( 8 )
      DOUBLE PRECISION CNTL( 4 )
      INTEGER ISAVE(17)
      DOUBLE PRECISION RSAVE(9)
      LOGICAL LSAVE(4)
      INTEGER I

      write(*,*) 'FLEXIBLE GMRES NOT AVAILABLE'      
      write(*,*) 'Source code for mi15 not given.'
      write(*,*) 'This subroutine is a dummy subroutine'
      write(*,*) 'to compile code in any case'
      write(*,*) 'Source code can be found'
      write(*,*) 'in http://www.hsl.rl.ac.uk/index.html'
      write(*,*) 'EXUTION WILL BE STOPPED'

      END

      SUBROUTINE MI15AD( IACT, N, M, W, LDW, Z, LDZ, LOCY, LOCZ, H, LDH,
     *                   RESID, ICNTL, CNTL, INFO, ISAVE, RSAVE, LSAVE )
      DOUBLE PRECISION RESID
      INTEGER          IACT, N, M, LDW, LDZ, LOCY, LOCZ, LDH
      DOUBLE PRECISION CNTL( 4 ), W( LDW, M + 7 ), H( LDH, M + 2 )
      DOUBLE PRECISION Z(LDZ, M)
      INTEGER          ICNTL( 8 ), INFO( 4 )
      INTEGER ISAVE(17)
      DOUBLE PRECISION RSAVE(9)
      LOGICAL LSAVE(4)

      write(*,*) 'FLEXIBLE GMRES NOT AVAILABLE'      
      write(*,*) 'Source code for mi15 not given.'
      write(*,*) 'This subroutine is a dummy subroutine'
      write(*,*) 'to compile code in any case'
      write(*,*) 'Source code can be found'
      write(*,*) 'in http://www.hsl.rl.ac.uk/index.html'
      write(*,*) 'EXUTION WILL BE STOPPED'

      END
