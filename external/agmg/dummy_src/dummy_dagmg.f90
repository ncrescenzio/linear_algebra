SUBROUTINE dagmg( n,a,ja,ia,f,x,ijob,iprint,nrest,iter,tol )
  IMPLICIT NONE
  INTEGER    :: n,ia(n+1),ja(*),ijob,iprint,nrest,iter
  REAL (kind(0.0d0)) :: a(*),f(n),x(n)
  REAL (kind(0.0d0)) :: tol

  Write(*,*) 'The source code for dagmg was not found.'
  Write(*,*) 'Check in http://www.agmg.eu/ to get the source file.'
  Write(*,*) 'This is a dummy subroutine written '
  Write(*,*) 'for making the library to work in any case'
  Write(*,*) 'Simulation will be stopped'
  
  stop
end SUBROUTINE dagmg
