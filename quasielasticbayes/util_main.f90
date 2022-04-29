      SUBROUTINE LUDCMP2(A,N,NP,INDX,D)
      integer N, NP, D
      real A(8, 8)
      integer INDX(8)
cf2py intent(in) :: N, NP
cf2py intent(in,out,copy) :: A
cf2py intent(in, out, copy) :: INDX, D                               !real parameters
	  CALL  LUDCMP(A,N,NP,INDX,D)
      RETURN
      END 
