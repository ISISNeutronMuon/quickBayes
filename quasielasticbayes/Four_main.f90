      SUBROUTINE four(DATA,N,NDIM,ISIGN,IFORM, y_out)
      INCLUDE 'mod_data.f90'
      INCLUDE 'mod_files.f90'
      INCLUDE 'options.f90'
      
	  INTEGER N(1), ISIGN, NDIM, IFORM
	  COMPLEX DATA(4098), y_out(4098)
      real x_in(m_d), y_in(m_d), e_in(m_d)
cf2py intent(in) :: DATA, N, ISIGN, NDIM, IFORM
cf2py intent(out) :: y_out                               !real parameters
      y_out = DATA
      OPEN(UNIT=11,FILE='C:\Users\BTR75544\work\f70.txt')
      WRITE(11,*) 'hi ', ISIGN	  
      DO KK = 1,8
      WRITE(11,*) y_out(KK), "   ", DATA(kk),"   ", KK
      END DO
      CLOSE(11)	  
	  CALL FOUR2(y_out,N,NDIM,ISIGN,IFORM)
      RETURN
      END
