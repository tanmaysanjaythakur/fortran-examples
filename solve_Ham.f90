! Author: tanmaysanjaythakur
! Description: This program solves for eigen values and vectors for complex
! hermitian matrix of order 4x4
! the documentation for cheevx can be found at http://www.netlib.org/lapack/explore-html/index.html

program Solve_Hamiltonian
      implicit none(type,external)

      complex    :: H(4,4), evec(4,4), dwork(1000)
      real*4     :: eval(4),temp, rwork(28),dr
      integer*4  :: i,j,info,lwork, iwork(20),ifail(4),di,M
      character  :: dc
      real       ::abstol
      complex,dimension(:),allocatable :: work
      parameter ( abstol=2*0.1175494351E-37 )  !for accurate eigenvalues abstol should be 2* underflow
      
      
      external   :: cheevx
      external   :: print_mat, print_sol
      
      !getting the optimal workspace
      lwork=-1
      call cheevx('V','A','U',4,H,4,dr,dr,di,di,abstol,M,eval,evec,4,dwork,lwork,rwork,iwork,ifail,info)
      lwork=int(dwork(1))
      allocate(work(lwork))
      !print *,'lwork is',lwork
   
      !some matrix to solve
      do i=1,4
        do j=1,4
            if( i==j) then
                H(i,j)=cmplx(1.0, 0.0)
            else
                H(i,j)=cmplx(0.0,0.0)
                H(j,i)=H(i,j)
            end if
        end do
      end do

      
      !Solving for eigenvalue and vectors for complex hermitian matrix
      call cheevx('V','A','U',4,H,4,dr,dr,di,di,abstol,M,eval,evec,4,work,lwork,rwork,iwork,ifail,info)
      
      !results
      if(info.gt.0) then
            print *,'Computation failed to converge'
      else
            print *,'number of evals found are',M
            print *,'The Hamiltonian is'
            call print_mat(m,4,H)
            print *,'The solution is'
            call print_sol(eval,evec,m,4)
      end if
      
      deallocate(work)

end program Solve_Hamiltonian

!some auxiallary routines for printing
subroutine print_mat(m,n,H)
      integer      :: m,n,i,j
      complex      :: H(m,n)
      
      do i = 1, m
         write(*,9998) ( H( i, j ), j = 1, N )
      end do
9998 format( 11(:,4X,'(',F5.2,',',F5.2,')') )
      return
end subroutine print_mat

subroutine print_sol(eval,evec,m,n)
      integer :: m,n,i,j
      complex :: evec(m,n)
      real*4  :: eval(m)
      
      do i=1,m
         write(*,9998) eval(i),evec(:,i)
      end do
9998  format(:,1X,F5.2,1X,':',1X,'{',4('(',F5.2,',',F5.2,')'),'}')
      return
end subroutine print_sol
!gfortran -L/usr/local/lib/ -o solve_ham solve_ham.f90 -lblas -llapack 
