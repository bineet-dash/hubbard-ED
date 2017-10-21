module array
real, dimension(:), allocatable :: evl_s
end module array

program io_test
use array
      integer :: n,i
      real y,T1,mu1,get_mu_s

      open (unit=10, file='eivals_size6_U10_sorted.txt', status='old', action='read')

      n=924

      allocate(evl_s(n))
read(5,*) T1
    do i=1,n
      read(10,*) y
      evl_s(i)=y
    end do
    mu1=get_mu_s(6,T1,n)
    print*,T1,mu1
end



!calculation of chemical potential for the system
!***************************************************************************
real function get_mu_s(fill,T,n)
 use array
!    use input
!    use global
	implicit none
	integer i,n,fill
	real T
  real f0, f, fL2, fR, mR, mL, rtmp,m_d
	mR = maxval(evl_s)       !right-side chemical potential
	fr=0.0d0
	do i=1,n
	fr=fr+(1.0d0/(exp((evl_s(i)-mR)/T)+1.0d0))
	end do
	 mL = minval(evl_s)       !left-side chemical potential
	 fL2=0.0d0
 	do i=1,n
 	fL2=fL2+(1.0d0/(exp((evl_s(i)-mL)/T)+1.0d0))
	end do
	m_d = 0.5d0*(mL+mR)    !middle chemical potential
	f=0.0d0
	do i=1,n
	f=f+(1.0d0/(exp((evl_s(i)-m_d)/T)+1.0d0))
	end do
	print*,f,fill

	do while(abs(f-fill).ge.1e-2)
	 m_d = 0.5d0*(mL+mR)
	f=0.0d0
	do i=1,n
	f=f+(1.0d0/(exp((evl_s(i)-m_d)/T)+1.0d0))
	end do
	 if(f.gt.fill)then
	  !if middle filling is above target, make it the new right bound.
	  mR = m_d
	  fR = f
	 elseif(f.lt.fill)then
	  !if middle filling is below target, make it the new left bound.
	  mL = m_d
	  fR = f
	 endif

   print*,f
	enddo

	!Return the middle value
	get_mu_s = m_d
	return
	end function get_mu_s
!****************************************************************************
