program isingmodel
implicit none

integer, parameter :: whatever = 10, numTs = 100, numBs = 100
integer :: iterations = 100000, antiloop = 1
integer :: i, j, k, i2, j2, k2, l, ib, it, randomindex
real, dimension(0:(whatever+1),0:(whatever+1),0:(whatever+1)) :: texas = 0
real*8 :: M = 0, M2 = 0, E = 0, E2 = 0, Esum = 0, E2sum = 0, Msum = 0, M2sum = 0, Etemp, Mtemp, r, Evar, Mvar, de
real*8 :: uB = 9.27401E-24, Jay = 2.0929E25, Kay = 1.3806E-23, B = 8, T = 1500, maxB = 15, maxT = 2500, total = 0
character :: magnettype
call init_random_seed()

!uncomment this section to work in eVs instead
!uB = 5.788E-5
!Kay = 8.617E-5
!Jay = 0.04

!REQUEST MATERIAL TYPE
!Assume a competent user that enters chars
print*,'Ferro? Antiferro? Para?'
print*,'(f/a/p)?'
read(*,*) magnettype
if (magnettype == 'a') then !If antiferro, allow for a "burn-in" run
	print*,'Antiferromagnet simulation running...'
	antiloop = 2
	Jay = -Jay
else if (magnettype == 'p') then
	print*,'Paramagnet simulation running...'
	Jay = 0.0
else
	print*,'Ferromagnet simulation running...'
end if

!Output files for E and M with variance
open(unit=100, file = "E.txt")
open(unit=200, file = "M.txt")

!ITERATE THROUGH APPLIED B
do ib = 1,numBs
B = float(ib)*maxB/float(numBs)

!ITERATE THROUGH TEMPERATURE
do it = 1,numTs
T = float(it)*maxT/float(numTs)

! Initialize spin state lattice to north
texas(1:whatever,1:whatever,1:whatever) = 1

!INITIAL ENERGY AND MAGNETIC MOMENT CALCULATION
E = 0.0
M = 0.0
do i = 1,whatever
	do j = 1,whatever
		do k = 1,whatever
		!Sum up each magnetic contribution
		M = M + texas(i,j,k)
		!Add each magnetic energy and neighbor interaction contribution
		E = E + (-uB)*texas(i,j,k)*(B+Jay*uB*(texas(i+1,j,k)+texas(i,j+1,k)&
		+ texas(i,j,k+1) + texas(i-1,j,k) + texas(i,j-1,k) &
		+ texas(i, j, k-1)))
		end do
	end do
end do

!BEGIN METROPOLIS-HASTINGS ALGORITHM
do l = 1,antiloop	!Performs a burn in for antiferro case

E2 = E**2.0
M2 = M**2.0
Msum = 0.0
Esum = 0.0
M2sum = 0.0
E2sum = 0.0

do i = 1,iterations	!Actual monte-carlo sweeps

	! Find random index and flip its spin value
	i2 = randomindex(whatever)
	j2 = randomindex(whatever)
	k2 = randomindex(whatever)
	texas(i2,j2,k2) = -1.0*texas(i2,j2,k2)
	
	!Recalculate new M and E values
	Mtemp = M + 2.0*texas(i2,j2,k2)
	Etemp = E + (-uB)*texas(i2,j2,k2)*(2.0*B + Jay*4.0*uB*(texas(i2+1,j2,k2) &
	+texas(i2,j2+1,k2) + texas(i2,j2,k2+1) + texas(i2-1,j2,k2) + texas(i2,j2-1,k2) &
	+ texas(i2, j2, k2-1)))
	total = float(i)
	de = Etemp - E
	
	!Call random number to decide acceptance/rejection
	call random_number(r)
	if (r < exp(-de/(Kay*T))) then
		E = Etemp
		M = Mtemp
		E2 = E**2
		M2 = M**2
	else
		texas(i2,j2,k2) = -1.0*texas(i2,j2,k2)
	end if
	Esum = Esum + E
	Msum = Msum + M
	E2sum = E2sum + E2
	M2sum = M2sum + M2
end do ! for MC sweep
end do ! for antiferro burn-in

!CALCULATE AVERAGES
Esum = Esum/total
Msum = Msum/total
E2sum = E2sum/total
M2sum = M2sum/total
Mvar = sqrt(abs(M2sum - Msum**2))
Evar = sqrt(abs(E2sum - Esum**2))

write(100,*) B, T, Esum, Evar
write(200,*) B, T, Msum, Mvar
total = 0

end do ! for T iterations
end do ! for B iterations

close(unit=100)
close(unit=200)
end program


!RANDOM SEED SUBROUTINE
subroutine init_random_seed()
     integer :: i, n, clock
     integer, dimension(:), allocatable :: seed

     CALL RANDOM_SEED(size = n)
     allocate(seed(n))

     CALL SYSTEM_CLOCK(COUNT=clock)

     seed = clock + 37 * (/ (i - 1, i = 1, n) /)
     CALL RANDOM_SEED(PUT = seed)

     deallocate(seed)
end subroutine

!RANDOM INDEX FUNCTION
integer function randomindex(length)
     real :: random
	call random_number(random)
     randomindex = floor(random*length + 1.0)
return
end
