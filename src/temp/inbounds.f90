logical function inbounds(nparams,proposalstate)
implicit none
integer, intent(in) :: nparams
real(8), intent(in) :: proposalstate(nparams,1)

integer :: j
real(8) :: bounds(2,nparams)
bounds = 0.0d0
bounds(1,:) = -100000.d0
bounds(2,:) = 100000.d0

! beta
bounds(1,1) = 0.0d0
bounds(2,1) = 1.d0

! phibar
bounds(1,2) = 0.0d0
bounds(2,2) = 2.d0

! gz
bounds(1,3) = -10.0d0
bounds(2,3) = 10.d0

! psil
bounds(1,4)= 0.0d0
bounds(2,4) = 1000.d0

! gamma
bounds(1,5) = 0.0d0
bounds(2,5) = 1.0d0

! sigmal
bounds(1,6) = 0.0d0
bounds(2,6) = 10000.d0

! phi
bounds(1,7) = 0.0d0
bounds(2,7) = 100000.0d0

! phiw
bounds(1,8) = 0.0d0
bounds(2,8) = 100000.0d0

! ep
bounds(1,9) = -100.0d0
bounds(2,9) = 100.0d0

! epw
bounds(1,10) = -100.0d0
bounds(2,10) = 100.0d0

! ap
bounds(1,11) = 0.d0
bounds(2,11) = 1.0d0

! aw
bounds(1,12) = 0.d0
bounds(2,12) = 1.0d0

! bw
bounds(1,13) = 0.0d0
bounds(2,13) = 1.0d0

! lamhp
bounds(1,14) = 0.0d0
bounds(2,14) = 100000.0d0

! alpha
bounds(1,15) = 0.0d0
bounds(2,15) = 1.0d0

! delta
bounds(1,16) = 0.0d0
bounds(2,16) = 1.0d0

! phii
bounds(1,17) = 0.0d0
bounds(2,17) = 100.0d0

! sigmaa
bounds(1,18) = 0.0d0
bounds(2,18) = 100.0d0

! gamrs
bounds(1,19) = -1.0d0
bounds(2,19) = 1.0d0

! gamdp
bounds(1,20) = 0.0d0
bounds(2,20) = 10.0d0

! gamdxhp
bounds(1,21) = -10.0d0
bounds(2,21) = 10.0d0

! gamdy
bounds(1,22) = 0.0d0
bounds(2,22) = 10.0d0

! shrgy
bounds(1,23) = 0.0d0
bounds(2,23) = 1.0d0

!sdevtech
bounds(1,24) = -0.0d0
bounds(2,24) = 100.0d0

bounds(1,25) = -1.0d0
bounds(2,25) = 1.0d0

bounds(1,26) = -0.0d0
bounds(2,26) = 100.0d0


bounds(1,27) = -1.0d0
bounds(2,27) = 1.0d0

bounds(1,28) = -0.0d0
bounds(2,28) = 100.0d0

bounds(1,29) = -1.0d0
bounds(2,29) = 1.0d0

bounds(1,30) = -0.0d0
bounds(2,30) = 100.0d0

bounds(1,31) = -1.0d0
bounds(2,31) = 1.0d0

bounds(1,32) = -0.0d0
bounds(2,32) = 100.0d0

bounds(1,33) = -1.0d0
bounds(2,33) = 1.0d0

bounds(1,34) = -0.0d0
bounds(2,34) = 100.0d0

bounds(1,35) = -1.0d0
bounds(2,35) = 1.0d0

bounds(1,36) = 0.0d0
bounds(2,36) = 1000.0d0

! measurementerrors
bounds(1,37) = -0.0d0
bounds(2,37) = 100.0d0
 
bounds(1,38) = -0.0d0
bounds(2,38) = 100.0d0
 
bounds(1,39) = -0.0d0
bounds(2,39) = 100.0d0

bounds(1,40) = -0.0d0
bounds(2,40) = 100.0d0
 
bounds(1,41) = -0.0d0
bounds(2,41) = 100.0d0
 
bounds(1,42) = -0.0d0
bounds(2,42) = 100.0d0

bounds(1,43) = -0.0d0
bounds(2,43) = 100.0d0



inbounds = .true.
do j = 1, nparams
    if (( proposalstate(j,1) < bounds(1,j)) .OR. (proposalstate(j,1) > bounds(2,j)   )  )then
        inbounds = .false.
        write(*,*) 'parameter ', j, ' was out of bounds'
        if (proposalstate(j,1) < bounds(1,j)) then
            write(*,*) proposalstate(j,1), ' < ' , bounds(1,j)
        else
            write(*,*) proposalstate(j,1), ' > ' , bounds(2,j)
        end if
        return
    end if
end do

end function inbounds
