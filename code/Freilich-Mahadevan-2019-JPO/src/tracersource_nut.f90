subroutine tracersource_nut(m,n,dtimel)
  !     -----------------------                                           
  USE header

  !     tracer sources                                                    
  !                                                                       
  integer  i,j,k,n,m,keuph
  REAL(kind=rc_kind) :: dtimel,fac,ka,kb,raa,la,kpar1,G,z
  !                                                                       
  !     consump is the uptake (or if negative, then addition) of tracer   

  !tau is the time scale for consump = 1day = 86400 seconds
  !tauinv = 2/(1.d0*86400.d0)    ! tauinv is in (per second)
  fac= dtimel/(1.d0*86400.d0*UL/LEN) ! dtime is non-dim time, TL=LEN/UL
  ! fac is non-dim

  ka = 5.d-2
  !kb = 3.d-2
  raa = 15.d-2
  !rbb = 15.d-3
  !la = 6.d-3
  kpar1 = 0.0461
  G = 2.d0

  keuph = 45 ! level of euphotic depth (130 meters)
  !
  !! Multiply by 1.d-3 to get this in MOLES per SECOND
  do k=keuph,NK
     do j=1,NJ
        do i=1,NI
			z = zc(i,j,k)*DL
            !Tr(1,i,j,k,n) = Tr(1,i,j,k,m)+  &
            !    fac*(-G*0.17*Tr(1,i,j,k,m)/(Tr(1,i,j,k,m)+ka)*Tr(2,i,j,k,m)+rbb*Tr(2,i,j,k,m))
            !Tr(2,i,j,k,n) = Tr(2,i,j,k,m)+  &
            !    fac*(G*0.17*Tr(1,i,j,k,m)/(Tr(1,i,j,k,m)+ka)*Tr(2,i,j,k,m)-rbb*Tr(2,i,j,k,m))
            Tr(1,i,j,k,n) = Tr(1,i,j,k,m)+fac*(-G*TANH(EXP(kpar1*z))*Tr(1,i,j,k,m)/(Tr(1,i,j,k,m)+ka)*Tr(2,i,j,k,m)+raa*Tr(2,i,j,k,m))
            Tr(2,i,j,k,n) = Tr(2,i,j,k,m)+fac*(G*TANH(EXP(kpar1*z))*Tr(1,i,j,k,m)/(Tr(1,i,j,k,m)+ka)*Tr(2,i,j,k,m)-raa*Tr(2,i,j,k,m))
            Tr(3,i,j,k,n) = Tr(3,i,j,k,m)
        end do
     end do
  end do
  
  !print *, G, ra, rb

  do k=1,keuph-1
    do j=1,NJ
        do i=1,NI
			z = zc(i,j,k)*DL
            !Tr(1,i,j,k,n) = Tr(1,i,j,k,m)+fac*(-G*TANH(EXP(kpar1*z))*Tr(1,i,j,k,m)/(Tr(1,i,j,k,m)+ka)*Tr(2,i,j,k,m)+raa*Tr(2,i,j,k,m))
            !Tr(2,i,j,k,n) = Tr(2,i,j,k,m)+fac*(G*TANH(EXP(kpar1*z))*Tr(1,i,j,k,m)/(Tr(1,i,j,k,m)+ka)*Tr(2,i,j,k,m)-raa*Tr(2,i,j,k,m))
            Tr(1,i,j,k,n) = Tr(1,i,j,k,m)+fac*(-G*TANH(EXP(kpar1*z))*Tr(1,i,j,k,m)/(Tr(1,i,j,k,m)+ka)*Tr(2,i,j,k,m)+raa*Tr(2,i,j,k,m))
            Tr(2,i,j,k,n) = Tr(2,i,j,k,m)+fac*(G*TANH(EXP(kpar1*z))*Tr(1,i,j,k,m)/(Tr(1,i,j,k,m)+ka)*Tr(2,i,j,k,m)-raa*Tr(2,i,j,k,m))
            Tr(3,i,j,k,n) = Tr(3,i,j,k,m)
        end do
     end do
  end do
  !
  do m=0,1
    do k=0,NK+1
        do j=0,NJ+1
            do it=1,ntr
                Tr(it,0,j,k,m)= Tr(it,NI,j,k,m)
                Tr(it,NI+1,j,k,m)= Tr(it,1,j,k,m)
            end do
        end do
    end do
  end do
  !Tr(:,:,:,:,0) = Tr(:,:,:,:,1)
  return
end subroutine tracersource_nut
