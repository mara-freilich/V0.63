subroutine tracersource(m,n,dtimel)
  !     -----------------------                                           
  USE header

  !     tracer sources                                                    
  !                                                                       
  integer  i,j,k,n,m,keuph
  REAL(kind=rc_kind) :: dtimel,fac,raa,rab,rac,rad,rae,raf,Ga,Gb,Gc,Gd,Ge,Gf,z,N0
  !                                                                       
  !     consump is the uptake (or if negative, then addition) of tracer   

  !tau is the time scale for consump = 1day = 86400 seconds
  !tauinv = 2/(1.d0*86400.d0)    ! tauinv is in (per second)
  fac= dtimel/(1.d0*86400.d0*UL/LEN) ! dtime is non-dim time, TL=LEN/UL
  ! fac is non-dim

  Ga = 3.d-2
  Gb = 3.d-1
  Gc = 3.d-0
  Gd = 3.d-2
  Ge = 3.d-1
  Gf = 3.d-0

  raa = 15.d-3
  rab = 15.d-3
  rac = 15.d-3
  rad = 15.d-1
  rae = 15.d-1
  raf = 15.d-1
  !
  !! Multiply by 1.d-3 to get this in MOLES per SECOND
  do k=1,NK
     do j=1,NJ
        do i=1,NI
			z = zc(i,j,k)*DL
			N0 = 2.5*(1-TANH((z+115)/30))
            Tr(1,i,j,k,n) = Tr(1,i,j,k,m)+fac*(-Ga*(Tr(1,i,j,k,m)-N0))
            Tr(2,i,j,k,n) = Tr(2,i,j,k,m)+fac*(Ga*(Tr(1,i,j,k,m)-N0)-raa*Tr(2,i,j,k,m))
            Tr(3,i,j,k,n) = Tr(3,i,j,k,m)+fac*(-Ga*(Tr(3,i,j,k,m)-N0))
            Tr(4,i,j,k,n) = Tr(4,i,j,k,m)+fac*(Ga*(Tr(3,i,j,k,m)-N0)-raa*Tr(4,i,j,k,m))
            Tr(5,i,j,k,n) = Tr(5,i,j,k,m)+fac*(-Ga*(Tr(5,i,j,k,m)-N0))
            Tr(6,i,j,k,n) = Tr(6,i,j,k,m)+fac*(Ga*(Tr(5,i,j,k,m)-N0)-raa*Tr(6,i,j,k,m))
            Tr(7,i,j,k,n) = Tr(7,i,j,k,m)+fac*(-Ga*(Tr(7,i,j,k,m)-N0))
            Tr(8,i,j,k,n) = Tr(8,i,j,k,m)+fac*(Ga*(Tr(7,i,j,k,m)-N0)-raa*Tr(8,i,j,k,m))
            Tr(9,i,j,k,n) = Tr(9,i,j,k,m)+fac*(-Ga*(Tr(9,i,j,k,m)-N0))
            Tr(10,i,j,k,n) = Tr(10,i,j,k,m)+fac*(Ga*(Tr(9,i,j,k,m)-N0)-raa*Tr(10,i,j,k,m))
            Tr(11,i,j,k,n) = Tr(11,i,j,k,m)+fac*(-Ga*(Tr(11,i,j,k,m)-N0))
            Tr(12,i,j,k,n) = Tr(12,i,j,k,m)+fac*(Ga*(Tr(11,i,j,k,m)-N0)-raa*Tr(12,i,j,k,m))
          end do
     end do
  end do
  
  !print *, G, ra, rb
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
end subroutine tracersource
