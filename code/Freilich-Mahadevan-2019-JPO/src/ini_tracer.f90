      subroutine tracerinit(stepl)
!     --------------------
!     initializes tracer fields 
!     
      USE header
      integer  i,j,k,l,n,stepl
      double precision z,mldensity,sig,trseed
!     
!     Tracer 1 is nitrate

!     if this is a continuation of a run, only trinit is to be initialized
!     Tr 1 is nutrient

!      phydepth= 30.d0
     trseed=1.d-1    ! Chl =0.08 mg/l, Chl2C=15
mldensity = 1025.8 ! nitracline density

!     TRACERS
!     =========

!     Initialze nutrient
      do j=0,NJ+1
         do i=0,NI+1
            do k=NK,1,-1
               !z = -zc(i,j,k)*DL
               sig = rho(i,j,k)
               Tr(1,i,j,k,0) = (sig-mldensity)*12.0
               Tr(2,i,j,k,0) = 1.d-5
               if (sig.lt.mldensity) Tr(1,i,j,k,0)=trseed
			   Tr(3,i,j,k,0) = Tr(1,i,j,k,0)
               !Tr(5,i,j,k,0) = Tr(1,i,j,k,0)
               !Tr(4,i,j,k,0) = Tr(2,i,j,k,0)
            end do
         end do
      end do


!     Periodic Boundaries (periodic e-w)
      do n=0,1
         do k=0,NK+1
            do j=0,NJ+1
               do it=1,ntr
                  Tr(it,0,j,k,n)= Tr(it,NI,j,k,n)
                  Tr(it,NI+1,j,k,n)= Tr(it,1,j,k,n)
               end do
            end do
         end do
      end do
      return
      end

