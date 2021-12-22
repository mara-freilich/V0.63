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
     trseed=1.d-3    ! Chl =0.08 mg/l, Chl2C=15
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
               Tr(2,i,j,k,0) = trseed*10
               !if (k.ge.40.and.k.lt.55) Tr(2,i,j,k,0)=2.d-2 ! initialize P on depth surfaces
               if (sig.lt.mldensity) Tr(1,i,j,k,0)=trseed
            end do
            Tr(3,i,j,45,0) = 1.d0
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

