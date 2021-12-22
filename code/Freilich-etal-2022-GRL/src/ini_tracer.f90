      subroutine tracerinit(stepl)
!     --------------------
!     initializes tracer fields 
!     
      USE header
	  implicit none
	  
#include "netcdf.inc"

      integer  i,j,k,l,n,stepl
      double precision z,mldensity,sig,trseed

!     TRACERS
!     =========
  integer :: idzSliceFile
  integer :: idtr

  REAL(kind=rc_kind) ::  rcode
  REAL(kind=rc_kind), dimension(0:NI+1,0:NJ+1, 0:NK+1,ntr)  :: Tr0

  !character (len = 550) :: zslice_data

  !integer start(4), count(4)

  !DATA start /1, 1, 1, 1/
  !count(1)= NI+2
  !count(2)= NJ+2
  !count(3)= NK+2
  !count(4)= ntr

  !idzSliceFile = ncopn('/Volumes/nandi/mara/BIOS24_timescales/adv_efeq3_full_20000.cdf', NCNOWRIT,rcode)

  !idtr = ncvid(idzSliceFile,'tr',rcode)
  !call ncvgt( idzSliceFile, idtr, start, count, Tr0, rcode )
  !call ncclos(idzSliceFile, rcode)

!     Periodic Boundaries (periodic e-w)

do n=0,1
	do k=0,NK+1
		do j=0,NJ+1
			do i=0,NI+1
				do it=1,ntr,2
					z = zc(i,j,k)*DL
					Tr(it,i,j,k,n)=2.5*(1-TANH((z+115)/30))
				end do
				do it=2,ntr,2
					Tr(it,i,j,k,n)=0.0
				end do
			end do
		end do
	end do
end do
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

