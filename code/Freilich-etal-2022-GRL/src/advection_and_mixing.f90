subroutine advection_and_mixing(m,n,dtimel,step) 
!     ---------------------------------------------                     
USE header

INTEGER :: i,j,k
INTEGER :: n,m,step,is_sink

REAL(kind=rc_kind) :: dtimel,w0_sink

REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1, 0:NK+1) :: var             ! var=(s,T,u,v,w,Tr(it,:,:,:)
REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1, 0:NK+1) :: uvarx    ! uvarx  is the source term, divergence of the advective fluxes
REAL(kind=rc_kind), dimension(    1:NI  ,1:NJ  , 1:NK  ) :: vardif                ! vardif is the source term from diabatic processes
INTEGER :: iv_compute_kz

iv_compute_kz=1

do it=1,ntr
	if(it==2) then; is_sink=1; w0_sink=2.d-2/86400. !4.d-1**3.5/86400.
		elseif(it==4) then; is_sink=1; w0_sink= 4.d-2/86400. !8.d-1**3.5/86400.
			elseif(it==6) then; is_sink=1; w0_sink= 0.7d-1/86400. !1.4**3.5/86400.
				else; is_sink=0; w0_sink=0.d0
	endif
    ! Flux terms initialized to 0
    uvarx = 0.d0; vardif = 0.d0
    var=Tr(it,:,:,:,m)
    ! ---------------------------------------------------------------
    ! computation of the advective fluxes, using QUICK scheme
    CALL advect(var,uvarx,is_sink,w0_sink)
    ! ---------------------------------------------------------------
    ! computation of the horizontal diabatic fluxes
    vardif=0.;
    call mixing_horizontal(var,vardif);
    !call mixing_isopycnal(var,vardif,1.);!PRINT*,"VISCOUS REDI";
    uvarx(1:NI,1:NJ,1:NK)=uvarx(1:NI,1:NJ,1:NK)-vardif(1:NI,1:NJ,1:NK)
    ! ---------------------------------------------------------------
    ! computation of the vertical diabatic fluxes
    vardif=0.;
    call mixing_vertical(var,vardif,m,step,iv_compute_kz);
    uvarx(1:NI,1:NJ,1:NK)=uvarx(1:NI,1:NJ,1:NK)-vardif(1:NI,1:NJ,1:NK)
    iv_compute_kz=0;
    ! ---------------------------------------------------------------
    ! final summation
    Tr(it,1:NI,1:NJ,1:NK,n) = Tr(it,1:NI,1:NJ,1:NK,0)-dtimel*Jacinv(1:NI,1:NJ,1:NK)*uvarx(1:NI,1:NJ,1:NK)
    if(it==2) Tr(it,1:NI,1:NJ,NK,n) = 1.d-3 !!!!!***** add to top b/c of sinking
  	if(it==4) Tr(it,1:NI,1:NJ,NK,n) = 1.d-3 !!!!!***** add to top b/c of sinking
    if(it==6) Tr(it,1:NI,1:NJ,NK,n) = 1.d-3 !!!!!***** add to top b/c of sinking
  	if(it==8) Tr(it,1:NI,1:NJ,NK,n) = 1.d-3 !!!!!***** add to top b/c of sinking
    if(it==10) Tr(it,1:NI,1:NJ,NK,n) = 1.d-3 !!!!!***** add to top b/c of sinking
  	if(it==12) Tr(it,1:NI,1:NJ,NK,n) = 1.d-3 !!!!!***** add to top b/c of sinking
enddo ! it loop

return 
                                                                        
END
