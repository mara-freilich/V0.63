  subroutine write_cdf(step,n)

  USE header,ONLY : NI,NJ,NK,ntr,nconsume,dirout,out1d_int,out2d_int,out3d_int,rc_kind,pvt,pv1,pv2,pv3
  IMPLICIT NONE

  INTEGER :: step,counter_2d,counter_3d,counter_1d,ksurf,islice,jslice,imooring,jmooring,n
  REAL :: sigma,z

! out1d_int  - frequency of 1d output
! out2d_int  - frequency of 2d output
! out3d_int  - frequency of 3d output

  ! 1D output
  if (mod(step,out1d_int).eq.0) then
    counter_1d= step/out1d_int +1 
    imooring=NI/2;jmooring=NJ/2;    call write_cdf_1D_mooring(imooring,jmooring,counter_1d)
  end if

  ! 2D output
  if (mod(step,out2d_int).eq.0) then
    counter_2d= step/out2d_int +1 

    call diag_pv(n);pvt=pv1+pv2+pv3;

 !  ksurf=NK/2;  call write_cdf_2D_sigma(ksurf,counter_2d,n)
 !   ksurf=9;  call write_cdf_2D_sigma(ksurf,counter_2d,n)
 !   ksurf=4;  call write_cdf_2D_sigma(ksurf,counter_2d,n)
 !      ksurf=1;  call write_cdf_2D_sigma(ksurf,counter_2d,n)
    ksurf=NK;    call write_cdf_2D_sigma(ksurf,counter_2d,n)    ! at 0.75 m depth
    ksurf=45;    call write_cdf_2D_sigma(ksurf,counter_2d,n)    ! at 125 m depth
    ksurf=44;    call write_cdf_2D_sigma(ksurf,counter_2d,n)    ! at 135 m depth
    ksurf=46;    call write_cdf_2D_sigma(ksurf,counter_2d,n)    ! at 117 m depth
 !   ksurf=58;    call write_cdf_2D_sigma(ksurf,counter_2d,n)    ! at 117 m depth

    islice=  2;call write_cdf_2D_x(islice,counter_2d,n)
 !   islice=  5;call write_cdf_2D_x(islice,counter_2d,n)
 !   islice=  8;call write_cdf_2D_x(islice,counter_2d,n)
 !   islice= 10;call write_cdf_2D_x(islice,counter_2d,n)
 !   islice= 12;call write_cdf_2D_x(islice,counter_2d,n)
 !   islice= 15;call write_cdf_2D_x(islice,counter_2d,n)
 !   islice= 18;call write_cdf_2D_x(islice,counter_2d,n)
 !   islice= 20;call write_cdf_2D_x(islice,counter_2d,n)
 !   islice= 22;call write_cdf_2D_x(islice,counter_2d,n)
 !   islice= 20;call write_cdf_2D_x(islice,counter_2d,n)
 !   islice= 40;call write_cdf_2D_x(islice,counter_2d,n)
  islice= NI/2;call write_cdf_2D_x(islice,counter_2d,n)
  islice= NI;call write_cdf_2D_x(islice,counter_2d,n)
 !  islice= 60;call write_cdf_2D_x(islice,counter_2d,n)
 !   islice= 80;call write_cdf_2D_x(islice,counter_2d,n)
 !  islice= 96;call write_cdf_2D_x(islic