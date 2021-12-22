!subroutine stprofile 
subroutine ini_st
  !     --------------------                                              
  USE header
!     initializes s as pot.density with a single vertical profile       
!     Across-front use the TANH profile with tightness as a measure of f
!     spread                                                            
!     Larger the factor, tighter the front and larger b_xx.             
!     tightness=4, gives the needed b_xx to make b_xx H/ (2f^2) < 1.    
  implicit none
  INTEGER, PARAMETER                 :: nps=355, nwork =50000
  real(kind=rc_kind) depth(nps),S_south(nps),S_north(nps),T_south(nps),T_north(nps), &
    &   svert1(64), svert2(64), svert3(64), tvert1(64), tvert2(64), tvert3(64)
!  REAL*8 :: bs(np),cs(np),ds(np), bt(np),ct(np),dt(np),                 &
!       xs(np),ys(np),xt(np),yt(np)
  !integer  i,j,k,it,iseed,npm1,n
  integer  i,j,k,iseed,npm1,n

  real(kind=rc_kind) :: dep(nps)
  real(kind=rc_kind) ::  bs1(nps),cs1(nps),ds1(nps),bT1(nps),cT1(nps),dT1(nps),  &
       sal(nps),temp(nps),z,seval,         &
       zmm,sbkgrnd,z1,z2,potdens
  !real(kind=rc_kind) :: slfac,dum,dscl,rdscl,yarg,ex2y,thy,ran3,              &
  real(kind=rc_kind) :: slfac,dscl,rdscl,yarg,ex2y,thy,ran3,              &
       perturb,slfacnew,zoffset,dz,bfsqbkgrnd,efold,dsalt(0:NJ+1),       &
       sd,am,da(NI),wiggles,amplitude
!     tightness = 10 represents a very tight front, =1 loose front      
!=      parameter (zoffset= 200.)  also used zoffset=50 ! if zoffset=0 m
!=       parameter (zoffset= 200.d0, tightness=0.03)                    
!=       parameter (zoffset= -10.d0)                                    
!       parameter (zoffset= -20.d0) 
  real(kind=rc_kind) :: dwork(nwork)
  parameter (zoffset= 0.d0)
  tightness=0.05d0

  data svert1 / 35.3109, 35.3373, 35.3618, 35.3947, 35.4478, 35.5256, &
	35.6040, 35.6846, 35.7621, 35.8754, 35.9776, 36.0322, 36.0835, &
	36.1579, 36.2321, 36.3038, 36.3674, 36.4245, 36.4798, 36.5248, &
	36.5699, 36.6067, 36.6365, 36.6602, 36.6772, 36.6932, 36.7074, &
	36.7177, 36.7242, 36.7289, 36.7300, 36.7306, 36.7301, 36.7282, &
	36.7251, 36.7248, 36.7267, 36.7277, 36.7277, 36.7317, 36.7447, &
	36.7613, 36.7891, 36.8210, 36.8671, 36.9189, 36.9993, 37.0623, &
	37.0915, 37.0943, 37.0914, 37.0921, 37.0923, 37.0923, 37.0923, &
	37.0929, 37.0933, 37.0933, 37.0933, 37.0933, 37.0933, 37.0935, &
	37.0939, 37.0939 /
  data svert2 / 35.3109, 35.3373, 35.3618, 35.3947, 35.4478, 35.5256, &
	35.6040, 35.6846, 35.7621, 35.8754, 35.9776, 36.0322, 36.0835, &
	36.1579, 36.2321, 36.3038, 36.3674, 36.4245, 36.4798, 36.5248, &
	36.5699, 36.6067, 36.6365, 36.6602, 36.6772, 36.6932, 36.7074, &
	36.7177, 36.7242, 36.7289, 36.7300, 36.7306, 36.7301, 36.7282, &
	36.7251, 36.7248, 36.7267, 36.7277, 36.7277, 36.7317, 36.7447, &
	36.7613, 36.7891, 36.8210, 36.8671, 36.9189, 36.9993, 37.0623, &
	37.0915, 37.0943, 37.0914, 37.0921, 37.0923, 37.0923, 37.0923, &
	37.0929, 37.0933, 37.0933, 37.0933, 37.0933, 37.0933, 37.0935, &
	37.0939, 37.0939 /
  data svert3 / 35.3109, 35.3373, 35.3618, 35.3947, 35.4478, 35.5256, &
	35.6040, 35.6846, 35.7621, 35.8754, 35.9776, 36.0322, 36.0835, &
	36.1579, 36.2321, 36.3038, 36.3674, 36.4245, 36.4798, 36.5248, &
	36.5699, 36.6067, 36.6365, 36.6602, 36.6772, 36.6932, 36.7074, &
	36.7177, 36.7242, 36.7289, 36.7300, 36.7306, 36.7301, 36.7282, &
	36.7251, 36.7248, 36.7267, 36.7277, 36.7277, 36.7317, 36.7447, &
	36.7613, 36.7891, 36.8210, 36.8671, 36.9189, 36.9993, 37.0623, &
	37.0915, 37.0943, 37.0914, 37.0921, 37.0923, 37.0923, 37.0923, &
	37.0929, 37.0933, 37.0933, 37.0933, 37.0933, 37.0933, 37.0935, &
	37.0939, 37.0939 /

  data tvert1 /  7.3703, 7.8303, 8.2207, 8.7004, 9.2721, 9.9915, &
    10.6424, 11.2810, 11.8923, 12.7481, 13.4400, 13.8004, &
    14.1611, 14.6329, 15.1081, 15.5436, 15.8973, 16.2098, &
    16.5541, 16.8621, 17.1734, 17.4402, 17.6718, 17.8752, &
    18.0548, 18.2407, 18.4245, 18.6050, 18.7657, 18.9134, &
    19.0235, 19.1510, 19.3157, 19.4822, 19.6502, 19.7986, &
    19.9312, 20.0542, 20.1679, 20.2985, 20.4560, 20.6427, &
    20.8997, 21.1466, 21.4604, 21.8244, 22.3599, 22.8786, &
    23.5023, 24.0171, 24.1167, 24.1234, 24.1254, 24.1259, &
    24.1273, 24.1280, 24.1290, 24.1303, 24.1310, 24.1314, &
    24.1311, 24.1297, 24.1256, 24.1256 /
  data tvert2 / 7.3703, 7.8303, 8.2207, 8.7004, 9.2721, 9.9915, &
	10.6424, 11.2810, 11.8923, 12.7481, 13.4400, 13.8004, &
	14.1611, 14.6329, 15.1081, 15.5436, 15.8973, 16.2098, &
	16.5211, 16.7790, 17.0417, 17.2614, 17.4474, 17.6066, &
	17.7434, 17.8878, 18.0315, 18.1731, 18.2961, 18.4072, &
	18.4820, 18.5752, 18.7067, 18.8410, 18.9779, 19.0961, &
	19.1994, 19.2941, 19.3803, 19.4844, 19.6161, 19.7779, &
	20.0107, 20.2342, 20.5253, 20.8673, 21.3815, 21.8796, &
	22.5023, 23.0171, 23.1167, 23.1234, 23.1254, 23.1259, &
	23.1273, 23.1280, 23.1290, 23.1303, 23.1310, 23.1314, &
	23.1311, 23.1297, 23.1256, 23.1256 /
  data tvert3 / 7.3703, 7.8303, 8.2207, 8.7004, 9.2721, 9.9915, &
    10.6424, 11.2810, 11.8923, 12.7481, 13.4400, 13.8004, &
    14.1611, 14.6329, 15.1081, 15.5436, 15.8973, 16.2098, &
    16.4881, 16.6959, 16.9100, 17.0826, 17.2230, 17.3380, &
    17.4320, 17.5349, 17.6385, 17.7412, 17.8265, 17.9010, &
    17.9405, 17.9994, 18.0977, 18.1998, 18.3056, 18.3936, &
    18.4676, 18.5340, 18.5927, 18.6703, 18.7762, 18.9131, &
    19.1217, 19.3218, 19.5902, 19.9102, 20.4031, 20.8806, &
    21.5023, 22.0171, 22.1167, 22.1234, 22.1254, 22.1259, &
    22.1273, 22.1280, 22.1290, 22.1303, 22.1310, 22.1314, &
    22.1311, 22.1297, 22.1256, 22.1256 /

  n=0
  it=1
  !mldepth = 80

  do k=1,NK
     sasouth(k)=svert1(k)
     Tsouth(k)=tvert1(k)
     write(60,*) k,sasouth(k),Tsouth(k)
  end do

  do k=1,NK
     smid(k)=svert2(k)
     Tmid(k)=tvert2(k)
     write(61,*) k,smid(k),Tmid(k)
  end do

  do k=1,NK
     snorth(k)=svert3(k)
     Tnorth(k)=tvert3(k)
     write(61,*) k,snorth(k),Tnorth(k)
  end do
 
  ! Perturb the FRONT                                                   
  ! Add random perturbation to front (white noise)
  iseed= 44294
  am=  0.d0    ! mean
  sd= 0.5d0   !standard deviaton
  call rannw (am, sd, iseed, da, NI, dwork, nwork)
  !     This generates a ran number array NI long, with mean 0 and sd=.1
  do j=0,NJ/2
    do i=0,NI+1
     !yfront= 0.5d0*(yc(NJ/2) +yc(0))+da(i)
     yfront= 96 + da(i)
     thy = 0.5 + 0.5*tanh(tightness*(yc(j)-yfront)*2.*PI) 
     !now thy goes from 0 to 1 ,  thy=0 at j=0

        do k=1,NK
            s(i,j,k,0)= sasouth(k) *(1.d0-thy) + smid(k)*thy
            T(i,j,k,0)= Tsouth(k) *(1.d0-thy) + Tmid(k)*thy
         end do
	 s(i,j,0,0)=s(i,j,1,0)
	 T(i,j,0,0)=T(i,j,1,0)
	 s(i,j,NK+1,0)=s(i,j,NK,0)
	 T(i,j,NK+1,0)=T(i,j,NK,0)
      end do
  end do
  !call evalrho(rho,0) 

    ! Perturb the FRONT                                                   
  ! Add random perturbation to front (white noise)
  iseed= 44294
  am=  0.d0    ! mean
  sd= 0.5d0   !standard deviaton
  call rannw (am, sd, iseed, da, NI, dwork, nwork)
  !     This generates a ran number array NI long, with mean 0 and sd=.1
  do j=(NJ/2+1),NJ+1
    do i=0,NI+1
     !yfront= 0.5d0*(yc(NJ+1) +yc(NJ/2+1))+da(i)
     yfront= 224 + da(i)
     thy = 0.5 + 0.5*tanh(tightness*(yc(j)-yfront)*2.*PI) 
     !now thy goes from 0 to 1 ,  thy=0 at j=NJ/2+1
        do k=1,NK
            s(i,j,k,0)= smid(k) *(1.d0-thy) + snorth(k)*thy
            T(i,j,k,0)= Tmid(k) *(1.d0-thy) + Tnorth(k)*thy
         end do
	 s(i,j,0,0)=s(i,j,1,0)
	 T(i,j,0,0)=T(i,j,1,0)
	 s(i,j,NK+1,0)=s(i,j,NK,0)
	 T(i,j,NK+1,0)=T(i,j,NK,0)
      end do
  end do
  !call evalrho(rho,0) 

  return 
END subroutine ini_st
