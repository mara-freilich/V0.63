!##################
!   size.h
!
! This include file contains the key dimensions of the model arrays.
!
!------------------

! ntr is the number of tracers
! nconsume is the number of columns of the consume array.
! Be careful, neither ntr nor nconsume can be zero.
INTEGER, PARAMETER :: ntr = 12, nconsume = 1


! NI, NJ, NK: dimensions of the model grid.
! ngrid, maxout, maxint,int1: key parameters for the computation of the nonhydrostatic pressure.

! INTEGER,PARAMETER                  :: NI=12,  NJ=12, NK=24,ngrid=3,maxout=6192,    maxint=3942,    int1=3456  
! INTEGER,PARAMETER                  :: NI=64,  NJ=64, NK=64,ngrid=6,maxout=333912,  maxint=299592,  int1=262144
! INTEGER,PARAMETER                  :: NI=128, NJ=128,NK=64,ngrid=6,maxout=1288296, maxint=1198368, int1=1048576
 INTEGER,PARAMETER                  :: NI=256, NJ=320,NK=64,ngrid=7,maxout=6313734, maxint=5991860, int1=5242880
! INTEGER,PARAMETER                  :: NI=128, NJ=192,NK=64,ngrid=7,maxout=1920804, maxint=1797558, int1= 1572864
! INTEGER,PARAMETER                  :: NI=24,  NJ=48, NK=24,ngrid=4,maxout=39992,   maxint=31590,   int1=127648
! INTEGER,PARAMETER                  :: NI=48,  NJ=24, NK=32,ngrid=4,maxout=13272,   maxint=9360,    int1=8192
! INTEGER,PARAMETER                  :: NI=48,  NJ=24, NK=32,ngrid=5,maxout=46416,   maxint=37448,   int1=32768
! INTEGER,PARAMETER                  :: NI=48,  NJ=32, NK=16,ngrid=4,maxout=36312,   maxint=28080,   int1=24576
! INTEGER,PARAMETER                  :: NI=48,  NJ=48, NK=24,ngrid=4,maxout=76352,   maxint=63180,   int1=55296
! INTEGER,PARAMETER                  :: NI=48,  NJ=48, NK=32,ngrid=4,maxout=99512,   maxint =84240,  int1=73728
! INTEGER,PARAMETER                  :: NI=48,  NJ=96, NK=24,ngrid=4,maxout=149072,  maxint=126360,  int1=110592
! INTEGER,PARAMETER                  :: NI=48,  NJ=96, NK=32,ngrid=5,maxout=194472,  maxint=168516,  int1=147456
! INTEGER,PARAMETER                  :: NI=96,  NJ=96, NK=24,ngrid=4,maxout=291092,  maxint=252720,  int1=221184
! INTEGER,PARAMETER                  :: NI=96,  NJ=96, NK=32,ngrid=4,maxout=379472,  maxint=336960,  int1=294912
! INTEGER,PARAMETER                 :: NI=96,  NJ=96, NK=64,ngrid=4,maxout=732992,  maxint=673920,  int1=589824
! INTEGER,PARAMETER                  :: NI=96,  NJ=192,NK=24,ngrid=4,maxout=575132,  maxint=505440,  int1=442368
                                                                                                    
! INTEGER,PARAMETER                  :: NI=96,  NJ=192,NK=32,ngrid=5,maxout=750240,  maxint=674064,  int1=589824
                                                
! INTEGER,PARAMETER                  :: NI=96,  NJ=192,NK=64,ngrid=6,maxout=1449264, maxint=1348164, int1=1179648
! INTEGER,PARAMETER                  :: NI=96,  NJ=384,NK=32,ngrid=5,maxout=1491264, maxint=1348128, int1=1179648
! INTEGER,PARAMETER                  :: NI=192, NJ=384,NK=32,ngrid=5,maxout=2946528, maxint=2696256, int1=2359296
! INTEGER,PARAMETER                  :: NI=192, NJ=192,NK=32,ngrid=5,maxout=1482336, maxint=1348128, int1=1179648
! INTEGER, PARAMETER                 :: NI=768, NJ=384,NK=64,ngrid=4,maxout=22553792,maxint=21565440,int1=18874368
! INTEGER, PARAMETER                 :: NI=1536,NJ=768,NK=64,ngrid=4,maxout=89804672,maxint=86261760,int1=75497472
! INTEGER, PARAMETER                 :: NI=384, NJ=384,NK=32,ngrid=4,maxout=5854352, maxint=5391360, int1=4718592
