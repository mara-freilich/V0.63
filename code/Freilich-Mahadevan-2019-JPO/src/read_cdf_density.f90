subroutine read_cdf_density(nstp)
  !subroutine readcdfcn(Tr,uf,vf,wf,xc,yc,zc,xf,yf,zf,nbegin)
  !     ------------------------------------------------------
  !     reads in the netCDF files written by the routine outcdf.f

  !use header, only : NI,NJ,NK,ntr,rc_kind
  use header

  implicit none

#include "netcdf.inc"
  ! include 'dims.f90'

  integer :: nstp
  integer :: idInFile,idzSliceFile
  integer :: idrho, ids, idT

  REAL(kind=rc_kind) ::  rcode

  character (len = 550) :: inname_data, zslice_data

  integer start(3), count(3), countuf(3), countvf(3), countwf(3)
  integer start2d(2), count2d(2)

  DATA start /1, 1, 1/
  DATA start2d /1, 1/
  DATA count2d /NI, NJ/
  
  count(1)= NI+2
  count(2)= NJ+2
  count(3)= NK+2

  countuf(1)= NI+1
  countuf(2)= NJ
  countuf(3)= NK

  countvf(1)= NI
  countvf(2)= NJ+1
  countvf(3)= NK

  countwf(1)= NI
  countwf(2)= NJ
  countwf(3)= NK+1

  WRITE(zslice_data,'("full_",I5.5,".cdf")') nstp
  idzSliceFile = ncopn(TRIM(dirout)//zslice_data, NCNOWRIT,rcode)

  idrho = ncvid(idzSliceFile,'rho',rcode)
  ids = ncvid(idzSliceFile,'s',rcode)
  idT = ncvid(idzSliceFile,'temp',rcode)
  call ncvgt( idzSliceFile, idrho, start, count, rho, rcode)
  call ncvgt( idzSliceFile, ids, start, count, s, rcode)
  call ncvgt( idzSliceFile, idT, start, count, T, rcode)
  call ncclos(idzSliceFile, rcode)

  return
end
