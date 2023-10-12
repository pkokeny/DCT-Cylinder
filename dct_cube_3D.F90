MODULE mesh
implicit none

contains
SUBROUTINE meshgrid(xgv, ygv, zgv, X, Y, Z)
  double precision,intent(in)   :: xgv(:), ygv(:), zgv(:)
  double precision,intent(out)  :: X(:,:,:), Y(:,:,:), Z(:,:,:)
  integer           :: sX, sY, sZ, i
  sX = size(xgv) ; sY = size(ygv) ; sZ = size(zgv)
 
    do i=1,sZ
      X(:,:,i) = spread( xgv, 1, sY )
      Y(:,:,i) = spread( ygv, 2, sX )
    enddo ! i

  do i=1,sX
    Z(i,:,:) = spread( zgv, 1, sY)
  enddo 
end SUBROUTINE meshgrid


SUBROUTINE meshgridX(xgv, ygv, zgv, X)
  double precision,intent(in)   :: xgv(:), ygv(:), zgv(:)
  double precision,intent(out)  :: X(:,:,:)
  integer           :: sX, sY, sZ, i
  sX = size(xgv) ; sY = size(ygv) ; sZ = size(zgv)
 
    do i=1,sZ
      X(:,:,i) = spread( xgv, 1, sY )
    enddo ! i
    
end SUBROUTINE meshgridX


SUBROUTINE meshgridY(xgv, ygv, zgv, Y)
  double precision,intent(in)   :: xgv(:), ygv(:), zgv(:)
  double precision,intent(out)  :: Y(:,:,:)
  integer           :: sX, sY, sZ, i
  sX = size(xgv) ; sY = size(ygv) ; sZ = size(zgv)
 
    do i=1,sZ
      Y(:,:,i) = spread( ygv, 2, sX )
    enddo ! i

end SUBROUTINE meshgridY

SUBROUTINE meshgridZ(xgv, ygv, zgv, Z)
  double precision,intent(in)   :: xgv(:), ygv(:), zgv(:)
  double precision,intent(out)  :: Z(:,:,:)
  integer           :: sX, sY, sZ, i
  sX = size(xgv) ; sY = size(ygv) ; sZ = size(zgv)
 
  do i=1,sX
    Z(i,:,:) = spread( zgv, 1, sY)
  enddo 
  
end SUBROUTINE meshgridZ

end module mesh


PROGRAM FFT_TEST

use, intrinsic :: iso_c_binding
use mesh
    
include 'fftw3.f03'

type(C_PTR) :: plan
integer :: dimxyz,radius,cropf,lor
integer :: partsp, npartFOVxyz
integer :: xc,yc,zc,t
double precision :: field, gamma, dX, TE, vf, pi, factor, phase
integer,dimension(3) :: highxyz
character(22) :: filenameCk
character(8) :: date
character(10) :: time
character(5) :: zone
integer,dimension(8) :: values

parameter (cropf = 85, lor = 4, dimxyz=((5*cropf)/2+1), radius=3, &
partsp=24, npartFOVxyz=((dimxyz-1)/partsp), lowxyz=((dimxyz-1)*2)/cropf, &
field=3.0D0, gamma=2.675D8, vf=6.7274306D-3, dX=(10.0D-6/64.0D0)/vf, pi=3.14159265D0)

double precision, dimension(dimxyz) :: Xgv,Ygv,Zgv
double precision, dimension(partsp) :: Xps,Yps,Zps
double precision, dimension(dimxyz,dimxyz,dimxyz) :: A, B, C
double precision, dimension(partsp,partsp,partsp) :: mag_unit, XX,YY,ZZ
complex*8, dimension(3,3,3) :: low_3D_k

highxyz(1)=dimxyz
highxyz(2)=dimxyz
highxyz(3)=dimxyz


!GREENS FUNCTION SECTION
do i = 1,dimxyz
    Xgv(i)=-1*(dimxyz)+i
    Ygv(i)=-1*(dimxyz)+i
    Zgv(i)=-1*(dimxyz)+i
end do

call meshgrid(Xgv, Ygv, Zgv, A, B, C)

!A-x^2+z^2
!B-y
!C-r^2
A=A**2+C**2
C=A+B**2
!Define Greens Function

!A-Greens Function
A=((2*B**2-A)/SQRT(C**5))/(4.0D0*pi)
B=0.0
!Put zero in middle
!where (C<=radius**2) A=0.0
where (C==0) A=0.0
C=0.0

print*, 'greens function defined'
call date_and_time(date,time,zone,values)
call date_and_time(DATE=date,ZONE=zone)
call date_and_time(TIME=time)
call date_and_time(VALUES=values)
print '(a,2x,a,2x,a)', date, time, zone
print '(8i5)', values   

!C-FFT of Greens Function
plan = fftw_plan_r2r(3,highxyz, A,C, & 
(/ FFTW_REDFT00, FFTW_REDFT00, FFTW_REDFT00 /),FFTW_ESTIMATE)
call fftw_execute_r2r(plan, A, C)
call fftw_destroy_plan(plan)
A=0.0

print*, 'greens function FFT'
call date_and_time(date,time,zone,values)
call date_and_time(DATE=date,ZONE=zone)
call date_and_time(TIME=time)
call date_and_time(VALUES=values)
print '(a,2x,a,2x,a)', date, time, zone
print '(8i5)', values   

do i = 1,partsp
    Xps(i)=-1*(partsp/2+1)+i
    Yps(i)=-1*(partsp/2+1)+i
    Zps(i)=-1*(partsp/2+1)+i
end do
    
mag_unit(:,:,:)=0.0
call meshgrid(Xps, Yps, Zps, XX, YY, ZZ)
where (XX**2+YY**2+ZZ**2<radius**2) mag_unit=1.0D0
!where (XX**2+YY**2+ZZ**2==0) mag_unit=1.0D0

!A-magnetization grid
do k=1,npartFOVxyz
    do j=1,npartFOVxyz
        do i=1,npartFOVxyz
            A(((i-1)*partsp+1):(i*partsp),((j-1)*partsp+1):(j*partsp), &
            ((k-1)*partsp+1):(k*partsp))=mag_unit(:,:,:)
        end do
    end do
end do

do i = 1,dimxyz
    Xgv(i)=i
    Ygv(i)=i
    Zgv(i)=i
end do

!B-meshgrids foe one dimension at a time, !C still used
call meshgridX(xgv, ygv, zgv, B)
where (B<=cropf) A=0.0
B=0.0

call meshgridY(xgv, ygv, zgv, B)
where (B<=cropf) A=0.0
B=0.0

call meshgridZ(xgv, ygv, zgv, B)
where (B<=cropf) A=0.0
B=0.0

A=A*field*dX

!B-FFT of magnetization
plan = fftw_plan_r2r(3,highxyz, A,B, & 
(/ FFTW_REDFT00, FFTW_REDFT00, FFTW_REDFT00 /),FFTW_ESTIMATE)
call fftw_execute_r2r(plan, A, B)
call fftw_destroy_plan(plan)
A=0.0

print*, 'magnetization FFT'
call date_and_time(date,time,zone,values)
call date_and_time(DATE=date,ZONE=zone)
call date_and_time(TIME=time)
call date_and_time(VALUES=values)
print '(a,2x,a,2x,a)', date, time, zone
print '(8i5)', values   

!C-Field data in Fourier Domain
C=B*C
B=0.0

!A-Field data in spatial domain
plan = fftw_plan_r2r(3,highxyz, C, A, & 
(/ FFTW_REDFT00, FFTW_REDFT00, FFTW_REDFT00 /),FFTW_ESTIMATE)
call fftw_execute_r2r(plan, C, A)
call fftw_destroy_plan(plan)
A=A/((2.0D0*dimxyz)**3.0D0)
C=0.0

print*, 'Field IFFT'
call date_and_time(date,time,zone,values)
call date_and_time(DATE=date,ZONE=zone)
call date_and_time(TIME=time)
call date_and_time(VALUES=values)
print '(a,2x,a,2x,a)', date, time, zone
print '(8i5)', values   

!open(unit = 1, file = "Field_DCT.dat", form="unformatted")
!write(1) A
!close(1)

do t=1,100,2

    factor=gamma*DBLE(t)*1.0D-3
    !B-real part of MR data
    !C-imaginary part of MR data
    do k=1,dimxyz
        do j=1,dimxyz
            do i=1,dimxyz
                phase=DBLE(A(i,j,k))*factor
                B(i,j,k)=COS(phase)
                C(i,j,k)=SIN(phase)
            end do
        end do
    end do
    
    !A-real part of k-space
    plan = fftw_plan_r2r(3,highxyz, B, B, & 
(/ FFTW_REDFT00, FFTW_REDFT00, FFTW_REDFT00 /),FFTW_ESTIMATE)
    call fftw_execute_r2r(plan, B, B)
    call fftw_destroy_plan(plan)

    !B-imaginary part of k-space
    plan = fftw_plan_r2r(3,highxyz, C, C, & 
(/ FFTW_REDFT00, FFTW_REDFT00, FFTW_REDFT00 /),FFTW_ESTIMATE)
    call fftw_execute_r2r(plan, C, C)
    call fftw_destroy_plan(plan)

    do k=1,3
        do j=1,3
            do i=1,3
                low_3D_k(i,j,k)=COMPLEX(B(i,j,k),C(i,j,k))
            end do
        end do
    end do

    WRITE(filenameCk,'("Complex_k_R03LT_",I2,".dat")') t
    open(unit = 1, file = filenameCk, form="unformatted")
    write(1) low_3D_k
    close(1)

    print*, t

end do

END PROGRAM FFT_TEST
