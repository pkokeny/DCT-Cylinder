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

SUBROUTINE meshgridXY(xgv, ygv, zgv, X, Y)
  double precision,intent(in)   :: xgv(:), ygv(:), zgv(:)
  double precision,intent(out)  :: X(:,:,:), Y(:,:,:)
  integer           :: sX, sY, sZ, i
  sX = size(xgv) ; sY = size(ygv) ; sZ = size(zgv)
 
    do i=1,sZ
      X(:,:,i) = spread( xgv, 1, sY )
      Y(:,:,i) = spread( ygv, 2, sX )
    enddo ! i

end SUBROUTINE meshgridXY


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
integer :: dimxy,dimz,radius,cropf,lor, s
integer :: partsp, npartFOVxy, npartFOVz, crop, field_only, manual_sum, npartFOV, greens
integer :: xc,yc,zc,t,FOM, n1, n2, n3, seed, field_sum, p, magnitude, partcube, gibbs
double precision :: field, gamma, dX, vf, pi, factor, phase, totalv
double precision,dimension(10) :: magscaleZ, magscaleXY
double precision,dimension(10) :: TE
integer,dimension(3) :: highxyz
character(100) :: filenameCk
character(8) :: date
character(10) :: time
character(5) :: zone
integer,dimension(8) :: values
real :: distance

parameter (cropf = 156, lor = 4, dimxy=10*cropf+1, dimz=30*cropf+1, radius=5*cropf, &
partsp=12, npartFOVxy=((dimxy-1)/partsp), npartFOVz=((dimz-1)/partsp), lowxy=(dimxy-1)/cropf+1, FOM=32, &
lowz=(dimz-1)/cropf+1, field=2.89D0, gamma=2.675D8, pi=3.14159265D0, &
npartFOV=npartFOVxy*npartFOVxy*npartFOVz, partcube=(lor+1)*2)

double precision, dimension(dimxy) :: Xgv,Ygv
double precision, dimension(dimz) ::  Zgv
double precision, dimension(partcube) :: Xps,Yps,Zps
double precision, dimension(dimxy,dimxy,dimz) :: A, B, C, D
double precision, dimension(dimxy,dimxy,cropf) :: field_test  
double precision, dimension(partcube,partcube,partcube) :: mag_unit, XX,YY,ZZ
complex*8, dimension(lowxy,lowxy,lowz) :: low_3D_k
integer,dimension(npartFOV) :: Xi,Yi,Zi
!complex*8, dimension(dimxy,dimxy,dimz) :: E

highxyz(1)=dimz
highxyz(2)=dimxy
highxyz(3)=dimxy
n1=dimz
n2=dimxy
n3=dimxy
!magscale(:)= (/ 0.96D0, 0.79D0, 0.64D0, 0.50D0, 0.42D0, 0.33D0, 0.28D0, 0.23D0, 0.19D0, 0.17D0, &
!0.15D0, 0.14D0, 0.14D0, 0.13D0, 0.13D0 /)
!magscale(:)= (/ 0.40D0, 0.38D0, 0.38D0, 0.41D0, 0.38D0, 0.44D0, 0.43D0 /)
!magscale(:)= (/ 0.92D0, 0.63D0, 0.43D0, 0.27D0, 0.18D0, 0.12D0, 0.10D0, 0.09D0, 0.08D0, 0.06D0, &
!0.03D0, 0.015D0, 0.027D0, 0.02D0, 0.02D0 /)
!magscale(:)= (/ 0.4966D0, 0.3460D0, 0.2372D0, 0.1664D0, 0.1243D0, 0.1143D0, 0.1000D0 /) !0.30ppm?
!magscale(:)= (/ 0.5090D0, 0.4114D0, 0.3325D0, 0.2687D0, 0.2172D0, 0.1755D0, 0.1419D0, &
!0.1147D0, 0.0927D0, 0.0749D0 /) !0.29ppm HeXie TE bulk
!magscaleZ(:)= (/ 0.2606D0, 0.1728D0, 0.1145D0, 0.0759D0, 0.0503D0, 0.0334D0, 0.0221D0, 0.0147D0, &
!0.015D0, 0.015D0 /) !0.56ppm HeXie bulk, estimated from 0.29
magscaleZ(:)= (/ 0.07D0, 0.03D0, 0.013D0, 0.01D0, 0.01D0, 0.01D0, 0.01D0, 0.01D0, 0.01D0, &
0.01D0 /) !1.11ppm HeXie bulk, estimated from 0.29
magscaleXY(:)= (/ 0.80D0, 0.39D0, 0.22D0, 0.18D0, 0.15D0, 0.13D0, 0.11D0, 0.13D0, 0.16D0, &
0.11D0  /) !1.11ppm HeXie bulk parallel
TE(:)= (/ 8.07D-3, 10.46D-3, 12.85D-3, 15.24D-3, 17.63D-3, 20.02D-3, 22.41D-3, 24.80D-3, 27.19D-3, 29.58D-3 /)

field_only=1
greens=0
field_sum=1
manual_sum=0
magnitude=1
crop=1
gibbs=1

do s=1,2,1

A=0.0D0
B=0.0D0
C=0.0D0
D=0.0D0

if (field_only==1) then
print*, 'field only'
!GREENS FUNCTION SECTION
do i = 1,dimxy
    Xgv(i)=-1*(dimxy)+i
    Ygv(i)=-1*(dimxy)+i
end do

do i = 1,dimz
    Zgv(i)=-1*(dimz)+i
end do


call meshgrid(Xgv, Ygv, Zgv, A, B, C)
Xgv(:)=0.0D0
Ygv(:)=0.0D0
Zgv(:)=0.0D0

!D-x^2+y^2
!C-z
!B-r^2
D=(A*A)+(B*B)
B=D+(C*C)

!Define Greens Function
if (greens==1) then
!A-Greens Function
A=((2.0D0*(C*C)-D)/SQRT(B*B*B*B*B))/(4.0D0*pi)
C=0.0D0
!Put zero in middle
!where (C<=lor**2) A=0.0D0
where (B==0) A=0.0D0
B=0.0D0

print*, 'greens function defined'
call date_and_time(date,time,zone,values)
call date_and_time(DATE=date,ZONE=zone)
call date_and_time(TIME=time)
call date_and_time(VALUES=values)
print '(a,2x,a,2x,a)', date, time, zone
print '(8i5)', values  

!C-FFT of Greens Function
plan = fftw_plan_r2r_3d(n1,n2,n3, A,C, & 
FFTW_REDFT00, FFTW_REDFT00, FFTW_REDFT00,FFTW_ESTIMATE)
call fftw_execute_r2r(plan, A, C)
call fftw_destroy_plan(plan)
A=0.0D0

open(unit = 1, file = "Greens_k_parallel_101030_156c.dat", form="unformatted")
write(1) C
close(1)

else
    
open(unit = 1, file = "Greens_k_parallel_101030_156c.dat", form="unformatted")
read(1) C
close(1)
A=0.0D0
B=0.0D0

end if

print*, 'greens function FFT'
call date_and_time(date,time,zone,values)
call date_and_time(DATE=date,ZONE=zone)
call date_and_time(TIME=time)
call date_and_time(VALUES=values)
print '(a,2x,a,2x,a)', date, time, zone
print '(8i5)', values 
 
!CONTINUOUS
!A(:,:,:)=1.0D0   
!mag_unit(:,:,:)=0.0
!call meshgrid(Xps, Yps, Zps, XX, YY, ZZ)
!where (XX**2+YY**2+ZZ**2<=radius**2) mag_unit=1.0D0
!where (XX**2+YY**2+ZZ**2==0) mag_unit=1.0D0

!QUASI-RANDOM ALTERNATE METHOD
!A-magnetization grid
!seed=values(5)+values(6)+values(7)+values(8)
!call srand(seed)
!do k=1,npartFOVz
!    do j=1,npartFOVxy
!        do i=1,npartFOVxy
!            mag_unit(:,:,:)=0.0D0
!            xc=NINT(RAND()*(FOM-1))+((partsp/2)-(FOM/2)+1)
!            yc=NINT(RAND()*(FOM-1))+((partsp/2)-(FOM/2)+1)
!            zc=NINT(RAND()*(FOM-1))+((partsp/2)-(FOM/2)+1)
!            mag_unit(xc,yc,zc)=1.0D0
!            A(((i-1)*partsp+1):(i*partsp),((j-1)*partsp+1):(j*partsp), &
!            ((k-1)*partsp+1):(k*partsp))=mag_unit(:,:,:)
!           !print *, xc, yc, zc
!        end do
!    end do
!end do

!RANDOM TRADITIONAL METHOD
do p = 1,partcube
    Xps(p)=-1*(partcube/2+1)+p
end do
call meshgrid(Xps, Xps, Xps, XX, YY, ZZ)
where (XX**2+YY**2+ZZ**2<=lor**2) mag_unit=1.0d0
seed=values(5)+values(6)+values(7)+values(8)
call srand(seed)

Xi(1)=NINT(RAND()*((dimxy-partcube)-2))+partcube/2+1
Yi(1)=NINT(RAND()*((dimxy-partcube)-2))+partcube/2+1
Zi(1)=NINT(RAND()*((dimz-partcube)-2))+partcube/2+1
A((Xi(1)-(partcube/2)):(Xi(1)+(partcube/2)-1),(Yi(1)-(partcube/2)):(Yi(1)+(partcube/2)-1), &
(Zi(1)-(partcube/2)):(Zi(1)+(partcube/2)-1))=mag_unit

do i=2,npartFOV
    Xi(i)=NINT(RAND()*((dimxy-partcube)-2))+partcube/2+1
    Yi(i)=NINT(RAND()*((dimxy-partcube)-2))+partcube/2+1
    Zi(i)=NINT(RAND()*((dimz-partcube)-2))+partcube/2+1
    print *, i, Xi(i), Yi(i), Zi(i)
    do j=1,(i-1)
        distance=sqrt(real((Xi(i)-Xi(j))**2+(Yi(i)-Yi(j))**2+(Zi(i)-Zi(j))**2))
        do while (distance <= real(lor))
            Xi(i)=NINT(RAND()*((dimxy-partcube)-2))+partcube/2+1
            Yi(i)=NINT(RAND()*((dimxy-partcube)-2))+partcube/2+1
            Zi(i)=NINT(RAND()*((dimz-partcube)-2))+partcube/2+1
            distance=sqrt(real((Xi(i)-Xi(j))**2+(Yi(i)-Yi(j))**2+(Zi(i)-Zi(j))**2))
            !print *, j, Xi(i), Yi(i), Zi(i), distance
        end do 
    !print *, j
    end do
    A((Xi(i)-(partcube/2)):(Xi(i)+(partcube/2)-1),(Yi(i)-(partcube/2)):(Yi(i)+(partcube/2)-1), &
    (Zi(i)-(partcube/2)):(Zi(i)+(partcube/2)-1))=mag_unit
end do

totalv=SUM(A)
vf=totalv/(DBLE(dimxy-1)*DBLE(dimxy-1)*DBLE(dimz-1))
!vf=1.0D0
dX=(1.11D-6)/vf

open (unit = 1, file = "par_values.txt")
write (1,*) "The total volume is ", totalv
write (1,*) "The volume fraction is ", vf
write (1,*) "The particle susceptibility is ", dX
close(1)

!QUASI-RANDOM TRADITIONAL METHOD
!seed=values(5)+values(6)+values(7)+values(8)
!call srand(seed)
!do k=1,npartFOVz
!    do j=1,npartFOVxy
!        do i=1,npartFOVxy
!            mag_unit(:,:,:)=0.0d0
!            xc=NINT(RAND()*(FOM-1))+((partsp/2)-(FOM/2)+1)
!            yc=NINT(RAND()*(FOM-1))+((partsp/2)-(FOM/2)+1)
!            zc=NINT(RAND()*(FOM-1))+((partsp/2)-(FOM/2)+1)
!            xc=xc-(partsp/2+1)
!            yc=yc-(partsp/2+1)
!            zc=zc-(partsp/2+1)
!            do p = 1,partsp
!                Xps(p)=-1*(partsp/2+1-xc)+p
!                Yps(p)=-1*(partsp/2+1-yc)+p
!                Zps(p)=-1*(partsp/2+1-zc)+p
!            end do
!            call meshgrid(Xps, Yps, Zps, XX, YY, ZZ)
!            where (XX**2+YY**2+ZZ**2<=lor**2) mag_unit=1.0d0
!            A(((i-1)*partsp+1):(i*partsp),((j-1)*partsp+1):(j*partsp), &
!            ((k-1)*partsp+1):(k*partsp))=mag_unit(:,:,:)
!        end do
!    end do
!end do

!MAGNITUDE
if (magnitude==1) then

B=1.0d0-A

open(unit = 1, file = "Magnitude_vf12x1r5p2197_rand_par_101030.dat", form="unformatted")
write(1) B
close(1)

B=0.0d0 

end if

!B-meshgrids for one dimension at a time, !C still used
!call meshgridXY(xgv, ygv, zgv, B)
where (D>=radius**2) A=0.0D0
D=0.0D0

do i = 1,dimxy
    Xgv(i)=i
    Ygv(i)=i
end do

do i = 1,dimz
    Zgv(i)=i
end do

call meshgridZ(xgv, ygv, zgv, D)
where (D<=(cropf*15)) A=0.0D0
D=0.0D0

Xgv(:)=0.0D0
Ygv(:)=0.0D0
Zgv(:)=0.0D0

A=A*field*dX

print*, 'Magnetization Grid Defined'
call date_and_time(date,time,zone,values)
call date_and_time(DATE=date,ZONE=zone)
call date_and_time(TIME=time)
call date_and_time(VALUES=values)
print '(a,2x,a,2x,a)', date, time, zone
print '(8i5)', values


!B-FFT of magnetization
plan = fftw_plan_r2r_3d(n1,n2,n3, A,B, & 
FFTW_REDFT00, FFTW_REDFT00, FFTW_REDFT00,FFTW_ESTIMATE)
call fftw_execute_r2r(plan, A, B)
call fftw_destroy_plan(plan)
A=0.0D0

print*, 'magnetization FFT'
call date_and_time(date,time,zone,values)
call date_and_time(DATE=date,ZONE=zone)
call date_and_time(TIME=time)
call date_and_time(VALUES=values)
print '(a,2x,a,2x,a)', date, time, zone
print '(8i5)', values   

!C-Field data in Fourier Domain
C=B*C
B=0.0D0

!A-Field data in spatial domain
plan = fftw_plan_r2r_3d(n1,n2,n3, C, A, & 
FFTW_REDFT00, FFTW_REDFT00, FFTW_REDFT00,FFTW_ESTIMATE)
call fftw_execute_r2r(plan, C, A)
call fftw_destroy_plan(plan)
A=A/(((2.0D0*dimxy)*(2.0D0*dimxy))*(2.0D0*dimz))
C=0.0D0

print*, 'Field IFFT'
call date_and_time(date,time,zone,values)
call date_and_time(DATE=date,ZONE=zone)
call date_and_time(TIME=time)
call date_and_time(VALUES=values)
print '(a,2x,a,2x,a)', date, time, zone
print '(8i5)', values

!field_test(:,:,:)=A(:,:,1:cropf)
!open(unit = 1, file = "Field_vf13x1r3p1728_rand3D_par.dat", form="unformatted")
!rite(1) field_test
!close(1) 

open(unit = 1, file = "Field_vf12x1r5p2197_rand_par_101030.dat", form="unformatted")
write(1) A
close(1) 

end if

if (field_sum==1) then
print*, 'field sum'

A=0.0d0
B=0.0d0
C=0.0d0
D=0.0d0

do i = 0,(dimxy-1)
    Xgv(i)=i
    Ygv(i)=i
end do

do i = 0,(dimz-1)
    Zgv(i)=i
end do

call meshgridXY(Xgv, Ygv, Zgv, A, B)
Xgv(:)=0.0d0
Ygv(:)=0.0d0
Zgv(:)=0.0d0

D=(A*A)+(B*B)
A=0.0d0
B=0.0d0

where (D<radius**2) C=((1.0d0/3.0d0)*(dX*vf)*field)
D=0.0d0
    
open(unit = 1, file = "Field_vf12x1r5p2197_rand_par_101030.dat", form="unformatted")
read(1) A
close(1)

open(unit = 1, file = "Field_vf12x1r5p2197_cont_par_101030.dat", form="unformatted")
read(1) B
close(1)

!B=B*(1.11D0/0.56D0)

!field_test(:,:,:)=A(:,:,1:168)
!open(unit = 1, file = "Field_vf14x29r3_quasi83D.dat", form="unformatted")
!write(1) field_test
!close(1) 

!field_test(:,:,:)=B(:,:,1:168)
!open(unit = 1, file = "Field_vf14x29r3_cont_quasi83D.dat", form="unformatted")
!write(1) field_test
!close(1) 

A=C-B+A
B=0.0d0

end if

if (crop==1) then

print*, 'crop section'
if (magnitude==1) then
D=0.0d0

open(unit = 1, file = "Magnitude_vf12x1r5p2197_rand_par_101030.dat", form="unformatted")
read(1) B
close(1) 

    do k=1,dimz
        do j=1,dimxy
            do i=1,dimxy
                D(i,j,k)=B(dimxy-(i-1),dimxy-(j-1),dimz-(k-1))
            end do
        end do
    end do

!field_test(:,:,:)=D(:,:,1:cropf)
!open(unit = 1, file = "Magnitude_vf21x03r3p1728_rand3D_par_flip.dat", form="unformatted")
!write(1) field_test
!close(1) 
B=0.0d0

end if
p=1
do t=1,10,1

if (gibbs==1) then
!!OUTSIDE MAGNITUDE
    do i = 1,dimxy
        Xgv(i)=i
        Ygv(i)=i
    end do

    do i = 1,dimz
        Zgv(i)=i
    end do

    call meshgridZ(Xgv, Ygv, Zgv, C)
    where (C>(cropf*10)) D=magscaleZ(p)
    C=0.0D0
    Xgv(:)=0.0d0
    Ygv(:)=0.0d0
    Zgv(:)=0.0d0    
    
    do i = 0,(dimxy-1)
        Xgv(i)=i
        Ygv(i)=i
    end do

    do i = 0,(dimz-1)
        Zgv(i)=i
    end do

    call meshgridXY(Xgv, Ygv, Zgv, B, C)
    Xgv(:)=0.0d0
    Ygv(:)=0.0d0
    Zgv(:)=0.0d0

    C=(B*B)+(C*C)
    B=0.0D0
    where (C>=radius**2) D=magscaleZ(p)/magscaleXY(p)
    !where (C>=radius**2) D=1.0d0
    C=0.0D0

end if
    !!CROP PART
    factor=gamma*TE(t)
    !B-real part of MR data
    !C-imaginary part of MR data
    do k=1,(dimz-1)
        do j=1,dimxy
            do i=1,dimxy
                phase=DBLE(A(i,j,k))*factor
                B(i,j,k)=D(i,j,k+1)*DCOS(phase)
                C(i,j,k)=D(i,j,k+1)*DSIN(phase)
            end do
        end do
    end do

k=dimz
    do j=1,dimxy
        do i=1,dimxy
            phase=DBLE(A(i,j,k))*factor
            B(i,j,k)=D(i,j,1)*DCOS(phase)
            C(i,j,k)=D(i,j,1)*DSIN(phase)
        end do
    end do


    !A-real part of k-space
    plan = fftw_plan_r2r_3d(n1,n2,n3, B, B, & 
FFTW_REDFT00, FFTW_REDFT00, FFTW_REDFT00,FFTW_ESTIMATE)
    call fftw_execute_r2r(plan, B, B)
    call fftw_destroy_plan(plan)

    !B-imaginary part of k-space
    plan = fftw_plan_r2r_3d(n1,n2,n3, C, C, & 
FFTW_REDFT00, FFTW_REDFT00, FFTW_REDFT00,FFTW_ESTIMATE)
    call fftw_execute_r2r(plan, C, C)
    call fftw_destroy_plan(plan)

    do k=1,lowz
        do j=1,lowxy
            do i=1,lowxy
                low_3D_k(i,j,k)=COMPLEX(B(i,j,k),C(i,j,k))
            end do
        end do
    end do

    WRITE(filenameCk,'("Complex_k_gibbsMATCH",I1,"_vf12x1r5p2197_rand_par",I2,".dat")') s, t
    open(unit = 1, file = filenameCk, form="unformatted")
    write(1) low_3D_k
    close(1)

    print*, t
    p=p+1
end do
end if

if (manual_sum==1) then
    
do t=1,6,1

    factor=gamma*TE(t)
    !B-complex
    do k=1,dimz
        do j=1,dimxy
            do i=1,dimxy
                phase=DBLE(A(i,j,k))*factor
!                E(i,j,k)=COMPLEX(COS(phase),SIN(phase))
            end do
        end do
    end do
    
    do k=1,lowz
        do j=1,lowxy
            do i=1,lowxy
!                low_3D_k(i,j,k)=SUM(E(((i-1)*cropf+1):(i*cropf),((j-1)*cropf+1):(j*cropf), &
!                ((k-1)*cropf+1):(k*cropf)))
            end do
        end do
    end do  

    WRITE(filenameCk,'("CmplxSum_vf17x56r5_FOV101020_",I2,".dat")') t
    open(unit = 1, file = filenameCk, form="unformatted")
    write(1) low_3D_k
    close(1)
    
end do            
end if

open(unit = 1, file = "Magnitude_vf12x1r5p2197_rand_par_101030.dat", form="unformatted")
close(1, status="delete")
open(unit = 1, file = "Field_vf12x1r5p2197_rand_par_101030.dat", form="unformatted")
close(1, status="delete")

end do

END PROGRAM FFT_TEST
