program freegas
implicit none

!==================================variables and parameters declaration=========

integer::i,j,k,n,ok,nn,tot,mm,zzz,col
real,allocatable,dimension(:)::pos,vel
INTEGER, DIMENSION (1) :: seed = 123456789
integer,parameter::m=1,niter=50000,outunit=50
real::avg_vx,avg_vy,avg_vz,vel1x,vel1y,vel1z,vel2x,vel2y,vel2z,num1,num2,num3,num11,num22,num33
real,parameter::lx=80.0,ly=80.0,lz=80.0,dt=0.005
character(len=30) ::fn
!===========================calling random seed=================================




!!=============Read the number of points========================================

print*,"give the number of particles"
read*,n


tot=n
col=0
!=================allocate the pos and vel matices===============================

allocate(pos(1:niter),stat=ok)

if(ok/=0)then
stop
end if
 

allocate(vel(1:niter),stat=ok)

if(ok/=0)then
stop
end if

!===================Initializing positions=========================================



CALL RANDOM_SEED ()

do nn=1,n


call random_number(num1)
call random_number(num2)
call random_number(num3)

pos(3*nn-2)=num1*lx
pos(3*nn-1)=num2*ly
pos(3*nn)  =num3*lz

 

end do

!=====================Initializing velocities=======================================
  

do nn=1,tot

call random_number(num11)
call random_number(num22)
call random_number(num33)

vel(3*nn-2)= (num11-0.5)*5.0
vel(3*nn-1)= (num22-0.5)*5.0
vel(3*nn)  = (num33-0.5)*5.0


end do

!============================average elocity=cm velocity for m=1=====================
avg_vx=0.0
avg_vy=0.0
avg_vz=0.0


do nn=1,tot

avg_vx=vel(3*nn-2)+avg_vx
avg_vy=vel(3*nn-1)+avg_vy
avg_vz=vel(3*nn)+avg_vz

end do

avg_vx=avg_vx/real(tot)
avg_vy=avg_vy/real(tot)
avg_vz=avg_vz/real(tot)


write(*,*)avg_vx,avg_vy,avg_vz
!===========================subtracting cm velocity===========================================================

do nn=1,tot

vel(3*nn-2)=vel(3*nn-2)-avg_vx
vel(3*nn-1)=vel(3*nn-1)-avg_vy
vel(3*nn)  =vel(3*nn)-avg_vz


end do

!=======================updating positions and velocities==========================================


do i=1,niter


  do nn=1,tot
    pos(3*nn-2)=pos(3*nn-2)+vel(3*nn-2)*dt
    pos(3*nn-1)=pos(3*nn-1)+vel(3*nn-1)*dt
    pos(3*nn)  =pos(3*nn)+vel(3*nn)*dt
  end do




!====================================================================================================
  do nn=1,tot
     do mm=1,tot
     
        if(mm>nn)then
          if(abs(pos(3*nn-2)-pos(3*mm-2))<=1 .and. abs(pos(3*nn-1)-pos(3*mm-1))<=1 .and. abs(pos(3*nn)-pos(3*mm))<=1)then
           vel1x=vel(3*nn-2)
           vel1y=vel(3*nn-1)
           vel1z=vel(3*nn)
           vel2x=vel(3*mm-2)
           vel2y=vel(3*mm-1)
           vel2z=vel(3*mm)
           vel(3*nn-2)=vel2x
           vel(3*nn-1)=vel2y
           vel(3*nn)=vel2z
           vel(3*mm-2)=vel1x
           vel(3*mm-1)=vel1y
           vel(3*mm)=vel1z
           col=col+1
          end if
        end if  
    end do
  end do 



  do nn=1,tot

     if(pos(3*nn-2)>=0.0 .and. pos(3*nn-2)<=80 .and. pos(3*nn)>=0 .and. pos(3*nn)<=80 &
     .and. (pos(3*nn-1)-80)>=0.01)then
     vel(3*nn-2)=vel(3*nn-2)
     vel(3*nn-1)=-vel(3*nn-1)
     vel(3*nn)=vel(3*nn)
     else if (int(pos(3*nn-2))>=0 .and. int(pos(3*nn-2))<=80 .and. int(pos(3*nn))>=0 .and. int(pos(3*nn))<=80 &
     .and.(pos(3*nn-1))<=-0.01)then
     vel(3*nn-2)=vel(3*nn-2)
     vel(3*nn-1)=-vel(3*nn-1)
     vel(3*nn)=vel(3*nn)
     else if (int(pos(3*nn-1))>=0 .and.int(pos(3*nn-1))<=80 .and. int(pos(3*nn))>=0 .and. int(pos(3*nn))<=80 &
     .and. (pos(3*nn-2)-80)>=0.01)then
     vel(3*nn-2)=-vel(3*nn-2)
     vel(3*nn-1)=vel(3*nn-1)
     vel(3*nn)=vel(3*nn)
     else if (int(pos(3*nn-1))>=0 .and. int(pos(3*nn-1))<=80 .and.int(pos(3*nn))>=0 .and.int(pos(3*nn))<=80 &
     .and.(pos(3*nn-2))<=-0.01)then
     vel(3*nn-2)=-vel(3*nn-2)
     vel(3*nn-1)=vel(3*nn-1)
     vel(3*nn)=vel(3*nn)
     else if (int(pos(3*nn-1))>=0 .and. int(pos(3*nn-1))<=80 .and.int(pos(3*nn-2))>=0 .and.int(pos(3*nn-2))<=80 &
     .and.(pos(3*nn)-80)>=0.01)then
     vel(3*nn-2)=vel(3*nn-2)
     vel(3*nn-1)=vel(3*nn-1)
     vel(3*nn)=-vel(3*nn) 
     else if (int(pos(3*nn-1))>=0 .and.int(pos(3*nn-1))<=80 .and.int(pos(3*nn-2))>=0 .and.int(pos(3*nn-2))<=80 &
     .and.(pos(3*nn))<=-0.01)then
     vel(3*nn-2)=vel(3*nn-2)
     vel(3*nn-1)=vel(3*nn-1)
     vel(3*nn)=-vel(3*nn)
     else
     vel(3*nn-2)=vel(3*nn-2)
     vel(3*nn-1)=vel(3*nn-1)
     vel(3*nn)=vel(3*nn)
     end if
  end do


   do nn=1,n
 
    pos(3*nn-2)=pos(3*nn-2)+vel(3*nn-2)*dt
    pos(3*nn-1)=pos(3*nn-1)+vel(3*nn-1)*dt
    pos(3*nn)  =pos(3*nn)+vel(3*nn)*dt
    
   end do
 
 if(mod(i,10)==0)then
  write(fn,fmt='(i0,a)') i,'.dat'
  open(unit=outunit,file=fn,form='formatted')
    do nn=1,tot
     write(outunit,*)pos(3*nn-2),pos(3*nn-1),pos(3*nn),vel(3*nn)
    end do
    close(outunit) 
  end if  
end do

print*,'col= ',col

end program
!======================subroutine============================================================

!====================================================================================================











