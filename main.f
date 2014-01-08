      program main
      implicit real*8 (a-h,o-z)
      integer,parameter::imax=5,jmax=5,nxy=25,nb=1
      dimension bxnb(nb),bynb(nb),txnb(nb),tynb(nb),qnxy(nxy),qnb(nb)
      real*8::dx=1.0,dy=1.0

      integer i,j
      dimension a(3,3),b(3),x(3)
      a = (/2,4,-2,2,7,4,3,7,5/)
      b = (/3,1,-7/)
      x = (/0,0,0/)
      
      do i=1,nxy
      qnxy(i) = i/10.
      enddo
      
      do i=1,nb
      bxnb(i) = 1.3
      bynb(i) = 1.7
      txnb(i) = 0.0
      tynb(i) = 0.1 
      qnb(i) = 10.4
      enddo
      call gauss(a,x,b,3)


cc      call solve(imax,jmax,dx,dy,qnxy,nxy,nb,bxnb,bynb,
cc     .    txnb,tynb,qnb)
      end program

      subroutine solve(imax,jmax,dx,dy,qnxy,nxy,nb,bxnb,bynb,
     .    txnb,tynb,qnb)
cc    imax: the num of nodes in x direction
cc    jmax: the num of nodes in y direction
cc    dxnxy: dimension in x direction of a cell 
cc    dynxy: dimension in y direction of a cell 
cc    qnxy: value stored in nodes
cc    nxy: total num of nodes 

      implicit real*8 (a-h,o-z)
      integer imax,jmax,nxy,nb
      dimension bxnb(nb),bynb(nb),txnb(nb),tynb(nb),qnxy(nxy),qnb(nb)
      real*8 dx,dy


      integer i,j,bxmin,bxmax,bymin,bymax,p2,p3,closetp1,closetp1tem
cc    closetp1 must be dimensioned for using minloc
      integer,dimension(4)::surrps1,surrps2
      dimension csurrpsx(4),csurrpsy(4),distob(4),cinter(3,3),
     . a(3),fai(3)
      real*8 bx,by,cosval,distobtem,costox,costoy
      real*8,parameter::cosminval=0.5,cosmaxval=0.86602
      do i=1, nb
      bx = bxnb(i)
      by = bynb(i)

      bxmin = floor(bx/dx)
      bxmax = bxmin + 1
      bymin = floor(by/dy)
      bymax = bymin + 1

      surrps1(1) = bymin*imax + bxmin + 1
      surrps1(2) = bymin*imax + bxmax + 1
      surrps1(3) = bymax*imax + bxmin + 1
      surrps1(4) = bymax*imax + bxmax + 1

      csurrpsx(1) = bxmin * dx
      csurrpsx(2) = bxmax * dx
      csurrpsx(3) = bxmin * dx
      csurrpsx(4) = bxmax * dx

      csurrpsy(1) = bymin * dy
      csurrpsy(2) = bymin * dy
      csurrpsy(3) = bymax * dy
      csurrpsy(4) = bymax * dy
cc    find the closet point  
      do j=1, 4
      tx2 = csurrpsx(j) - bx
      ty2 = csurrpsy(j) - by
      distobtem = sqrt(tx2**2 + ty2**2)

cc    caculate cosine of angle between two vectors
      cosval = (txnb(i)*tx2 + tynb(i)*ty2)/
     . (distobtem*sqrt(txnb(i)**2+tynb(i)**2))

      if(cosval .lt. 0) then
cc    if node locate inner solid then set dis to maximum          
          distob(j) = sqrt(dx**2+dy**2) + 1
      else 
          distob(j) = distobtem
      endif
      enddo
      closetp1tem = minloc(distob,1)
      closetp1 = surrps1(closetp1tem)
      tx2 = csurrpsx(closetp1tem) - bx
      ty2 = csurrpsy(closetp1tem) - by
      costox = abs(tx2/sqrt(tx2**2+ty2**2))
cc    0.86602 represent 30degree, 0.5 represent 60degree      
      if(tx2.gt.0 .and. ty2.gt.0) then
          if(costox .gt. cosmaxval) then
              p2 = closetp1 + 1
              p3 = closetp1 + imax + 1
          else if(costox.gt.cosminval .and. costox.lt.cosmaxval) then
              p2 = closetp1 + 1
              p3 = closetp1 + imax
          else
              p2 = closetp1 + imax
              p3 = closetp1 + imax + 1
          endif
      else if(tx2.lt.0 .and. ty2.gt.0) then
          if(costox .gt. cosmaxval) then
              p2 = closetp1 - 1
              p3 = closetp1 + imax - 1
          else if(costox.gt.cosminval .and. costox.lt.cosmaxval) then
              p2 = closetp1 - 1
              p3 = closetp1 + imax
          else
              p2 = closetp1 + imax
              p3 = closetp1 + imax - 1
          endif
      else if(tx2.lt.0 .and. ty2.lt.0) then
          if(costox .gt. cosmaxval) then
              p2 = closetp1 - 1
              p3 = closetp1 - imax - 1
          else if(costox.gt.cosminval .and. costox.lt.cosmaxval) then
              p2 = closetp1 - 1
              p3 = closetp1 - imax
          else
              p2 = closetp1 - imax
              p3 = closetp1 - imax - 1
          endif
      else if(tx2.gt.0 .and. ty2.lt.0) then
          if(costox .gt. cosmaxval) then
              p2 = closetp1 + 1
              p3 = closetp1 - imax + 1
          else if(costox.gt.cosminval .and. costox.lt.cosmaxval) then
              p2 = closetp1 + 1
              p3 = closetp1 - imax
          else
              p2 = closetp1 - imax
              p3 = closetp1 - imax + 1
          endif
      endif


      do k=1,3
      cinter(k,1) = 1
      enddo
      cinter(1,2) = bx
      cinter(1,3) = by
      cinter(2,2) = (p2%imax - 1)*dx
      cinter(2,3) = floor(p2/imax)*dy
      cinter(3,2) = (p3%imax - 1)*dx
      cinter(3,3) = floor(p3/imax)*dy
      fai(1) = qnb(i)
      fai(2) = qnxy(p2)
      fai(3) = qnxy(p3)






      
      




      write(*,100) closetp1,p2,p3
      enddo

          
  100 format('result is',3I10)   
      end subroutine

      subroutine gauss(a,x,fai,n)
      implicit real*8(a-h,o-z)
      dimension fai(n),a(n,n),x(n)
      integer n
      real*8 c

      do i=1,n-1
      do j=i+1,n
      c = a(j,i)/a(i,i)
      do k=i,n
      a(j,k) = a(j,k) - c*a(i,k)
      fai(j) = fai(j) - c*fai(i)
      enddo
      enddo
      enddo

      do i=n,1
      c=fai(i)
      do j=i+1,n
      c = c - a(i,j)*x(j)
      enddo
      x(i)=c/a(i,i)
      write(*,*) fai(i)
      enddo

      end

