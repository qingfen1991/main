      program main
      implicit real*8 (a-h,o-z)
      integer imax,jmax,nxy,nb
      parameter(imax=5,jmax=5,nxy=25,nb=1)
      dimension bxnb(nb),bynb(nb),txnb(nb),tynb(nb),qnxy(nxy),qnb(nb),
     .          in_flag(nxy) 
      real*8 dx,dy
      dx = 1.0
      dy = 1.0
cc    test case  
      do i=1,nxy
          qnxy(i) = i/10.
          in_flag(i) = 0
      enddo
      in_flag(1) = 1
      in_flag(2) = 1
      in_flag(6) = 1

      
      do i=1,nb
          bxnb(i) = 0.1
          bynb(i) = 0.1
          txnb(i) = 0.0
          tynb(i) = 0.1 
          qnb(i) = 1.15
      enddo

      write(*,*) qnb(1)
      call solve(imax,jmax,dx,dy,qnxy,nxy,in_flag,nb,bxnb,bynb,
     .    txnb,tynb,qnb)
      write(*,*) qnb(1)
      end program

      subroutine solve(imax,jmax,dx,dy,qnxy,nxy,in_flag,nb,bxnb,bynb,
     .    txnb,tynb,qnb)
cc    imax: the num of nodes in x direction
cc    jmax: the num of nodes in y direction
cc    dxnxy: dimension in x direction of a cell 
cc    dynxy: dimension in y direction of a cell 
cc    qnxy: value stored in nodes
cc    nxy: total num of nodes 

      implicit real*8 (a-h,o-z)
      integer imax,jmax,nxy,nb
      dimension bxnb(nb),bynb(nb),txnb(nb),tynb(nb),qnxy(nxy),qnb(nb),
     .in_flag(nxy)
      real*8 dx,dy


      integer i,j,bxmin,bxmax,bymin,bymax,p2,p3,closetp1,closetp1tem,
     . surrps(4)
      dimension csurrpsx(4),csurrpsy(4),distob(4),cinter(3,3),fai(3),
     . x(3) 
      real*8 bx,by,cosval,distobtem,costox,costoy,cosminval,cosmaxval
      parameter(cosminval=0.5,cosmaxval=0.86602)
      do i=1, nb
      bx = bxnb(i)
      by = bynb(i)

      bxmin = floor(bx/dx)
      bxmax = bxmin + 1
      bymin = floor(by/dy)
      bymax = bymin + 1

      surrps(1) = bymin*imax + bxmin + 1
      surrps(2) = bymin*imax + bxmax + 1
      surrps(3) = bymax*imax + bxmin + 1
      surrps(4) = bymax*imax + bxmax + 1

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
cc--      cosval = (txnb(i)*tx2 + tynb(i)*ty2)/
cc--     . (distobtem*sqrt(txnb(i)**2+tynb(i)**2))

      if(in_flag(surrps(j)) .eq. 1) then
cc    if node locate inner solid then set dis to maximum          
          distob(j) = sqrt(dx**2+dy**2) + 1
      else 
          distob(j) = distobtem
      endif
      enddo
      closetp1tem = minloc(distob,1)
      closetp1 = surrps(closetp1tem)

cc    ==================take care of boundary==============  
      if(closetp1/imax.eq.0 .or. closetp1/imax.eq.(jmax-1) 
     ..or. mod(closetp1,imax).eq.0 .or. mod(closetp1,imax).eq.1) then
cc    no interpolation on the boundary      
        qnb(i) = qnxy(closetp1)
      else
cc    start interpolating          

          tx2 = csurrpsx(closetp1tem) - bx
          ty2 = csurrpsy(closetp1tem) - by
          costox = abs(tx2/sqrt(tx2**2+ty2**2))
          if(tx2.gt.0 .and. ty2.gt.0) then
              if(costox .gt. cosmaxval) then
                  p2 = closetp1 + 1
                  p3 = closetp1 + imax + 1
              else if(costox.gt.cosminval .and. costox.lt.cosmaxval)then 
                  p2 = closetp1 + 1
                  p3 = closetp1 + imax
              else
                  p2 = closetp1 + imax
                  p3 = closetp1 + imax +1
              endif
          else if(tx2.lt.0 .and. ty2.gt.0) then
              if(costox .gt. cosmaxval) then
                  p2 = closetp1 - 1
                  p3 = closetp1 + imax - 1
              else if(costox.gt.cosminval .and. costox.lt.cosmaxval)then
                  p2 = closetp1 - 1
                  p3 = closetp1 + imax
              else
                  p2 = closetp1 + imax
                  p3 = closetp1 + imax -1
              endif
          else if(tx2.lt.0 .and. ty2.lt.0) then
              if(costox .gt. cosmaxval) then
                  p2 = closetp1 - 1
                  p3 = closetp1 - imax - 1
              else if(costox.gt.cosminval .and. costox.lt.cosmaxval)then
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
              else if(costox.gt.cosminval .and. costox.lt.cosmaxval)then
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
          if(mod(p2,imax) .eq. 0) then
              cinter(2,2) = (mod(p2,imax) - 1)*dx
          else
              cinter(2,2) = (imax - 1)*dx
          endif
          cinter(2,3) = floor(real(p2/imax))*dy
cc        take care of the boundary     
          if(mod(p3,imax) .eq. 0) then
              cinter(3,2) = (mod(p3,imax) - 1)*dx
          else 
              cinter(3,2) = (imax - 1)*dx
          endif
          cinter(3,3) = floor(real(p3/imax))*dy
          
          fai(1) = qnb(i)
          fai(2) = qnxy(p2)
          fai(3) = qnxy(p3)
          call gauss(cinter,x,fai,3)
          qnb(i) = x(1) + x(2)*(mod(closetp1,imax)-1)*dx + 
     .             x(3)*floor(real(closetp1/imax))*dy
      endif
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
      enddo
      fai(j) = fai(j) - c*fai(i)
      enddo
      enddo

      do i=n,1,-1
      c=fai(i)
      do j=i+1,n
      c = c - a(i,j)*x(j)
      enddo
      x(i)=c/a(i,i)
      enddo

      end
