
subroutine bisection2(a,n,x, ifind)
  integer :: n, i
  integer, intent(out) :: ifind
  real    :: a(n), x
  integer :: i1, i2, imid
  i1=1
  i2=n
  if(x>=a(n)) then
     ifind=n
     return
  endif
  if(x<=a(1)) then
     ifind=1
     return
  endif
  do while (i2-i1>1)
     imid=(i1+i2)/2
     if(a(imid)>x) then
        i2=imid
     else 
        if (a(imid)<x) then
           i1=imid
        else
           ifind=imid
           return
        endif
     endif
  end do
  
  ifind=i1
  
end subroutine bisection2
