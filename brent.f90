module brent_method
    implicit none
    
    real(kind = 8), parameter :: CGOLD = 0.381966013D0
    
    contains

    
    !******************************************************************
    !@brief: find the minimum of a function by using parabolic interpolation
    !@param x1 y1 x2 y2 x3 y3: three points
    !@return: p, q, where the new evaluated point(x4, the minmium point) equals to x3 + p/q
    
    !The x of minimum for a paraobolic interpolation by given three points is 
    !x = (x1^2*(y2-y3)+x2^2*(y3-y1)+x3^2*(y1-y2))/(2*(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2)))
    !which is invariant under the permutation of the three points.
    
    !the default value x3< x1, x2
    subroutine parabolic_interpolation(x1, y1, x2, y2, x3, y3, p, q)
        real(kind=8), intent(in):: x1, x2, x3, y1, y2, y3
        real(kind=8), intent(out):: p, q
        
        p = (x2 - x3)**2*(y3 - y1) + (x1 - x3)**2*(y2 - y3)
        q = 2D0*((x2 - x3)*(y3 - y1) + (x1 - x3)*(y2 - y3))
        
    end subroutine parabolic_interpolation
    
    !******************************************************************
    !@param a, b: boundary of the interval
    !@param x, fx: least value of f
    !@param w, fw: second least value of f
    !@param v, fv: the old value of w
    !@param dx: the variation of x in each iteration
    !@param dxold: the variation of x in the 2nd previous iteration
    !@return u: last evaluated point 
    
      
    !fx < fw < fv
    subroutine brent(a, b, v, w, x, fv, fw, fx, u, dx, dxold)
        real(kind=8), intent(inout):: a, b, v, w, x, fv, fw, fx, dx, dxold
        
        real(kind=8):: u, fu         
        real(kind=8):: p, q, deltax, xm
        integer:: qsafe = 0, parasafe = 0
        
        xm = 0.5D0*(a + b)
        
        
        !parabolic interpolation
        call parabolic_interpolation(v, fv, w, fw, x, fx, p, q)
        
        if (q /= 0) qsafe = 1
        
        if (qsafe) then
            deltax = p/q
        else
            deltax = 0D0
        end if
        
        !q /= 0
        !x + deltax is in the interval
        !deltax is less than half of the 2nd previous dx
        if (qsafe .and. (a < x + deltax) .and. (x + deltax < b) .and. (abs(deltax) < 0.5D0*abs(dxold))) parasafe = 1
        
        !here we keep the 2nd previous dx -> dxold
        if (parasafe) then
            dxold = dx
            dx = deltax
        else
            if (x >= xm) then
                dxold = a - x
            else
                dxold = b - x
            end if
            dx = CGOLD*dxold
        end if
        
        u = x + dx
        
        !update v, w, x
        !need a external function to evaluate the function value of f(u)
        
    end subroutine brent
    
    subroutine BrentMethod(a0, b0, xguess, epsilon)
        real(kind=8), intent(in):: a0, b0, xguess
    end subroutine BrentMethod
    
end module brent_method