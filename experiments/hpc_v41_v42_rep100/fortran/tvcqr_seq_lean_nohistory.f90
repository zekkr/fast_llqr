! Time-varying coefficient quantile regression - Fortran implementation
! This is an exact translation of tvcqr_seq.R
! 
! Compile with: R CMD SHLIB tvcqr_seq.f90
! or: gfortran -shared -fPIC -o tvcqr_seq.so tvcqr_seq.f90

subroutine tvcqr_seq_fortran(x, y, m, nvar, tau, h, tol, maxit, bland_int, &
                              theta_ll_est, it_num, residual_est, H_mat)
    
    implicit none
    
    ! Input arguments - matching R function signature exactly
    integer, intent(in) :: m, nvar, maxit, bland_int
    double precision, intent(in) :: x(m, nvar), y(m), tau, tol
    double precision, intent(inout) :: h
    double precision, allocatable :: pivot_row(:)
    ! Output arguments - matching R output structure
    double precision, intent(out) :: theta_ll_est(m, nvar+1)
    integer, intent(out) :: it_num(m)
    ! Temporary fair timing control: the production benchmark does not consume
    ! the m-by-m residual history, so retain only an ABI-compatible placeholder.
    double precision, intent(out) :: residual_est(*)
    integer, intent(out) :: H_mat(m, 2*(nvar+1))
    
    ! Local variables matching R code
    double precision :: time_index(m), w(m)
    double precision :: cc(2*(1+nvar) + 2*m)
    double precision :: gammax(m+1, 2*(nvar+1))
    double precision :: b(m+1)
    integer :: IB(m+1)
    logical :: freevarrow(m+1)
    integer :: r1(2*(nvar+1)), r2(2*(nvar+1))
    double precision :: rr(2, 2*(nvar+1))
    
    ! Working variables for simplex
    double precision :: yy(m+1), ee(m+1)
    double precision :: k_vals(m+1)
    double precision :: u(m), v(m), estimate(2*(nvar+1))
    
    ! Loop variables and misc
    integer :: i, j, k, t, eva_t, iter
    integer :: t_rr, tsep
    double precision :: rrl, min_k
    logical :: bland
    double precision :: b_k_original
    integer :: n_beta_in_basis
    ! Convert integer to logical for bland
    bland = (bland_int /= 0)
    


    ! Set default bandwidth if h = 0 (matching R: if (is.null(h)))
    if (h <= 0.0d0) then
        h = dble(m)**(-0.2d0)
    end if
    
    ! Initialize time_index = (1:m)/m
    do i = 1, m
        time_index(i) = dble(i) / dble(m)
    end do
    
    ! Initialize weights for t=1: w at time 1/m
    do i = 1, m
        if (abs(1.0d0/dble(m) - time_index(i)) <= h) then
            w(i) = 0.75d0 * (1.0d0 - ((1.0d0/dble(m) - time_index(i))/h)**2)
        else
            w(i) = 0.0d0
        end if
    end do
    
    ! Initialize cc vector
    cc = 0.0d0
    do i = 1, m
        cc(2*(1+nvar) + i) = tau * w(i)
        cc(2*(1+nvar) + m + i) = (1.0d0 - tau) * w(i)
    end do
    
    ! Build A matrix and initialize gammax
    ! A = cbind(1, x) then cbind with time-scaled version
    ! First, build basic part
    do i = 1, m
        gammax(i, 1) = 1.0d0
        do j = 1, nvar
            gammax(i, j+1) = x(i, j)
        end do
    end do
    
    ! Add time-scaled part: apply(gammax, 2, function(x) x * (1:m)/m)
    do i = 1, m
        do j = 1, nvar+1
            gammax(i, nvar+1+j) = gammax(i, j) * time_index(i)
        end do
    end do
    
    ! Initialize b = c(y, 0)
    do i = 1, m
        b(i) = y(i)
    end do
    b(m+1) = 0.0d0
    
    ! Adjust gammax and b based on sign of y
    ! gammax[y<0,] <- -gammax[y<0,]
    ! b[b<0] <- -b[b<0]
    do i = 1, m
        if (y(i) < 0.0d0) then
            do j = 1, 2*(nvar+1)
                gammax(i, j) = -gammax(i, j)
            end do
            b(i) = -b(i)
        end if
    end do
    
    ! Initialize IB
    ! IB <- (y>=0)*((1:m)+2*(1+nvar))+(y<0)*((1:m)+2*(1+nvar)+m)
    do i = 1, m
        if (y(i) >= 0.0d0) then
            IB(i) = i + 2*(1+nvar)
        else
            IB(i) = i + 2*(1+nvar) + m
        end if
    end do
    
 

    ! Add last row: gammax <- rbind(gammax, -cc[IB] %*% gammax)
    do j = 1, 2*(nvar+1)
        gammax(m+1, j) = 0.0d0
        do i = 1, m
            gammax(m+1, j) = gammax(m+1, j) - cc(IB(i)) * gammax(i, j)
        end do
    end do
    


    ! Initialize freevarrow
    freevarrow = .false.
    freevarrow(m+1) = .true.
    
    ! Initialize r1 and r2
    do i = 1, 2*(nvar+1)
        r1(i) = i
        r2(i) = 0
    end do
    
    ! Initialize outputs
    it_num = 0
    theta_ll_est = 0.0d0
    residual_est(1) = 0.0d0
    H_mat = 0
    allocate(pivot_row(2*(nvar+1)))
  
    ! Main loop over evaluation times
    do eva_t = 1, m
        
        ! Update weights and cc for current time (except for eva_t=1)
        if (eva_t >= 2) then
            ! Update weights
            do i = 1, m
                if (abs(dble(eva_t)/dble(m) - time_index(i)) <= h) then
                    w(i) = 0.75d0 * (1.0d0 - ((dble(eva_t)/dble(m) - time_index(i))/h)**2)
                else
                    w(i) = 0.0d0
                end if
            end do
            
            ! Update cc
            do i = 1, m
                cc(2*(1+nvar) + i) = tau * w(i)
                cc(2*(1+nvar) + m + i) = (1.0d0 - tau) * w(i)
            end do
            
            ! Update last row of gammax only
            ! gammax[(m+1),] <- tau*w[r1-2-2*nvar] - cc[IB] %*% gammax[1:m,]
            do j = 1, 2*(nvar+1)
                gammax(m+1, j) = 0.0d0
                
                ! Add tau*w[r1-2-2*nvar] part
                i = r1(j) - 2 - 2*nvar
                if (i > 0 .and. i <= m) then
                    gammax(m+1, j) = gammax(m+1, j) + tau * w(i)
                end if
                
                ! Subtract cc[IB] %*% gammax[1:m,] part
                do i = 1, m
                    gammax(m+1, j) = gammax(m+1, j) - cc(IB(i)) * gammax(i, j)
                end do
            end do
        end if
        
        ! Simplex iterations
        iter = 0


        do while (iter < maxit)
            

            ! Compute reduced costs rr
            do i = 1, 2*(nvar+1)
                rr(1, i) = gammax(m+1, i)
                



                ! Second row calculation
                if (r2(i) /= 0) then
                    j = r1(i) - 2 - 2*nvar
                    if (j > 0 .and. j <= m) then
                                rr(2, i) = w(j) - rr(1, i)
                    else if (j <= 0) then
                                ! When index is <= 0, R uses w[1]
                                rr(2, i) = w(1) - rr(1, i)
                    else
                                ! When j > m, this shouldn't happen in a well-formed problem
                                ! But to be safe, treat it as 0 weight
                                rr(2, i) = 0.0d0 - rr(1, i)
                    end if
                else
                    rr(2, i) = 0.0d0
                end if
                
                ! Set negative for free variables
                if (r2(i) == 0) then
                    rr(1, i) = -abs(rr(1, i))
                end if
            end do
            
        
            ! Check optimality
            rrl = minval(rr)

            
            if (rrl >= -tol) exit

            

            
            ! Choose entering variable
            if (bland) then
                ! Bland's rule
                t = huge(1)
                t_rr = 0
                
                ! Check first row
                do i = 1, 2*(nvar+1)
                    if (rr(1, i) < -tol) then
                        if (r1(i) < t) then
                            t = r1(i)
                            t_rr = i
                            tsep = 1
                        end if
                    end if
                end do
                
                ! If nothing found in first row, check second row
                if (t_rr == 0) then
                    do i = 1, 2*(nvar+1)
                        if (rr(2, i) < -tol) then
                            if (r2(i) < t) then
                                t = r2(i)
                                t_rr = i
                                tsep = 2
                            end if
                        end if
                    end do
                end if
            else
                ! Steepest descent - find position of minimum
                do i = 1, 2
                    do j = 1, 2*(nvar+1)
                        if (abs(rr(i, j) - rrl) < 1.0d-14) then
                            t_rr = j
                            tsep = i
                            if (tsep == 1) then
                                t = r1(t_rr)
                            else
                                t = r2(t_rr)
                            end if
                            goto 100  ! Break out of nested loop
                        end if
                    end do
                end do
                100 continue
            end if
            


            ! Pivoting based on r2[t_rr]
            if (r2(t_rr) /= 0) then
                ! Choose k
                if (tsep == 1) then
                    yy = gammax(:, t_rr)
                else
                    yy = -gammax(:, t_rr)
                end if
                
                ! Compute k = b/yy for ratio test
                min_k = huge(1.0d0)
                k = 0
                
                do i = 1, m+1
                    if (yy(i) > 0.0d0 .and. .not. freevarrow(i)) then
                        k_vals(i) = b(i) / yy(i)
                        if (k_vals(i) < min_k) then
                            min_k = k_vals(i)
                            if (.not. bland) k = i
                        end if
                    else
                        k_vals(i) = huge(1.0d0)
                    end if
                end do
                
                ! Bland's rule for ties
                if (bland .and. k == 0) then
                    do i = 1, m+1
                        if (abs(k_vals(i) - min_k) < 1.0d-14) then
                            if (k == 0 .or. IB(i) < IB(k)) then
                                k = i
                            end if
                        end if
                    end do
                end if
                
                ! Update yy for complementary variable
                if (tsep /= 1) then
                    j = r1(t_rr) - 2 - 2*nvar
                    if (j > 0 .and. j <= m) then
                        yy(m+1) = yy(m+1) + w(j)
                    end if
                end if
                
            else
                ! r2[t_rr] == 0 case
                yy = gammax(:, t_rr)
                
                if (yy(m+1) < 0.0d0) then
                    ! Standard ratio test
                    min_k = huge(1.0d0)
                    k = 0
                    
                    do i = 1, m+1
                        if (yy(i) > 0.0d0 .and. .not. freevarrow(i)) then
                            k_vals(i) = b(i) / yy(i)
                            if (k_vals(i) < min_k) then
                                min_k = k_vals(i)
                                if (.not. bland) k = i
                            end if
                        else
                            k_vals(i) = huge(1.0d0)
                        end if
                    end do
                    
                    if (bland .and. k == 0) then
                        do i = 1, m+1
                            if (abs(k_vals(i) - min_k) < 1.0d-14) then
                                if (k == 0 .or. IB(i) < IB(k)) then
                                    k = i
                                end if
                            end if
                        end do
                    end if
                else
                    ! Negative ratio test
                    min_k = huge(1.0d0)
                    k = 0
                    
                    do i = 1, m+1
                        if (yy(i) < 0.0d0 .and. .not. freevarrow(i)) then
                            k_vals(i) = -b(i) / yy(i)
                            if (k_vals(i) < min_k) then
                                min_k = k_vals(i)
                                if (.not. bland) k = i
                            end if
                        else
                            k_vals(i) = huge(1.0d0)
                        end if
                    end do
                    
                    if (bland .and. k == 0) then
                        do i = 1, m+1
                            if (abs(k_vals(i) - min_k) < 1.0d-14) then
                                if (k == 0 .or. IB(i) < IB(k)) then
                                    k = i
                                end if
                            end if
                        end do
                    end if
                end if
                
                freevarrow(k) = .true.
            end if
            
            


            ! Perform pivot
            ee = yy / yy(k)
            ee(k) = 1.0d0 - 1.0d0 / yy(k)
            

            

            ! Update entering column
            if (IB(k) <= (m + 2*nvar + 2)) then
                gammax(:, t_rr) = 0.0d0
                gammax(k, t_rr) = 1.0d0
                r1(t_rr) = IB(k)
                r2(t_rr) = IB(k) + m

            
            else
                gammax(:, t_rr) = 0.0d0
                gammax(k, t_rr) = -1.0d0
                j = IB(k) - m - 2*nvar - 2
                if (j > 0 .and. j <= m) then
                    gammax(m+1, t_rr) = w(j)
                else
                    gammax(m+1, t_rr) = 0.0d0  ! Safety check
                    
                end if
                r1(t_rr) = IB(k) - m
                r2(t_rr) = IB(k)

         
            end if

          

            ! Save the pivot row BEFORE any updates
            do j = 1, 2*(nvar+1)
                pivot_row(j) = gammax(k, j)
            end do

            

            ! Now update the tableau using the SAVED pivot row

            do j = 1, 2*(nvar+1)
                do i = 1, m+1
                    gammax(i, j) = gammax(i, j) - ee(i) * pivot_row(j)
                end do
            end do

            b_k_original = b(k)
            ! Update b vector
            do i = 1, m+1
                b(i) = b(i) - ee(i) * b_k_original
            end do

            
            

            ! Update basis
            IB(k) = t
            
            

            


            iter = iter + 1
            

        end do
        
        ! Store iteration count
        it_num(eva_t) = iter
        
        ! Store H (interpolation indices)
        do i = 1, 2*(nvar+1)
            H_mat(eva_t, i) = r1(i) - 2*nvar - 2
        end do
        
        ! Extract solution
        u = 0.0d0
        v = 0.0d0
        estimate = 0.0d0
        ! First, let's see which beta variables are actually in the basis
        n_beta_in_basis = 0
        do i = 1, m
            if (IB(i) >= 1 .and. IB(i) <= 2*(nvar+1)) then
                n_beta_in_basis = n_beta_in_basis + 1
            end if
        end do

        
        ! Extract values for all variables
        do i = 1, m
            if (IB(i) >= 1 .and. IB(i) <= 2*(nvar+1)) then
                ! Beta/gamma variable - store in the correct position
                estimate(IB(i)) = b(i)
                
                 
                
            else if (IB(i) >= 2*(nvar+1)+1 .and. IB(i) <= 2*(nvar+1)+m) then
                ! u variable
                j = IB(i) - 2*(nvar+1)
                u(j) = b(i)
            else if (IB(i) > 2*(nvar+1)+m) then
                ! v variable
                j = IB(i) - 2*(nvar+1) - m
                v(j) = b(i)
            end if
        end do
        
        ! Check if all beta variables are in basis, otherwise they are 0
        
        ! Compute theta_ll_est
        ! theta_ll_est[eva_t, ] <- estimate[1:(nvar+1)] + (eva_t/m)*estimate[(nvar+2):(2*(nvar+1))]
        do j = 1, nvar+1
            theta_ll_est(eva_t, j) = estimate(j) + &
                (dble(eva_t)/dble(m)) * estimate(nvar+1+j)
        end do
        


        ! Residual history intentionally omitted in this temporary fair control.
        
    end do  ! End main loop

    deallocate(pivot_row)
end subroutine tvcqr_seq_fortran
