! Local Linear Quantile Regression - Sequential Algorithm (Fortran)
! OPTIMIZED VERSION - Eliminates redundant multiplications in nested loops
!
! Key optimization: Precomputes coefficient vector cc (matching TVCQR and R code)
! instead of redundantly computing tau*w and (1-tau)*w inside nested loops.
! This eliminates (nvar+1) * m redundant multiplications and conditionals per eval point!
!
! Compile with: gfortran -shared -fPIC -O3 -march=native -funroll-loops -ffast-math -o llqr_seq.so llqr_seq.f90

subroutine llqr_seq_fortran(x, y, z, m, nvar, rounds, tau, h, tol, maxit, &
                            case_int, bland_int, ll_est, d_ll_est, it_num, residual_est, H_mat)
    
    implicit none
    
    ! Input arguments
    integer, intent(in) :: m, nvar, rounds, maxit, case_int, bland_int
    double precision, intent(in) :: x(m), y(m), z(rounds), tau, tol
    double precision, intent(inout) :: h
    
    ! Output arguments
    double precision, intent(out) :: ll_est(rounds)
    double precision, intent(out) :: d_ll_est(rounds)
    integer, intent(out) :: it_num(rounds)
    double precision, intent(out) :: residual_est(*)
    integer, intent(out) :: H_mat(rounds, nvar+1)
    
    ! Local variables - matching tvcqr_seq.f90 structure
    double precision :: A(m, nvar+1)     ! Design matrix [1, x]
    double precision :: w(m)              ! Kernel weights
    double precision :: eva_z(m)          ! z - x for kernel
    double precision :: cc(nvar+1+2*m)    ! Precomputed coefficient vector (OPTIMIZATION!)
    double precision :: gammax(m+1, nvar+1)
    double precision :: b(m+1)
    integer :: IB(m+1)
    logical :: freevarrow(m+1)
    integer :: r1(nvar+1), r2(nvar+1)
    double precision :: rr(2, nvar+1)

    ! Sorting variables (CRITICAL: R code sorts z!)
    double precision :: z_sorted(rounds)
    integer :: z_order(rounds)
    double precision :: ll_est_sorted(rounds), d_ll_est_sorted(rounds)
    integer :: it_num_sorted(rounds)
    integer :: H_mat_sorted(rounds, nvar+1)

    ! Working variables
    double precision :: yy(m+1), ee(m+1), k_vals(m+1)
    double precision :: u(m), v(m), estimate(nvar+1)
    integer :: i, j, k, rd, iter, t_rr, tsep, t
    integer :: idx_r1_minus_offset
    double precision :: rrl, min_k, pi, u_val
    logical :: bland, z_is_sorted
    
    ! Constants
    pi = 4.0d0 * atan(1.0d0)
    bland = (bland_int /= 0)
    
    ! Set default bandwidth
    if (h <= 0.0d0) then
        h = dble(m)**(-0.2d0)
    end if

    ! ============================================================
    ! CRITICAL: Sort z array (matching R behavior!)
    ! R code: original_order <- order(z); z <- z[original_order]
    ! ============================================================
    do i = 1, rounds
        z_sorted(i) = z(i)
        z_order(i) = i
    end do

    z_is_sorted = .true.
    do i = 2, rounds
        if (z_sorted(i) < z_sorted(i-1)) then
            z_is_sorted = .false.
            exit
        end if
    end do
    if (.not. z_is_sorted) then
        do i = 1, rounds-1
            do j = i+1, rounds
                if (z_sorted(i) > z_sorted(j)) then
                    min_k = z_sorted(i)
                    z_sorted(i) = z_sorted(j)
                    z_sorted(j) = min_k
                    k = z_order(i)
                    z_order(i) = z_order(j)
                    z_order(j) = k
                end if
            end do
        end do
    end if

    ! ============================================================
    ! Build design matrix A = [1, x] (constant across eval points)
    ! ============================================================
    do i = 1, m
        A(i, 1) = 1.0d0
        A(i, 2) = x(i)
    end do

    ! ============================================================
    ! Initialize ALL arrays to avoid uninitialized memory issues
    ! ============================================================
    cc = 0.0d0
    gammax = 0.0d0
    b = 0.0d0
    IB = 0
    freevarrow = .false.
    r1 = 0
    r2 = 0
    rr = 0.0d0
    yy = 0.0d0
    ee = 0.0d0
    k_vals = 0.0d0
    estimate = 0.0d0
    u = 0.0d0
    v = 0.0d0

    ! ============================================================
    ! Main loop over evaluation points (USING SORTED z!)
    ! ============================================================
    do rd = 1, rounds

        ! Compute kernel weights for current SORTED z
        do i = 1, m
            eva_z(i) = z_sorted(rd) - x(i)
            u_val = eva_z(i) / h
            if (case_int == 2) then
                if (abs(u_val) <= 1.0d0) then
                    w(i) = 0.75d0 * (1.0d0 - u_val * u_val)
                else
                    w(i) = 0.0d0
                end if
            else
                w(i) = exp(-0.5d0 * u_val**2) / sqrt(2.0d0 * pi)
            end if
        end do

        ! ============================================================
        ! OPTIMIZATION: Precompute coefficient vector cc
        ! This eliminates redundant multiplications in nested loops!
        ! Matching TVCQR pattern and R code: cc <- c(rep(0, nvar+1), tau*w, (1-tau)*w)
        ! ============================================================
        do i = 1, m
            cc(nvar + 1 + i) = tau * w(i)           ! u_i coefficients
            cc(nvar + 1 + m + i) = (1.0d0 - tau) * w(i)  ! v_i coefficients
        end do

! DEBUG: Print kernel weights for first eval point
        if (rd == 1) then
!            write(*,*) '=== FORTRAN: Kernel weights ==='
            do i = 1, m
!                write(*,'(A,I2,A,F10.6)') 'w(', i, ') = ', w(i)
            end do
!            write(*,*) ''
        end if

        ! CRITICAL FIX: Only initialize state for rd==1 (warm-start for rd>=2)
        if (rd == 1) then
            ! Initialize gammax from design matrix
            do i = 1, m
                do j = 1, nvar+1
                    gammax(i, j) = A(i, j)
                end do
            end do

            ! Flip signs for negative y values
            do i = 1, m
                if (y(i) < 0.0d0) then
                    gammax(i, :) = -gammax(i, :)
                end if
            end do

            ! Initialize b vector
            do i = 1, m
                b(i) = abs(y(i))
            end do
            b(m+1) = 0.0d0

            ! Initialize basis IB
            do i = 1, m
                if (y(i) >= 0.0d0) then
                    IB(i) = (nvar + 1) + i  ! u_i
                else
                    IB(i) = (nvar + 1) + m + i  ! v_i
                end if
            end do
            IB(m+1) = nvar + 1

! DEBUG: Print initial basis
!            write(*,*) '=== FORTRAN: Initial basis ==='
            do i = 1, m
!                write(*,'(A,I2,A,I4)') 'IB(', i, ') = ', IB(i)
            end do
!            write(*,*) ''

            ! Initialize freevarrow
            do i = 1, m
                freevarrow(i) = .false.
            end do
            freevarrow(m+1) = .true.

            ! Initialize r1, r2 (non-basic variables)
            do i = 1, nvar+1
                r1(i) = i
                r2(i) = 0
            end do

            rr = 0.0d0

            ! First evaluation point: -cc[IB] %*% gammax
            ! OPTIMIZED: Use precomputed cc array - no conditionals, no redundant multiplications!
            do j = 1, nvar+1
                gammax(m+1, j) = 0.0d0
                do i = 1, m
                    gammax(m+1, j) = gammax(m+1, j) - cc(IB(i)) * gammax(i, j)
                end do
            end do
        else
            ! Subsequent evaluation points (rd >= 2): tau*w[r1-1-nvar] - cc[IB] %*% gammax[1:m,]
            ! Only update objective row - everything else carries over (warm-start)
            ! OPTIMIZED: Use precomputed cc array - no conditionals, no redundant multiplications!
            do j = 1, nvar+1
                ! First term: tau*w[r1-1-nvar]
                idx_r1_minus_offset = r1(j) - 1 - nvar
                if (idx_r1_minus_offset > 0) then
                    gammax(m+1, j) = tau * w(idx_r1_minus_offset)
                else
                    gammax(m+1, j) = tau * w(1)
                end if

                ! Second term: - cc[IB] %*% gammax[1:m,]
                ! OPTIMIZED: Simple lookup, no conditionals!
                do i = 1, m
                    gammax(m+1, j) = gammax(m+1, j) - cc(IB(i)) * gammax(i, j)
                end do
            end do
        end if
        
        ! DEBUG: Print objective row
        if (rd == 1) then
!            write(*,*) '=== FORTRAN: Objective row (last row of gammax) ==='
            do j = 1, nvar+1
!                write(*,'(A,I2,A,F12.6)') 'gammax(m+1,', j, ') = ', gammax(m+1, j)
            end do
!            write(*,*) ''
        end if

        ! ============================================================
        ! Simplex iterations
        ! ============================================================
        iter = 0
        do while (iter < maxit)
            
            ! Step 2: Compute reduced costs
            ! Following the exact pattern from tvcqr_seq.f90 and llqr R code
            do i = 1, nvar+1
                rr(1, i) = gammax(m+1, i)
                
                ! Compute index: r1-1-nvar for LLQR (vs r1-2-2*nvar for TVCQR)
                
                ! Match R logic: (w[(r1-1-nvar) * (r1-1-nvar>0) + (r1-1-nvar<=0)] - rr[1,]) * (r2!=0)
                if (r2(i) /= 0) then

                    idx_r1_minus_offset = r1(i) - 1 - nvar

                    if (idx_r1_minus_offset > 0) then
                        rr(2, i) = w(idx_r1_minus_offset) - rr(1, i)
                    else
                        rr(2, i) = w(1) - rr(1, i)
                    end if
                else
                    rr(2, i) = 0.0d0
                    rr(1, i) = -abs(rr(1, i))
                end if
            end do
            
            ! Special handling for rd==1 (matching R code exactly)
            if (rd == 1) then
                do i = 1, nvar+1
                    if (r2(i) == 0) then
                        rr(2, i) = 0.0d0
                        rr(1, i) = -abs(rr(1, i))
                    end if
                end do
            end if
            
            ! DEBUG: Print reduced costs for first iteration
            if (rd == 1 .and. iter == 0) then
!                write(*,*) '=== FORTRAN: First iteration reduced costs ==='
!                write(*,*) 'r1 array:', r1(1:nvar+1)
!                write(*,*) 'r2 array:', r2(1:nvar+1)
                do i = 1, nvar+1
!                    write(*,'(A,I2,A,F12.6)') 'rr(1,', i, ') = ', rr(1,i)
!                    write(*,'(A,I2,A,F12.6)') 'rr(2,', i, ') = ', rr(2,i)
                end do
!                write(*,'(A,F12.6)') 'Minimum reduced cost: ', minval(rr)
!                write(*,'(A,L)') 'Is optimal (rrl >= -tol)? ', (minval(rr) >= -tol)
!                write(*,*) ''
            end if

            ! Check optimality
            rrl = minval(rr)
            if (rrl >= -tol) then
                exit  ! Optimal solution found
            end if
            
            ! Step 3: Select entering variable
            if (bland) then
                ! Bland's rule
                t_rr = 0
                tsep = 1
                do i = 1, nvar+1
                    if (rr(1, i) < -tol) then
                        if (t_rr == 0 .or. r1(i) < r1(t_rr)) then
                            t_rr = i
                        end if
                    end if
                end do
                if (t_rr == 0) then
                    tsep = 2
                    do i = 1, nvar+1
                        if (rr(2, i) < -tol) then
                            if (t_rr == 0 .or. r2(i) < r2(t_rr)) then
                                t_rr = i
                            end if
                        end if
                    end do
                else
                    tsep = 1
                end if

                ! Get the actual entering variable value
                if (tsep == 1) then
                    t = r1(t_rr)
                else
                    t = r2(t_rr)
                end if


            else
                ! Find the first (i,j) such that rr(i,j) == rrl within tolerance
                tsep = 1
                t_rr = 1
                do j = 1, nvar+1         ! scan columns first to mimic R's which(..., arr.ind=TRUE)[1,]
                  do i = 1, 2
                    if (abs(rr(i, j) - rrl) <= tol) then
                      tsep = i
                      t_rr = j
                      exit
                    end if
                  end do
                  if (abs(rr(tsep, t_rr) - rrl) <= tol) exit
                end do
                if (tsep == 1) then
                  t = r1(t_rr)
                else
                  t = r2(t_rr)
                end if

            end if
            

! DEBUG: Print which variable is entering
            if (rd == 1 .and. iter < 10) then
!                write(*,'(A,I2,A)') '--- Starting iteration ', iter+1, ' ---'
!                write(*,'(A,I2,A,I2)') '  Entering variable index in r1: t_rr = ', t_rr, &
!                                       '  which is variable ', r1(t_rr)
!                write(*,'(A,F12.6)') '  Its reduced cost rr(1,t_rr) = ', rr(1, t_rr)
            end if

            ! Get pivot column
            !do i = 1, m+1
            !    yy(i) = gammax(i, t_rr)
            !end do
            
            ! Step 4-5: Determine pivot column and leaving variable
            k = 0
            min_k = huge(1.0d0)

            if (r2(t_rr) /= 0) then
                ! Case 1: Entering variable is a residual variable (u or v)
                ! Standard ratio test with appropriate sign for yy
                if (tsep == 1) then
                    do i = 1, m+1
                        yy(i) = gammax(i, t_rr)
                    end do
                else
                    do i = 1, m+1
                        yy(i) = -gammax(i, t_rr)
                    end do
                end if
                
                ! Find leaving variable with standard ratio test
                do i = 1, m+1
                    if (yy(i) > tol .and. .not. freevarrow(i)) then
                        k_vals(i) = b(i) / yy(i)
                        if (k_vals(i) < min_k - tol .or. &
                            (abs(k_vals(i) - min_k) < tol .and. (k == 0 .or. IB(i) < IB(k)))) then
                            min_k = k_vals(i)
                            k = i
                        end if
                    end if
                end do
                
                ! Step 6': Adjust last element of yy if needed
                if (tsep /= 1) then
                    yy(m+1) = yy(m+1) + w(r1(t_rr) - 1 - nvar)
                end if
                
            else
                ! Case 2: Entering variable is a FREE coefficient variable (beta)
                do i = 1, m+1
                    yy(i) = gammax(i, t_rr)
                end do
                
                if (yy(m+1) < 0.0d0) then
                    ! Beta can increase from 0 - use positive yy values
                    do i = 1, m+1
                        if (yy(i) > tol .and. .not. freevarrow(i)) then
                            k_vals(i) = b(i) / yy(i)
                            if (k_vals(i) < min_k - tol .or. &
                                (abs(k_vals(i) - min_k) < tol .and. (k == 0 .or. IB(i) < IB(k)))) then
                                min_k = k_vals(i)
                                k = i
                            end if
                        end if
                    end do
                else
                    ! Beta can decrease from 0 - use negative yy values
                    do i = 1, m+1
                        if (yy(i) < -tol .and. .not. freevarrow(i)) then
                            k_vals(i) = -b(i) / yy(i)
                            if (k_vals(i) < min_k - tol .or. &
                                (abs(k_vals(i) - min_k) < tol .and. (k == 0 .or. IB(i) < IB(k)))) then
                                min_k = k_vals(i)
                                k = i
                            end if
                        end if
                    end do
                end if
                
                ! Mark leaving variable as free (special for coefficient variables)
                if (k > 0) then
                    freevarrow(k) = .true.
                end if
            end if
            
            if (k == 0) exit
            
            ! Mark leaving variable as free
            !freevarrow(k) = .true.
            
! DEBUG: Print which variable is leaving
!            if (rd == 1 .and. iter < 10) then
!                write(*,'(A,I2)') '  Leaving row: k = ', k
!                write(*,'(A,I4)') '  Variable leaving basis: IB(k) = ', IB(k)
!                write(*,'(A,F12.6)') '  Value at leaving position: b(k) = ', b(k)
!                write(*,'(A,I2)') '  r2[t_rr] = ', r2(t_rr)
!            end if


            ! Step 6: Pivoting
            do i = 1, m+1
                if (i == k) then
                    ee(i) = 1.0d0 - 1.0d0 / yy(k)
                else
                    ee(i) = yy(i) / yy(k)
                end if
            end do

            ! DEBUG: Print ee values for iteration 4
!            if (rd == 1 .and. iter == 3) then
!                write(*,*) '=== ITERATION 4 DEBUG ==='
!                write(*,'(A,F15.10)') '  yy[k] = ', yy(k)
!                write(*,'(A,F15.10)') '  ee[k] = ', ee(k)
!                write(*,'(A,F15.10)') '  ee[m+1] = ', ee(m+1)
!                write(*,'(A,F15.10,A,F15.10)') '  gammax[k,:] BEFORE pivot update = ', &
!                    gammax(k,1), ', ', gammax(k,2)
!                write(*,'(A,F15.10,A,F15.10)') '  gammax[m+1,:] BEFORE pivot update = ', &
!                    gammax(m+1,1), ', ', gammax(m+1,2)
!            end if

            ! Update pivot column and r1, r2
            ! Following tvcqr_seq.f90 pattern exactly
            if (IB(k) <= (nvar + m + 1)) then
                ! u_i is leaving
                gammax(:, t_rr) = 0.0d0
                gammax(k, t_rr) = 1.0d0
                r1(t_rr) = IB(k)
                r2(t_rr) = IB(k) + m
            else
                ! v_i is leaving
                gammax(:, t_rr) = 0.0d0
                gammax(k, t_rr) = -1.0d0
                gammax(m+1, t_rr) = w(IB(k) - m - nvar - 1)
                r1(t_rr) = IB(k) - m
                r2(t_rr) = IB(k)
            end if

            ! DEBUG: Print after pivot column update
            if (rd == 1 .and. iter == 3) then
!                write(*,'(A,F15.10,A,F15.10)') '  gammax[k,:] AFTER pivot update = ', &
!                    gammax(k,1), ', ', gammax(k,2)
!                write(*,'(A,F15.10,A,F15.10)') '  gammax[m+1,:] AFTER pivot update = ', &
!                    gammax(m+1,1), ', ', gammax(m+1,2)
            end if

            ! CRITICAL FIX: Save gammax[k,:] AFTER pivot column update
            ! This must be done BEFORE row operations to match R's tcrossprod behavior
            do j = 1, nvar+1
                estimate(j) = gammax(k, j)
            end do

            ! Update all columns (including the pivot column) and b
            ! Use SAVED values (estimate) to avoid in-place modification issues
            do j = 1, nvar+1
              do i = 1, m+1
                gammax(i, j) = gammax(i, j) - ee(i) * estimate(j)
              end do
            end do

            !do j = 1, nvar+1
            !  if (j /= t_rr) then
            !    do i = 1, m+1
            !      gammax(i, j) = gammax(i, j) - ee(i) * gammax(k, j)
            !    end do
            !  end if
            !end do

            ! CRITICAL FIX: Save b[k] before updating b vector
            min_k = b(k)
            do i = 1, m+1
              b(i) = b(i) - ee(i) * min_k
            end do

            IB(k) = t !r1(t_rr)

            ! DEBUG: Print after row operations
            if (rd == 1 .and. iter == 3) then
!                write(*,'(A,F15.10,A,F15.10)') '  gammax[k,:] AFTER row operations = ', &
!                    gammax(k,1), ', ', gammax(k,2)
!                write(*,'(A,F15.10,A,F15.10)') '  gammax[m+1,:] AFTER row operations = ', &
!                    gammax(m+1,1), ', ', gammax(m+1,2)
!                write(*,'(A,I2,A,I2)') '  r1 = [', r1(1), ', ', r1(2), ']'
!                write(*,'(A,I2,A,I2)') '  r2 = [', r2(1), ', ', r2(2), ']'
!                write(*,*) '=== END ITERATION 4 DEBUG ==='
            end if


! DEBUG: Print result of pivot
!            if (rd == 1 .and. iter < 10) then
!                write(*,'(A,I4)') '  After pivot: IB(k) now = ', IB(k)
!                write(*,'(A,I2)') '  Iteration ', iter+1, ' complete'
!                write(*,*) ''
!            end if




            iter = iter + 1
            
        end do
        
        ! ============================================================
        ! Extract solution
        ! ============================================================
        it_num_sorted(rd) = iter
        

! DEBUG: Print final basis before extraction
        if (rd == 1) then
!            write(*,*) '=== FORTRAN: Final basis after simplex ==='
!            write(*,*) 'Basis variables (IB array):'
            do i = 1, m
                if (IB(i) <= nvar + 1) then
!                    write(*,'(A,I2,A,I2,A,I2,A,F12.6)') '  IB(', i, ') = ', IB(i), &
!                                                   ' (coefficient var)  b(', i, ') = ', b(i)
                else
!                    write(*,'(A,I2,A,I2,A,I2,A,F12.6)') '  IB(', i, ') = ', IB(i), &
!                                                   ' (residual var)  b(', i, ') = ', b(i)
                end if
            end do
!            write(*,*) ''
        end if



        u = 0.0d0
        v = 0.0d0
        estimate = 0.0d0
        
        ! Extract from basis (matching R pattern)
        do i = 1, m
            if (IB(i) > nvar + 1 .and. IB(i) <= nvar + 1 + m) then
                u(IB(i) - nvar - 1) = b(i)
            else if (IB(i) > nvar + 1 + m) then
                v(IB(i) - nvar - 1 - m) = b(i)
            else if (IB(i) >= 1 .and. IB(i) <= nvar + 1) then
                estimate(IB(i)) = b(i)
            end if
        end do
        
        ! Transform back: β₀ = α₀ + α₁*z, β₁ = α₁ (using SORTED z)
        ! Store in SORTED arrays first
        ll_est_sorted(rd) = estimate(1) + estimate(2) * z_sorted(rd)
        d_ll_est_sorted(rd) = estimate(2)

! DEBUG: Print final solution
        if (rd == 1) then
!            write(*,*) '=== FORTRAN: Final solution ==='
!            write(*,'(A,I4)') 'Total iterations: ', iter
            do i = 1, nvar+1
!                write(*,'(A,I2,A,F12.6)') 'estimate(', i, ') = ', estimate(i)
            end do
!            write(*,'(A,F12.6)') 'll_est = estimate(1) + estimate(2)*z = ', ll_est_sorted(rd)
!            write(*,'(A,F12.6)') 'd_ll_est = estimate(2) = ', d_ll_est_sorted(rd)
!            write(*,*) ''
        end if

        ! Store H matrix in SORTED order
        do i = 1, nvar+1
            j = r1(i) - nvar - 1
            if (j >= 1 .and. j <= m) then
                H_mat_sorted(rd, i) = j
            else
                H_mat_sorted(rd, i) = 0
            end if
        end do

    end do

    ! ============================================================
    ! Unsort results back to original order (matching R track_order=TRUE)
    ! R code: if (track_order) ll_est <- ll_est[order(original_order)]
    ! ============================================================
    do rd = 1, rounds
        k = z_order(rd)  ! Original position of this sorted element
        ll_est(k) = ll_est_sorted(rd)
        d_ll_est(k) = d_ll_est_sorted(rd)
        it_num(k) = it_num_sorted(rd)
        do i = 1, nvar+1
            H_mat(k, i) = H_mat_sorted(rd, i)
        end do
    end do

end subroutine llqr_seq_fortran
