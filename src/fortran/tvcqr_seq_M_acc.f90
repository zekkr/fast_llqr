! Time-varying coefficient quantile regression - Sequential preprocessing
! Faithful translation of tvcqr_seq_ppro.R with focus on correctness for eva_t <= 4
!
! Compile with: R CMD SHLIB tvcqr_seq_corrected.f90
! or: gfortran -shared -fPIC -o tvcqr_seq_M_acc.so tvcqr_seq_M_acc.f90 -llapack -lblas

subroutine tvcqr_seq_ppro_fortran(x, y, m, nvar, tau, h, h_factor, tol, maxit, &
                                   bland_int, Mm_factor, eps, store_residual_int, &
                                   theta_ll_est, it_num, residual_est, &
                                   M_out, n_sub, H_seq, ierr)
    
    implicit none
    
    ! Input arguments
    integer, intent(in) :: m, nvar, maxit, bland_int, store_residual_int
    double precision, intent(in) :: x(m, nvar), y(m), tau, tol, h_factor
    double precision, intent(in) :: Mm_factor, eps
    double precision, intent(inout) :: h
    
    ! Output arguments
    double precision, intent(out) :: theta_ll_est(m, nvar+1)
    integer, intent(out) :: it_num(m)
    double precision, intent(out) :: residual_est(m, m)
    double precision, intent(out) :: M_out
    integer, intent(out) :: n_sub(m)
    integer, intent(out) :: H_seq(m, 2*(nvar+1))
    integer, intent(out) :: ierr
    
    ! Local variables
    double precision :: x_norms(m), mm_thresh, mmm_thresh, M_threshold
    logical :: sl(m), sh(m), not_jl_or_jh(m)
    integer :: idx_not_jl_or_jh(m)
    double precision :: temp_check
    integer :: min_subsample_size, n_potential_S
    double precision :: residual_scale
    double precision :: abs_r(m)  ! NEW - for median of absolute residuals
    
    integer :: ms, ms_org
    double precision :: ws(m+3)
    
    double precision :: glob_wx(2*(nvar+1)), glob_wy
    double precision :: ghib_wx(2*(nvar+1)), ghib_wy
    double precision :: wsl(m), wsh(m)
    
    !double precision :: gammaxs_temp(m+3, 2*(nvar+1))
    !double precision :: bs_temp(m+3)
    !double precision :: gammaxs(m+3, 2*(nvar+1))
    !double precision :: bs(m+3)
    double precision :: xhinv(2*(nvar+1), 2*(nvar+1))
    double precision :: Pxhbarxhinv(m, 2*(nvar+1))
    double precision :: lambda(m+3)

    double precision :: best_pivot_value
    integer :: best_k
    double precision :: ratio
    
    ! Storage for reuse in eva_t >= 3
    double precision :: gammaxs_pos(m, 2*(nvar+1))
    double precision :: gammaxs_neg(m, 2*(nvar+1))
    double precision :: bs_pos(m), bs_neg(m)
    integer :: idx_Hbar_pos(m), idx_Hbar_neg(m)
    integer :: n_Hbar_pos, n_Hbar_neg
    
    double precision :: time_index(m), w(m)
    double precision :: cc(2*(1+nvar) + 2*m)
    !double precision :: A(m, 2*(nvar+1))
    !double precision :: gammax(m+1, 2*(nvar+1))
    double precision :: b(m+1)
    integer :: IB(m+1), IBs(m+3)
    logical :: freevarrow(m+3)
    integer :: r1(2*(nvar+1)), r2(2*(nvar+1))
    double precision :: rr(2, 2*(nvar+1))
    
    double precision :: yy(m+3), ee(m+3), k_vals(m+3)
    double precision :: u(m), v(m), estimate(2*(nvar+1))
    double precision :: r(m), r_prev(m)
    double precision :: pivot_row(2*(nvar+1))
    
    integer :: i, j, k, t, eva_t, iter
    integer :: t_rr, tsep
    double precision :: rrl, min_k, temp_sum
    logical :: bland, store_residual
    double precision :: b_k_original
    logical :: not_optimal, not_new_sl_sh, debug_active
    integer :: bad_signs, n_sl, n_sh, n_sure_signs
    integer :: idpos(m), idneg(m), n_idpos, n_idneg
    integer :: H_indices(2*(nvar+1)), Hbar_indices(m+2)
    integer :: u_in_IBs(m), v_in_IBs(m)
    integer :: idx
    
    ! Variables for eva_t >= 3
    integer :: idx_Hbar_pos2(m), idx_Hbar_neg2(m)
    logical :: matched_rows_pos(m), matched_rows_neg(m)
    integer :: rows_pos(m), rows_neg(m)
    integer :: n_matched_pos, n_matched_neg
    integer :: id_gammaxs_Hbar(m)
    double precision :: temp_vec(2*(nvar+1))

    integer :: preprocessing_attempts
    logical :: valid_H
    double precision :: temp_vec1(2*(nvar+1))
    integer :: ii, kk
    logical :: force_full_sample, accept_subsample
    integer :: empty_pivot_count, max_empty_pivot_retries
    double precision :: res_tol

    
    integer :: n_hbar_rows        ! For the critical lambda/Pxhbarxhinv dimension
    integer :: idx_count, idx_loop
    double precision :: weight_sum

    double precision, allocatable :: A(:,:)           ! Size: m × 2(nvar+1)
    double precision, allocatable :: gammax(:,:)      ! Size: (m+1) × 2(nvar+1)
    double precision, allocatable :: gammaxs_temp(:,:)! Size: (m+3) × 2(nvar+1)
    double precision, allocatable :: bs_temp(:)       ! Size: m+3
    double precision, allocatable :: gammaxs(:,:)     ! Size: (m+3) × 2(nvar+1)
    double precision, allocatable :: bs(:)            ! Size: m+3
    logical :: unbounded_detected
    logical :: simplex_converged
    integer :: total_simplex_iterations
    integer :: total_preprocessing_loops
    integer :: max_iter_at_any_t



 

    ierr = 0
    total_simplex_iterations = 0
    total_preprocessing_loops = 0
    max_iter_at_any_t = 0
    



    ! Start of executable code
    
    ! Convert integer to logical for bland
    bland = (bland_int /= 0)
    store_residual = (store_residual_int /= 0)
    max_empty_pivot_retries = 3
    res_tol = 1.0d-8
    
    ! Allocate the big arrays
    allocate(A(m, 2*(nvar+1)))
    allocate(gammax(m+1, 2*(nvar+1)))
    allocate(gammaxs_temp(m+3, 2*(nvar+1)))
    allocate(bs_temp(m+3))
    allocate(gammaxs(m+3, 2*(nvar+1)))
    allocate(bs(m+3))

    ! Set default bandwidth if h = 0
    if (h <= 0.0d0) then
        h = dble(m)**(-0.2d0) * h_factor
    end if
    
    ! Calculate x norms
    do i = 1, m
        x_norms(i) = 0.0d0
        do j = 1, nvar
            x_norms(i) = x_norms(i) + x(i,j)**2
        end do
        x_norms(i) = sqrt(x_norms(i))
    end do
    
    ! Calculate mm threshold - matching R: mm <- log(m)^{4} * h^2 *max(x.norms)
    mm_thresh = log(dble(m))**4 * h**2 * maxval(x_norms)
    ! Add a minimum threshold to prevent degeneracy
    !mm_thresh = max(mm_thresh, 1.0d0)  ! Ensure mm_thresh is at least 1.0
    ! Initialize time_index = (1:m)/m
    do i = 1, m
        time_index(i) = dble(i) / dble(m)
    end do
    
    ! Build A matrix
    do i = 1, m
        A(i, 1) = 1.0d0
        do j = 1, nvar
            A(i, j+1) = x(i, j)
        end do
        do j = 1, nvar+1
            A(i, nvar+1+j) = A(i, j) * time_index(i)
        end do
    end do
    
    ! Initialize outputs
    it_num = 0
    theta_ll_est = 0.0d0
    if (store_residual) then
        residual_est = 0.0d0
    else
        residual_est(1, 1) = 0.0d0
    end if
    n_sub = 0
    H_seq = 0
    M_out = 0.0d0
    
    ! ============================================
    ! EVA_T = 1: Standard simplex (no preprocessing)
    ! ============================================
    
    ! Calculate weights for t=1
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
    
    ! Build gammax matrix
    do i = 1, m
        do j = 1, 2*(nvar+1)
            gammax(i, j) = A(i, j)
        end do
    end do
    
    ! Initialize b
    do i = 1, m
        b(i) = y(i)
    end do
    b(m+1) = 0.0d0
    
    ! Adjust gammax and b based on sign of y - matching R logic
    do i = 1, m
        if (y(i) < 0.0d0) then
            do j = 1, 2*(nvar+1)
                gammax(i, j) = -gammax(i, j)
            end do
            b(i) = -b(i)
        end if
    end do
    
    ! Initialize IB - matching R: IB <- (y >= 0) * ((1:m) + 2 * (1 + nvar)) + (y < 0) * ((1:m) + 2 * (1 + nvar) + m)
    do i = 1, m
        if (y(i) >= 0.0d0) then
            IB(i) = i + 2*(1+nvar)
        else
            IB(i) = i + 2*(1+nvar) + m
        end if
    end do
    
    ! Add last row to gammax
    do j = 1, 2*(nvar+1)
        gammax(m+1, j) = 0.0d0
        do i = 1, m
            gammax(m+1, j) = gammax(m+1, j) - cc(IB(i)) * gammax(i, j)
        end do
    end do
    
    ! Initialize freevarrow
    freevarrow(1:m+1) = .false.
    freevarrow(m+1) = .true.
    
    ! Initialize r1 and r2
    do i = 1, 2*(nvar+1)
        r1(i) = i
        r2(i) = 0
    end do
    
    n_sub(1) = m
    
    ! Simplex iterations for eva_t = 1
    iter = 0




    do while (iter < maxit)


! Compute reduced costs
        do i = 1, 2*(nvar+1)
            rr(1, i) = gammax(m+1, i)
            
            ! Calculate the index for weight lookup
            j = r1(i) - 2 - 2*nvar
            
            ! R logic: (r1 - 2 - 2 * nvar) * (r1 - 2 - 2 * nvar > 0) + (r1 - 2 - 2 * nvar <= 0)
            if (j > 0) then
                ! Use j when positive
                if (j <= m) then
                    rr(2, i) = w(j) - rr(1, i)
                else
                    rr(2, i) = w(1) - rr(1, i)  ! Bounds check
                end if
            else
                ! When j <= 0, the R expression evaluates to 1
                rr(2, i) = w(1) - rr(1, i)
            end if
            
            ! Apply the (r2 != 0) condition from R
            if (r2(i) == 0) then
                rr(2, i) = 0.0d0
            end if
            
            ! Handle the sign flip for r2 == 0
            if (r2(i) == 0) then
                rr(1, i) = -abs(rr(1, i))
            end if
        end do
        
        ! Check optimality
        rrl = minval(rr)
        if (rrl >= -tol) exit
        
        ! Choose entering variable
        if (bland) then
            t = huge(1)
            t_rr = 0
            
            do i = 1, 2*(nvar+1)
                if (rr(1, i) < -tol) then
                    if (r1(i) < t) then
                        t = r1(i)
                        t_rr = i
                        tsep = 1
                    end if
                end if
            end do
            
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
                        goto 100
                    end if
                end do
            end do
            100 continue
        end if
        
        ! Choose leaving variable
        if (r2(t_rr) /= 0) then
            if (tsep == 1) then
                yy(1:m+1) = gammax(:, t_rr)
            else
                yy(1:m+1) = -gammax(:, t_rr)
            end if
            
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
            
            if (tsep /= 1) then
                j = r1(t_rr) - 2 - 2*nvar
                if (j > 0 .and. j <= m) then
                    yy(m+1) = yy(m+1) + w(j)
                end if
            end if
            
        else
            yy(1:m+1) = gammax(:, t_rr)
            
            if (yy(m+1) < 0.0d0) then
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
        ee(1:m+1) = yy(1:m+1) / yy(k)
        ee(k) = 1.0d0 - 1.0d0 / yy(k)
        
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
                gammax(m+1, t_rr) = 0.0d0
            end if
            r1(t_rr) = IB(k) - m
            r2(t_rr) = IB(k)
        end if
        
        do j = 1, 2*(nvar+1)
            pivot_row(j) = gammax(k, j)
        end do
        
        do j = 1, 2*(nvar+1)
            do i = 1, m+1
                gammax(i, j) = gammax(i, j) - ee(i) * pivot_row(j)
            end do
        end do
        
        b_k_original = b(k)
        do i = 1, m+1
            b(i) = b(i) - ee(i) * b_k_original
        end do
        
        IB(k) = t
        
        iter = iter + 1
    end do
    
    it_num(1) = iter
    
    ! Extract solution for t=1
    u = 0.0d0
    v = 0.0d0
    estimate = 0.0d0
    



    do i = 1, m
        if (IB(i) >= 1 .and. IB(i) <= 2*(nvar+1)) then
            estimate(IB(i)) = b(i)
        else if (IB(i) >= 2*(nvar+1)+1 .and. IB(i) <= 2*(nvar+1)+m) then
            j = IB(i) - 2*(nvar+1)
            u(j) = b(i)
        else if (IB(i) > 2*(nvar+1)+m) then
            j = IB(i) - 2*(nvar+1) - m
            v(j) = b(i)
        end if
    end do
    
    do j = 1, nvar+1
        theta_ll_est(1, j) = estimate(j) + (1.0d0/dble(m)) * estimate(nvar+1+j)
    end do
    
    do i = 1, m
        r_prev(i) = u(i) - v(i)
        if (store_residual) then
            residual_est(1, i) = r_prev(i)
        end if
    end do


    
    ! Extract H observations - only those that correspond to actual observations
    do i = 1, 2*(nvar+1)
        j = r1(i) - 2 - 2*nvar
        if (j > 0 .and. j <= m) then
            H_seq(1, i) = j
        else
            ! This r1 element doesn't correspond to an observation
            ! This can happen when a coefficient variable is in the basis
            H_seq(1, i) = 0  ! Mark as invalid
        end if
    end do
    
    ! Save X(h)^{-1} for next iteration (stored in first 2*(nvar+1) rows of gammax)
    do i = 1, 2*(nvar+1)
        do j = 1, 2*(nvar+1)
            gammaxs_temp(i, j) = gammax(i, j)
        end do
    end do
    
    ! ============================================
    ! EVA_T = 2 to 4: Preprocessing iterations
    ! ============================================
    
    do eva_t = 2, m






        debug_active = .false.
        not_optimal = .true.
        not_new_sl_sh = .true.
        force_full_sample = .false.
        empty_pivot_count = 0
        






        if (debug_active) then
            write(6, *) ''
            write(6, *) '========================================='
            write(6, *) 'FORTRAN: eva_t =', eva_t
            write(6, *) '========================================='
        end if


        ! Initialize counters for this time point
        !n_Hbar_pos = 0
        !n_Hbar_neg = 0

        ! Update weights for current t
        do i = 1, m
            if (abs(dble(eva_t)/dble(m) - time_index(i)) <= h) then
                w(i) = 0.75d0 * (1.0d0 - ((dble(eva_t)/dble(m) - time_index(i))/h)**2)
            else
                w(i) = 0.0d0
            end if
        end do


 

        
        mmm_thresh = mm_thresh
        
        j = 0
        iter = 0  ! Initialize iteration counter here, outside preprocessing loop

        preprocessing_attempts = 0
        do while (not_optimal)
            unbounded_detected = .false.  ! Initialize for each attempt
            simplex_converged = .false.

            preprocessing_attempts = preprocessing_attempts + 1






            total_preprocessing_loops = total_preprocessing_loops + 1
            ! Add safety valve to prevent infinite preprocessing
            if (preprocessing_attempts > max(10, 2*(nvar+1))) then
                force_full_sample = .true.
                not_new_sl_sh = .true.
            end if
            ! Get residuals from previous time
            do i = 1, m
                r(i) = r_prev(i)
            end do
            



            if (debug_active) then
                ! Show residual statistics that will be used for partitioning
                write(6, '(A)', advance='no') 'FORTRAN: Residual stats from t-1: '
                write(6, '(A,F10.6)', advance='no') 'min=', minval(r)
                write(6, '(A,F10.6)', advance='no') ', max=', maxval(r)
                write(6, '(A,F10.6)', advance='no') ', mean=', sum(r)/dble(m)
                write(6, '(A,F10.6)') ', median=', median_value(r, m)  ! You'll need a median function
                
                ! This is crucial - show the exact residuals that will be used
                write(6, '(A)', advance='no') 'FORTRAN: First 10 residuals: '
                do i = 1, min(10, m)
                    write(6, '(F8.4,1X)', advance='no') r(i)
                end do
                write(6, *)
            end if

            ! Calculate threshold and partition observations
            if (not_new_sl_sh) then
                ! Create array of absolute residuals
                do i = 1, m
                    abs_r(i) = abs(r(i))
                end do
                ! Fix 4: Scale threshold to residual magnitude (KEEP THIS!)
                residual_scale = median_value(abs_r, m)
                M_threshold = max(Mm_factor * mmm_thresh * log(log(dble(m))), 0.1d0 * residual_scale)

                
                ! NEW: Ensure minimum subsample size
                min_subsample_size = max(5 * (nvar + 1), ceiling(0.2d0 * dble(m)))

                ! Count potential observations in S
                n_potential_S = 0
                do i = 1, m
                    if (abs(r(i)) <= M_threshold) then
                        n_potential_S = n_potential_S + 1
                    end if
                end do

                ! If too few observations would remain, increase M
                do while (n_potential_S < min_subsample_size .and. M_threshold < maxval(abs(r)))
                    M_threshold = M_threshold * 1.5d0
                    n_potential_S = 0
                    do i = 1, m
                        if (abs(r(i)) <= M_threshold) then
                            n_potential_S = n_potential_S + 1
                        end if
                    end do
                end do
                
                
                do i = 1, m
                    sl(i) = r(i) < -M_threshold
                    sh(i) = r(i) > M_threshold
                end do




            end if






            
            ! Count observations in each partition
            if (force_full_sample) then
                sl = .false.
                sh = .false.
            end if
            n_sl = 0
            n_sh = 0
            ms = 0
            do i = 1, m
                not_jl_or_jh(i) = .not. (sl(i) .or. sh(i))
                if (not_jl_or_jh(i)) then
                    ms = ms + 1
                    idx_not_jl_or_jh(ms) = i
                end if
                if (sl(i)) n_sl = n_sl + 1
                if (sh(i)) n_sh = n_sh + 1
            end do
            



            ! Check if too many observations are classified as sure-sign
            ! After counting ms (uncertain observations)
            if (ms < 2*(nvar+1)) then
                write(6, *) 'WARNING: Only', ms, 'uncertain observations at eva_t=', eva_t
                mmm_thresh = 2.0d0 * mmm_thresh
                not_new_sl_sh = .true.
                if (preprocessing_attempts >= max_empty_pivot_retries) then
                    force_full_sample = .true.
                end if
                cycle
            end if


            ! ===== INSERT NEW CODE HERE =====
            ! CRITICAL: Ensure H observations are ALWAYS in the subsample
            ! This is a theoretical requirement of the algorithm
            if (eva_t > 1) then
                do k = 1, 2*(nvar+1)
                    if (H_seq(eva_t-1, k) > 0) then
                        idx = H_seq(eva_t-1, k)
                        if (sl(idx) .or. sh(idx)) then
                            
                            sl(idx) = .false.
                            sh(idx) = .false.
                            not_jl_or_jh(idx) = .true.
                        end if
                    end if
                end do
                

                ! IMPORTANT: Rebuild not_jl_or_jh array after forcing
                do i = 1, m
                    not_jl_or_jh(i) = .not. (sl(i) .or. sh(i))
                end do

                ! Recount and rebuild index array
                ms = 0
                do i = 1, m
                    if (not_jl_or_jh(i)) then
                        ms = ms + 1
                        idx_not_jl_or_jh(ms) = i
                    end if
                end do
                

                if (debug_active) then
                    write(6, *) 'FORTRAN: Preprocessing partition:'
                    write(6, '(A,F10.6)') '  M_threshold = ', M_threshold
                    write(6, '(A,I4,A,I4,A,I4)') '  n_sl = ', n_sl, ', n_sh = ', n_sh, ', ms = ', ms
                    
        
                end if


                !k = 0
                !do i = 1, m
                !    if (not_jl_or_jh(i)) then
                !        k = k + 1
                !        idx_not_jl_or_jh(k) = i
                !    end if
                !end do
            end if
            ! ===== END OF NEW CODE =====

     


            ms_org = ms  ! Store the original subsample size before adding aggregated obs

            
            ! Extract subsample weights
            do i = 1, ms
                ws(i) = w(idx_not_jl_or_jh(i))
            end do

            
            ! Initialize for eva_t = 2
            if (eva_t == 2) then
                do i = 1, ms
                    do j = 1, 2*(nvar+1)
                        gammaxs_temp(i, j) = A(idx_not_jl_or_jh(i), j)
                    end do
                    bs_temp(i) = y(idx_not_jl_or_jh(i))
                end do




            end if
            
            ! Add aggregated observations if JL is not empty
            if (n_sl > 0) then
                glob_wx = 0.0d0
                glob_wy = 0.0d0
                
                do i = 1, m
                    if (sl(i)) then
                        wsl(i) = w(i)
                        do j = 1, 2*(nvar+1)
                            glob_wx(j) = glob_wx(j) + A(i, j) * wsl(i)
                        end do
                        glob_wy = glob_wy + y(i) * wsl(i)
                    end if
                end do

                
                

                
                do j = 1, 2*(nvar+1)
                    gammaxs_temp(m + 1, j) = glob_wx(j)
                end do
                bs_temp(m + 1) = glob_wy
                ms = ms + 1
                ws(ms) = 1.0d0
              



            end if
            
            ! Add aggregated observations if JH is not empty
            if (n_sh > 0) then
                ghib_wx = 0.0d0
                ghib_wy = 0.0d0
                
                do i = 1, m
                    if (sh(i)) then
                        wsh(i) = w(i)
                        do j = 1, 2*(nvar+1)
                            ghib_wx(j) = ghib_wx(j) + A(i, j) * wsh(i)
                        end do
                        ghib_wy = ghib_wy + y(i) * wsh(i)
                    end if
                end do
                
                

                do j = 1, 2*(nvar+1)
                    gammaxs_temp(m + 2, j) = ghib_wx(j)
                end do
                bs_temp(m + 2) = ghib_wy
                ms = ms + 1
                ws(ms) = 1.0d0



            end if
            
        

            ! Find positive and negative residuals in subsample
            n_idpos = 0
            n_idneg = 0
            do i = 1, ms_org
                if (r(idx_not_jl_or_jh(i)) > 0.0d0) then
                    n_idpos = n_idpos + 1
                    idpos(n_idpos) = i
                else if (r(idx_not_jl_or_jh(i)) < 0.0d0) then
                    n_idneg = n_idneg + 1
                    idneg(n_idneg) = i
                end if
            end do
            
            ! Map H indices to subsample
            do i = 1, 2*(nvar+1)

                H_indices(i) = 0
                if (H_seq(eva_t-1, i) > 0 .and. H_seq(eva_t-1, i) <= m) then
                    ! Valid H observation - try to find it in subsample
                    do j = 1, ms_org
                        if (H_seq(eva_t-1, i) == idx_not_jl_or_jh(j)) then
                            H_indices(i) = j
                            exit
                        end if
                    end do
                end if
            end do

           
 


            if (debug_active) then
                write(6, *) 'FORTRAN: H observation mapping:'
                write(6, '(A)', advance='no') '  H_indices (subsample) = '
                do i = 1, 2*(nvar+1)
                    write(6, '(I3,1X)', advance='no') H_indices(i)
                end do
                write(6, *)
                
                ! Show the actual observations these map to
                write(6, '(A)', advance='no') '  H observations (original) = '
                do i = 1, 2*(nvar+1)
                    if (H_indices(i) > 0 .and. H_indices(i) <= ms_org) then
                        write(6, '(I3,1X)', advance='no') idx_not_jl_or_jh(H_indices(i))
                    else
                        write(6, '(A4,1X)', advance='no') 'NA'
                    end if
                end do
                write(6, *)
            end if
            
            ! Check if any H_indices are missing (important for degenerate cases)
            k = 0  ! Count valid H observations
            do i = 1, 2*(nvar+1)
                if (H_seq(eva_t-1, i) > 0 .and. H_seq(eva_t-1, i) <= m) then
                    if (H_indices(i) == 0) then
                        write(6, *) 'WARNING: H observation', H_seq(eva_t-1, i), 'not in subsample at eva_t=', eva_t
                    else
                        k = k + 1
                    end if
                end if
            end do

            ! Check if we have enough valid H observations
            if (k < 2*(nvar+1)) then
                write(6, *) 'WARNING: Only', k, 'valid H observations out of', 2*(nvar+1), 'at eva_t=', eva_t
                write(6, *) 'Problem may be degenerate or ill-conditioned'
                
                ! For severely degenerate cases, we might need special handling
                if (k < nvar+1) then
                    write(6, *) 'CRITICAL: Too few H observations for a unique solution'
                    ! You might want to use a regularized solution or
                    ! fall back to previous time's estimate here
                end if
            end if

            ! CRITICAL: Reset freevarrow for EVERY preprocessing attempt
            ! This must be done before building the tableau
            freevarrow(1:ms+1) = .false.
            do i = 1, 2*(nvar+1)
                freevarrow(i) = .true.
            end do
            freevarrow(ms+1) = .true.  ! objective row

            
            ! Remove H from idpos and idneg
            do i = 1, 2*(nvar+1)
                if (H_indices(i) > 0) then
                    do j = 1, n_idpos
                        if (idpos(j) == H_indices(i)) then
                            idpos(j) = 0
                        end if
                    end do
                    do j = 1, n_idneg
                        if (idneg(j) == H_indices(i)) then
                            idneg(j) = 0
                        end if
                    end do
                end if
            end do
            




            ! Compact idpos and idneg
            k = 0
            do i = 1, n_idpos
                if (idpos(i) /= 0) then
                    k = k + 1
                    idpos(k) = idpos(i)
                end if
            end do
            n_idpos = k
            
            k = 0
            do i = 1, n_idneg
                if (idneg(i) /= 0) then
                    k = k + 1
                    idneg(k) = idneg(i)
                end if
            end do
            n_idneg = k
 

            ! Check for near-degeneracy
            if (n_idpos + n_idneg <= 1) then
                write(6, *) 'DEGENERATE case at eva_t=', eva_t, ': n_idpos=', n_idpos, ', n_idneg=', n_idneg
                
                ! For degenerate cases, the solution is determined by the H observations
                ! We need to compute xhinv and use it directly
                
                if (eva_t == 2) then
                    ! Build the matrix from H observations
                    do i = 1, 2*(nvar+1)
                        if (H_indices(i) > 0) then
                            do j = 1, 2*(nvar+1)
                                xhinv(i, j) = A(idx_not_jl_or_jh(H_indices(i)), j)
                            end do
                        else
                            ! This shouldn't happen if preprocessing is working correctly
                            write(6, *) 'FATAL: Missing H observation in degenerate case'
                            return
                        end if
                    end do
                    
                    ! Compute the inverse
                    call matrix_inverse_2p(xhinv, xhinv, 2*(nvar+1), ierr)
                    if (ierr /= 0) then
                        write(6,*) 'ERROR: X(h) is singular in degenerate case at eva_t=', eva_t
                        return
                    end if
                else
                    ! For eva_t > 2, xhinv should already be in gammaxs_temp
                    do i = 1, 2*(nvar+1)
                        do j = 1, 2*(nvar+1)
                            xhinv(i, j) = gammaxs_temp(i, j)
                        end do
                    end do
                end if
                
                ! Compute the estimate using xhinv
                do i = 1, 2*(nvar+1)




                    estimate(i) = 0.0d0
                    do j = 1, 2*(nvar+1)
                        if (H_indices(j) > 0) then
                            estimate(i) = estimate(i) + xhinv(i, j) * y(idx_not_jl_or_jh(H_indices(j)))
                        end if
                    end do
                end do
                
         

                ! Store xhinv in gammaxs_temp for next iteration
                do i = 1, 2*(nvar+1)
                    do j = 1, 2*(nvar+1)
                        gammaxs_temp(i, j) = xhinv(i, j)
                    end do
                end do
                


                ! For degenerate case, H observations should be preserved from previous iteration
                ! This matches R's behavior where H <- r1 - 2 - 2 * nvar
                do i = 1, 2*(nvar+1)
                    if (eva_t > 1) then
                        ! In degenerate case, keep previous H observations
                        H_seq(eva_t, i) = H_seq(eva_t-1, i)
                    else
                        ! For eva_t=1, this shouldn't happen in degenerate case
                        H_seq(eva_t, i) = i  ! Default to first observations
                    end if
                end do


                
                ! Skip the tableau setup and simplex iterations
                iter = 0
                goto 300  ! Jump to residual calculation
            end if
            
            ! Set up u_in_IBs and v_in_IBs
            do i = 1, n_idpos
                u_in_IBs(i) = idpos(i) + 2*(nvar + 1)
            end do
            do i = 1, n_idneg
                v_in_IBs(i) = idneg(i) + 2*(nvar + 1) + ms
            end do
            
            ! Set up r1 and r2 for subsample
            do i = 1, 2*(nvar+1)
                r1(i) = H_indices(i) + 2*(nvar + 1)
                r2(i) = r1(i) + ms
            end do
            





            ! Now handle the initialization based on eva_t
            if (eva_t == 2) then
                ! Initialize IBs and freevarrow based on four cases
                if (n_sl > 0 .and. n_sh > 0) then
                    ! Case 1: Both JL and JH non-empty
                    do i = 1, 2*(nvar+1)
                        IBs(i) = i
                    end do
                    do i = 1, n_idpos
                        IBs(2*(nvar+1) + i) = u_in_IBs(i)
                    end do
                    do i = 1, n_idneg
                        IBs(2*(nvar+1) + n_idpos + i) = v_in_IBs(i)
                    end do
                    IBs(2*(nvar+1) + n_idpos + n_idneg + 1) = 2*(nvar+1) + 2*ms - 1
                    IBs(2*(nvar+1) + n_idpos + n_idneg + 2) = 2*(nvar+1) + ms
                    
                    do i = 1, n_idpos
                        Hbar_indices(i) = idpos(i)
                    end do
                    do i = 1, n_idneg
                        Hbar_indices(n_idpos + i) = idneg(i)
                    end do
                    Hbar_indices(n_idpos + n_idneg + 1) = ms - 1
                    Hbar_indices(n_idpos + n_idneg + 2) = ms
                    
                    do i = 1, ms-2
                        do j = 1, 2*(nvar+1)
                            gammaxs(i, j) = gammaxs_temp(i, j)
                        end do
                        bs(i) = bs_temp(i)
                    end do
                    do j = 1, 2*(nvar+1)
                        gammaxs(ms-1, j) = gammaxs_temp(m+1, j)
                        gammaxs(ms, j) = gammaxs_temp(m+2, j)
                    end do
                    bs(ms-1) = bs_temp(m+1)
                    bs(ms) = bs_temp(m+2)
                    
                    do i = 1, n_idpos
                        lambda(i) = tau * ws(idpos(i))
                    end do
                    do i = 1, n_idneg
                        lambda(n_idpos + i) = (1.0d0 - tau) * ws(idneg(i))
                    end do
                    lambda(n_idpos + n_idneg + 1) = 1.0d0 - tau
                    lambda(n_idpos + n_idneg + 2) = tau

                    
                    freevarrow(1:ms+1) = .false.
                    do i = 1, 2*(nvar+1)
                        freevarrow(i) = .true.
                    end do
                    freevarrow(ms-1) = .true.
                    freevarrow(ms) = .true.
                    freevarrow(ms+1) = .true.
                    




                else if (n_sl > 0) then
                    ! Case 2: Only JL non-empty
                    do i = 1, 2*(nvar+1)
                        IBs(i) = i
                    end do
                    do i = 1, n_idpos
                        IBs(2*(nvar+1) + i) = u_in_IBs(i)
                    end do
                    do i = 1, n_idneg
                        IBs(2*(nvar+1) + n_idpos + i) = v_in_IBs(i)
                    end do
                    IBs(2*(nvar+1) + n_idpos + n_idneg + 1) = 2*(nvar+1) + 2*ms
                    
                    do i = 1, n_idpos
                        Hbar_indices(i) = idpos(i)
                    end do
                    do i = 1, n_idneg
                        Hbar_indices(n_idpos + i) = idneg(i)
                    end do
                    Hbar_indices(n_idpos + n_idneg + 1) = ms
                    
                    do i = 1, ms-1
                        do j = 1, 2*(nvar+1)
                            gammaxs(i, j) = gammaxs_temp(i, j)
                        end do
                        bs(i) = bs_temp(i)
                    end do
                    do j = 1, 2*(nvar+1)
                        gammaxs(ms, j) = gammaxs_temp(m+1, j)
                    end do
                    bs(ms) = bs_temp(m+1)
                    
                    do i = 1, n_idpos
                        lambda(i) = tau * ws(idpos(i))
                    end do
                    do i = 1, n_idneg
                        lambda(n_idpos + i) = (1.0d0 - tau) * ws(idneg(i))
                    end do
                    lambda(n_idpos + n_idneg + 1) = 1.0d0 - tau
                    
                    freevarrow(1:ms+1) = .false.
                    do i = 1, 2*(nvar+1)
                        freevarrow(i) = .true.
                    end do
                    freevarrow(ms) = .true.
                    freevarrow(ms+1) = .true.
                    
                else if (n_sh > 0) then
                    ! Case 3: Only JH non-empty
                    do i = 1, 2*(nvar+1)
                        IBs(i) = i
                    end do
                    do i = 1, n_idpos
                        IBs(2*(nvar+1) + i) = u_in_IBs(i)
                    end do
                    do i = 1, n_idneg
                        IBs(2*(nvar+1) + n_idpos + i) = v_in_IBs(i)
                    end do
                    IBs(2*(nvar+1) + n_idpos + n_idneg + 1) = 2*(nvar+1) + ms
                    
                    do i = 1, n_idpos
                        Hbar_indices(i) = idpos(i)
                    end do
                    do i = 1, n_idneg
                        Hbar_indices(n_idpos + i) = idneg(i)
                    end do
                    Hbar_indices(n_idpos + n_idneg + 1) = ms
                    
                    do i = 1, ms-1
                        do j = 1, 2*(nvar+1)
                            gammaxs(i, j) = gammaxs_temp(i, j)
                        end do
                        bs(i) = bs_temp(i)
                    end do
                    do j = 1, 2*(nvar+1)
                        gammaxs(ms, j) = gammaxs_temp(m+2, j)
                    end do
                    bs(ms) = bs_temp(m+2)
                    
                    do i = 1, n_idpos
                        lambda(i) = tau * ws(idpos(i))
                    end do
                    do i = 1, n_idneg
                        lambda(n_idpos + i) = (1.0d0 - tau) * ws(idneg(i))
                    end do
                    lambda(n_idpos + n_idneg + 1) = tau
                    
                    freevarrow(1:ms+1) = .false.
                    do i = 1, 2*(nvar+1)
                        freevarrow(i) = .true.
                    end do
                    freevarrow(ms) = .true.
                    freevarrow(ms+1) = .true.


    
                    
                else
                    ! Case 4: Both JL and JH empty
                    do i = 1, 2*(nvar+1)
                        IBs(i) = i
                    end do
                    do i = 1, n_idpos
                        IBs(2*(nvar+1) + i) = u_in_IBs(i)
                    end do
                    do i = 1, n_idneg
                        IBs(2*(nvar+1) + n_idpos + i) = v_in_IBs(i)
                    end do
                    
                    do i = 1, n_idpos
                        Hbar_indices(i) = idpos(i)
                    end do
                    do i = 1, n_idneg
                        Hbar_indices(n_idpos + i) = idneg(i)
                    end do
                    
                    do i = 1, ms
                        do j = 1, 2*(nvar+1)
                            gammaxs(i, j) = gammaxs_temp(i, j)
                        end do
                        bs(i) = bs_temp(i)
                    end do
                    
                    do i = 1, n_idpos
                        lambda(i) = tau * ws(idpos(i))
                    end do
                    do i = 1, n_idneg
                        lambda(n_idpos + i) = (1.0d0 - tau) * ws(idneg(i))
                    end do
                    
                    freevarrow(1:ms+1) = .false.
                    do i = 1, 2*(nvar+1)
                        freevarrow(i) = .true.
                    end do
                    freevarrow(ms+1) = .true.
                end if
                
                ! Compute initial tableau for t=2
                call matrix_inverse_2p(gammaxs(H_indices, :), xhinv, 2*(nvar+1), ierr)
                
                if (ierr /= 0) then
                    write(6,*) 'X(h) is singular at eva_t = ', eva_t, ', info = ', ierr
                    return
                end if
                
                ! Compute Pxhbarxhinv = P @ gammaxs[Hbar, :] @ xhinv
                ! For positive residuals (u variables), P has +1
                do i = 1, n_idpos
                    do j = 1, 2*(nvar+1)
                        Pxhbarxhinv(i, j) = 0.0d0
                        do k = 1, 2*(nvar+1)
                            Pxhbarxhinv(i, j) = Pxhbarxhinv(i, j) + gammaxs(Hbar_indices(i), k) * xhinv(k, j)
                        end do
                    end do
                end do
                
                ! For negative residuals (v variables), P has -1
                do i = 1, n_idneg
                    do j = 1, 2*(nvar+1)
                        Pxhbarxhinv(n_idpos + i, j) = 0.0d0
                        do k = 1, 2*(nvar+1)
                            Pxhbarxhinv(n_idpos + i, j) = Pxhbarxhinv(n_idpos + i, j) - &
                                                         gammaxs(Hbar_indices(n_idpos + i), k) * xhinv(k, j)
                        end do
                    end do
                end do
                
                ! Handle JL and JH based on case
                if (n_sl > 0 .and. n_sh > 0) then
                    ! For v_L (JL aggregated), P has -1
                    do j = 1, 2*(nvar+1)
                        Pxhbarxhinv(n_idpos + n_idneg + 1, j) = 0.0d0
                        do k = 1, 2*(nvar+1)
                            Pxhbarxhinv(n_idpos + n_idneg + 1, j) = Pxhbarxhinv(n_idpos + n_idneg + 1, j) - &
                                                                   gammaxs(ms-1, k) * xhinv(k, j)
                        end do
                    end do
                    
                    ! For u_H (JH aggregated), P has +1
                    do j = 1, 2*(nvar+1)
                        Pxhbarxhinv(n_idpos + n_idneg + 2, j) = 0.0d0
                        do k = 1, 2*(nvar+1)
                            Pxhbarxhinv(n_idpos + n_idneg + 2, j) = Pxhbarxhinv(n_idpos + n_idneg + 2, j) + &
                                                                   gammaxs(ms, k) * xhinv(k, j)
                        end do
                    end do
                else if (n_sl > 0) then
                    ! Only v_L
                    do j = 1, 2*(nvar+1)
                        Pxhbarxhinv(n_idpos + n_idneg + 1, j) = 0.0d0
                        do k = 1, 2*(nvar+1)
                            Pxhbarxhinv(n_idpos + n_idneg + 1, j) = Pxhbarxhinv(n_idpos + n_idneg + 1, j) - &
                                                                   gammaxs(ms, k) * xhinv(k, j)
                        end do
                    end do
                else if (n_sh > 0) then
                    ! Only u_H
                    do j = 1, 2*(nvar+1)
                        Pxhbarxhinv(n_idpos + n_idneg + 1, j) = 0.0d0
                        do k = 1, 2*(nvar+1)
                            Pxhbarxhinv(n_idpos + n_idneg + 1, j) = Pxhbarxhinv(n_idpos + n_idneg + 1, j) + &
                                                                   gammaxs(ms, k) * xhinv(k, j)
                        end do
                    end do
                end if
                
                ! Build the tableau
                ! First rows are xhinv
                do i = 1, 2*(nvar+1)
                    do j = 1, 2*(nvar+1)
                        gammaxs(i, j) = xhinv(i, j)
                    end do
                end do



                
                ! Next rows are -Pxhbarxhinv
                k = n_idpos + n_idneg
                if (n_sl > 0) k = k + 1
                if (n_sh > 0) k = k + 1
                
                do i = 1, k
                    do j = 1, 2*(nvar+1)
                        gammaxs(2*(nvar+1) + i, j) = -Pxhbarxhinv(i, j)
                    end do
                end do
                
                ! Last row is the objective function row
                do j = 1, 2*(nvar+1)
                    if (H_indices(j) < 1 .or. H_indices(j) > ms) then
                        write(6, *) 'ERROR: Invalid H_indices(', j, ')=', H_indices(j), ' at eva_t=', eva_t
                        write(6, *) '  Valid range is 1 to', ms
                        gammaxs(ms + 1, j) = 0.0d0  ! Safety fallback
                    else
                        gammaxs(ms + 1, j) = tau * ws(H_indices(j))
                    end if
                    do i = 1, k
                        gammaxs(ms + 1, j) = gammaxs(ms + 1, j) + lambda(i) * Pxhbarxhinv(i, j)
                    end do
                end do
                
                ! Build bs vector
                ! First part: xhinv @ bs[H]


                do i = 1, 2*(nvar+1)
                    bs(i) = 0.0d0
                    do j = 1, 2*(nvar+1)
                        bs(i) = bs(i) + xhinv(i, j) * bs_temp(H_indices(j))
                    end do
                end do
                
                ! Second part: -Pxhbarxhinv @ bs[H] + P @ bs[Hbar]
                ! For positive residuals
                do i = 1, n_idpos
                    bs(2*(nvar+1) + i) = bs_temp(Hbar_indices(i))
                    do j = 1, 2*(nvar+1)
                        bs(2*(nvar+1) + i) = bs(2*(nvar+1) + i) - Pxhbarxhinv(i, j) * bs_temp(H_indices(j))
                    end do
                end do
                
                ! For negative residuals
                do i = 1, n_idneg
                    bs(2*(nvar+1) + n_idpos + i) = -bs_temp(Hbar_indices(n_idpos + i))
                    do j = 1, 2*(nvar+1)
                        bs(2*(nvar+1) + n_idpos + i) = bs(2*(nvar+1) + n_idpos + i) - &
                                                      Pxhbarxhinv(n_idpos + i, j) * bs_temp(H_indices(j))
                    end do
                end do
                
                ! Handle JL and JH
                if (n_sl > 0 .and. n_sh > 0) then
                    bs(ms-1) = -bs_temp(m+1)
                    do j = 1, 2*(nvar+1)
                        bs(ms-1) = bs(ms-1) - Pxhbarxhinv(n_idpos + n_idneg + 1, j) * bs_temp(H_indices(j))
                    end do
                    
                    bs(ms) = bs_temp(m+2)
                    do j = 1, 2*(nvar+1)
                        bs(ms) = bs(ms) - Pxhbarxhinv(n_idpos + n_idneg + 2, j) * bs_temp(H_indices(j))
                    end do
                else if (n_sl > 0) then
                    bs(ms) = -bs_temp(m+1)
                    do j = 1, 2*(nvar+1)
                        bs(ms) = bs(ms) - Pxhbarxhinv(n_idpos + n_idneg + 1, j) * bs_temp(H_indices(j))
                    end do
                else if (n_sh > 0) then
                    bs(ms) = bs_temp(m+2)
                    do j = 1, 2*(nvar+1)
                        bs(ms) = bs(ms) - Pxhbarxhinv(n_idpos + n_idneg + 1, j) * bs_temp(H_indices(j))
                    end do
                end if
                
                bs(ms+1) = 0.0d0



                
            else
                ! eva_t >= 3: Reuse computations from previous iteration
           




                ! xhinv is already stored in gammaxs_temp from previous iteration



                do i = 1, 2*(nvar+1)
                    do j = 1, 2*(nvar+1)
                        xhinv(i, j) = gammaxs_temp(i, j)
                    end do
                end do
                



                ! Map current idpos to original indices
                do i = 1, n_idpos
                    idx_Hbar_pos2(i) = idx_not_jl_or_jh(idpos(i))
                end do
                
                ! Check which positive residuals can be reused
                n_matched_pos = 0
                do i = 1, n_idpos
                    matched_rows_pos(i) = .false.
                    do j = 1, n_Hbar_pos
                        if (idx_Hbar_pos2(i) == idx_Hbar_pos(j)) then
                            matched_rows_pos(i) = .true.
                            n_matched_pos = n_matched_pos + 1
                            rows_pos(n_matched_pos) = j
                            exit
                        end if
                    end do
                end do

       


                
                ! Initialize positive residual rows in gammaxs_temp
                ! First, copy the rows that can be reused
                k = 0
                do i = 1, n_idpos
                    if (matched_rows_pos(i)) then
                        k = k + 1
                        do j = 1, 2*(nvar+1)
                            gammaxs_temp(2*(nvar+1) + i, j) = gammaxs_pos(rows_pos(k), j)
                        end do
                        bs_temp(2*(nvar+1) + i) = bs_pos(rows_pos(k))
                    end if
                end do
                
                ! For new positive residuals, compute fresh
                do i = 1, n_idpos
                    if (.not. matched_rows_pos(i)) then
                        ! Compute -P @ X(hbar) @ X(h)^{-1}
                        do j = 1, 2*(nvar+1)
                            temp_vec(j) = 0.0d0
                            do k = 1, 2*(nvar+1)
                                temp_vec(j) = temp_vec(j) + A(idx_Hbar_pos2(i), k) * xhinv(k, j)
                            end do
                        end do
                        
                        do j = 1, 2*(nvar+1)
                            gammaxs_temp(2*(nvar+1) + i, j) = -temp_vec(j)
                        end do
                        
                        ! Compute -P @ X(hbar) @ X(h)^{-1} @ y(h) + P @ y(hbar)
                        bs_temp(2*(nvar+1) + i) = y(idx_Hbar_pos2(i))
                        do j = 1, 2*(nvar+1)
                            bs_temp(2*(nvar+1) + i) = bs_temp(2*(nvar+1) + i) - &
                                                    temp_vec(j) * y(idx_not_jl_or_jh(H_indices(j)))
                        end do
                    end if
                end do
                
                ! Similarly for negative residuals
                do i = 1, n_idneg
                    idx_Hbar_neg2(i) = idx_not_jl_or_jh(idneg(i))
                end do
                
                n_matched_neg = 0
                do i = 1, n_idneg
                    matched_rows_neg(i) = .false.
                    do j = 1, n_Hbar_neg
                        if (idx_Hbar_neg2(i) == idx_Hbar_neg(j)) then
                            matched_rows_neg(i) = .true.
                            n_matched_neg = n_matched_neg + 1
                            rows_neg(n_matched_neg) = j
                            exit
                        end if
                    end do
                end do
                
                








                ! Copy reusable negative residual rows
                k = 0
                do i = 1, n_idneg
                    if (matched_rows_neg(i)) then
                        k = k + 1

        

                        do j = 1, 2*(nvar+1)
                            gammaxs_temp(2*(nvar+1) + n_idpos + i, j) = gammaxs_neg(rows_neg(k), j)
                        end do
                        bs_temp(2*(nvar+1) + n_idpos + i) = bs_neg(rows_neg(k))
                    end if
                end do
                



                ! Compute fresh for new negative residuals
                do i = 1, n_idneg
                    if (.not. matched_rows_neg(i)) then
                        ! Compute P @ X(hbar) @ X(h)^{-1} (note: P is -1 for v variables)
                        do j = 1, 2*(nvar+1)
                            temp_vec(j) = 0.0d0
                            do k = 1, 2*(nvar+1)
                                temp_vec(j) = temp_vec(j) - A(idx_Hbar_neg2(i), k) * xhinv(k, j)
                            end do
                        end do
                        
                        do j = 1, 2*(nvar+1)
                            gammaxs_temp(2*(nvar+1) + n_idpos + i, j) = -temp_vec(j)
                        end do
                        
                        bs_temp(2*(nvar+1) + n_idpos + i) = -y(idx_Hbar_neg2(i))
                        do j = 1, 2*(nvar+1)
                            bs_temp(2*(nvar+1) + n_idpos + i) = bs_temp(2*(nvar+1) + n_idpos + i) - &
                                                              temp_vec(j) * y(idx_not_jl_or_jh(H_indices(j)))
                        end do
                    end if
                end do
                


                ! Handle JL and JH aggregated rows
                if (n_sl > 0 .and. n_sh > 0) then



                    ! For JL: compute X_L^T @ X(h)^{-1}
                    do j = 1, 2*(nvar+1)
                        temp_vec(j) = 0.0d0
                        do k = 1, 2*(nvar+1)
                            temp_vec(j) = temp_vec(j) + gammaxs_temp(m+1, k) * xhinv(k, j)
                        end do
                        
                    end do

                    do j = 1, 2*(nvar+1)
                        gammaxs_temp(m+1, j) = temp_vec(j)
                    end do

                  

                    ! Compute corresponding bs entry
                    temp_sum = 0.0d0
                    do j = 1, 2*(nvar+1)
                        temp_sum = temp_sum + gammaxs_temp(m+1, j) * y(idx_not_jl_or_jh(H_indices(j)))
                    end do
                    bs_temp(m+1) = temp_sum - bs_temp(m+1)
                    

                    
                    ! For JH: compute -X_H^T @ X(h)^{-1}
                    do j = 1, 2*(nvar+1)
                        temp_vec(j) = 0.0d0
                        do k = 1, 2*(nvar+1)
                            temp_vec(j) = temp_vec(j) - gammaxs_temp(m+2, k) * xhinv(k, j)
                        end do
                        !gammaxs_temp(m+2, j) = temp_vec(j)
                    end do
                    
                    ! Now copy the result back AFTER computing all elements
                    do j = 1, 2*(nvar+1)
                        gammaxs_temp(m+2, j) = temp_vec(j)
                    end do

                    temp_sum = 0.0d0
                    do j = 1, 2*(nvar+1)
                        temp_sum = temp_sum + gammaxs_temp(m+2, j) * y(idx_not_jl_or_jh(H_indices(j)))
                    end do
                    bs_temp(m+2) = temp_sum + bs_temp(m+2)
                    
                    ! Copy to gammaxs
                    do i = 1, ms-2
                        do j = 1, 2*(nvar+1)
                            gammaxs(i, j) = gammaxs_temp(i, j)
                        end do
                        bs(i) = bs_temp(i)
                    end do
                    do j = 1, 2*(nvar+1)
                        gammaxs(ms-1, j) = gammaxs_temp(m+1, j)
                        gammaxs(ms, j) = gammaxs_temp(m+2, j)
                    end do
                    bs(ms-1) = bs_temp(m+1)
                    bs(ms) = bs_temp(m+2)
                    
                    do i = 1, n_idpos
                        lambda(i) = tau * ws(idpos(i))
                    end do
                    do i = 1, n_idneg
                        lambda(n_idpos + i) = (1.0d0 - tau) * ws(idneg(i))
                    end do
                    lambda(n_idpos + n_idneg + 1) = 1.0d0 - tau
                    lambda(n_idpos + n_idneg + 2) = tau
                    
                else if (n_sl > 0) then
                    ! Similar for only JL case
                    do j = 1, 2*(nvar+1)
                        temp_vec(j) = 0.0d0
                        do k = 1, 2*(nvar+1)
                            temp_vec(j) = temp_vec(j) + gammaxs_temp(m+1, k) * xhinv(k, j)
                        end do
                        
                    end do
                    
                    do j = 1, 2*(nvar+1)
                        gammaxs_temp(m+1, j) = temp_vec(j)
                    end do





                    
                    temp_sum = 0.0d0
                    do j = 1, 2*(nvar+1)
                        temp_sum = temp_sum + gammaxs_temp(m+1, j) * y(idx_not_jl_or_jh(H_indices(j)))
                    end do
                    bs_temp(m+1) = temp_sum - bs_temp(m+1)
                    
                    do i = 1, ms-1
                        do j = 1, 2*(nvar+1)
                            gammaxs(i, j) = gammaxs_temp(i, j)
                        end do
                        bs(i) = bs_temp(i)
                    end do
                    do j = 1, 2*(nvar+1)
                        gammaxs(ms, j) = gammaxs_temp(m+1, j)
                    end do
                    bs(ms) = bs_temp(m+1)
                    



                    do i = 1, n_idpos
                        lambda(i) = tau * ws(idpos(i))
                    end do
                    do i = 1, n_idneg
                        lambda(n_idpos + i) = (1.0d0 - tau) * ws(idneg(i))
                    end do
                    lambda(n_idpos + n_idneg + 1) = 1.0d0 - tau
                    


                else if (n_sh > 0) then





                    ! Similar for only JH case
                    do j = 1, 2*(nvar+1)
                        temp_vec(j) = 0.0d0
                        do k = 1, 2*(nvar+1)
                            temp_vec(j) = temp_vec(j) - gammaxs_temp(m+2, k) * xhinv(k, j)
                        end do
                        !gammaxs_temp(m+2, j) = temp_vec(j)
                    end do

                    ! Now copy the result back AFTER computing all elements
                    do j = 1, 2*(nvar+1)
                        gammaxs_temp(m+2, j) = temp_vec(j)
                    end do




                    
                    temp_sum = 0.0d0
                    do j = 1, 2*(nvar+1)
                        temp_sum = temp_sum + gammaxs_temp(m+2, j) * y(idx_not_jl_or_jh(H_indices(j)))
                    end do
                    bs_temp(m+2) = temp_sum + bs_temp(m+2)
                    
                    do i = 1, ms-1
                        do j = 1, 2*(nvar+1)
                            gammaxs(i, j) = gammaxs_temp(i, j)
                        end do
                        bs(i) = bs_temp(i)
                    end do
                    do j = 1, 2*(nvar+1)
                        gammaxs(ms, j) = gammaxs_temp(m+2, j)
                    end do
                    bs(ms) = bs_temp(m+2)
                    
                    do i = 1, n_idpos
                        lambda(i) = tau * ws(idpos(i))
                    end do
                    do i = 1, n_idneg
                        lambda(n_idpos + i) = (1.0d0 - tau) * ws(idneg(i))
                    end do
                    lambda(n_idpos + n_idneg + 1) = tau
                    
                else
                    ! Both empty
                    do i = 1, ms
                        do j = 1, 2*(nvar+1)
                            gammaxs(i, j) = gammaxs_temp(i, j)
                        end do
                        bs(i) = bs_temp(i)
                    end do
                    
                    do i = 1, n_idpos
                        lambda(i) = tau * ws(idpos(i))
                    end do
                    do i = 1, n_idneg
                        lambda(n_idpos + i) = (1.0d0 - tau) * ws(idneg(i))
                    end do
                end if
                
                ! Set up IBs and freevarrow (same logic as eva_t = 2)
                if (n_sl > 0 .and. n_sh > 0) then
                    do i = 1, 2*(nvar+1)
                        IBs(i) = i
                    end do
                    do i = 1, n_idpos
                        IBs(2*(nvar+1) + i) = u_in_IBs(i)
                    end do
                    do i = 1, n_idneg
                        IBs(2*(nvar+1) + n_idpos + i) = v_in_IBs(i)
                    end do
                    IBs(2*(nvar+1) + n_idpos + n_idneg + 1) = 2*(nvar+1) + 2*ms - 1
                    IBs(2*(nvar+1) + n_idpos + n_idneg + 2) = 2*(nvar+1) + ms
                    
                    do i = 1, n_idpos
                        Hbar_indices(i) = idpos(i)
                    end do
                    do i = 1, n_idneg
                        Hbar_indices(n_idpos + i) = idneg(i)
                    end do
                    Hbar_indices(n_idpos + n_idneg + 1) = ms - 1
                    Hbar_indices(n_idpos + n_idneg + 2) = ms
                    
                    freevarrow(1:ms+1) = .false.
                    do i = 1, 2*(nvar+1)
                        freevarrow(i) = .true.
                    end do
                    freevarrow(ms-1) = .true.
                    freevarrow(ms) = .true.
                    freevarrow(ms+1) = .true.




                    
                else if (n_sl > 0) then
                    do i = 1, 2*(nvar+1)
                        IBs(i) = i
                    end do
                    do i = 1, n_idpos
                        IBs(2*(nvar+1) + i) = u_in_IBs(i)
                    end do
                    do i = 1, n_idneg
                        IBs(2*(nvar+1) + n_idpos + i) = v_in_IBs(i)
                    end do
                    IBs(2*(nvar+1) + n_idpos + n_idneg + 1) = 2*(nvar+1) + 2*ms
                    
                    do i = 1, n_idpos
                        Hbar_indices(i) = idpos(i)
                    end do
                    do i = 1, n_idneg
                        Hbar_indices(n_idpos + i) = idneg(i)
                    end do
                    Hbar_indices(n_idpos + n_idneg + 1) = ms
                    
                    freevarrow(1:ms+1) = .false.
                    do i = 1, 2*(nvar+1)
                        freevarrow(i) = .true.
                    end do
                    freevarrow(ms) = .true.
                    freevarrow(ms+1) = .true.
                    
                else if (n_sh > 0) then
                    do i = 1, 2*(nvar+1)
                        IBs(i) = i
                    end do
                    do i = 1, n_idpos
                        IBs(2*(nvar+1) + i) = u_in_IBs(i)
                    end do
                    do i = 1, n_idneg
                        IBs(2*(nvar+1) + n_idpos + i) = v_in_IBs(i)
                    end do
                    IBs(2*(nvar+1) + n_idpos + n_idneg + 1) = 2*(nvar+1) + ms
                    
                    do i = 1, n_idpos
                        Hbar_indices(i) = idpos(i)
                    end do
                    do i = 1, n_idneg
                        Hbar_indices(n_idpos + i) = idneg(i)
                    end do
                    Hbar_indices(n_idpos + n_idneg + 1) = ms
                    
                    freevarrow(1:ms+1) = .false.
                    do i = 1, 2*(nvar+1)
                        freevarrow(i) = .true.
                    end do
                    freevarrow(ms) = .true.
                    freevarrow(ms+1) = .true.
                    
                else
                    do i = 1, 2*(nvar+1)
                        IBs(i) = i
                    end do
                    do i = 1, n_idpos
                        IBs(2*(nvar+1) + i) = u_in_IBs(i)
                    end do
                    do i = 1, n_idneg
                        IBs(2*(nvar+1) + n_idpos + i) = v_in_IBs(i)
                    end do
                    
                    do i = 1, n_idpos
                        Hbar_indices(i) = idpos(i)
                    end do
                    do i = 1, n_idneg
                        Hbar_indices(n_idpos + i) = idneg(i)
                    end do
                    
                    freevarrow(1:ms+1) = .false.
                    do i = 1, 2*(nvar+1)
                        freevarrow(i) = .true.
                    end do
                    freevarrow(ms+1) = .true.
                end if
                
                ! Extract Pxhbarxhinv from the stored -Pxhbarxhinv in gammaxs
                ! R: Pxhbarxhinv <- - gammaxs[(2*(nvar+1)+1):ms,]
                n_hbar_rows = ms - 2*(nvar+1)  ! This is the number of Hbar rows




                do i = 1, n_hbar_rows  ! Use the descriptive name
                    do j = 1, 2*(nvar+1)
                        Pxhbarxhinv(i, j) = -gammaxs(2*(nvar+1) + i, j)
                    end do
                end do

                





                ! Build the last row
                do j = 1, 2*(nvar+1)

                    ! Add bounds checking for H_indices
                    if (H_indices(j) > 0 .and. H_indices(j) <= ms_org) then
                        gammaxs(ms + 1, j) = tau * ws(H_indices(j))
                    else
                        write(6, *) 'ERROR: Invalid H_indices(', j, ')=', H_indices(j), ' at eva_t=', eva_t
                        write(6, *) '  Valid range is 1 to', ms_org
                        ! This should not happen if H observations are correctly forced into subsample
                        stop 'Invalid H_indices in objective row calculation'
                    end if



   


                    do i = 1, n_hbar_rows  ! Use the descriptive name
                        gammaxs(ms + 1, j) = gammaxs(ms + 1, j) + lambda(i) * Pxhbarxhinv(i, j)
                    end do
                end do

  






                bs(ms+1) = 0.0d0
                


            end if  ! End of eva_t == 2 vs eva_t >= 3
            

   

            ! Check if we have enough observations for a well-posed problem
            if (ms < 2*(nvar+1)) then
                write(6, *) 'WARNING: Degenerate problem at eva_t=', eva_t
                write(6, *) 'Only', ms, 'rows for', 2*(nvar+1), 'parameters'
                write(6, *) 'Skipping simplex and using H-observation solution'
                
                ! Use the H observations directly
                do i = 1, 2*(nvar+1)
                    estimate(i) = 0.0d0
                    do j = 1, 2*(nvar+1)
                        if (H_indices(j) > 0) then
                            estimate(i) = estimate(i) + xhinv(i, j) * y(idx_not_jl_or_jh(H_indices(j)))
                        end if
                    end do
                end do
                
                ! Skip simplex entirely
                iter = 0
                goto 300  ! Jump to residual calculation
            end if



            ! Invariant checks before simplex
            if (ms < 2*(nvar+1)) then
                write(6, *) 'ERROR: ms too small at eva_t=', eva_t, ', ms=', ms
                write(6, *) 'This will cause degenerate simplex behavior'
            end if

            ! Check that all beta columns are basic
            if (.not. all(freevarrow(1:2*(nvar+1)))) then
                write(6, *) 'ERROR: Beta not basic at eva_t=', eva_t
                do i = 1, 2*(nvar+1)
                    if (.not. freevarrow(i)) then
                        write(6, *) '  Beta column', i, 'is not basic'
                    end if
                end do
            end if





            ! Simplex iterations for the reduced problem
            !iter = 0



            do while (iter < maxit)
                total_preprocessing_loops = total_preprocessing_loops + 1
                ! Add a safety check to prevent infinite preprocessing:
                if (preprocessing_attempts > 100) then
                    write(6, *) 'ERROR: Preprocessing stuck in infinite loop at eva_t=', eva_t
                    write(6, *) 'M_threshold:', M_threshold
                    write(6, *) 'bad_signs:', bad_signs
                    write(6, *) 'ms:', ms
                    stop 'Preprocessing not converging'
                end if
                ! Step 2: Compute reduced costs - CRITICAL SECTION
                ! This must match R: rr[2, ] <- (ws[r1 - 2 - 2 * nvar] - rr[1, ])
                do i = 1, 2*(nvar+1)
                      rr(1, i) = gammaxs(ms + 1, i)
                      


                      ! Get the subsample index
                      j = r1(i) - 2 - 2*nvar
                      
                      ! Use subsample weights ws, matching R's behavior
                      if (j > 0 .and. j <= ms) then
                          rr(2, i) = ws(j) - rr(1, i)
                      else if (j <= 0) then
                          ! R would use ws[1] when index is non-positive
                          rr(2, i) = ws(1) - rr(1, i)
                      else
                        ! j > ms - use ws(1) as fallback, matching R's behavior
                          rr(2, i) = ws(1) - rr(1, i)
                      end if
                end do




                if (eva_t == 1) then
                    do i = 1, 2*(nvar+1)
                        if (r2(i) == 0) then
                            rr(2, i) = 0.0d0
                            rr(1, i) = -abs(rr(1, i))
                        end if
                    end do
                end if





                ! Check optimality
                rrl = minval(rr)
                

                ! Verify the mask is working correctly
                if (rrl >= 0.0d0) then
                    do i = 1, 2*(nvar+1)
                        if (r2(i) == 0 .and. rr(2, i) /= 0.0d0) then
                            write(6, *) 'ERROR: r2 mask failed at eva_t=', eva_t
                            write(6, *) 'Column', i, 'has r2=0 but rr(2,i)=', rr(2, i)
                            stop 'Mask implementation error'
                        end if
                    end do
                end if

                ! Find location of minimum
                call minloc2d(rr, rrl, tsep, t_rr)




                ! ADD SAFETY VALVE HERE:
                if (iter > 50000) then



                    write(6, *) 'SAFETY: Forcing exit from simplex at eva_t=', eva_t
                    write(6, *) 'This prevents memory exhaustion from excessive iterations'
                    exit  ! Break out of the simplex loop
                end if


                




                if (rrl >= -tol) then

  


                    simplex_converged = .true.  ! Mark as converged
                    exit
                end if




                
                ! Additional check: if we've been iterating too long with little progress,
                ! the problem might be degenerate
                if (iter > 100 .and. abs(rrl) < 1.0d-6) then
                    write(6, *) 'WARNING: Simplex making little progress at eva_t=', eva_t
                    write(6, *) 'rrl=', rrl, 'iter=', iter
                    ! Force exit to avoid infinite loop
                    exit
                end if



                ! Step 3: Choose entering variable
                if (bland) then
                    t = huge(1)
                    t_rr = 0
                    
                    ! Check first row of reduced costs
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
                    ! Steepest descent rule
                    do j = 1, 2*(nvar+1)
                        do i = 1, 2
                            if (abs(rr(i, j) - rrl) < 1.0d-14) then
                                t_rr = j
                                tsep = i
                                if (tsep == 1) then
                                    t = r1(t_rr)
                                else
                                    t = r2(t_rr)
                                end if
                                goto 200
                            end if
                        end do
                    end do
                    200 continue
                end if
                




                ! Step 4: Get pivot column
                if (tsep == 1) then
                    yy(1:ms+1) = gammaxs(:, t_rr)
                else
                    yy(1:ms+1) = -gammaxs(:, t_rr)
                end if
                

 

                ! Step 5: Choose leaving variable (ratio test)
                min_k = huge(1.0d0)
                k = 0
                
                
                ! First, find all valid ratios
                do i = 1, ms
                    if (yy(i) > 0.0d0 .and. .not. freevarrow(i)) then
                        k_vals(i) = bs(i) / yy(i)
                        if (k_vals(i) < min_k) then  !k_vals(i) >= -tol .and.
                            min_k = k_vals(i)
                            if (.not. bland) then
                                k = i
                            end if
                        end if
                    else
                        k_vals(i) = huge(1.0d0)
                    end if
                end do





                ! Check if we found any valid pivots
                if (min_k >= huge(1.0d0)) then

     

 



                    ! No valid pivots found - problem is unbounded
                    write(6, *) 'The problem is unbounded, doubling m at time', eva_t
                    

     
                    empty_pivot_count = empty_pivot_count + 1
                    if (empty_pivot_count >= 2) then
                        ierr = 31
                        return
                    end if
                    mmm_thresh = 2.0d0 * mmm_thresh
                    not_new_sl_sh = .true.
                    if (empty_pivot_count >= max_empty_pivot_retries) then
                        force_full_sample = .true.
                    end if

                    unbounded_detected = .true.  ! SET THE FLAG
                    ! Exit the simplex loop to restart preprocessing
                    exit
                end if




                ! For Bland rule, choose the one with smallest index
                if (bland .and. min_k < huge(1.0d0)) then
                    k = 0
                    do i = 1, ms
                        if (abs(k_vals(i) - min_k) < 1.0d-14) then
                            if (k == 0 .or. IBs(i) < IBs(k)) then
                                k = i
                            end if
                        end if
                    end do
                end if




                ! Check if we found a valid pivot
                if (k == 0 .or. min_k >= huge(1.0d0)) then
                    write(6, *) 'WARNING: No valid leaving variable at eva_t=', eva_t, ', iter=', iter
                    write(6, *) 'DETAILED DEBUG INFO:'
                    write(6, *) '  Entering variable t=', t, ', t_rr=', t_rr, ', tsep=', tsep
                    write(6, *) '  Total rows ms=', ms

                    ! Count and show non-free variables
                    k = 0
                    do i = 1, ms
                        if (.not. freevarrow(i)) k = k + 1
                    end do
                    write(6, *) '  Number of non-free variables:', k

                    ! Show details for non-free variables with positive yy
                    write(6, *) '  Non-free vars with yy > 0:'
                    do i = 1, ms
                        if (.not. freevarrow(i) .and. yy(i) > 0.0d0) then
                            write(6, '(A,I3,A,F12.6,A,F12.6,A,F12.6,A,I6)') &
                                '    Row ', i, ': yy=', yy(i), ', bs=', bs(i), &
                                ', ratio=', bs(i)/yy(i), ', IBs=', IBs(i)
                        end if
                    end do

                    ! Check for negative ratios
                    write(6, *) '  Non-free vars with yy > 0 but bs < 0 (infeasible):'
                    k = 0
                    do i = 1, ms
                        if (.not. freevarrow(i) .and. yy(i) > 0.0d0 .and. bs(i) < 0.0d0) then
                            k = k + 1
                            if (k <= 5) then  ! Show first 5
                                write(6, '(A,I3,A,F12.6,A,F12.6)') &
                                    '    Row ', i, ': yy=', yy(i), ', bs=', bs(i)
                            end if
                        end if
                    end do
                    if (k > 5) write(6, *) '    ... and', k-5, 'more'
                    ! Also print the current objective row
                    write(6, *) '  Current objective row (first 8):', (gammaxs(ms+1, j), j=1, 2*(nvar+1))

                    write(6, *) 'This indicates an unbounded or degenerate problem'
                    
                    ! Try to diagnose the issue
                    write(6, *) 'Diagnostic info:'
                    write(6, *) '  Number of positive yy:', count(yy(1:ms) > tol)
                    write(6, *) '  Number of non-free vars:', count(.not. freevarrow(1:ms))
                    write(6, *) '  Min yy:', minval(yy(1:ms)), 'Max yy:', maxval(yy(1:ms))
                    empty_pivot_count = empty_pivot_count + 1
                    mmm_thresh = 2.0d0 * mmm_thresh
                    not_new_sl_sh = .true.
                    if (empty_pivot_count >= max_empty_pivot_retries) then
                        force_full_sample = .true.
                    end if
                    unbounded_detected = .true.
                    exit
                end if
                
                ! Pivoting step 6': Adjust for entering v_i
                if (tsep /= 1) then
                    ! Use subsample weight exactly as R does
                    j = r1(t_rr) - 2 - 2*nvar
                    if (j >= 1) then
                        yy(ms + 1) = yy(ms + 1) + ws(j)
                    end if
                end if
                

                ! Check if we actually found a valid pivot
                if (k == 0) then
                    write(6, *) 'ERROR: No leaving variable found (k=0) at eva_t=', eva_t, ', iter=', iter
                    write(6, *) 'min_k=', min_k
                    write(6, *) 'This means no feasible pivot exists'
                    empty_pivot_count = empty_pivot_count + 1
                    mmm_thresh = 2.0d0 * mmm_thresh
                    not_new_sl_sh = .true.
                    if (empty_pivot_count >= max_empty_pivot_retries) then
                        force_full_sample = .true.
                    end if
                    unbounded_detected = .true.
                    exit
                end if


                ! Perform pivot

                ! Replace the entire "Check for near-zero pivot before division" block with:
                if (abs(yy(k)) < 1.0d-10) then
                    write(6, *) 'WARNING: Near-zero pivot element at eva_t=', eva_t, ', iter=', iter
                    write(6, *) 'Original pivot value:', yy(k), 'at position k=', k
                    
                    ! Instead of giving up, let's try to find a better pivot
                    ! This is what R's solver does implicitly
                    
                    ! Initialize search for best alternative pivot
                    best_pivot_value = 0.0d0
                    best_k = 0
                    
                    ! Look through all possible pivot candidates
                    do i = 1, ms
                        ! Skip variables that must stay in the basis (freevarrow)
                        if (.not. freevarrow(i)) then
                            ! We need a non-zero element in the pivot column
                            if (abs(yy(i)) > best_pivot_value) then
                                ! But we also need to maintain feasibility
                                ! Check if pivoting on this element keeps solution non-negative
                                
                                if (yy(i) > 1.0d-12) then
                                    ! For positive pivot element
                                    if (bs(i) >= -1.0d-10) then  ! bs(i) is essentially non-negative
                                        ratio = bs(i) / yy(i)
                                        if (ratio >= -1.0d-10) then  ! Would give non-negative result
                                            best_pivot_value = abs(yy(i))
                                            best_k = i
                                        end if
                                    end if
                                else if (yy(i) < -1.0d-12) then
                                    ! For negative pivot element
                                    if (bs(i) <= 1.0d-10) then  ! bs(i) is essentially non-positive
                                        ratio = bs(i) / yy(i)
                                        if (ratio >= -1.0d-10) then  ! Would give non-negative result
                                            best_pivot_value = abs(yy(i))
                                            best_k = i
                                        end if
                                    end if
                                end if
                            end if
                        end if
                    end do
                    
                    ! Did we find a viable alternative?
                    if (best_k > 0 .and. best_pivot_value > 1.0d-12) then
                        write(6, *) 'Found alternative pivot with value:', best_pivot_value, 'at position:', best_k
                        k = best_k  ! Use the alternative pivot
                        ! Continue with the normal pivot operation
                    else
                        write(6, *) 'No viable pivot found - problem is truly degenerate'
                        write(6, *) 'Terminating simplex and using current best solution'
                        
                        ! Extract whatever solution we have so far
                        do i = 1, 2*(nvar+1)
                            estimate(i) = bs(i)
                        end do
                        
                        ! Mark this as a numerical failure
                        iter = maxit  ! This will trigger any fallback handling
                        exit
                    end if
                end if

                ee(1:ms+1) = yy(1:ms+1) / yy(k)
                ee(k) = 1.0d0 - 1.0d0 / yy(k)
                




                ! Update basis representation
                if (IBs(k) <= (ms + 2*nvar + 2)) then
                    gammaxs(:, t_rr) = 0.0d0
                    gammaxs(k, t_rr) = 1.0d0
                    r1(t_rr) = IBs(k)
                    r2(t_rr) = IBs(k) + ms
                else
                    gammaxs(:, t_rr) = 0.0d0
                    gammaxs(k, t_rr) = -1.0d0
                    j = IBs(k) - ms - 2*nvar - 2
                    if (j >= 1 .and. j <= ms) then
                        gammaxs(ms + 1, t_rr) = ws(j)
                    else
                        gammaxs(ms + 1, t_rr) = 0.0d0
                    end if
                    r1(t_rr) = IBs(k) - ms
                    r2(t_rr) = IBs(k)
                end if
                
                ! Update tableau
                do j = 1, 2*(nvar+1)
                    pivot_row(j) = gammaxs(k, j)
                end do
                

 


                do j = 1, 2*(nvar+1)
                    do i = 1, ms+1
                        gammaxs(i, j) = gammaxs(i, j) - ee(i) * pivot_row(j)
                    end do
                end do
                
                b_k_original = bs(k)
                do i = 1, ms+1
                    bs(i) = bs(i) - ee(i) * b_k_original
                end do
                
                IBs(k) = t
                iter = iter + 1



            end do
            



            ! ADD ITERATION TRACKING HERE:
            total_simplex_iterations = total_simplex_iterations + iter
            if (iter > max_iter_at_any_t) then
                max_iter_at_any_t = iter
            end if

            ! Log concerning patterns:
            if (iter > 10000) then
                write(6, *) 'WARNING: Excessive iterations at eva_t=', eva_t
                write(6, *) 'Iterations:', iter
                write(6, *) 'Subsample size ms:', ms
                write(6, *) 'Consider the following:'
                write(6, *) '- The problem may be degenerate'
                write(6, *) '- The tolerance may be too tight'
                write(6, *) '- There may be cycling in the simplex'
            end if





            ! Always extract the current solution from the tableau
            ! This ensures we have an estimate even if unbounded
            do i = 1, 2*(nvar+1)
               ! if (i <= ms) then
                    estimate(i) = bs(i)



                !else
                   ! estimate(i) = 0.0d0  ! Safety fallback
              !  end if
            end do


   


            if (iter == maxit) then
                write(6, *) 'WARNING: Max iterations reached at eva_t =', eva_t
                write(6, *) 'Simplex failed to converge. Increasing threshold / full-sample fallback.'
                empty_pivot_count = empty_pivot_count + 1
                mmm_thresh = 2.0d0 * mmm_thresh
                not_new_sl_sh = .true.
                if (empty_pivot_count >= max_empty_pivot_retries) then
                    force_full_sample = .true.
                end if
                cycle
            else
                ! Normal case - simplex converged
                do i = 1, 2*(nvar+1)
                    estimate(i) = bs(i)
                end do
            end if

 

            300 continue
            
            ! Check signs of residuals using FULL sample
            do i = 1, m
                r(i) = y(i)
                do j = 1, 2*(nvar+1)
                    r(i) = r(i) - A(i, j) * estimate(j)
                end do
            end do




            ! Count bad signs regardless of simplex convergence
            ! This matches R's behavior
            bad_signs = 0
            n_sure_signs = n_sl + n_sh  ! Total sure-sign observations





            ! Only check bad signs if we didn't hit max iterations
            ! (but DO check them if unbounded was detected)
            if (iter < maxit) then
                do i = 1, m
                    if ((r(i) < 0.0d0) .and. sh(i)) then
                        bad_signs = bad_signs + 1


                    end if
                    if ((r(i) > 0.0d0) .and. sl(i)) then
                        bad_signs = bad_signs + 1

                    end if
                end do


            else
                ! If we hit max iterations, force exit without bad signs check
                not_optimal = .false.
            end if




            ! Handle bad signs
            if (bad_signs > 0) then

                if (bad_signs > int(0.1d0 * dble(ms))) then
                    mmm_thresh = 2.0d0 * mmm_thresh
                    not_new_sl_sh = .true.
                else
                    ! Fix bad signs ONLY when bad_signs <= 0.1*ms
                    do i = 1, m
                        if ((r(i) < 0.0d0) .and. sh(i)) then
                            sh(i) = .false.
                        end if
                        if ((r(i) > 0.0d0) .and. sl(i)) then
                            sl(i) = .false.
                        end if
                    end do
                    not_new_sl_sh = .false.







                end if

            else
                ! No bad signs - we've reached optimality
                do i = 1, 2*(nvar+1)
                    j = r1(i) - 2 - 2*nvar
                    if (j > 0 .and. j <= ms_org) then
                        H_indices(i) = idx_not_jl_or_jh(j)
                    else
                        H_indices(i) = 0
                    end if
                end do
                accept_subsample = certify_tvcqr_candidate(H_indices, r, A, m, 2*(nvar+1), res_tol)
                if (.not. accept_subsample) then
                    if ((.not. any(sl)) .and. (.not. any(sh)) .and. (ms >= m)) then
                        ierr = 1
                        return
                    end if
                    mmm_thresh = 2.0d0 * mmm_thresh
                    not_new_sl_sh = .true.
                    cycle
                end if
                not_optimal = .false.
                





                ! Extract H for next iteration
                do i = 1, 2*(nvar+1)
                    if (H_indices(i) > 0) then
                        H_seq(eva_t, i) = H_indices(i)
                    else
                        ! This matches R's behavior - r1 should always give valid indices
                        ! If not, there's a bug in the simplex algorithm
                        write(6, *) 'ERROR: Invalid certified H index at eva_t=', eva_t, ', i=', i
                        write(6, *) 'r1(i)=', r1(i), ', ms_org=', ms_org
                        H_seq(eva_t, i) = 0
                    end if




                end do
                


               


                ! If not the last time point, save information for reuse
                if (eva_t < m) then



  


  
                    ! Calculate ms_org for saving
                    if (n_sl > 0 .and. n_sh > 0) then
                        ms_org = ms - 2
                    else if (n_sl > 0 .or. n_sh > 0) then
                        ms_org = ms - 1
                    else
                        ms_org = ms
                    end if
                    
                    ! Check if we came from degenerate case
                    !if (iter > 0) then
                        ! Normal case - save xhinv from gammaxs
                        do i = 1, 2*(nvar+1)
                            do j = 1, 2*(nvar+1)
                                gammaxs_temp(i, j) = gammaxs(i, j)
                            end do
                            bs_temp(i) = bs(i)
                        end do
                    !end if
                    ! If iter == 0 (degenerate case), xhinv is already in gammaxs_temp
                    
                    ! Initialize counters
                    n_Hbar_pos = 0
                    n_Hbar_neg = 0
                    
                    ! Only process Hbar information if there are non-aggregated observations
                    if (ms_org > 0 .and. iter > 0) then






                        ! Process the IBs to extract Hbar information
                        do idx = 2*(nvar+1) + 1, ms_org
                            i = idx - 2*(nvar+1)! idx = 2*(nvar+1) + i
                            if (IBs(idx) > 2*(nvar+1)+ms) then
                                ! This is a v variable
                                n_Hbar_neg = n_Hbar_neg + 1


               



                                id_gammaxs_Hbar(i) = IBs(idx) - 2*(nvar+1) - ms
                                idx_Hbar_neg(n_Hbar_neg) = idx_not_jl_or_jh(id_gammaxs_Hbar(i))
                                do j = 1, 2*(nvar+1)
                                    gammaxs_neg(n_Hbar_neg, j) = gammaxs(idx, j)
                                end do
                                bs_neg(n_Hbar_neg) = bs(idx)
                            else
                                ! This is a u variable
                                n_Hbar_pos = n_Hbar_pos + 1
                                id_gammaxs_Hbar(i) = IBs(idx) - 2*(nvar+1)
                                idx_Hbar_pos(n_Hbar_pos) = idx_not_jl_or_jh(id_gammaxs_Hbar(i))
                                do j = 1, 2*(nvar+1)
                                    gammaxs_pos(n_Hbar_pos, j) = gammaxs(idx, j)
                                end do
                                bs_pos(n_Hbar_pos) = bs(idx)
                            end if
                        end do

     


                    end if

                end if
                
            end if
            
        end do  ! End while not_optimal
        
        ! Store results for this time point
        it_num(eva_t) = iter
        n_sub(eva_t) = ms
        M_out = M_threshold
        

        ! CRITICAL FIX: Always recompute residuals with final estimate
        ! This ensures residual_est is consistent with the estimate
        do i = 1, m
            r(i) = y(i)
            do j = 1, 2*(nvar+1)
                r(i) = r(i) - A(i, j) * estimate(j)
            end do
        end do




        ! Compute theta_ll_est for current eva_t
        do j = 1, nvar+1
            theta_ll_est(eva_t, j) = estimate(j) + (dble(eva_t)/dble(m)) * estimate(nvar+1+j)
        end do

        ! DEBUG: Show estimates at problem times

        if (debug_active) then
            write(6, *) 'FORTRAN: Raw estimate vector:'
            write(6, '(A)', advance='no') '  estimate = '
            do j = 1, 2*(nvar+1)
                write(6, '(F10.6,1X)', advance='no') estimate(j)
            end do
            write(6, *)
        end if





        do i = 1, m
            r_prev(i) = r(i)
            if (store_residual) then
                residual_est(eva_t, i) = r(i)
            end if
        end do


        
    end do  ! End of eva_t loop
    
    
    deallocate(A)
    deallocate(gammax)
    deallocate(gammaxs_temp)
    deallocate(bs_temp)
    deallocate(gammaxs)
    deallocate(bs)

contains

    logical function certify_tvcqr_candidate(H_idx, r_vec, A_mat, m_loc, p, res_tol_loc)
        implicit none
        integer, intent(in) :: m_loc, p
        integer, intent(in) :: H_idx(p)
        double precision, intent(in) :: r_vec(m_loc), A_mat(m_loc, p), res_tol_loc
        integer :: i, j, ierr_local
        double precision :: AH(p, p), AH_inv(p, p)

        certify_tvcqr_candidate = .true.

        do i = 1, p
            if (H_idx(i) < 1 .or. H_idx(i) > m_loc) then
                certify_tvcqr_candidate = .false.
                return
            end if
            do j = i + 1, p
                if (H_idx(i) == H_idx(j)) then
                    certify_tvcqr_candidate = .false.
                    return
                end if
            end do
            if (abs(r_vec(H_idx(i))) > res_tol_loc) then
                certify_tvcqr_candidate = .false.
                return
            end if
        end do

        do i = 1, p
            do j = 1, p
                AH(i, j) = A_mat(H_idx(i), j)
            end do
        end do
        call matrix_inverse_2p(AH, AH_inv, p, ierr_local)
        if (ierr_local /= 0) then
            certify_tvcqr_candidate = .false.
            return
        end if
    end function certify_tvcqr_candidate

    ! Helper subroutine to save algorithm state
    subroutine save_checkpoint(eva_t, filename_prefix)
        integer, intent(in) :: eva_t
        character(len=*), intent(in) :: filename_prefix
        character(len=100) :: filename
        integer :: unit_num
        double precision :: res_min, res_max, res_mean

        ! Save coefficient estimates
        write(filename, '(A,A,I4.4,A)') trim(filename_prefix), '_coef_t', eva_t, '.dat'
        open(newunit=unit_num, file=filename, status='replace')
        write(unit_num, '(A,I4)') '# Coefficients at eva_t = ', eva_t
        do i = 1, eva_t
            write(unit_num, '(I5,20F12.6)') i, (theta_ll_est(i,j), j=1,nvar+1)
        end do
        close(unit_num)
        
        ! Save residual information
        write(filename, '(A,A,I4.4,A)') trim(filename_prefix), '_resid_t', eva_t, '.dat'
        open(newunit=unit_num, file=filename, status='replace')
        write(unit_num, '(A,I4)') '# Residual summary at eva_t = ', eva_t
        do i = 1, eva_t

            res_min = minval(residual_est(i, 1:m))
            res_max = maxval(residual_est(i, 1:m))
            res_mean = sum(residual_est(i, 1:m)) / dble(m)
            write(unit_num, '(I5,3F12.6,2I8)') i, res_min, res_max, res_mean, &
                it_num(i), n_sub(i)
        end do
        close(unit_num)
        
        write(6, '(A,A)') 'Checkpoint saved with prefix: ', trim(filename_prefix)
    end subroutine save_checkpoint

    ! Helper function to create comma-separated list of indices where logical array is true
    function int_list(logical_array, true_val) result(str)
        logical, intent(in) :: logical_array(:)
        integer, intent(in) :: true_val
        character(len=200) :: str
        integer :: i, pos
        
        str = ''
        pos = 1
        do i = 1, size(logical_array)
            if (logical_array(i) .eqv. (true_val == 0)) then
                if (pos > 1) then
                    write(str(pos:pos), '(A1)') ','
                    pos = pos + 1
                end if
                write(str(pos:), '(I0)') i
                pos = pos + len_trim(str(pos:)) + 1
            end if
        end do
    end function int_list

    ! Helper to find location of minimum in 2D array
    subroutine minloc2d(arr, min_val, row, col)
        double precision, intent(in) :: arr(:,:)
        double precision, intent(in) :: min_val
        integer, intent(out) :: row, col
        integer :: i, j
        
        do i = 1, size(arr, 1)
            do j = 1, size(arr, 2)
                if (abs(arr(i,j) - min_val) < 1.0d-14) then
                    row = i
                    col = j
                    return
                end if
            end do
        end do
        row = 0
        col = 0
    end subroutine minloc2d
    ! Matrix inversion subroutine using LAPACK
    subroutine matrix_inverse_2p(A_in, A_inv, n, ierr)
        implicit none
        integer, intent(in) :: n
        double precision, intent(in) :: A_in(n, n)
        double precision, intent(out) :: A_inv(n, n)
        integer, intent(out) :: ierr
        
        ! LAPACK workspace variables
        integer :: lwork, info
        integer :: ipiv(n)
        double precision, allocatable :: work(:)
        
        ierr = 0
        
        ! Copy input matrix (LAPACK overwrites)
        A_inv = A_in
        
        ! LU factorization
        call DGETRF(n, n, A_inv, n, ipiv, info)
        if (info /= 0) then
            ierr = info
            return
        end if
        
        ! Query optimal workspace
        allocate(work(1))
        lwork = -1
        call DGETRI(n, A_inv, n, ipiv, work, lwork, info)
        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork))
        
        ! Compute inverse
        call DGETRI(n, A_inv, n, ipiv, work, lwork, info)
        if (info /= 0) then
            ierr = info
            return
        end if
        
        deallocate(work)
        
    end subroutine matrix_inverse_2p


    recursive subroutine quicksort_real(arr, left, right)
        implicit none
        double precision, intent(inout) :: arr(:)
        integer, intent(in) :: left, right
        integer :: i, j
        double precision :: pivot, temp

        if (left >= right) return
        if (right - left <= 16) then
            call insertion_sort_real(arr, left, right)
            return
        end if

        pivot = arr((left + right) / 2)
        i = left
        j = right

        do
            do while (arr(i) < pivot)
                i = i + 1
            end do
            do while (arr(j) > pivot)
                j = j - 1
            end do
            if (i <= j) then
                temp = arr(i)
                arr(i) = arr(j)
                arr(j) = temp
                i = i + 1
                j = j - 1
            end if
            if (i > j) exit
        end do

        if (left < j) call quicksort_real(arr, left, j)
        if (i < right) call quicksort_real(arr, i, right)
    end subroutine quicksort_real

    subroutine insertion_sort_real(arr, left, right)
        implicit none
        double precision, intent(inout) :: arr(:)
        integer, intent(in) :: left, right
        integer :: i, j
        double precision :: value

        do i = left + 1, right
            value = arr(i)
            j = i - 1
            do while (j >= left .and. arr(j) > value)
                arr(j + 1) = arr(j)
                j = j - 1
            end do
            arr(j + 1) = value
        end do
    end subroutine insertion_sort_real

    double precision function median_value(arr, n)
        implicit none
        integer, intent(in) :: n
        double precision, intent(in) :: arr(n)
        double precision :: sorted(n)

        sorted = arr
        call quicksort_real(sorted, 1, n)

        if (mod(n, 2) == 0) then
            median_value = 0.5d0 * (sorted(n/2) + sorted(n/2 + 1))
        else
            median_value = sorted((n+1)/2)
        end if
    end function median_value

    
end subroutine tvcqr_seq_ppro_fortran
