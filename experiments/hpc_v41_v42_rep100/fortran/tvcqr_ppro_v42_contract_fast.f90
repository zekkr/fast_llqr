! Time-varying coefficient quantile regression - Sequential preprocessing
! Faithful translation of tvcqr_seq_ppro.R with focus on correctness for eva_t <= 4
! Temporary v42 contract-aligned ablation: v39 plus fit-wide inverse-workspace
! reuse, fused retained-row/aggregate accumulation, and exact contiguous-support
! loop simplification. The first point retains the original full-m solve.
! The accepted residual image is materialized only on the current/next exact
! positive-weight windows and exceptional previous/candidate H rows. Dense
! consumers and independent/full-active recovery retain the original full-m
! path. Retained rows are merged in original observation-ID order from the
! exact active window and at most q forced-H exceptions, avoiding full-m mask
! clears and packing scans. Screening, literal verification, H zeroing, and all
! repair/recovery decisions are unchanged.
!
! Compile with: R CMD SHLIB tvcqr_seq_corrected.f90
! or: gfortran -shared -fPIC -o tvcqr_seq_M_acc.so tvcqr_seq_M_acc.f90 -llapack -lblas

subroutine tvcqr_seq_ppro_fortran(x, y, m, nvar, tau, h, h_factor, tol, maxit, &
                                   bland_int, Mm_factor, eps, store_residual_int, debug_int, &
                                   theta_ll_est, beta_full_est, it_num, residual_est, &
                                   M_out, first_n_sub, repair_count, final_n_sub, &
                                   init_mode, init_trigger, H_seq, same_h_refit_attempted, &
                                   same_h_refit_recovered, ierr, failed_eval, min_subsample_size_in, &
                                   always_same_h_refit_int, threshold_lower_bound_int, &
                                   threshold_scale_mode_int)
    
    implicit none
    
    ! Input arguments
    integer, intent(in) :: m, nvar, maxit, bland_int, store_residual_int, debug_int
    integer, intent(in) :: min_subsample_size_in, always_same_h_refit_int
    integer, intent(in) :: threshold_lower_bound_int, threshold_scale_mode_int
    double precision, intent(in) :: x(m, nvar), y(m), tau, tol, h_factor
    double precision, intent(in) :: Mm_factor, eps
    double precision, intent(inout) :: h
    
    ! Output arguments
    double precision, intent(out) :: theta_ll_est(m, nvar+1)
    double precision, intent(out) :: beta_full_est(m, 2*(nvar+1))
    integer, intent(out) :: it_num(m)
    double precision, intent(out) :: residual_est(*)
    double precision, intent(out) :: M_out
    integer, intent(out) :: first_n_sub(m)
    integer, intent(out) :: repair_count(m)
    integer, intent(out) :: final_n_sub(m)
    integer, intent(out) :: init_mode(m)
    integer, intent(out) :: init_trigger(m)
    integer, intent(out) :: H_seq(m, 2*(nvar+1))
    integer, intent(out) :: same_h_refit_attempted(m)
    integer, intent(out) :: same_h_refit_recovered(m)
    integer, intent(out) :: ierr
    integer, intent(out) :: failed_eval
    
    ! Local variables
    double precision :: x_norms(m), mm_thresh, mmm_thresh, M_threshold, threshold_scale
    logical :: sl(m), sh(m), not_jl_or_jh(m)
    logical :: has_sl_agg, has_sh_agg
    integer :: idx_not_jl_or_jh(m)
    integer :: forced_h_ids(2*(nvar+1)), n_forced_h, insert_id, insert_pos
    double precision :: temp_check
    integer :: min_subsample_size, n_potential_S
    integer :: n_active, target_min
    double precision :: residual_scale, pivot_tol
    double precision :: abs_r(m)  ! NEW - for median of absolute residuals
    
    integer :: ms, ms_org
    double precision :: ws(m+3)
    
    double precision :: glob_wx(2*(nvar+1)), glob_wy
    double precision :: ghib_wx(2*(nvar+1)), ghib_wy
    
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
    double precision :: estimate_refit(2*(nvar+1)), r_refit(m)
    double precision :: r(m), r_prev(m)
    double precision :: beta_prev(2*(nvar+1))
    double precision :: pivot_row(2*(nvar+1))
    
    integer :: i, j, k, t, eva_t, iter
    integer :: active_lo, active_hi, previous_active_lo, previous_active_hi
    integer :: next_active_lo, next_active_hi
    integer :: candidate_residual_lo, candidate_residual_hi
    integer :: r_prev_valid_lo, r_prev_valid_hi, r_prev_extra_count
    integer :: r_prev_extra(4*(nvar+1))
    integer :: probe_idx
    integer :: t_rr, tsep
    double precision :: rrl, min_k, temp_sum, eval_time, probe_weight
    logical :: bland, store_residual, debug_requested, always_same_h_refit, threshold_lower_bound
    double precision :: b_k_original
    logical :: not_optimal, not_new_sl_sh, debug_active, same_h_ok, refit_used
    logical :: lazy_residual_mode, candidate_residual_full, r_prev_full
    logical :: next_window_ok, same_h_solve_ok
    logical :: duplicate_h
    logical :: first_n_sub_recorded
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
    integer :: inv_info
    logical :: valid_H
    double precision :: temp_vec1(2*(nvar+1))
    integer :: ii, kk
    logical :: force_full_sample, accept_subsample, no_pivot_flag
    logical :: use_independent_init, previous_H_valid, shifted_init_success
    logical :: full_recovery_success
    integer :: independent_trigger
    integer :: empty_pivot_count, max_empty_pivot_retries
    double precision :: res_tol

    
    integer :: n_hbar_rows        ! For the critical lambda/Pxhbarxhinv dimension
    integer :: idx_count, idx_loop
    double precision :: weight_sum, sum_w_sl, sum_w_sh

    double precision, allocatable :: A(:,:)           ! Size: m × 2(nvar+1)
    double precision, allocatable :: gammax(:,:)      ! Size: (m+1) × 2(nvar+1)
    double precision, allocatable :: gammaxs_temp(:,:)! Size: (m+3) × 2(nvar+1)
    double precision, allocatable :: bs_temp(:)       ! Size: m+3
    double precision, allocatable :: gammaxs(:,:)     ! Size: (m+3) × 2(nvar+1)
    double precision, allocatable :: bs(:)            ! Size: m+3
    logical :: unbounded_detected
    logical :: simplex_converged
    integer :: iter_attempt
    integer :: total_simplex_iterations
    integer :: total_preprocessing_loops
    integer :: max_iter_at_any_t
    integer :: inverse_lwork, inverse_query_info
    integer :: inverse_query_piv(2*(nvar+1))
    double precision :: inverse_query_work(1)
    double precision, allocatable :: inverse_work(:)



 

    ierr = 0
    failed_eval = 0
    total_simplex_iterations = 0
    total_preprocessing_loops = 0
    max_iter_at_any_t = 0
    active_lo = 0
    active_hi = 0
    next_active_lo = 1
    next_active_hi = 0
    candidate_residual_lo = 1
    candidate_residual_hi = 0
    r_prev_valid_lo = 1
    r_prev_valid_hi = 0
    r_prev_extra_count = 0
    r_prev_full = .true.
    



    ! Start of executable code
    
    ! Convert integer to logical for bland
    bland = (bland_int /= 0)
    store_residual = (store_residual_int /= 0)
    debug_requested = (debug_int /= 0)
    always_same_h_refit = (always_same_h_refit_int /= 0)
    threshold_lower_bound = (threshold_lower_bound_int /= 0)
    ! The first exact lazy path deliberately excludes every interface mode that
    ! consumes or returns a dense residual vector, plus the rescue-only refit
    ! mode whose alternate acceptance path remains byte-for-byte dense.
    lazy_residual_mode = always_same_h_refit .and. (.not. threshold_lower_bound) .and. &
                         (.not. debug_requested) .and. (.not. store_residual)
    pivot_tol = max(10.0d0 * tol, 1.0d-12)
    max_empty_pivot_retries = 3
    res_tol = 1.0d-6
    
    ! Allocate the big arrays
    allocate(A(m, 2*(nvar+1)))
    allocate(gammax(m+1, 2*(nvar+1)))
    allocate(gammaxs_temp(m+3, 2*(nvar+1)))
    allocate(bs_temp(m+3))
    allocate(gammaxs(m+3, 2*(nvar+1)))
    allocate(bs(m+3))

    ! All inversions have the same q-by-q shape. Query LAPACK once and reuse
    ! its preferred workspace for the complete fit.
    xhinv = 0.0d0
    do i = 1, 2*(nvar+1)
        inverse_query_piv(i) = i
    end do
    inverse_lwork = -1
    call DGETRI(2*(nvar+1), xhinv, 2*(nvar+1), inverse_query_piv, &
                inverse_query_work, inverse_lwork, inverse_query_info)
    if (inverse_query_info /= 0) then
        ierr = 2
        failed_eval = 1
        return
    end if
    inverse_lwork = max(2*(nvar+1), int(inverse_query_work(1)))
    allocate(inverse_work(inverse_lwork))

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
    beta_full_est = 0.0d0
    if (store_residual) then
        do i = 1, m*m
            residual_est(i) = 0.0d0
        end do
    end if
    first_n_sub = 0
    repair_count = 0
    final_n_sub = 0
    init_mode = 0
    init_trigger = 0
    H_seq = 0
    same_h_refit_attempted = 0
    same_h_refit_recovered = 0
    M_out = 0.0d0
    if (threshold_scale_mode_int == 1) then
        threshold_scale = log(log(dble(m)))
    else if (threshold_scale_mode_int == 2) then
        threshold_scale = log(dble(m))
    else
        ierr = 5
        failed_eval = 1
        return
    end if
    if ((.not. threshold_lower_bound) .and. Mm_factor <= 0.0d0) then
        ierr = 5
        failed_eval = 1
        return
    end if

    ! ============================================
    ! EVA_T = 1: Standard simplex (no preprocessing)
    ! ============================================
    
    ! Calculate weights for t=1 and record the exact positive-weight window.
    n_active = 0
    do i = 1, m
        if (abs(1.0d0/dble(m) - time_index(i)) <= h) then
            w(i) = 0.75d0 * (1.0d0 - ((1.0d0/dble(m) - time_index(i))/h)**2)
        else
            w(i) = 0.0d0
        end if
        if (w(i) > 0.0d0) then
            n_active = n_active + 1
            if (active_lo == 0) active_lo = i
            active_hi = i
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
    
    first_n_sub(1) = m
    repair_count(1) = 0
    final_n_sub(1) = m
    
    ! Simplex iterations for eva_t = 1
    iter = 0
    simplex_converged = .false.




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
        if (rrl >= -tol) then
            simplex_converged = .true.
            exit
        end if
        
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
            do j = 1, 2*(nvar+1)
                do i = 1, 2
                    if (abs(rr(i, j) - rrl) <= 1.0d-14) then
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
                if (yy(i) > pivot_tol .and. .not. freevarrow(i)) then
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

            if (k == 0 .or. min_k >= huge(1.0d0)) then
                call set_failure(2, 1)
                return
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
                    if (yy(i) > pivot_tol .and. .not. freevarrow(i)) then
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
                    if (yy(i) < -pivot_tol .and. .not. freevarrow(i)) then
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

            if (k == 0 .or. min_k >= huge(1.0d0)) then
                call set_failure(2, 1)
                return
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
    
    if (.not. simplex_converged) then
        call set_failure(2, 1)
        return
    end if

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
    do j = 1, 2*(nvar+1)
        beta_full_est(1, j) = estimate(j)
    end do
    
    do i = 1, m
        r_prev(i) = u(i) - v(i)
        if (store_residual) then
            residual_est((i - 1) * m + 1) = r_prev(i)
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
    do i = 1, 2*(nvar+1)
        if (H_seq(1, i) < 1 .or. H_seq(1, i) > m) then
            call set_failure(3, 1)
            return
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






        debug_active = debug_requested
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

        ! The time grid is increasing and the Epanechnikov support is compact,
        ! so both positive-weight endpoints move monotonically to the right.
        ! Probe endpoints with the exact production weight expression, clear
        ! observations leaving on the left, and recompute every weight inside
        ! the new window. Values to the right have never entered and remain zero.
        previous_active_lo = active_lo
        previous_active_hi = active_hi
        eval_time = dble(eva_t) / dble(m)

        active_lo = max(1, previous_active_lo)
        do while (active_lo <= m)
            if (abs(eval_time - time_index(active_lo)) <= h) then
                probe_weight = 0.75d0 * (1.0d0 - ((eval_time - time_index(active_lo))/h)**2)
            else
                probe_weight = 0.0d0
            end if
            if (probe_weight > 0.0d0) exit
            active_lo = active_lo + 1
        end do
        if (active_lo > m) then
            call set_failure(5, eva_t)
            return
        end if

        active_hi = max(previous_active_hi, active_lo - 1)
        do while (active_hi < m)
            probe_idx = active_hi + 1
            if (abs(eval_time - time_index(probe_idx)) <= h) then
                probe_weight = 0.75d0 * (1.0d0 - ((eval_time - time_index(probe_idx))/h)**2)
            else
                probe_weight = 0.0d0
            end if
            if (probe_weight <= 0.0d0) exit
            active_hi = probe_idx
        end do

        if (previous_active_lo > 0 .and. active_lo > previous_active_lo) then
            do i = previous_active_lo, min(active_lo - 1, previous_active_hi)
                w(i) = 0.0d0
            end do
        end if

        do i = active_lo, active_hi
            if (abs(eval_time - time_index(i)) <= h) then
                w(i) = 0.75d0 * (1.0d0 - ((eval_time - time_index(i))/h)**2)
            else
                w(i) = 0.0d0
            end if
        end do
        ! Epanechnikov support on the ordered regular grid has no interior holes:
        ! every index between the exact first/last positive endpoints is active.
        n_active = active_hi - active_lo + 1

        ! Read-only lookahead for residual ownership. This evaluates the same
        ! stored-time kernel expression as the live weight update but never
        ! mutates w, active, or the current endpoints.
        if (lazy_residual_mode .and. eva_t < m) then
            call find_positive_window(eva_t + 1, active_lo, active_hi, &
                                      next_active_lo, next_active_hi, next_window_ok)
            if (.not. next_window_ok) then
                ! This cannot occur on the regular time grid because eva_t+1 is
                ! itself a strictly positive anchor. Preserve correctness by
                ! promoting this call to the dense residual path.
                lazy_residual_mode = .false.
                next_active_lo = 1
                next_active_hi = 0
            end if
        else
            next_active_lo = 1
            next_active_hi = 0
            next_window_ok = .true.
        end if

        candidate_residual_lo = active_lo
        candidate_residual_hi = active_hi
        if (next_active_hi >= next_active_lo) then
            candidate_residual_lo = min(candidate_residual_lo, next_active_lo)
            candidate_residual_hi = max(candidate_residual_hi, next_active_hi)
        end if


 

        
        mmm_thresh = mm_thresh
        
        j = 0
        iter = 0  ! Initialize iteration counter here, outside preprocessing loop
        repair_count(eva_t) = 0
        first_n_sub_recorded = .false.
        use_independent_init = .false.
        independent_trigger = 0
        init_mode(eva_t) = 1
        init_trigger(eva_t) = 0

        preprocessing_attempts = 0
        do while (not_optimal)
            unbounded_detected = .false.  ! Initialize for each attempt
            simplex_converged = .false.
            no_pivot_flag = .false.
            iter_attempt = 0

            preprocessing_attempts = preprocessing_attempts + 1






            total_preprocessing_loops = total_preprocessing_loops + 1
            ! Add safety valve to prevent infinite preprocessing
            if (preprocessing_attempts > max(10, 2*(nvar+1))) then
                force_full_sample = .true.
                not_new_sl_sh = .true.
            end if
            ! Get only the accepted residual entries that this preparation can
            ! consume. A transition to any dense mode first reconstructs the
            ! complete previous accepted image from its saved beta and H.
            if (lazy_residual_mode) then
                if (.not. previous_range_covered(active_lo, active_hi)) then
                    call promote_previous_residual(eva_t - 1)
                end if
                do i = active_lo, active_hi
                    r(i) = r_prev(i)
                end do
            else
                if (.not. r_prev_full) call promote_previous_residual(eva_t - 1)
                do i = 1, m
                    r(i) = r_prev(i)
                end do
            end if

            if (.not. use_independent_init) then
                call validate_previous_H(H_seq(eva_t-1, :), previous_H_valid)
                if (.not. previous_H_valid) then
                    use_independent_init = .true.
                    independent_trigger = 1
                    init_mode(eva_t) = 2
                    init_trigger(eva_t) = independent_trigger
                end if
            end if

            ! Ordinary retained-row sign initialization also reads previous-H
            ! residuals. They may lie outside both active windows, so copy them
            ! explicitly without widening the contiguous residual span.
            if (lazy_residual_mode .and. (.not. use_independent_init)) then
                do k = 1, 2*(nvar+1)
                    idx = H_seq(eva_t-1, k)
                    if (.not. previous_id_covered(idx)) then
                        call promote_previous_residual(eva_t - 1)
                    end if
                    r(idx) = r_prev(idx)
                end do
            end if
            



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

            ! Calculate threshold and partition observations. sl/sh are only
            ! recomputed when requested; all dependent masks/counts are rebuilt
            ! below on every attempt.
            if (force_full_sample) then
                do i = active_lo, active_hi
                    sl(i) = .false.
                    sh(i) = .false.
                end do
            else if (not_new_sl_sh) then
                if (threshold_lower_bound) then
                    ! The median is used only by the optional lower-bound mode.
                    do i = 1, m
                        abs_r(i) = abs(r(i))
                    end do
                    residual_scale = median_value(abs_r, m)
                    M_threshold = max(Mm_factor * mmm_thresh * threshold_scale, 0.1d0 * residual_scale)
                else
                    M_threshold = Mm_factor * mmm_thresh * threshold_scale
                end if

                ! NEW: Ensure minimum subsample size
                if (min_subsample_size_in >= 0) then
                    min_subsample_size = min_subsample_size_in
                else
                    min_subsample_size = max(5 * (2 * (nvar + 1)), ceiling(0.2d0 * dble(m)))
                end if

                target_min = min(min_subsample_size, n_active)

                ! Count active potential observations in S
                n_potential_S = 0
                min_k = 0.0d0
                do i = active_lo, active_hi
                    if (abs(r(i)) > min_k) min_k = abs(r(i))
                    if (abs(r(i)) <= M_threshold) then
                        n_potential_S = n_potential_S + 1
                    end if
                end do

                ! If too few active observations would remain, increase M.
                do while (n_potential_S < target_min .and. M_threshold < min_k)
                    M_threshold = M_threshold * 1.5d0
                    n_potential_S = 0
                    do i = active_lo, active_hi
                        if (abs(r(i)) <= M_threshold) then
                            n_potential_S = n_potential_S + 1
                        end if
                    end do
                end do
                
                
                do i = active_lo, active_hi
                    sl(i) = r(i) < -M_threshold
                    sh(i) = r(i) > M_threshold
                end do




            end if






            ! Rebuild current active signs. Values outside the exact active
            ! window are never read by the sparse packer below.
            if (force_full_sample) then
                do i = active_lo, active_hi
                    sl(i) = .false.
                    sh(i) = .false.
                end do
            end if

            ! Always force previous H observations into the retained subsample.
            ! Zero-weight previous-H rows are basis padding, not active screened rows.
            if (.not. use_independent_init) then
                do k = 1, 2*(nvar+1)
                    idx = H_seq(eva_t-1, k)
                    sl(idx) = .false.
                    sh(idx) = .false.
                end do
            end if

            ! Recompute all counts and retained indices after force-H. The
            ! active window is already in original-ID order. Only previous-H
            ! observations outside that window need a small sorted exception
            ! list; merging these two ordered sources exactly reproduces the
            ! former full 1:m packing order.
            n_forced_h = 0
            if (.not. use_independent_init) then
                do k = 1, 2*(nvar+1)
                    idx = H_seq(eva_t-1, k)
                    if (idx < active_lo .or. idx > active_hi) then
                        duplicate_h = .false.
                        do j = 1, n_forced_h
                            if (forced_h_ids(j) == idx) duplicate_h = .true.
                        end do
                        if (.not. duplicate_h) then
                            n_forced_h = n_forced_h + 1
                            forced_h_ids(n_forced_h) = idx
                        end if
                    end if
                end do
                do k = 2, n_forced_h
                    insert_id = forced_h_ids(k)
                    insert_pos = k - 1
                    do while (insert_pos >= 1)
                        if (forced_h_ids(insert_pos) <= insert_id) exit
                        forced_h_ids(insert_pos + 1) = forced_h_ids(insert_pos)
                        insert_pos = insert_pos - 1
                    end do
                    forced_h_ids(insert_pos + 1) = insert_id
                end do
            end if

            n_sl = 0
            n_sh = 0
            sum_w_sl = 0.0d0
            sum_w_sh = 0.0d0
            glob_wx = 0.0d0
            glob_wy = 0.0d0
            ghib_wx = 0.0d0
            ghib_wy = 0.0d0
            ms_org = 0
            k = 1
            do i = active_lo, active_hi
                do while (k <= n_forced_h)
                    if (forced_h_ids(k) >= i) exit
                    ms_org = ms_org + 1
                    idx_not_jl_or_jh(ms_org) = forced_h_ids(k)
                    k = k + 1
                end do
                if (k <= n_forced_h) then
                    if (forced_h_ids(k) == i) then
                        ms_org = ms_org + 1
                        idx_not_jl_or_jh(ms_org) = i
                        k = k + 1
                    else if (.not. (sl(i) .or. sh(i))) then
                        ms_org = ms_org + 1
                        idx_not_jl_or_jh(ms_org) = i
                    end if
                else if (.not. (sl(i) .or. sh(i))) then
                    ms_org = ms_org + 1
                    idx_not_jl_or_jh(ms_org) = i
                end if
                if (sl(i)) then
                    n_sl = n_sl + 1
                    sum_w_sl = sum_w_sl + w(i)
                    do j = 1, 2*(nvar+1)
                        glob_wx(j) = glob_wx(j) + A(i, j) * w(i)
                    end do
                    glob_wy = glob_wy + y(i) * w(i)
                else if (sh(i)) then
                    n_sh = n_sh + 1
                    sum_w_sh = sum_w_sh + w(i)
                    do j = 1, 2*(nvar+1)
                        ghib_wx(j) = ghib_wx(j) + A(i, j) * w(i)
                    end do
                    ghib_wy = ghib_wy + y(i) * w(i)
                end if
            end do
            do while (k <= n_forced_h)
                ms_org = ms_org + 1
                idx_not_jl_or_jh(ms_org) = forced_h_ids(k)
                k = k + 1
            end do
            has_sl_agg = (sum_w_sl > 0.0d0)
            has_sh_agg = (sum_w_sh > 0.0d0)
            ms = ms_org

            if (ms_org < 2*(nvar+1)) then
                if (debug_active) write(6, *) 'WARNING: Only', ms_org, 'uncertain observations at eva_t=', eva_t
                if (force_full_sample) then
                    call set_failure(5, eva_t)
                    return
                end if
                call record_repair()
                mmm_thresh = 2.0d0 * mmm_thresh
                not_new_sl_sh = .true.
                if (preprocessing_attempts >= max_empty_pivot_retries) then
                    force_full_sample = .true.
                end if
                cycle
            end if

            if (debug_active) then
                write(6, *) 'FORTRAN: Preprocessing partition:'
                write(6, '(A,F10.6)') '  M_threshold = ', M_threshold
                write(6, '(A,I4,A,I4,A,I4)') '  n_sl = ', n_sl, ', n_sh = ', n_sh, ', ms = ', ms
            end if

     


            ms_org = ms  ! Store the original subsample size before adding aggregated obs

            
            ! Extract subsample weights
            do i = 1, ms
                ws(i) = w(idx_not_jl_or_jh(i))
            end do

            
            ! Initialize for eva_t = 2
            if (eva_t == 2 .or. use_independent_init) then
                do i = 1, ms
                    do j = 1, 2*(nvar+1)
                        gammaxs_temp(i, j) = A(idx_not_jl_or_jh(i), j)
                    end do
                    bs_temp(i) = y(idx_not_jl_or_jh(i))
                end do




            end if
            
            ! Add aggregated observations if JL is not empty
            if (has_sl_agg) then
                do j = 1, 2*(nvar+1)
                    gammaxs_temp(m + 1, j) = glob_wx(j)
                end do
                bs_temp(m + 1) = glob_wy
                ms = ms + 1
                ws(ms) = 1.0d0
              



            end if
            
            ! Add aggregated observations if JH is not empty
            if (has_sh_agg) then
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
                        if (debug_active) then
                            write(6, *) 'WARNING: H observation', H_seq(eva_t-1, i), &
                                'not in subsample at eva_t=', eva_t
                        end if
                    else
                        k = k + 1
                    end if
                end if
            end do

            ! Check if we have enough valid H observations
            if ((.not. use_independent_init) .and. k < 2*(nvar+1)) then
                if (debug_active) then
                    write(6, *) 'WARNING: Only', k, 'valid H observations out of', 2*(nvar+1), &
                        'at eva_t=', eva_t
                end if
                use_independent_init = .true.
                independent_trigger = 2
                init_mode(eva_t) = 2
                init_trigger(eva_t) = independent_trigger
            end if

            if (.not. use_independent_init) then
                call matrix_inverse_2p(A(H_seq(eva_t-1, :), :), xhinv, 2*(nvar+1), inv_info)
                if (inv_info /= 0) then
                    use_independent_init = .true.
                    independent_trigger = 3
                    init_mode(eva_t) = 2
                    init_trigger(eva_t) = independent_trigger
                end if
            end if

            if (use_independent_init) then
                do i = 1, ms_org
                    do j = 1, 2*(nvar+1)
                        gammaxs_temp(i, j) = A(idx_not_jl_or_jh(i), j)
                    end do
                    bs_temp(i) = y(idx_not_jl_or_jh(i))
                end do
                do i = 1, 2*(nvar+1)
                    beta_prev(i) = beta_full_est(eva_t-1, i)
                end do
                if (.not. first_n_sub_recorded) then
                    first_n_sub(eva_t) = ms
                    first_n_sub_recorded = .true.
                end if
                call run_fresh_initialization(beta_prev, ms_org, ms, shifted_init_success)
                if (.not. shifted_init_success) then
                    init_mode(eva_t) = 3
                    call run_full_active_recovery(full_recovery_success)
                    if (.not. full_recovery_success) then
                        call set_failure(2, eva_t)
                        return
                    end if
                end if
                simplex_converged = .true.
                no_pivot_flag = .false.
                goto 24430
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
                if (has_sl_agg .and. has_sh_agg) then
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
                    




                else if (has_sl_agg) then
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
                    
                else if (has_sh_agg) then
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
                call matrix_inverse_2p(gammaxs(H_indices, :), xhinv, 2*(nvar+1), inv_info)
                
                if (inv_info /= 0) then
                    if (debug_active) write(6,*) 'X(h) is singular at eva_t = ', eva_t, ', info = ', inv_info
                    if (force_full_sample) then
                        call set_failure(4, eva_t)
                        return
                    end if
                    empty_pivot_count = empty_pivot_count + 1
                    call record_repair()
                    mmm_thresh = 2.0d0 * mmm_thresh
                    not_new_sl_sh = .true.
                    if (empty_pivot_count >= max_empty_pivot_retries) then
                        force_full_sample = .true.
                    end if
                    cycle
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
                if (has_sl_agg .and. has_sh_agg) then
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
                else if (has_sl_agg) then
                    ! Only v_L
                    do j = 1, 2*(nvar+1)
                        Pxhbarxhinv(n_idpos + n_idneg + 1, j) = 0.0d0
                        do k = 1, 2*(nvar+1)
                            Pxhbarxhinv(n_idpos + n_idneg + 1, j) = Pxhbarxhinv(n_idpos + n_idneg + 1, j) - &
                                                                   gammaxs(ms, k) * xhinv(k, j)
                        end do
                    end do
                else if (has_sh_agg) then
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
                if (has_sl_agg) k = k + 1
                if (has_sh_agg) k = k + 1
                
                do i = 1, k
                    do j = 1, 2*(nvar+1)
                        gammaxs(2*(nvar+1) + i, j) = -Pxhbarxhinv(i, j)
                    end do
                end do
                
                ! Last row is the objective function row
                do j = 1, 2*(nvar+1)
                    if (H_indices(j) < 1 .or. H_indices(j) > ms) then
                        if (debug_active) then
                            write(6, *) 'ERROR: Invalid H_indices(', j, ')=', H_indices(j), ' at eva_t=', eva_t
                            write(6, *) '  Valid range is 1 to', ms
                        end if
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
                if (has_sl_agg .and. has_sh_agg) then
                    bs(ms-1) = -bs_temp(m+1)
                    do j = 1, 2*(nvar+1)
                        bs(ms-1) = bs(ms-1) - Pxhbarxhinv(n_idpos + n_idneg + 1, j) * bs_temp(H_indices(j))
                    end do
                    
                    bs(ms) = bs_temp(m+2)
                    do j = 1, 2*(nvar+1)
                        bs(ms) = bs(ms) - Pxhbarxhinv(n_idpos + n_idneg + 2, j) * bs_temp(H_indices(j))
                    end do
                else if (has_sl_agg) then
                    bs(ms) = -bs_temp(m+1)
                    do j = 1, 2*(nvar+1)
                        bs(ms) = bs(ms) - Pxhbarxhinv(n_idpos + n_idneg + 1, j) * bs_temp(H_indices(j))
                    end do
                else if (has_sh_agg) then
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
                if (has_sl_agg .and. has_sh_agg) then



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
                    
                else if (has_sl_agg) then
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
                    


                else if (has_sh_agg) then





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
                if (has_sl_agg .and. has_sh_agg) then
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




                    
                else if (has_sl_agg) then
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
                    
                else if (has_sh_agg) then
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
                        if (debug_active) then
                            write(6, *) 'ERROR: Invalid H_indices(', j, ')=', H_indices(j), ' at eva_t=', eva_t
                            write(6, *) '  Valid range is 1 to', ms_org
                        end if
                        ! This should not happen if H observations are correctly forced into subsample
                        call set_failure(3, eva_t)
                        return
                    end if



   


                    do i = 1, n_hbar_rows  ! Use the descriptive name
                        gammaxs(ms + 1, j) = gammaxs(ms + 1, j) + lambda(i) * Pxhbarxhinv(i, j)
                    end do
                end do

  






                bs(ms+1) = 0.0d0
                


            end if  ! End of eva_t == 2 vs eva_t >= 3
            

   

            ! Invariant checks before simplex
            if (ms < 2*(nvar+1)) then
                if (debug_active) write(6, *) 'WARNING: Reduced problem too small at eva_t=', eva_t, ', ms=', ms
                if (force_full_sample) then
                    call set_failure(5, eva_t)
                    return
                end if
                empty_pivot_count = empty_pivot_count + 1
                call record_repair()
                mmm_thresh = 2.0d0 * mmm_thresh
                not_new_sl_sh = .true.
                if (empty_pivot_count >= max_empty_pivot_retries) then
                    force_full_sample = .true.
                end if
                cycle
            end if

            ! Check that all beta columns are basic
            if (.not. all(freevarrow(1:2*(nvar+1)))) then
                if (debug_active) write(6, *) 'ERROR: Beta not basic at eva_t=', eva_t
                do i = 1, 2*(nvar+1)
                    if (.not. freevarrow(i)) then
                        if (debug_active) write(6, *) '  Beta column', i, 'is not basic'
                    end if
                end do
            end if





            ! Simplex iterations for the reduced problem
            !iter = 0

            if (.not. first_n_sub_recorded) then
                first_n_sub(eva_t) = ms
                first_n_sub_recorded = .true.
            end if


            do while (iter < maxit)
                total_preprocessing_loops = total_preprocessing_loops + 1
                ! Add a safety check to prevent infinite preprocessing:
                if (preprocessing_attempts > 100) then
                    if (debug_active) then
                        write(6, *) 'ERROR: Preprocessing stuck in infinite loop at eva_t=', eva_t
                        write(6, *) 'M_threshold:', M_threshold
                        write(6, *) 'bad_signs:', bad_signs
                        write(6, *) 'ms:', ms
                    end if
                    call set_failure(2, eva_t)
                    return
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
                            if (debug_active) then
                                write(6, *) 'ERROR: r2 mask failed at eva_t=', eva_t
                                write(6, *) 'Column', i, 'has r2=0 but rr(2,i)=', rr(2, i)
                            end if
                            call set_failure(5, eva_t)
                            return
                        end if
                    end do
                end if

                ! Find location of minimum
                call minloc2d(rr, rrl, tsep, t_rr)




                ! ADD SAFETY VALVE HERE:
                if (iter > 50000) then



                    if (debug_active) then
                        write(6, *) 'SAFETY: Forcing exit from simplex at eva_t=', eva_t
                        write(6, *) 'This prevents memory exhaustion from excessive iterations'
                    end if
                    no_pivot_flag = .true.
                    exit  ! Break out of the simplex loop
                end if


                




                if (rrl >= -tol) then

  


                    simplex_converged = .true.  ! Mark as converged
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
                    yy(1:ms+1) = gammaxs(1:ms+1, t_rr)
                else
                    yy(1:ms+1) = -gammaxs(1:ms+1, t_rr)
                end if
                

 

                ! Step 5: Choose leaving variable (ratio test)
                min_k = huge(1.0d0)
                k = 0
                
                
                ! First, find all valid ratios
                do i = 1, ms
                    if (yy(i) > pivot_tol .and. .not. freevarrow(i)) then
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
                    if (debug_active) write(6, *) 'The problem is unbounded, doubling m at time', eva_t
                    

     
                    no_pivot_flag = .true.
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
                    if (debug_active) then
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
                        write(6, *) '  Non-free vars with yy > pivot_tol:'
                        do i = 1, ms
                            if (.not. freevarrow(i) .and. yy(i) > pivot_tol) then
                                write(6, '(A,I3,A,F12.6,A,F12.6,A,F12.6,A,I6)') &
                                    '    Row ', i, ': yy=', yy(i), ', bs=', bs(i), &
                                    ', ratio=', bs(i)/yy(i), ', IBs=', IBs(i)
                            end if
                        end do

                        ! Check for negative ratios
                        write(6, *) '  Non-free vars with yy > pivot_tol but bs < 0 (infeasible):'
                        k = 0
                        do i = 1, ms
                            if (.not. freevarrow(i) .and. yy(i) > pivot_tol .and. bs(i) < 0.0d0) then
                                k = k + 1
                                if (k <= 5) then
                                    write(6, '(A,I3,A,F12.6,A,F12.6)') &
                                        '    Row ', i, ': yy=', yy(i), ', bs=', bs(i)
                                end if
                            end if
                        end do
                        if (k > 5) write(6, *) '    ... and', k-5, 'more'
                        write(6, *) '  Current objective row (first 8):', (gammaxs(ms+1, j), j=1, 2*(nvar+1))

                        write(6, *) 'This indicates an unbounded or degenerate problem'
                        write(6, *) 'Diagnostic info:'
                        write(6, *) '  Number of positive yy:', count(yy(1:ms) > pivot_tol)
                        write(6, *) '  Number of non-free vars:', count(.not. freevarrow(1:ms))
                        write(6, *) '  Min yy:', minval(yy(1:ms)), 'Max yy:', maxval(yy(1:ms))
                    end if
                    no_pivot_flag = .true.
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
                    if (debug_active) then
                        write(6, *) 'ERROR: No leaving variable found (k=0) at eva_t=', eva_t, ', iter=', iter
                        write(6, *) 'min_k=', min_k
                        write(6, *) 'This means no feasible pivot exists'
                    end if
                    no_pivot_flag = .true.
                    unbounded_detected = .true.
                    exit
                end if


                ! Perform pivot only when the selected pivot is safely away from zero.
                if (abs(yy(k)) <= pivot_tol) then
                    if (debug_active) then
                        write(6, *) 'WARNING: Near-zero pivot element at eva_t=', eva_t, ', iter=', iter
                        write(6, *) 'Original pivot value:', yy(k), 'at position k=', k
                        write(6, *) 'Terminating simplex attempt without accepting candidate'
                    end if
                    no_pivot_flag = .true.
                    exit
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
                iter_attempt = iter_attempt + 1



            end do
            



            ! ADD ITERATION TRACKING HERE:
            total_simplex_iterations = total_simplex_iterations + iter_attempt
            if (iter > max_iter_at_any_t) then
                max_iter_at_any_t = iter
            end if

            ! Log concerning patterns:
            if (iter > 10000) then
                if (debug_active) then
                    write(6, *) 'WARNING: Excessive iterations at eva_t=', eva_t
                    write(6, *) 'Iterations:', iter
                    write(6, *) 'Subsample size ms:', ms
                    write(6, *) 'Consider the following:'
                    write(6, *) '- The problem may be degenerate'
                    write(6, *) '- The tolerance may be too tight'
                    write(6, *) '- There may be cycling in the simplex'
                end if
            end if
            if (no_pivot_flag .or. (.not. simplex_converged)) then
                if (debug_active) write(6, *) 'WARNING: Reduced simplex did not converge at eva_t =', eva_t
                if (force_full_sample) then
                    call set_failure(2, eva_t)
                    return
                end if
                empty_pivot_count = empty_pivot_count + 1
                call record_repair()
                mmm_thresh = 2.0d0 * mmm_thresh
                not_new_sl_sh = .true.
                if (empty_pivot_count >= max_empty_pivot_retries .or. iter >= maxit) then
                    force_full_sample = .true.
                end if
                cycle
            end if

24430       continue
            do i = 1, 2*(nvar+1)
                estimate(i) = bs(i)
            end do

            ! Extract the final H from the converged reduced simplex basis.
            do i = 1, 2*(nvar+1)
                j = r1(i) - 2 - 2*nvar
                if (j > 0 .and. j <= ms_org) then
                    H_indices(i) = idx_not_jl_or_jh(j)
                else
                    H_indices(i) = 0
                end if
            end do

            ! Independent initialization and full-active recovery stay dense in
            ! v38. The ordinary default path materializes exactly the rows read
            ! by current verification/caching and next-point screening.
            candidate_residual_full = (.not. lazy_residual_mode) .or. use_independent_init
            refit_used = .false.
            if (always_same_h_refit .and. tvcqr_H_basis_valid(H_indices, A, m, 2*(nvar+1))) then
                same_h_refit_attempted(eva_t) = 1
                if (candidate_residual_full) then
                    call same_h_refit_tvcqr(H_indices, A, y, m, 2*(nvar+1), res_tol, &
                                            estimate_refit, r_refit, same_h_ok)
                else
                    ! Keep the same dummy-array arithmetic context as v37's
                    ! dense same-H helper while shortening only its row domain.
                    call same_h_refit_lazy_tvcqr(H_indices, H_seq(eva_t-1, :), A, y, &
                                                 m, 2*(nvar+1), res_tol, estimate_refit, &
                                                 r, same_h_ok)
                end if
                if (same_h_ok) then
                    do i = 1, 2*(nvar+1)
                        estimate(i) = estimate_refit(i)
                    end do
                    if (candidate_residual_full) then
                        do i = 1, m
                            r(i) = r_refit(i)
                        end do
                    end if
                    refit_used = .true.
                end if
            end if

            ! If same-H refit failed or was not attempted, form the candidate
            ! image from the tableau estimate. The dense branch is the original
            ! full-m recurrence; the lazy branch uses the same recurrence and
            ! coefficient order on its exact required domain.
            if (.not. refit_used) then
                if (candidate_residual_full) then
                    do i = 1, m
                        r(i) = y(i)
                        do j = 1, 2*(nvar+1)
                            r(i) = r(i) - A(i, j) * estimate(j)
                        end do
                    end do
                else
                    call materialize_candidate_residual(estimate, H_seq(eva_t-1, :), &
                                                        H_indices, .false.)
                end if
            end if




            ! Literally inspect every omitted current-active observation.
            bad_signs = 0
            n_sure_signs = n_sl + n_sh  ! Total sure-sign observations

            do i = active_lo, active_hi
                if ((r(i) <= 0.0d0) .and. sh(i)) then
                    bad_signs = bad_signs + 1
                end if
                if ((r(i) >= 0.0d0) .and. sl(i)) then
                    bad_signs = bad_signs + 1
                end if
            end do




            ! Handle bad signs
            if (bad_signs > 0) then
                call handle_current_bad_signs()

            else
                ! No bad signs - we've reached optimality
                accept_subsample = certify_tvcqr_candidate(H_indices, r, A, m, 2*(nvar+1), res_tol)
                if ((.not. accept_subsample) .and. (.not. always_same_h_refit)) then
                    if (tvcqr_H_basis_valid(H_indices, A, m, 2*(nvar+1))) then
                        same_h_refit_attempted(eva_t) = 1
                        call same_h_refit_tvcqr(H_indices, A, y, m, 2*(nvar+1), res_tol, &
                                                estimate_refit, r_refit, same_h_ok)
                        if (same_h_ok) then
                            do i = 1, 2*(nvar+1)
                                estimate(i) = estimate_refit(i)
                            end do
                            do i = 1, m
                                r(i) = r_refit(i)
                            end do
                            refit_used = .true.
                            bad_signs = 0
                            do i = active_lo, active_hi
                                if ((r(i) <= 0.0d0) .and. sh(i)) then
                                    bad_signs = bad_signs + 1
                                end if
                                if ((r(i) >= 0.0d0) .and. sl(i)) then
                                    bad_signs = bad_signs + 1
                                end if
                            end do
                            if (bad_signs == 0) then
                                accept_subsample = certify_tvcqr_candidate(H_indices, r, A, m, &
                                                                           2*(nvar+1), res_tol)
                            end if
                        end if
                    end if
                end if
                if (bad_signs > 0) then
                    call handle_current_bad_signs()
                    cycle
                end if
                if (.not. accept_subsample) then
                    if (((n_sl == 0) .and. (n_sh == 0)) .or. force_full_sample) then
                        call set_failure(1, eva_t)
                        return
                    end if
                    call record_repair()
                    mmm_thresh = 2.0d0 * mmm_thresh
                    not_new_sl_sh = .true.
                    cycle
                end if
                not_optimal = .false.
                if (refit_used) then
                    same_h_refit_recovered(eva_t) = 1
                end if
                do i = 1, 2*(nvar+1)
                    bs(i) = estimate(i)
                end do
                do i = 1, 2*(nvar+1)
                    r(H_indices(i)) = 0.0d0
                end do
                





                ! Extract H for next iteration
                do i = 1, 2*(nvar+1)
                    if (H_indices(i) > 0) then
                        H_seq(eva_t, i) = H_indices(i)
                    else
                        ! This matches R's behavior - r1 should always give valid indices
                        ! If not, there's a bug in the simplex algorithm
                        if (debug_active) then
                            write(6, *) 'ERROR: Invalid certified H index at eva_t=', eva_t, ', i=', i
                            write(6, *) 'r1(i)=', r1(i), ', ms_org=', ms_org
                        end if
                        call set_failure(3, eva_t)
                        return
                    end if




                end do
                


               


                ! If not the last time point, save information for reuse
                if (eva_t < m) then



  


  
                    ! Calculate ms_org for saving
                    if (has_sl_agg .and. has_sh_agg) then
                        ms_org = ms - 2
                    else if (has_sl_agg .or. has_sh_agg) then
                        ms_org = ms - 1
                    else
                        ms_org = ms
                    end if
                    
                    do i = 1, 2*(nvar+1)
                        do j = 1, 2*(nvar+1)
                            gammaxs_temp(i, j) = gammaxs(i, j)
                        end do
                        bs_temp(i) = bs(i)
                    end do
                    
                    ! Initialize counters
                    n_Hbar_pos = 0
                    n_Hbar_neg = 0
                    
                    ! Only process Hbar information if there are non-aggregated observations
                    if (ms_org > 0) then






                        ! Process the IBs to extract Hbar information
                        do idx = 2*(nvar+1) + 1, ms_org
                            i = idx - 2*(nvar+1)! idx = 2*(nvar+1) + i
                            if (IBs(idx) > 2*(nvar+1)+ms) then
                                ! This is a v variable
                                n_Hbar_neg = n_Hbar_neg + 1


               



                                id_gammaxs_Hbar(i) = IBs(idx) - 2*(nvar+1) - ms
                                if (id_gammaxs_Hbar(i) < 1 .or. id_gammaxs_Hbar(i) > ms_org) then
                                    call set_failure(5, eva_t)
                                    return
                                end if
                                idx_Hbar_neg(n_Hbar_neg) = idx_not_jl_or_jh(id_gammaxs_Hbar(i))
                                if (idx_Hbar_neg(n_Hbar_neg) < 1 .or. idx_Hbar_neg(n_Hbar_neg) > m) then
                                    call set_failure(5, eva_t)
                                    return
                                end if
                                do j = 1, 2*(nvar+1)
                                    gammaxs_neg(n_Hbar_neg, j) = gammaxs(idx, j)
                                end do
                                bs_neg(n_Hbar_neg) = max(-r(idx_Hbar_neg(n_Hbar_neg)), 0.0d0)
                            else
                                ! This is a u variable
                                n_Hbar_pos = n_Hbar_pos + 1
                                id_gammaxs_Hbar(i) = IBs(idx) - 2*(nvar+1)
                                if (id_gammaxs_Hbar(i) < 1 .or. id_gammaxs_Hbar(i) > ms_org) then
                                    call set_failure(5, eva_t)
                                    return
                                end if
                                idx_Hbar_pos(n_Hbar_pos) = idx_not_jl_or_jh(id_gammaxs_Hbar(i))
                                if (idx_Hbar_pos(n_Hbar_pos) < 1 .or. idx_Hbar_pos(n_Hbar_pos) > m) then
                                    call set_failure(5, eva_t)
                                    return
                                end if
                                do j = 1, 2*(nvar+1)
                                    gammaxs_pos(n_Hbar_pos, j) = gammaxs(idx, j)
                                end do
                                bs_pos(n_Hbar_pos) = max(r(idx_Hbar_pos(n_Hbar_pos)), 0.0d0)
                            end if
                        end do

     


                    end if

                end if
                
            end if
            
        end do  ! End while not_optimal
        
        ! Store results for this time point
        it_num(eva_t) = iter
        final_n_sub(eva_t) = ms
        M_out = M_threshold
        

        ! Compute theta_ll_est for current eva_t
        do j = 1, nvar+1
            theta_ll_est(eva_t, j) = estimate(j) + (dble(eva_t)/dble(m)) * estimate(nvar+1+j)
        end do
        do j = 1, 2*(nvar+1)
            beta_full_est(eva_t, j) = estimate(j)
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





        if (candidate_residual_full) then
            do i = 1, m
                r_prev(i) = r(i)
                if (store_residual) then
                    residual_est((i - 1) * m + eva_t) = r(i)
                end if
            end do
            r_prev_full = .true.
            r_prev_valid_lo = 1
            r_prev_valid_hi = m
            r_prev_extra_count = 0
        else
            ! Commit only the candidate's proven valid domain. The old accepted
            ! image remains untouched elsewhere and is never read without a
            ! descriptor check/promotion.
            do i = candidate_residual_lo, candidate_residual_hi
                r_prev(i) = r(i)
            end do
            do k = 1, 2*(nvar+1)
                idx = H_seq(eva_t, k)
                r_prev(idx) = r(idx)
            end do
            r_prev_full = .false.
            r_prev_valid_lo = candidate_residual_lo
            r_prev_valid_hi = candidate_residual_hi
            r_prev_extra_count = 0
            do k = 1, 2*(nvar+1)
                call append_previous_extra(H_seq(eva_t, k))
            end do
        end if


        
    end do  ! End of eva_t loop
    
    
    deallocate(A)
    deallocate(gammax)
    deallocate(gammaxs_temp)
    deallocate(bs_temp)
    deallocate(gammaxs)
    deallocate(bs)
    deallocate(inverse_work)

contains

    subroutine find_positive_window(eval_idx, lo_hint, hi_hint, lo_out, hi_out, found)
        implicit none
        integer, intent(in) :: eval_idx, lo_hint, hi_hint
        integer, intent(out) :: lo_out, hi_out
        logical, intent(out) :: found
        integer :: wi
        double precision :: eval_loc, weight_loc

        found = .false.
        lo_out = 1
        hi_out = 0
        if (eval_idx < 1 .or. eval_idx > m) return

        eval_loc = dble(eval_idx) / dble(m)
        lo_out = max(1, lo_hint)
        do while (lo_out <= m)
            if (abs(eval_loc - time_index(lo_out)) <= h) then
                weight_loc = 0.75d0 * (1.0d0 - ((eval_loc - time_index(lo_out))/h)**2)
            else
                weight_loc = 0.0d0
            end if
            if (weight_loc > 0.0d0) exit
            lo_out = lo_out + 1
        end do
        if (lo_out > m) return

        hi_out = max(hi_hint, lo_out - 1)
        do while (hi_out < m)
            wi = hi_out + 1
            if (abs(eval_loc - time_index(wi)) <= h) then
                weight_loc = 0.75d0 * (1.0d0 - ((eval_loc - time_index(wi))/h)**2)
            else
                weight_loc = 0.0d0
            end if
            if (weight_loc <= 0.0d0) exit
            hi_out = wi
        end do
        found = (hi_out >= lo_out)
    end subroutine find_positive_window

    logical function previous_range_covered(lo_idx, hi_idx)
        implicit none
        integer, intent(in) :: lo_idx, hi_idx

        if (hi_idx < lo_idx) then
            previous_range_covered = .true.
        else if (r_prev_full) then
            previous_range_covered = .true.
        else
            previous_range_covered = (lo_idx >= r_prev_valid_lo .and. &
                                      hi_idx <= r_prev_valid_hi)
        end if
    end function previous_range_covered

    logical function previous_id_covered(row_idx)
        implicit none
        integer, intent(in) :: row_idx
        integer :: ei

        previous_id_covered = .false.
        if (row_idx < 1 .or. row_idx > m) return
        if (r_prev_full) then
            previous_id_covered = .true.
            return
        end if
        if (row_idx >= r_prev_valid_lo .and. row_idx <= r_prev_valid_hi) then
            previous_id_covered = .true.
            return
        end if
        do ei = 1, r_prev_extra_count
            if (r_prev_extra(ei) == row_idx) then
                previous_id_covered = .true.
                return
            end if
        end do
    end function previous_id_covered

    subroutine promote_previous_residual(previous_eval)
        implicit none
        integer, intent(in) :: previous_eval
        integer :: pi, pj, ph

        if (r_prev_full) return
        do pi = 1, m
            r_prev(pi) = y(pi)
            do pj = 1, 2*(nvar+1)
                r_prev(pi) = r_prev(pi) - A(pi, pj) * beta_full_est(previous_eval, pj)
            end do
        end do
        do pj = 1, 2*(nvar+1)
            ph = H_seq(previous_eval, pj)
            if (ph >= 1 .and. ph <= m) r_prev(ph) = 0.0d0
        end do
        r_prev_full = .true.
        r_prev_valid_lo = 1
        r_prev_valid_hi = m
        r_prev_extra_count = 0
    end subroutine promote_previous_residual

    subroutine materialize_candidate_residual(beta_owner, H_previous, H_candidate, make_full)
        implicit none
        double precision, intent(in) :: beta_owner(2*(nvar+1))
        integer, intent(in) :: H_previous(2*(nvar+1)), H_candidate(2*(nvar+1))
        logical, intent(in) :: make_full
        integer :: ri, rj, rk, rh
        logical :: already_materialized

        if (make_full) then
            do ri = 1, m
                r(ri) = y(ri)
                do rj = 1, 2*(nvar+1)
                    r(ri) = r(ri) - A(ri, rj) * beta_owner(rj)
                end do
            end do
            return
        end if

        do ri = candidate_residual_lo, candidate_residual_hi
            r(ri) = y(ri)
            do rj = 1, 2*(nvar+1)
                r(ri) = r(ri) - A(ri, rj) * beta_owner(rj)
            end do
        end do

        do ri = 1, 2*(nvar+1)
            rh = H_previous(ri)
            if (rh >= 1 .and. rh <= m) then
                if (rh < candidate_residual_lo .or. rh > candidate_residual_hi) then
                    r(rh) = y(rh)
                    do rj = 1, 2*(nvar+1)
                        r(rh) = r(rh) - A(rh, rj) * beta_owner(rj)
                    end do
                end if
            end if
        end do

        do ri = 1, 2*(nvar+1)
            rh = H_candidate(ri)
            if (rh >= 1 .and. rh <= m) then
                if (rh < candidate_residual_lo .or. rh > candidate_residual_hi) then
                    already_materialized = .false.
                    do rk = 1, 2*(nvar+1)
                        if (H_previous(rk) == rh) then
                            already_materialized = .true.
                            exit
                        end if
                    end do
                    if (.not. already_materialized) then
                        r(rh) = y(rh)
                        do rj = 1, 2*(nvar+1)
                            r(rh) = r(rh) - A(rh, rj) * beta_owner(rj)
                        end do
                    end if
                end if
            end if
        end do
    end subroutine materialize_candidate_residual

    subroutine append_previous_extra(row_idx)
        implicit none
        integer, intent(in) :: row_idx
        integer :: ei

        if (row_idx < 1 .or. row_idx > m) return
        if (row_idx >= r_prev_valid_lo .and. row_idx <= r_prev_valid_hi) return
        do ei = 1, r_prev_extra_count
            if (r_prev_extra(ei) == row_idx) return
        end do
        if (r_prev_extra_count < 4*(nvar+1)) then
            r_prev_extra_count = r_prev_extra_count + 1
            r_prev_extra(r_prev_extra_count) = row_idx
        end if
    end subroutine append_previous_extra

    subroutine validate_previous_H(H_idx, valid)
        implicit none
        integer, intent(in) :: H_idx(2*(nvar+1))
        logical, intent(out) :: valid
        integer :: vi, vj

        valid = .true.
        do vi = 1, 2*(nvar+1)
            if (H_idx(vi) < 1 .or. H_idx(vi) > m) then
                valid = .false.
                return
            end if
            do vj = vi + 1, 2*(nvar+1)
                if (H_idx(vi) == H_idx(vj)) then
                    valid = .false.
                    return
                end if
            end do
        end do
    end subroutine validate_previous_H

    subroutine reorder_fresh_tableau(theta_offset, n_individual, n_total, success)
        implicit none
        double precision, intent(in) :: theta_offset(2*(nvar+1))
        integer, intent(in) :: n_individual, n_total
        logical, intent(out) :: success
        double precision :: gx_tmp(m+3, 2*(nvar+1)), bv_tmp(m+3)
        integer :: IB_tmp(m+3)
        logical :: fvr_tmp(m+3), used_row(m+2)
        integer :: ci, cj, src, dest, qdim

        success = .false.
        qdim = 2*(nvar+1)
        used_row = .false.

        do ci = 1, qdim
            src = 0
            do cj = 1, n_total
                if (IBs(cj) == ci) then
                    src = cj
                    exit
                end if
            end do
            if (src < 1 .or. src > n_individual) return
            used_row(src) = .true.
            do cj = 1, qdim
                gx_tmp(ci, cj) = gammaxs(src, cj)
            end do
            bv_tmp(ci) = bs(src) + theta_offset(ci)
            IB_tmp(ci) = IBs(src)
            fvr_tmp(ci) = freevarrow(src)
        end do

        dest = qdim
        do src = 1, n_individual
            if (.not. used_row(src)) then
                dest = dest + 1
                do cj = 1, qdim
                    gx_tmp(dest, cj) = gammaxs(src, cj)
                end do
                bv_tmp(dest) = bs(src)
                IB_tmp(dest) = IBs(src)
                fvr_tmp(dest) = freevarrow(src)
            end if
        end do
        do src = n_individual + 1, n_total
            dest = dest + 1
            do cj = 1, qdim
                gx_tmp(dest, cj) = gammaxs(src, cj)
            end do
            bv_tmp(dest) = bs(src)
            IB_tmp(dest) = IBs(src)
            fvr_tmp(dest) = freevarrow(src)
        end do
        if (dest /= n_total) return

        do cj = 1, qdim
            gx_tmp(n_total + 1, cj) = gammaxs(n_total + 1, cj)
        end do
        bv_tmp(n_total + 1) = bs(n_total + 1)
        IB_tmp(n_total + 1) = IBs(n_total + 1)
        fvr_tmp(n_total + 1) = freevarrow(n_total + 1)

        do ci = 1, n_total + 1
            do cj = 1, qdim
                gammaxs(ci, cj) = gx_tmp(ci, cj)
            end do
            bs(ci) = bv_tmp(ci)
            IBs(ci) = IB_tmp(ci)
            freevarrow(ci) = fvr_tmp(ci)
        end do
        success = .true.
    end subroutine reorder_fresh_tableau

    subroutine run_fresh_initialization(theta_offset, n_individual, n_total, success)
        implicit none
        double precision, intent(in) :: theta_offset(2*(nvar+1))
        integer, intent(in) :: n_individual, n_total
        logical, intent(out) :: success
        double precision :: shifted_rhs, delta_est(2*(nvar+1))
        integer :: si, sj, src, qdim, reduced_idx, iter_fresh, remaining
        logical :: no_pivot_fresh, converged_fresh, reorder_ok

        success = .false.
        qdim = 2*(nvar+1)
        if (n_individual < qdim .or. n_total < n_individual .or. n_total > m + 2) return

        do si = 1, n_total
            if (si <= n_individual) then
                src = si
            else if (has_sl_agg .and. si == n_individual + 1) then
                src = m + 1
            else
                src = m + 2
            end if

            shifted_rhs = bs_temp(src)
            do sj = 1, qdim
                shifted_rhs = shifted_rhs - gammaxs_temp(src, sj) * theta_offset(sj)
            end do

            if (src == m + 1) then
                if (shifted_rhs > 0.0d0) return
                do sj = 1, qdim
                    gammaxs(si, sj) = -gammaxs_temp(src, sj)
                end do
                bs(si) = -shifted_rhs
                IBs(si) = qdim + n_total + si
            else if (src == m + 2) then
                if (shifted_rhs < 0.0d0) return
                do sj = 1, qdim
                    gammaxs(si, sj) = gammaxs_temp(src, sj)
                end do
                bs(si) = shifted_rhs
                IBs(si) = qdim + si
            else if (shifted_rhs < 0.0d0) then
                do sj = 1, qdim
                    gammaxs(si, sj) = -gammaxs_temp(src, sj)
                end do
                bs(si) = -shifted_rhs
                IBs(si) = qdim + n_total + si
            else
                do sj = 1, qdim
                    gammaxs(si, sj) = gammaxs_temp(src, sj)
                end do
                bs(si) = shifted_rhs
                IBs(si) = qdim + si
            end if
            freevarrow(si) = (si > n_individual)
        end do

        IBs(n_total + 1) = 0
        freevarrow(n_total + 1) = .true.
        bs(n_total + 1) = 0.0d0
        do sj = 1, qdim
            gammaxs(n_total + 1, sj) = 0.0d0
            do si = 1, n_total
                if (IBs(si) > qdim .and. IBs(si) <= qdim + n_total) then
                    gammaxs(n_total + 1, sj) = gammaxs(n_total + 1, sj) - &
                        tau * ws(si) * gammaxs(si, sj)
                else if (IBs(si) > qdim + n_total) then
                    gammaxs(n_total + 1, sj) = gammaxs(n_total + 1, sj) - &
                        (1.0d0 - tau) * ws(si) * gammaxs(si, sj)
                end if
            end do
        end do

        do si = 1, qdim
            r1(si) = si
            r2(si) = 0
        end do
        rr = 0.0d0
        remaining = maxit - iter
        if (remaining <= 0) return
        call run_simplex_full_tvcqr(gammaxs, bs, IBs, freevarrow, r1, r2, rr, ws, &
                                    m + 3, n_total, qdim, tol, remaining, bland, &
                                    iter_fresh, no_pivot_fresh, converged_fresh)
        iter = iter + iter_fresh
        total_simplex_iterations = total_simplex_iterations + iter_fresh
        if (no_pivot_fresh .or. (.not. converged_fresh)) return

        delta_est = 0.0d0
        do si = 1, n_total
            if (IBs(si) >= 1 .and. IBs(si) <= qdim) then
                delta_est(IBs(si)) = bs(si)
            end if
        end do
        do si = 1, qdim
            estimate(si) = theta_offset(si) + delta_est(si)
            reduced_idx = r1(si) - qdim
            if (reduced_idx < 1 .or. reduced_idx > n_individual) return
            H_indices(si) = idx_not_jl_or_jh(reduced_idx)
        end do
        if (.not. tvcqr_H_basis_valid(H_indices, A, m, qdim)) return

        call reorder_fresh_tableau(theta_offset, n_individual, n_total, reorder_ok)
        if (.not. reorder_ok) return
        success = .true.
    end subroutine run_fresh_initialization

    subroutine run_full_active_recovery(success)
        implicit none
        logical, intent(out) :: success
        double precision :: zero_offset(2*(nvar+1))
        integer :: fi, fj

        success = .false.
        zero_offset = 0.0d0
        sl = .false.
        sh = .false.
        has_sl_agg = .false.
        has_sh_agg = .false.
        n_sl = 0
        n_sh = 0
        ms_org = m
        ms = m
        do fi = 1, m
            idx_not_jl_or_jh(fi) = fi
            not_jl_or_jh(fi) = .true.
            ws(fi) = w(fi)
            do fj = 1, 2*(nvar+1)
                gammaxs_temp(fi, fj) = A(fi, fj)
            end do
            bs_temp(fi) = y(fi)
        end do
        call run_fresh_initialization(zero_offset, m, m, success)
    end subroutine run_full_active_recovery

    subroutine run_simplex_full_tvcqr(gx, bv, IBv, fvr, r1v, r2v, rrv, wv, &
                                      ldgx, mv, p, tl, mxit, bld, iters, no_pivot, converged)
        implicit none
        integer, intent(in) :: ldgx, mv, p, mxit
        double precision, intent(inout) :: gx(ldgx, p), bv(mv+1)
        integer, intent(inout) :: IBv(mv+1), r1v(p), r2v(p)
        logical, intent(inout) :: fvr(mv+1)
        double precision, intent(inout) :: rrv(2, p)
        double precision, intent(in) :: wv(mv), tl
        logical, intent(in) :: bld
        integer, intent(out) :: iters
        logical, intent(out) :: no_pivot, converged
        double precision :: yyv(mv+1), eev(mv+1), kval(mv+1)
        double precision :: rrlv, min_kv, pivot_val
        integer :: ii, jj, kk, enter_col, enter_side, enter_var, idx_offset

        iters = 0
        no_pivot = .false.
        converged = .false.
        do while (iters < mxit)
            do ii = 1, p
                rrv(1, ii) = gx(mv+1, ii)
                if (r2v(ii) /= 0) then
                    idx_offset = r1v(ii) - p
                    if (idx_offset >= 1 .and. idx_offset <= mv) then
                        rrv(2, ii) = wv(idx_offset) - rrv(1, ii)
                    else
                        rrv(2, ii) = -rrv(1, ii)
                    end if
                else
                    rrv(2, ii) = 0.0d0
                    rrv(1, ii) = -abs(rrv(1, ii))
                end if
            end do

            rrlv = minval(rrv)
            if (rrlv >= -tl) then
                converged = .true.
                exit
            end if

            enter_col = 0
            enter_side = 0
            enter_var = huge(1)
            if (bld) then
                do jj = 1, p
                    if (rrv(1, jj) < -tl .and. r1v(jj) < enter_var) then
                        enter_var = r1v(jj)
                        enter_col = jj
                        enter_side = 1
                    end if
                end do
                if (enter_col == 0) then
                    do jj = 1, p
                        if (rrv(2, jj) < -tl .and. r2v(jj) < enter_var) then
                            enter_var = r2v(jj)
                            enter_col = jj
                            enter_side = 2
                        end if
                    end do
                end if
            else
                do jj = 1, p
                    do ii = 1, 2
                        if (abs(rrv(ii, jj) - rrlv) < tl) then
                            enter_col = jj
                            enter_side = ii
                            if (ii == 1) then
                                enter_var = r1v(jj)
                            else
                                enter_var = r2v(jj)
                            end if
                            exit
                        end if
                    end do
                    if (enter_col /= 0) exit
                end do
            end if
            if (enter_col == 0) then
                no_pivot = .true.
                exit
            end if

            if (r2v(enter_col) /= 0) then
                if (enter_side == 1) then
                    yyv = gx(1:mv+1, enter_col)
                else
                    yyv = -gx(1:mv+1, enter_col)
                end if
                min_kv = huge(1.0d0)
                kk = 0
                do ii = 1, mv+1
                    if (yyv(ii) > tl .and. .not. fvr(ii)) then
                        kval(ii) = bv(ii) / yyv(ii)
                        if (kval(ii) < min_kv - tl) then
                            min_kv = kval(ii)
                            kk = ii
                        else if (abs(kval(ii) - min_kv) < tl .and. bld) then
                            if (kk == 0 .or. IBv(ii) < IBv(kk)) kk = ii
                        end if
                    end if
                end do
                if (kk == 0) then
                    no_pivot = .true.
                    exit
                end if
                if (enter_side == 2) then
                    idx_offset = r1v(enter_col) - p
                    if (idx_offset >= 1 .and. idx_offset <= mv) then
                        yyv(mv+1) = yyv(mv+1) + wv(idx_offset)
                    end if
                end if
            else
                yyv = gx(1:mv+1, enter_col)
                min_kv = huge(1.0d0)
                kk = 0
                if (yyv(mv+1) < 0.0d0) then
                    do ii = 1, mv+1
                        if (yyv(ii) > tl .and. .not. fvr(ii)) then
                            kval(ii) = bv(ii) / yyv(ii)
                            if (kval(ii) < min_kv - tl) then
                                min_kv = kval(ii)
                                kk = ii
                            else if (abs(kval(ii) - min_kv) < tl .and. bld) then
                                if (kk == 0 .or. IBv(ii) < IBv(kk)) kk = ii
                            end if
                        end if
                    end do
                else
                    do ii = 1, mv+1
                        if (yyv(ii) < -tl .and. .not. fvr(ii)) then
                            kval(ii) = -bv(ii) / yyv(ii)
                            if (kval(ii) < min_kv - tl) then
                                min_kv = kval(ii)
                                kk = ii
                            else if (abs(kval(ii) - min_kv) < tl .and. bld) then
                                if (kk == 0 .or. IBv(ii) < IBv(kk)) kk = ii
                            end if
                        end if
                    end do
                end if
                if (kk == 0) then
                    no_pivot = .true.
                    exit
                end if
                fvr(kk) = .true.
            end if

            do ii = 1, mv+1
                if (ii == kk) then
                    eev(ii) = 1.0d0 - 1.0d0 / yyv(kk)
                else
                    eev(ii) = yyv(ii) / yyv(kk)
                end if
            end do

            if (IBv(kk) <= p + mv) then
                gx(1:mv+1, enter_col) = 0.0d0
                gx(kk, enter_col) = 1.0d0
                r1v(enter_col) = IBv(kk)
                r2v(enter_col) = IBv(kk) + mv
            else
                gx(1:mv+1, enter_col) = 0.0d0
                gx(kk, enter_col) = -1.0d0
                idx_offset = IBv(kk) - p - mv
                if (idx_offset >= 1 .and. idx_offset <= mv) then
                    gx(mv+1, enter_col) = wv(idx_offset)
                end if
                r1v(enter_col) = IBv(kk) - mv
                r2v(enter_col) = IBv(kk)
            end if

            do jj = 1, p
                pivot_val = gx(kk, jj)
                do ii = 1, mv+1
                    gx(ii, jj) = gx(ii, jj) - eev(ii) * pivot_val
                end do
            end do
            pivot_val = bv(kk)
            do ii = 1, mv+1
                bv(ii) = bv(ii) - eev(ii) * pivot_val
            end do
            IBv(kk) = enter_var
            iters = iters + 1
        end do
    end subroutine run_simplex_full_tvcqr

    subroutine record_repair()
        implicit none

        if (eva_t >= 1 .and. eva_t <= m) then
            repair_count(eva_t) = repair_count(eva_t) + 1
        end if
    end subroutine record_repair

    subroutine handle_current_bad_signs()
        implicit none
        integer :: bi

        call record_repair()
        if (bad_signs > int(0.1d0 * dble(ms))) then
            mmm_thresh = 2.0d0 * mmm_thresh
            not_new_sl_sh = .true.
        else
            do bi = active_lo, active_hi
                if ((r(bi) <= 0.0d0) .and. sh(bi)) then
                    sh(bi) = .false.
                end if
                if ((r(bi) >= 0.0d0) .and. sl(bi)) then
                    sl(bi) = .false.
                end if
            end do
            not_new_sl_sh = .false.
        end if
    end subroutine handle_current_bad_signs

    subroutine set_failure(code, eval_idx)
        implicit none
        integer, intent(in) :: code, eval_idx

        ierr = code
        if (eval_idx >= 1 .and. eval_idx <= m) then
            failed_eval = eval_idx
            H_seq(eval_idx, :) = 0
        else
            failed_eval = 1
        end if
    end subroutine set_failure

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

    logical function tvcqr_H_basis_valid(H_idx, A_mat, m_loc, p)
        implicit none
        integer, intent(in) :: m_loc, p
        integer, intent(in) :: H_idx(p)
        double precision, intent(in) :: A_mat(m_loc, p)
        integer :: i, j, ierr_local
        double precision :: AH(p, p), AH_inv(p, p)

        tvcqr_H_basis_valid = .true.

        do i = 1, p
            if (H_idx(i) < 1 .or. H_idx(i) > m_loc) then
                tvcqr_H_basis_valid = .false.
                return
            end if
            do j = i + 1, p
                if (H_idx(i) == H_idx(j)) then
                    tvcqr_H_basis_valid = .false.
                    return
                end if
            end do
        end do

        do i = 1, p
            do j = 1, p
                AH(i, j) = A_mat(H_idx(i), j)
            end do
        end do
        call matrix_inverse_2p(AH, AH_inv, p, ierr_local)
        if (ierr_local /= 0) then
            tvcqr_H_basis_valid = .false.
        end if
    end function tvcqr_H_basis_valid

    subroutine same_h_refit_beta_tvcqr(H_idx, A_mat, y_vec, m_loc, p, &
                                       estimate_out, solved)
        implicit none
        integer, intent(in) :: m_loc, p
        integer, intent(in) :: H_idx(p)
        double precision, intent(in) :: A_mat(m_loc, p), y_vec(m_loc)
        double precision, intent(out) :: estimate_out(p)
        logical, intent(out) :: solved
        integer :: si, sj, ierr_local
        double precision :: AH(p, p), AH_inv(p, p), yH(p)

        solved = .false.
        estimate_out = 0.0d0

        ! This repeats the same validity check used inside v37's dense helper;
        ! the caller has already made the same outer check as v37.
        if (.not. tvcqr_H_basis_valid(H_idx, A_mat, m_loc, p)) return

        do si = 1, p
            yH(si) = y_vec(H_idx(si))
            do sj = 1, p
                AH(si, sj) = A_mat(H_idx(si), sj)
            end do
        end do
        call matrix_inverse_2p(AH, AH_inv, p, ierr_local)
        if (ierr_local /= 0) return

        do si = 1, p
            estimate_out(si) = 0.0d0
            do sj = 1, p
                estimate_out(si) = estimate_out(si) + AH_inv(si, sj) * yH(sj)
            end do
        end do
        solved = .true.
    end subroutine same_h_refit_beta_tvcqr

    subroutine same_h_refit_lazy_tvcqr(H_idx, H_previous, A_mat, y_vec, m_loc, p, &
                                       res_tol_loc, estimate_out, r_out, recovered)
        implicit none
        integer, intent(in) :: m_loc, p
        integer, intent(in) :: H_idx(p), H_previous(p)
        double precision, intent(in) :: A_mat(m_loc, p), y_vec(m_loc), res_tol_loc
        double precision, intent(out) :: estimate_out(p)
        double precision, intent(inout) :: r_out(m_loc)
        logical, intent(out) :: recovered
        integer :: si, sj, sk, sh_idx, ierr_local
        logical :: already_materialized
        double precision :: AH(p, p), AH_inv(p, p), yH(p)

        recovered = .false.
        estimate_out = 0.0d0

        if (.not. tvcqr_H_basis_valid(H_idx, A_mat, m_loc, p)) return

        do si = 1, p
            yH(si) = y_vec(H_idx(si))
            do sj = 1, p
                AH(si, sj) = A_mat(H_idx(si), sj)
            end do
        end do
        call matrix_inverse_2p(AH, AH_inv, p, ierr_local)
        if (ierr_local /= 0) return

        do si = 1, p
            estimate_out(si) = 0.0d0
            do sj = 1, p
                estimate_out(si) = estimate_out(si) + AH_inv(si, sj) * yH(sj)
            end do
        end do

        do si = candidate_residual_lo, candidate_residual_hi
            r_out(si) = y_vec(si)
            do sj = 1, p
                r_out(si) = r_out(si) - A_mat(si, sj) * estimate_out(sj)
            end do
        end do

        do si = 1, p
            sh_idx = H_previous(si)
            if (sh_idx >= 1 .and. sh_idx <= m_loc) then
                if (sh_idx < candidate_residual_lo .or. sh_idx > candidate_residual_hi) then
                    r_out(sh_idx) = y_vec(sh_idx)
                    do sj = 1, p
                        r_out(sh_idx) = r_out(sh_idx) - A_mat(sh_idx, sj) * estimate_out(sj)
                    end do
                end if
            end if
        end do
        do si = 1, p
            sh_idx = H_idx(si)
            if (sh_idx >= 1 .and. sh_idx <= m_loc) then
                if (sh_idx < candidate_residual_lo .or. sh_idx > candidate_residual_hi) then
                    already_materialized = .false.
                    do sk = 1, p
                        if (H_previous(sk) == sh_idx) then
                            already_materialized = .true.
                            exit
                        end if
                    end do
                    if (.not. already_materialized) then
                        r_out(sh_idx) = y_vec(sh_idx)
                        do sj = 1, p
                            r_out(sh_idx) = r_out(sh_idx) - A_mat(sh_idx, sj) * estimate_out(sj)
                        end do
                    end if
                end if
            end if
        end do

        recovered = certify_tvcqr_candidate(H_idx, r_out, A_mat, m_loc, p, res_tol_loc)
    end subroutine same_h_refit_lazy_tvcqr

    subroutine same_h_refit_tvcqr(H_idx, A_mat, y_vec, m_loc, p, res_tol_loc, &
                                  estimate_out, r_out, recovered)
        implicit none
        integer, intent(in) :: m_loc, p
        integer, intent(in) :: H_idx(p)
        double precision, intent(in) :: A_mat(m_loc, p), y_vec(m_loc), res_tol_loc
        double precision, intent(out) :: estimate_out(p), r_out(m_loc)
        logical, intent(out) :: recovered
        integer :: i, j, ierr_local
        double precision :: AH(p, p), AH_inv(p, p), yH(p)

        recovered = .false.
        estimate_out = 0.0d0
        r_out = 0.0d0

        if (.not. tvcqr_H_basis_valid(H_idx, A_mat, m_loc, p)) then
            return
        end if

        do i = 1, p
            yH(i) = y_vec(H_idx(i))
            do j = 1, p
                AH(i, j) = A_mat(H_idx(i), j)
            end do
        end do
        call matrix_inverse_2p(AH, AH_inv, p, ierr_local)
        if (ierr_local /= 0) then
            return
        end if

        do i = 1, p
            estimate_out(i) = 0.0d0
            do j = 1, p
                estimate_out(i) = estimate_out(i) + AH_inv(i, j) * yH(j)
            end do
        end do

        do i = 1, m_loc
            r_out(i) = y_vec(i)
            do j = 1, p
                r_out(i) = r_out(i) - A_mat(i, j) * estimate_out(j)
            end do
        end do

        recovered = certify_tvcqr_candidate(H_idx, r_out, A_mat, m_loc, p, res_tol_loc)
    end subroutine same_h_refit_tvcqr

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

            res_min = residual_est(i)
            res_max = residual_est(i)
            res_mean = 0.0d0
            do j = 1, m
                res_min = min(res_min, residual_est((j - 1) * m + i))
                res_max = max(res_max, residual_est((j - 1) * m + i))
                res_mean = res_mean + residual_est((j - 1) * m + i)
            end do
            res_mean = res_mean / dble(m)
            write(unit_num, '(I5,3F12.6,2I8)') i, res_min, res_max, res_mean, &
                it_num(i), final_n_sub(i)
        end do
        close(unit_num)
        
        if (debug_active) write(6, '(A,A)') 'Checkpoint saved with prefix: ', trim(filename_prefix)
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
        
        do j = 1, size(arr, 2)
            do i = 1, size(arr, 1)
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
        integer :: info
        integer :: ipiv(n)
        
        ierr = 0
        
        ! Copy input matrix (LAPACK overwrites)
        A_inv = A_in
        
        ! LU factorization
        call DGETRF(n, n, A_inv, n, ipiv, info)
        if (info /= 0) then
            ierr = info
            return
        end if
        
        ! Compute inverse
        call DGETRI(n, A_inv, n, ipiv, inverse_work, inverse_lwork, info)
        if (info /= 0) then
            ierr = info
            return
        end if
        
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
            do while (j >= left)
                if (arr(j) <= value) exit
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
