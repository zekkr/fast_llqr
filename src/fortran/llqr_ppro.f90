! Local Linear Quantile Regression - PPRO Algorithm (M Threshold Warm Start)
! Fortran implementation of llqr_tau_seq_ppro.R
!
! Compile with: gfortran -shared -fPIC -O3 -march=native -funroll-loops -ffast-math -o llqr_ppro.so llqr_ppro.f90

! Helper function: 2x2 matrix inverse
subroutine inv22(mat, inv_mat, success)
    implicit none
    double precision, intent(in) :: mat(2, 2)
    double precision, intent(out) :: inv_mat(2, 2)
    logical, intent(out) :: success
    double precision :: det

    det = mat(1,1) * mat(2,2) - mat(1,2) * mat(2,1)

    if (abs(det) < 1.0d-15) then
        success = .false.
        inv_mat = 0.0d0
        return
    end if

    ! Match R's inv22: [d, -b; -c, a] / det
    ! For mat = [[a,b],[c,d]], inv = [[d,-b],[-c,a]] / det
    inv_mat(1,1) = mat(2,2) / det    ! d
    inv_mat(1,2) = -mat(1,2) / det   ! -b
    inv_mat(2,1) = -mat(2,1) / det   ! -c
    inv_mat(2,2) = mat(1,1) / det    ! a
    success = .true.
end subroutine inv22

! Helper function: compute mean of absolute values - O(n) instead of median's O(n log n)
! Helper function: median of absolute values
! IMPORTANT: Must match R's median(abs(r)) exactly
function median_abs(arr, n) result(median_val)
    ! NOTE: Despite the name, this computes MEAN(abs(arr)) for performance
    ! Sorting for true median is too expensive for large n (O(n log n) vs O(n))
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: arr(n)
    double precision :: median_val
    integer :: i
    double precision :: sum_abs

    ! Compute mean of absolute values (fast O(n) operation)
    sum_abs = 0.0d0
    do i = 1, n
        sum_abs = sum_abs + abs(arr(i))
    end do
    median_val = sum_abs / dble(n)
end function median_abs

! Helper function: max of array
function max_array(arr, n) result(maxval)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: arr(n)
    double precision :: maxval
    integer :: i

    maxval = arr(1)
    do i = 2, n
        if (arr(i) > maxval) maxval = arr(i)
    end do
end function max_array

! Main PPRO subroutine
subroutine llqr_ppro_fortran(x, y, z, m, nvar, rounds, tau, h, tol, maxit, &
                             Mm_factor, bland_int, ll_est, d_ll_est, it_num, &
                             residual_est, H_mat)

    implicit none

    ! Input arguments
    integer, intent(in) :: m, nvar, rounds, maxit, bland_int
    double precision, intent(in) :: x(m), y(m), z(rounds), tau, tol, Mm_factor
    double precision, intent(inout) :: h

    ! Output arguments
    double precision, intent(out) :: ll_est(rounds)
    double precision, intent(out) :: d_ll_est(rounds)
    integer, intent(out) :: it_num(rounds)
    double precision, intent(out) :: residual_est(rounds, m)
    integer, intent(out) :: H_mat(rounds, nvar+1)

    ! Local variables for full problem (round 1)
    double precision :: A(m, nvar+1)     ! Design matrix [1, x]
    double precision :: w(m)              ! Kernel weights
    double precision :: eva_z(m)          ! z - x for kernel
    double precision :: gammax(m+1, nvar+1)
    double precision :: b(m+1)
    integer :: IB(m+1)
    logical :: freevarrow(m+1)
    integer :: r1(nvar+1), r2(nvar+1)
    double precision :: rr(2, nvar+1)

    ! Working variables
    double precision :: yy(m+1), ee(m+1), k_vals(m+1)
    double precision :: u(m), v(m), estimate(nvar+1)
    integer :: i, j, k, rd, iter, t_rr, tsep, t, jj
    integer :: idx_r1_minus_offset
    double precision :: rrl, min_k, pi, pivot_row_value
    logical :: bland

    ! PPRO-specific variables
    double precision :: mm, mmm, M_threshold, residual_scale
    double precision :: x_norms(m)
    logical :: sl(m), sh(m), not_jl_or_jh(m)
    integer :: ms  ! subsample size
    integer :: H_prev(nvar+1)
    integer :: n_bad_signs
    logical :: not_optimal, not_new_sl_sh  ! Bad signs loop control

    ! Subsample arrays (will be allocated dynamically in concept, but use max size)
    double precision :: gammaxs(m+2, nvar+1)  ! max possible size
    double precision :: bs(m+2)
    double precision :: ws(m+1)
    integer :: IBs(m+2)
    logical :: freevarrows(m+3)
    integer :: r1s(nvar+1), r2s(nvar+1)

    ! Temporary storage
    double precision :: gammaxs_temp(m+2, nvar+1)
    double precision :: bs_temp(m+2)
    integer :: idx_not_jl_or_jh(m)
    integer :: n_subsample, n_potential_S
    integer :: min_subsample_size
    double precision :: r(m)

    ! Z-sorting variables (CRITICAL FIX: match R's z-sorting behavior)
    double precision :: z_sorted(rounds)
    integer :: z_order(rounds)
    double precision :: ll_est_sorted(rounds), d_ll_est_sorted(rounds)
    integer :: it_num_sorted(rounds)
    double precision :: residual_est_sorted(rounds, m)
    integer :: H_mat_sorted(rounds, nvar+1)

    ! Warm start variables (for rd==2)
    double precision :: xh(nvar+1, nvar+1)
    double precision :: xhinv(nvar+1, nvar+1)
    logical :: inv_success
    integer :: H_subsample(nvar+1)
    integer :: idpos(m), idneg(m)
    integer :: n_idpos, n_idneg, n_hbar
    integer :: Hbar(m)
    double precision :: P(m)
    integer :: u_in_IBs(m), v_in_IBs(m)
    double precision :: Pxhbar(m, nvar+1)
    double precision :: Pxhbarxhinv(m, nvar+1)
    double precision :: bs_hbar(m)
    double precision :: lambda(m)
    double precision :: obj_row(nvar+1)
    integer :: r_idx
    logical :: is_in_H
    double precision :: u_subsample(m), v_subsample(m)
    ! Temporary arrays for simplex with correct dimensions
    double precision :: gammaxs_simplex(m, nvar+1)
    double precision :: bs_simplex(m)

    ! rd>2 incremental update variables (matching R's approach)
    double precision :: xhinv_stored(nvar+1, nvar+1)
    double precision :: bs_stored(nvar+1)
    double precision :: gammaxs_pos(m, nvar+1)  ! Stored positive Hbar rows
    double precision :: bs_pos(m)
    integer :: idx_Hbar_pos(m), n_pos_prev
    double precision :: gammaxs_neg(m, nvar+1)  ! Stored negative Hbar rows
    double precision :: bs_neg(m)
    integer :: idx_Hbar_neg(m), n_neg_prev
    integer :: idx_Hbar_pos2(m), idx_Hbar_neg2(m)
    logical :: matched_rows_pos(m), matched_rows_neg(m)
    integer :: n_matched_pos, n_matched_neg, n_unmatched_pos, n_unmatched_neg
    integer :: matched_idx_pos(m), matched_idx_neg(m)
    integer :: unmatched_idx_pos(m), unmatched_idx_neg(m)
    integer :: row_map_pos(m), row_map_neg(m)
    double precision :: Pxhbarxhinv_new(m, nvar+1)
    double precision :: temp_vec(nvar+1)
    integer :: ms_org, curr_idx, prev_idx, match_pos
    integer :: n_xhinv_hbar, ms_current, n_lambda  ! For rd>2 tableau assembly

    ! Helper function declarations
    double precision :: median_abs, max_array

    ! Constants
    pi = 4.0d0 * atan(1.0d0)
    bland = (bland_int /= 0)

    ! Set default bandwidth
    if (h <= 0.0d0) then
        h = dble(m)**(-0.2d0)
    end if

    ! ============================================================
    ! CRITICAL FIX: Sort z array (matching R behavior!)
    ! R code: original_order <- order(z); z <- z[original_order]
    ! This is ESSENTIAL for PPRO's warm-start to work correctly!
    ! ============================================================
    do i = 1, rounds
        z_sorted(i) = z(i)
        z_order(i) = i
    end do

    ! Simple bubble sort (sufficient for typical evaluation points)
    do i = 1, rounds-1
        do j = i+1, rounds
            if (z_sorted(i) > z_sorted(j)) then
                ! Swap z values
                min_k = z_sorted(i)
                z_sorted(i) = z_sorted(j)
                z_sorted(j) = min_k
                ! Swap order indices
                k = z_order(i)
                z_order(i) = z_order(j)
                z_order(j) = k
            end if
        end do
    end do

    ! ============================================================
    ! Build design matrix A = [1, x] (constant across eval points)
    ! ============================================================
    do i = 1, m
        A(i, 1) = 1.0d0
        A(i, 2) = x(i)
    end do

    ! Compute mm = max(x_norms) for M threshold
    ! NOTE: R computes x.norms from raw x, NOT from design matrix A
    ! R code: x.norms <- apply(x, 1, function(row) sqrt(sum(row^2)))
    ! For univariate x, this is just abs(x)
    do i = 1, m
        x_norms(i) = abs(x(i))
    end do
    mm = max_array(x_norms, m)

    ! ============================================================
    ! ROUND 1: Standard simplex (identical to seq)
    ! ============================================================
    rd = 1

    ! Compute kernel weights for first z (USING SORTED Z!)
    do i = 1, m
        eva_z(i) = z_sorted(rd) - x(i)
        w(i) = exp(-0.5d0 * (eva_z(i)/h)**2) / sqrt(2.0d0 * pi)
    end do

    ! Initialize gammax from design matrix
    do i = 1, m
        do j = 1, nvar+1
            gammax(i, j) = A(i, j)
        end do
    end do

    ! Flip signs for negative y
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
            IB(i) = i + nvar + 1  ! u_i
        else
            IB(i) = i + nvar + 1 + m  ! v_i
        end if
    end do
    ! Objective row (m+1) is protected by freevarrow, so IB(m+1) doesn't matter
    ! Set to 0 to avoid conflicts with actual variable indices
    IB(m+1) = 0

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

    ! Compute objective row for round 1
    do j = 1, nvar+1
        gammax(m+1, j) = 0.0d0
        do i = 1, m
            if (IB(i) > nvar + 1 .and. IB(i) <= nvar + 1 + m) then
                ! u_i is basic: use tau * w[i]
                idx_r1_minus_offset = IB(i) - nvar - 1
                if (idx_r1_minus_offset >= 1 .and. idx_r1_minus_offset <= m) then
                    gammax(m+1, j) = gammax(m+1, j) - tau * w(idx_r1_minus_offset) * gammax(i, j)
                end if
            else if (IB(i) > nvar + 1 + m) then
                ! v_i is basic: use (1-tau) * w[i]
                idx_r1_minus_offset = IB(i) - nvar - 1 - m
                if (idx_r1_minus_offset >= 1 .and. idx_r1_minus_offset <= m) then
                    gammax(m+1, j) = gammax(m+1, j) - (1.0d0 - tau) * w(idx_r1_minus_offset) * gammax(i, j)
                end if
            end if
        end do
    end do

    ! Run simplex for round 1 (using PPRO simplex which also works for cold start)
    call run_simplex_ppro(gammax, b, IB, freevarrow, r1, r2, rr, w, m+1, m, nvar, &
                          tau, tol, maxit, bland, iter)

    ! Extract solution for round 1
    call extract_solution(gammax, b, IB, m, nvar, estimate, u, v, r1)

    ! Store results in SORTED arrays (will unsort at the end)
    ll_est_sorted(rd) = estimate(1) + estimate(2) * z_sorted(rd)
    d_ll_est_sorted(rd) = estimate(2)
    it_num_sorted(rd) = iter
    residual_est_sorted(rd, :) = u - v
    H_mat_sorted(rd, :) = r1 - 1 - nvar

    ! BUG FIX: Do NOT sort H_mat for rd=1!
    ! R code keeps H in simplex order: H <- r1 - 1 - nvar (no sort)
    ! The order affects xhinv = solve(A[H,]) and ws[H] in subsequent rounds
    ! Sorting was causing rd=2 to use different H indices than R


    ! ============================================================
    ! ROUNDS 2+: PPRO with M threshold
    ! ============================================================
    do rd = 2, rounds
        ! Initialize bad signs loop control
        not_optimal = .true.
        not_new_sl_sh = .true.
        mmm = mm

        ! Compute kernel weights for current z (USING SORTED Z!)
        do i = 1, m
            eva_z(i) = z_sorted(rd) - x(i)
            w(i) = exp(-0.5d0 * (eva_z(i)/h)**2) / sqrt(2.0d0 * pi)
        end do

        ! BAD SIGNS OUTER LOOP: Keep trying until solution is good
        do while (not_optimal)

            ! ========================================================
            ! Step 1: Compute M threshold from previous residuals
            ! ========================================================
            if (not_new_sl_sh) then
                r = residual_est_sorted(rd-1, :)
                residual_scale = median_abs(r, m)
                M_threshold = max(Mm_factor * mmm * log(log(dble(m))), 0.1d0 * residual_scale)

                ! R's minimum subsample size check (lines 528-533 in R code)
                ! R uses constant 20, NOT min(100, m/2)
                min_subsample_size = 20

                ! OPTIMIZATION: Compute max_r once instead of in loop
                min_k = maxval(abs(r))
                n_subsample = count(abs(r) <= M_threshold)
                do while (n_subsample < min_subsample_size .and. M_threshold < min_k)
                    M_threshold = M_threshold * 1.5d0
                    n_subsample = count(abs(r) <= M_threshold)
                end do

                ! Classify observations
                do i = 1, m
                    sl(i) = r(i) < -M_threshold
                    sh(i) = r(i) > M_threshold
                    not_jl_or_jh(i) = .not. (sl(i) .or. sh(i))
                end do


                ! Force H observations into subsample (S)
                H_prev = H_mat_sorted(rd-1, :)
                do i = 1, nvar+1
                    if (H_prev(i) > 0 .and. H_prev(i) <= m) then
                        if (sl(H_prev(i)) .or. sh(H_prev(i))) then
                            sl(H_prev(i)) = .false.
                            sh(H_prev(i)) = .false.
                            not_jl_or_jh(H_prev(i)) = .true.
                        end if
                    end if
                end do
            end if

        ! Count subsample size
        ms = count(not_jl_or_jh)
        n_subsample = ms

        ! Get indices of subsample observations
        j = 0
        do i = 1, m
            if (not_jl_or_jh(i)) then
                j = j + 1
                idx_not_jl_or_jh(j) = i
            end if
        end do


        ! ========================================================
        ! Step 2: Build subsample problem
        ! ========================================================
        ! NOTE: Subsample observation rows are built ONLY for rd==2 (R line 567)
        ! For rd>2, only aggregates are built; subsample rows are NOT rebuilt!
        ms = n_subsample

        if (rd == 2) then
            ! Copy subsample observations (ONLY for rd==2!)
            do i = 1, ms
                do j = 1, nvar+1
                    gammaxs_temp(i, j) = A(idx_not_jl_or_jh(i), j)
                end do
                bs_temp(i) = y(idx_not_jl_or_jh(i))
            end do
        end if

        ! Copy weights for all rd>=2
        do i = 1, n_subsample
            ws(i) = w(idx_not_jl_or_jh(i))
        end do

        ! Build aggregates for both rd==2 and rd>2 (R lines 572-591)
        if (rd >= 2) then
            ! Aggregate sl observations if any
            ! NOTE: sl is a logical array of size m, so loop over ALL m observations
            ! Store at index m+1 to match R's storage at n+1
            ! R: glob.wx <- colSums(gammaxsl * wsl) where wsl = w[sl]
            ! R: glob.wy <- sum(y[sl] * wsl)
            if (any(sl)) then
                do j = 1, nvar+1
                    gammaxs_temp(m+1, j) = 0.0d0
                    do i = 1, m
                        if (sl(i)) then
                            gammaxs_temp(m+1, j) = gammaxs_temp(m+1, j) + A(i, j) * w(i)
                        end if
                    end do
                end do
                bs_temp(m+1) = 0.0d0
                do i = 1, m
                    if (sl(i)) then
                        bs_temp(m+1) = bs_temp(m+1) + y(i) * w(i)
                    end if
                end do
                ms = ms + 1
                ws(ms) = 1.0d0  ! Weight for sl aggregate is 1
            end if

            ! Aggregate sh observations if any
            ! NOTE: sh is a logical array of size m, so loop over ALL m observations
            ! Store at index m+2 to match R's storage at n+2
            ! R: glob.wx <- colSums(gammaxsh * wsh) where wsh = w[sh]
            ! R: glob.wy <- sum(y[sh] * wsh)
            if (any(sh)) then
                do j = 1, nvar+1
                    gammaxs_temp(m+2, j) = 0.0d0
                    do i = 1, m
                        if (sh(i)) then
                            gammaxs_temp(m+2, j) = gammaxs_temp(m+2, j) + A(i, j) * w(i)
                        end if
                    end do
                end do
                bs_temp(m+2) = 0.0d0
                do i = 1, m
                    if (sh(i)) then
                        bs_temp(m+2) = bs_temp(m+2) + y(i) * w(i)
                    end if
                end do
                ms = ms + 1
                ws(ms) = 1.0d0  ! Weight for sh aggregate is 1
            end if

        end if

        ! ========================================================
        ! Step 3: Inverse-based warm start
        ! ========================================================
        if (rd == 2) then
            ! Map H from previous round to subsample coordinates
            do i = 1, nvar+1
                H_subsample(i) = 0
                do j = 1, n_subsample
                    if (idx_not_jl_or_jh(j) == H_prev(i)) then
                        H_subsample(i) = j
                        exit
                    end if
                end do
                if (H_subsample(i) == 0) then
                    return
                end if
            end do

            ! Extract xh = gammaxs[H,]
            do i = 1, nvar+1
                do j = 1, nvar+1
                    xh(i, j) = gammaxs_temp(H_subsample(i), j)
                end do
            end do

            ! Compute inverse
            call inv22(xh, xhinv, inv_success)

            if (.not. inv_success) then
                return
            end if

            ! Identify positive and negative residuals in subsample (excluding H)
            n_idpos = 0
            n_idneg = 0
            do i = 1, n_subsample
                ! Check if this index is in H
                is_in_H = .false.
                do j = 1, nvar+1
                    if (H_subsample(j) == i) then
                        is_in_H = .true.
                        exit
                    end if
                end do

                if (.not. is_in_H) then
                    ! Get residual from previous round for this subsample observation
                    ! BUG FIX: Use residual_est_sorted (populated after each round) not residual_est (only at end)
                    r_idx = idx_not_jl_or_jh(i)
                    if (residual_est_sorted(rd-1, r_idx) > 0.0d0) then
                        n_idpos = n_idpos + 1
                        idpos(n_idpos) = i
                    else if (residual_est_sorted(rd-1, r_idx) < 0.0d0) then
                        n_idneg = n_idneg + 1
                        idneg(n_idneg) = i
                    else
                    end if
                end if
            end do



            ! Build Hbar (concatenate idpos, idneg, and aggregate indices)
            ! R code: Hbar <- c(idpos,idneg,ms-1,ms) when both sl and sh exist
            n_hbar = n_idpos + n_idneg
            do i = 1, n_idpos
                Hbar(i) = idpos(i)
            end do
            do i = 1, n_idneg
                Hbar(n_idpos + i) = idneg(i)
            end do

            ! Add aggregate row indices if they exist
            ! NOTE: In Fortran, Hbar contains indices into gammaxs_temp for accessing original rows
            ! Aggregates are stored at m+1 (sl) and m+2 (sh) in gammaxs_temp
            ! The Fortran code accesses gammaxs_temp(Hbar(i), :) and stores result in gammaxs(nvar+1+i, :)
            ! So Hbar should contain temp array indices, not final tableau indices
            if (any(sl) .and. any(sh)) then
                ! Both sl and sh exist
                Hbar(n_hbar + 1) = m + 1  ! sl aggregate in temp
                Hbar(n_hbar + 2) = m + 2  ! sh aggregate in temp
                n_hbar = n_hbar + 2
            else if (any(sl)) then
                ! Only sl exists
                Hbar(n_hbar + 1) = m + 1  ! sl aggregate in temp
                n_hbar = n_hbar + 1
            else if (any(sh)) then
                ! Only sh exists
                Hbar(n_hbar + 1) = m + 2  ! sh aggregate in temp
                n_hbar = n_hbar + 1
            end if

            ! Build P vector (1 for idpos, -1 for idneg, sign for aggregates)
            do i = 1, n_idpos
                P(i) = 1.0d0
            end do
            do i = 1, n_idneg
                P(n_idpos + i) = -1.0d0
            end do

            ! Add P values for aggregates
            if (any(sl) .and. any(sh)) then
                P(n_idpos + n_idneg + 1) = -1.0d0  ! sl aggregate gets -1
                P(n_idpos + n_idneg + 2) = 1.0d0   ! sh aggregate gets +1
            else if (any(sl)) then
                P(n_idpos + n_idneg + 1) = 1.0d0   ! Only sl: gets +1? Check R code
            else if (any(sh)) then
                P(n_idpos + n_idneg + 1) = -1.0d0  ! Only sh: gets -1? Check R code
            end if

!DEBUG             write(*,*) '=== FORTRAN Round 2 Warm Start ==='
!DEBUG             write(*,*) 'H_prev:', H_prev
!DEBUG             write(*,*) 'H_subsample:', H_subsample(1:nvar+1)
!DEBUG             write(*,*) 'idpos:', idpos(1:n_idpos)
!DEBUG             write(*,*) 'idneg:', idneg(1:n_idneg)
!DEBUG             write(*,*) 'Hbar:', Hbar(1:n_hbar)
!DEBUG             write(*,*) 'P:', P(1:n_hbar)
!DEBUG             write(*,*) ''
!DEBUG             write(*,*) 'xhinv:'
!DEBUG             do i = 1, nvar+1
!DEBUG                 write(*,'(2F12.6)') xhinv(i, :)
!DEBUG             end do
!DEBUG             write(*,*) ''


            ! Build u_in_IBs and v_in_IBs
            ! R code: u.in.IBs <- idpos + (nvar + 1)
            !         v.in.IBs <- idneg + (nvar + 1) + ms
            ! FIXED: Use ms (with aggregates), not n_subsample!
            do i = 1, n_idpos
                u_in_IBs(i) = idpos(i) + nvar + 1
            end do
            do i = 1, n_idneg
                v_in_IBs(i) = idneg(i) + nvar + 1 + ms
            end do


            ! CRITICAL: Initialize IBs array to avoid garbage values
            ! This prevents segfault from uninitialized memory access
            ! Note: freevarrows is initialized separately below (lines 706-714)
            do i = 1, m+2
                IBs(i) = nvar + 1  ! Safe sentinel value
            end do

            ! Build IBs (initial basis)
            ! R code: IBs <- c(1:(nvar + 1), u.in.IBs, v.in.IBs, nvar+1+2*ms-1, nvar+1+ms)
            ! First nvar+1 entries: 1:(nvar+1) - beta coefficients are basic
            do i = 1, nvar+1
                IBs(i) = i
            end do
            ! Next: u_in_IBs
            do i = 1, n_idpos
                IBs(nvar + 1 + i) = u_in_IBs(i)
            end do
            ! Next: v_in_IBs
            do i = 1, n_idneg
                IBs(nvar + 1 + n_idpos + i) = v_in_IBs(i)
            end do
            ! Next: Add aggregate variable indices if they exist
            ! R code: if both sl and sh exist, add nvar+1+2*ms-1 (v_L for sl) and nvar+1+ms (u_H for sh)
            if (any(sl) .and. any(sh)) then
                IBs(nvar + 1 + n_idpos + n_idneg + 1) = nvar + 1 + 2*ms - 1  ! v_L for sl aggregate
                IBs(nvar + 1 + n_idpos + n_idneg + 2) = nvar + 1 + ms        ! u_H for sh aggregate
            else if (any(sl)) then
                ! Only sl: R code uses nvar+1+2*ms for v_L
                IBs(nvar + 1 + n_idpos + n_idneg + 1) = nvar + 1 + 2*ms
            else if (any(sh)) then
                ! Only sh: R code uses nvar+1+ms for u_H
                IBs(nvar + 1 + n_idpos + n_idneg + 1) = nvar + 1 + ms
            end if
            ! REMOVED: Don't add dummy for objective row - R doesn't do this
            ! IBs should have exactly ms elements, matching R


            ! Initialize freevarrows
            ! R code: freevarrow <- c(rep(TRUE,nvar + 1), rep(FALSE,length(u.in.IBs)),
            !                         rep(FALSE,length(v.in.IBs)), TRUE, TRUE, TRUE)
            ! First nvar+1 are TRUE (beta coefficients are free variables)
            do i = 1, nvar+1
                freevarrows(i) = .true.
            end do
            ! Next n_hbar are FALSE (u/v variables for Hbar)
            do i = nvar+2, nvar+1+n_hbar
                freevarrows(i) = .false.
            end do
            ! BUG FIX: Aggregate variables and objective row must be TRUE!
            ! After Hbar variables, set remaining entries to TRUE
            ! For sl & sh case: v_L, u_H, and objective row are TRUE
            ! For sl only: v_L and objective row are TRUE
            ! For sh only: u_H and objective row are TRUE
            if (any(sl) .and. any(sh)) then
                ! Set v_L, u_H, and objective row to TRUE
                freevarrows(nvar + 1 + n_hbar + 1) = .true.  ! v_L
                freevarrows(nvar + 1 + n_hbar + 2) = .true.  ! u_H
                freevarrows(ms + 1) = .true.  ! objective row
            else if (any(sl)) then
                ! Set v_L and objective row to TRUE
                freevarrows(nvar + 1 + n_hbar + 1) = .true.  ! v_L
                freevarrows(ms + 1) = .true.  ! objective row
            else if (any(sh)) then
                ! Set u_H and objective row to TRUE
                freevarrows(nvar + 1 + n_hbar + 1) = .true.  ! u_H
                freevarrows(ms + 1) = .true.  ! objective row
            else
                ! No aggregates, just objective row
                freevarrows(ms + 1) = .true.  ! objective row
            end if

            ! Initialize r1s, r2s for subsample BEFORE debug output
            ! Non-basic variables are H observations
            ! NOTE: r2 uses ms (subsample size WITH aggregates), not n_subsample
            do i = 1, nvar+1
                r1s(i) = H_subsample(i) + nvar + 1
                r2s(i) = r1s(i) + ms
            end do

            ! Build warm start tableau
            ! First, copy non-H subsample rows from gammaxs_temp to gammaxs
            ! R code: gammaxs <- gammaxs.temp[c(1:(ms-2), n+1, n+2),]
            ! Rows 1:(ms-2) are non-H subsample observations
            ! NOTE: In Fortran, subsample rows are in positions 1:n_subsample of gammaxs_temp,
            ! but aggregates are at m+1 and m+2, not at n_subsample+1 and n_subsample+2!
            do i = 1, n_subsample
                do j = 1, nvar+1
                    gammaxs(i, j) = gammaxs_temp(i, j)
                end do
            end do

            ! Then append aggregate rows from m+1 and m+2 if they exist
            if (any(sl)) then
                do j = 1, nvar+1
                    gammaxs(n_subsample + 1, j) = gammaxs_temp(m+1, j)
                end do
            end if
            if (any(sh)) then
                k = n_subsample + 1
                if (any(sl)) k = k + 1  ! If sl exists, sh goes to next position
                do j = 1, nvar+1
                    gammaxs(k, j) = gammaxs_temp(m+2, j)
                end do
            end if

            ! Then overwrite first nvar+1 rows with xhinv
            do i = 1, nvar+1
                do j = 1, nvar+1
                    gammaxs(i, j) = xhinv(i, j)
                end do
            end do

            ! Compute bs[1:(nvar+1)] = xhinv %*% bs_temp[H]
            do i = 1, nvar+1
                bs(i) = 0.0d0
                do j = 1, nvar+1
                    bs(i) = bs(i) + xhinv(i, j) * bs_temp(H_subsample(j))
                end do
            end do

!DEBUG             write(*,*) 'IBs:', IBs(1:nvar+1+n_hbar+1)
!DEBUG             write(*,*) 'r1s:', r1s
!DEBUG             write(*,*) 'r2s:', r2s
!DEBUG             write(*,*) ''
!DEBUG             write(*,*) 'bs[1:', nvar+1, ']:', bs(1:nvar+1)
!DEBUG             write(*,*) ''


            ! Hbar rows: -P * gammaxs[Hbar,] * xhinv
            if (n_hbar > 0) then
                do i = 1, nvar+1
                end do

                ! Compute Pxhbarxhinv and bs row by row, storing each immediately
                ! IMPORTANT: R computes Pxhbarxhinv WITHOUT leading negative sign!
                ! R code:
                !   Pxhbar <- gammaxs[Hbar, ] * P  (element-wise, NO negative!)
                !   Pxhbarxhinv <- Pxhbar %*% xhinv
                !   gammaxs <- rbind(xhinv, - Pxhbarxhinv)  (negative added when storing to tableau)
                ! For objective row: obj_row <- tau * ws[H] + t(lambda) %*% Pxhbarxhinv (uses Pxhbarxhinv WITHOUT negative)
                do i = 1, n_hbar
                    ! Step 1: Compute P[i] * gammaxs_ORIGINAL[Hbar[i],]
                    ! FIXED: Removed incorrect leading negative sign
                    ! IMPORTANT: Must use gammaxs_temp (original design matrix), NOT gammaxs!
                    ! Because gammaxs[1:2,] has been overwritten with xhinv
                    do j = 1, nvar+1
                        Pxhbar(i, j) = P(i) * gammaxs_temp(Hbar(i), j)
                    end do

                    ! Step 2: Multiply by xhinv to get Pxhbarxhinv[i,]
                    do j = 1, nvar+1
                        Pxhbarxhinv(i, j) = 0.0d0
                        do k = 1, nvar+1
                            Pxhbarxhinv(i, j) = Pxhbarxhinv(i, j) + Pxhbar(i, k) * xhinv(k, j)
                        end do
                    end do

                    ! Step 3: IMMEDIATELY store into gammaxs[nvar+1+i,]
                    ! R stores -Pxhbarxhinv in the tableau: gammaxs <- rbind(xhinv, - Pxhbarxhinv)
                    do j = 1, nvar+1
                        gammaxs(nvar + 1 + i, j) = -Pxhbarxhinv(i, j)
                    end do

                    ! Step 4: Compute bs[nvar+1+i]
                    ! R formula: bs <- c(xhinv %*% bs[H], - Pxhbarxhinv %*% bs[H] + bs[Hbar] * P)
                    ! FIXED: Use Pxhbarxhinv, not Pxhbar!
                    ! bs[nvar+1+i] = -Pxhbarxhinv[i,] @ bs_temp[H] + bs_temp[Hbar[i]] * P[i]
                    bs_hbar(i) = 0.0d0
                    do j = 1, nvar+1
                        bs_hbar(i) = bs_hbar(i) - Pxhbarxhinv(i, j) * bs_temp(H_subsample(j))
                    end do
                    bs_hbar(i) = bs_hbar(i) + bs_temp(Hbar(i)) * P(i)

                    ! Step 5: IMMEDIATELY store bs[nvar+1+i]
                    bs(nvar + 1 + i) = bs_hbar(i)
                end do

                ! Compute objective row
                ! lambda = [tau * ws[idpos], (1-tau) * ws[idneg], (1-tau), tau]
                ! where the last two elements are for sl and sh aggregates
                do i = 1, n_idpos
                    lambda(i) = tau * ws(idpos(i))
                end do
                do i = 1, n_idneg
                    lambda(n_idpos + i) = (1.0d0 - tau) * ws(idneg(i))
                end do
                ! Add aggregate elements to lambda
                ! NOTE: Hbar is [idpos, idneg, m+1, m+2] so aggregates are at positions n_idpos+n_idneg+1 and +2
                if (any(sl)) then
                    lambda(n_idpos + n_idneg + 1) = 1.0d0 - tau  ! sl aggregate: 1-tau
                end if
                if (any(sh)) then
                    k = n_idpos + n_idneg + 1
                    if (any(sl)) k = k + 1
                    lambda(k) = tau  ! sh aggregate: tau
                end if

                ! obj_row = tau * ws[H] + t(lambda) %*% Pxhbarxhinv
                ! R code: obj_row <- tau * ws[H] + crossprod(lambda, Pxhbarxhinv)
                do j = 1, nvar+1
                    ! First term: tau * ws[H[j]] (element-wise, not matrix multiplication!)
                    obj_row(j) = tau * ws(H_subsample(j))
                    ! Second term: sum_i lambda[i] * Pxhbarxhinv[i, j]
                    do i = 1, n_hbar
                        obj_row(j) = obj_row(j) + lambda(i) * Pxhbarxhinv(i, j)
                    end do
                end do

                ! Store objective row at ms+1 (last row of tableau)
                do j = 1, nvar+1
                    gammaxs(ms + 1, j) = obj_row(j)
                end do


!DEBUG                 write(*,*) '-Pxhbarxhinv:'
!DEBUG                 do i = 1, min(n_hbar, 4)
!DEBUG                     write(*,'(2F12.6)') gammaxs(nvar+1+i, :)
!DEBUG                 end do
!DEBUG                 write(*,*) 'bs[3:', nvar+1+n_hbar, ']:', bs(nvar+2:min(nvar+1+n_hbar, nvar+5))
!DEBUG                 write(*,*) ''
!DEBUG                 write(*,*) 'lambda:', lambda(1:min(n_hbar, 8))
!DEBUG                 write(*,*) 'ws[H]:', ws(H_subsample(1:nvar+1))
!DEBUG                 write(*,*) 'Objective row:', obj_row
!DEBUG                 write(*,*) 'bs[objective] =', bs(nvar + 1 + n_hbar + 1)
!DEBUG                 write(*,*) ''
            end if

            ! Set bs for objective row at ms+1 (last row of tableau)
            bs(ms + 1) = 0.0d0


            ! Run PPRO-specific simplex on subsample
            ! NOTE: Pass ms (subsample size WITH aggregates), not n_subsample

            ! FIX: Copy to correctly-sized array to avoid column-major dimension mismatch
            do i = 1, ms+1
                do j = 1, nvar+1
                    gammaxs_simplex(i, j) = gammaxs(i, j)
                end do
                bs_simplex(i) = bs(i)
            end do

            ! Pass m (actual leading dimension), not ms+1
            ! The array gammaxs_simplex has dimension (m, nvar+1), not (ms+1, nvar+1)
            call run_simplex_ppro(gammaxs_simplex, bs_simplex, IBs, freevarrows, r1s, r2s, rr, ws, &
                                  m, ms, nvar, tau, tol, maxit, bland, iter)


            ! Extract solution from subsample
            ! NOTE: Use ms for extraction too
            ! Use bs_simplex which has the correct values after simplex
            call extract_solution(gammaxs_simplex, bs_simplex, IBs, ms, nvar, estimate, &
                                  u_subsample, v_subsample, r1s)


            ! Compute ll_est and d_ll_est (STORE IN SORTED ARRAYS!)
            ll_est_sorted(rd) = estimate(1) + estimate(2) * z_sorted(rd)
            d_ll_est_sorted(rd) = estimate(2)
            it_num_sorted(rd) = iter

            ! Compute FULL residuals for all observations (needed for next point's subsample selection)
            ! R code line 836: r <- y - A %*% estimate
            do i = 1, m
                r(i) = y(i) - (A(i, 1) * estimate(1) + A(i, 2) * estimate(2))
            end do

            ! Map H back to full indices
            ! R code line 857: H <- r1 - 1 - nvar
            ! R code line 858: H <- idx_not_jl_or_jh[H]
            ! No conditional fallback - always map through idx_not_jl_or_jh
            do i = 1, nvar+1
                H_mat_sorted(rd, i) = idx_not_jl_or_jh(r1s(i) - nvar - 1)
            end do

            ! DEBUG: Print H values and solution for rd=2
!DEBUG            if (rd == 2) then
!DEBUG                write(*,*) '=== H extraction at rd=2 (successful iteration) ==='
!DEBUG                write(*,*) 'r1s:', r1s
!DEBUG                write(*,*) 'n_subsample:', n_subsample
!DEBUG                write(*,*) 'H_mat_sorted[2,]:', H_mat_sorted(rd, :)
!DEBUG                write(*,*) 'estimate:', estimate
!DEBUG                write(*,*) 'll_est:', ll_est_sorted(rd)
!DEBUG                write(*,*) 'd_ll_est:', d_ll_est_sorted(rd)
!DEBUG                write(*,*) 'iter:', iter
!DEBUG                write(*,*) 'mmm:', mmm
!DEBUG                write(*,*) 'bs_simplex[1:3]:', bs_simplex(1:3)
!DEBUG                write(*,*) 'gammaxs_simplex[1,]:', gammaxs_simplex(1, :)
!DEBUG                write(*,*) 'gammaxs_simplex[2,]:', gammaxs_simplex(2, :)
!DEBUG            end if

            ! BUG FIX: Do NOT sort H_mat! R code does not sort it.
            ! The order returned by the simplex (r1 indices) must be preserved
            ! for warm-start to work correctly.
            ! if (H_mat_sorted(rd, 1) > H_mat_sorted(rd, 2)) then
            !     j = H_mat_sorted(rd, 1)
            !     H_mat_sorted(rd, 1) = H_mat_sorted(rd, 2)
            !     H_mat_sorted(rd, 2) = j
            ! end if


            ! Set H residuals to 0 (R code line 856: r[H] <- 0)
            do i = 1, nvar+1
                if (H_mat_sorted(rd, i) > 0 .and. H_mat_sorted(rd, i) <= m) then
                    r(H_mat_sorted(rd, i)) = 0.0d0
                end if
            end do

            ! Store full residuals (needed for next point's subsample selection)
            residual_est_sorted(rd, :) = r

            ! Check for bad signs (R's fixup logic)
            if (count(sl) > 0 .or. count(sh) > 0) then
                ! Check bad signs
                n_bad_signs = 0
                do i = 1, m
                    if ((sh(i) .and. r(i) < 0.0d0) .or. (sl(i) .and. r(i) > 0.0d0)) then
                        n_bad_signs = n_bad_signs + 1
                    end if
                end do

                if (n_bad_signs > 0) then
                    ! R logic: if bad.signs > 0.1 * ms, double mmm; otherwise remove bad obs
                    if (dble(n_bad_signs) > 0.1d0 * dble(ms)) then
                        ! Too many bad signs: double M and retry
                        mmm = mmm * 2.0d0
                        not_new_sl_sh = .true.
                        ! Continue while loop - will rebuild with larger M
                    else
                        ! Few bad signs: remove them from sl/sh and retry
                        do i = 1, m
                            if (sh(i) .and. r(i) < 0.0d0) sh(i) = .false.
                            if (sl(i) .and. r(i) > 0.0d0) sl(i) = .false.
                        end do
                        not_new_sl_sh = .false.
                        ! Continue while loop - will rebuild with adjusted sl/sh
                    end if
                else
                    ! No bad signs: success!
                    ! Store xhinv and Hbar information for next round (rd>2 will reuse)
                    ! R code lines 860-876: gammaxs.temp[1:(nvar+1),] <- gammaxs[1:(nvar+1),]
                    do i = 1, nvar+1
                        do j = 1, nvar+1
                            xhinv_stored(i, j) = gammaxs_simplex(i, j)
                        end do
                        bs_stored(i) = bs_simplex(i)
                    end do

                    ! Store Hbar rows separated by sign (for rd>2 incremental update)
                    ! R: ms.org <- ms - any(sl) - any(sh)
                    ms_org = ms
                    if (any(sl)) ms_org = ms_org - 1
                    if (any(sh)) ms_org = ms_org - 1

                    ! DEBUG: Print for rd=2
!DEBUG                    if (rd == 2) then
!DEBUG                        write(*,*) '=== STORAGE at rd=2 (with sl/sh) ==='
!DEBUG                        write(*,*) 'ms:', ms
!DEBUG                        write(*,*) 'ms_org:', ms_org
!DEBUG                        write(*,*) 'any(sl):', any(sl)
!DEBUG                        write(*,*) 'any(sh):', any(sh)
!DEBUG                        write(*,*) 'IBs[nvar+2]:', IBs(nvar+2)
!DEBUG                        write(*,*) 'IBs[nvar+3]:', IBs(nvar+3)
!DEBUG                        write(*,*) 'IBs[nvar+4]:', IBs(nvar+4)
!DEBUG                    end if

                    ! Extract and separate Hbar rows by sign (R lines 864-875)
                    n_pos_prev = 0
                    n_neg_prev = 0
                    do i = nvar+2, ms_org
                        ! Determine if this row is u (positive) or v (negative)
                        ! R: id_gammaxs_Hbar <- IBs[(nvar+2):ms.org]
                        ! R: p_Hbar <- ifelse(id_gammaxs_Hbar > nvar+1+ms, -1, 1)
                        if (IBs(i) > nvar + 1 + ms) then
                            ! This is a v variable (negative residual)
                            n_neg_prev = n_neg_prev + 1
                            ! Get original data index
                            ! R: id_gammaxs_Hbar <- ifelse(..., id - nvar - 1 - ms, id - nvar - 1)
                            curr_idx = IBs(i) - nvar - 1 - ms
                            idx_Hbar_neg(n_neg_prev) = idx_not_jl_or_jh(curr_idx)
                            ! Store the tableau row
                            do j = 1, nvar+1
                                gammaxs_neg(n_neg_prev, j) = gammaxs_simplex(i, j)
                            end do
                            bs_neg(n_neg_prev) = bs_simplex(i)
                        else
                            ! This is a u variable (positive residual)
                            n_pos_prev = n_pos_prev + 1
                            ! Get original data index
                            curr_idx = IBs(i) - nvar - 1
                            idx_Hbar_pos(n_pos_prev) = idx_not_jl_or_jh(curr_idx)

                            ! DEBUG: Print first few stored values for rd=2
!DEBUG                            if (rd == 2 .and. n_pos_prev <= 5) then
!DEBUG                                write(*,'(A,I2,A,I5,A,I5,A,I5,A,I5)') '  Storing pos row ', n_pos_prev, &
!DEBUG                                    ': IBs(', i, ')=', IBs(i), ', curr_idx=', curr_idx, &
!DEBUG                                    ', idx_Hbar_pos=', idx_Hbar_pos(n_pos_prev)
!DEBUG                            end if

                            ! Store the tableau row
                            do j = 1, nvar+1
                                gammaxs_pos(n_pos_prev, j) = gammaxs_simplex(i, j)
                            end do
                            bs_pos(n_pos_prev) = bs_simplex(i)
                        end if
                    end do

                    ! DEBUG: Print summary for rd=2
!DEBUG                    if (rd == 2) then
!DEBUG                        write(*,*) 'Stored n_pos_prev:', n_pos_prev
!DEBUG                        write(*,*) 'Stored n_neg_prev:', n_neg_prev
!DEBUG                        if (n_pos_prev > 0) write(*,*) 'idx_Hbar_pos (first 5):', idx_Hbar_pos(1:min(5, n_pos_prev))
!DEBUG                        if (n_neg_prev > 0) write(*,*) 'idx_Hbar_neg (first 5):', idx_Hbar_neg(1:min(5, n_neg_prev))
!DEBUG                    end if

                    not_optimal = .false.
                end if
            else
                ! No sl or sh: automatically good
                ! Store xhinv and Hbar information for next round
                do i = 1, nvar+1
                    do j = 1, nvar+1
                        xhinv_stored(i, j) = gammaxs_simplex(i, j)
                    end do
                    bs_stored(i) = bs_simplex(i)
                end do

                ! Store Hbar rows separated by sign
                ! R line 863: ms.org <- ms - any(sl) - any(sh)
                ms_org = ms
                if (any(sl)) ms_org = ms_org - 1
                if (any(sh)) ms_org = ms_org - 1
                n_pos_prev = 0
                n_neg_prev = 0

                ! DEBUG: Print for rd=2
!DEBUG                if (rd == 2) then
!DEBUG                    write(*,*) '=== STORAGE at rd=2 ==='
!DEBUG                    write(*,*) 'ms:', ms
!DEBUG                    write(*,*) 'ms_org:', ms_org
!DEBUG                    write(*,*) 'any(sl):', any(sl)
!DEBUG                    write(*,*) 'any(sh):', any(sh)
!DEBUG                    write(*,*) 'IBs[nvar+2]:', IBs(nvar+2)
!DEBUG                    write(*,*) 'IBs[nvar+3]:', IBs(nvar+3)
!DEBUG                    write(*,*) 'IBs[nvar+4]:', IBs(nvar+4)
!DEBUG                end if

                do i = nvar+2, ms_org
                    if (IBs(i) > nvar + 1 + ms) then
                        n_neg_prev = n_neg_prev + 1
                        curr_idx = IBs(i) - nvar - 1 - ms
                        idx_Hbar_neg(n_neg_prev) = idx_not_jl_or_jh(curr_idx)
                        do j = 1, nvar+1
                            gammaxs_neg(n_neg_prev, j) = gammaxs_simplex(i, j)
                        end do
                        bs_neg(n_neg_prev) = bs_simplex(i)
                    else
                        n_pos_prev = n_pos_prev + 1
                        curr_idx = IBs(i) - nvar - 1
                        idx_Hbar_pos(n_pos_prev) = idx_not_jl_or_jh(curr_idx)

                        ! DEBUG: Print first few stored values for rd=2
!DEBUG                        if (rd == 2 .and. n_pos_prev <= 5) then
!DEBUG                            write(*,'(A,I2,A,I5,A,I5,A,I5)') '  Storing pos row ', n_pos_prev, &
!DEBUG                                ': IBs(', i, ')=', IBs(i), ', curr_idx=', curr_idx, &
!DEBUG                                ', idx_Hbar_pos=', idx_Hbar_pos(n_pos_prev)
!DEBUG                        end if

                        do j = 1, nvar+1
                            gammaxs_pos(n_pos_prev, j) = gammaxs_simplex(i, j)
                        end do
                        bs_pos(n_pos_prev) = bs_simplex(i)
                    end if
                end do

                ! DEBUG: Print summary for rd=2
!DEBUG                if (rd == 2) then
!DEBUG                    write(*,*) 'Stored n_pos_prev:', n_pos_prev
!DEBUG                    write(*,*) 'Stored n_neg_prev:', n_neg_prev
!DEBUG                    if (n_pos_prev > 0) write(*,*) 'idx_Hbar_pos (first 5):', idx_Hbar_pos(1:min(5, n_pos_prev))
!DEBUG                    if (n_neg_prev > 0) write(*,*) 'idx_Hbar_neg (first 5):', idx_Hbar_neg(1:min(5, n_neg_prev))
!DEBUG                end if

                not_optimal = .false.
            end if

        else if (rd > 2) then
            ! ========================================================
            ! Step 3 (rd>2): Incremental warm start (R lines 664-732)
            ! ========================================================
            ! Reuse xhinv from previous round (R line 666)
            ! R: xhinv <- gammaxs.temp[1:(nvar + 1),]
            do i = 1, nvar+1
                do j = 1, nvar+1
                    xhinv(i, j) = xhinv_stored(i, j)
                end do
            end do

            ! Map H from previous round to subsample coordinates (same as rd==2)
            do i = 1, nvar+1
                H_subsample(i) = 0
                do j = 1, n_subsample
                    if (idx_not_jl_or_jh(j) == H_prev(i)) then
                        H_subsample(i) = j
                        exit
                    end if
                end do
                if (H_subsample(i) == 0) then
                    return
                end if
            end do

            ! DEBUG: Print for rd=3
!DEBUG            if (rd == 3) then
!DEBUG                write(*,*) '=== FORTRAN DEBUG: rd=3 ==='
!DEBUG                write(*,*) 'H_prev:', H_prev
!DEBUG                write(*,*) 'H_subsample:', H_subsample
!DEBUG                write(*,*) 'n_subsample:', n_subsample
!DEBUG                write(*,*) 'ms (with aggregates):', ms
!DEBUG                write(*,*) 'n_pos_prev:', n_pos_prev
!DEBUG                write(*,*) 'n_neg_prev:', n_neg_prev
!DEBUG            end if

            ! Identify positive and negative residuals in subsample (same as rd==2)
            n_idpos = 0
            n_idneg = 0
            do i = 1, n_subsample
                ! Check if this index is in H
                is_in_H = .false.
                do j = 1, nvar+1
                    if (H_subsample(j) == i) then
                        is_in_H = .true.
                        exit
                    end if
                end do

                if (.not. is_in_H) then
                    r_idx = idx_not_jl_or_jh(i)
                    if (residual_est_sorted(rd-1, r_idx) > 0.0d0) then
                        n_idpos = n_idpos + 1
                        idpos(n_idpos) = i
                    else if (residual_est_sorted(rd-1, r_idx) < 0.0d0) then
                        n_idneg = n_idneg + 1
                        idneg(n_idneg) = i
                    end if
                end if
            end do

            ! Build Hbar (same as rd==2)
            n_hbar = n_idpos + n_idneg
            do i = 1, n_idpos
                Hbar(i) = idpos(i)
            end do
            do i = 1, n_idneg
                Hbar(n_idpos + i) = idneg(i)
            end do
            if (any(sl) .and. any(sh)) then
                Hbar(n_hbar + 1) = m + 1
                Hbar(n_hbar + 2) = m + 2
                n_hbar = n_hbar + 2
            else if (any(sl)) then
                Hbar(n_hbar + 1) = m + 1
                n_hbar = n_hbar + 1
            else if (any(sh)) then
                Hbar(n_hbar + 1) = m + 2
                n_hbar = n_hbar + 1
            end if

            ! Reuse stored basis block exactly (top nvar+1 rows)
            ! R: gammaxs.temp[1:(nvar+1),] <- gammaxs.temp[1:(nvar+1),]
            !    bs.temp[1:(nvar+1)] <- bs.temp[1:(nvar+1)]
            ! NOTE: For rd>2, these are the previous LP's basis rows - DO NOT recompute!
            do i = 1, nvar+1
                do j = 1, nvar+1
                    gammaxs_temp(i, j) = xhinv_stored(i, j)
                end do
                bs_temp(i) = bs_stored(i)
            end do

            ! DEBUG: Print xhinv and bs for rd=3
!DEBUG            if (rd == 3) then
!DEBUG                write(*,*) 'Step 1: xhinv reuse'
!DEBUG                write(*,*) '  xhinv[1,]:', xhinv(1, :)
!DEBUG                write(*,*) '  xhinv[2,]:', xhinv(2, :)
!DEBUG                write(*,*) '  bs_temp[1:2]:', bs_temp(1:2)
!DEBUG                write(*,*) '  y[H]:', y(idx_not_jl_or_jh(H_subsample(1))), y(idx_not_jl_or_jh(H_subsample(2)))
!DEBUG            end if

            ! ========================================================
            ! Incremental update of positive Hbar rows (R lines 668-682)
            ! ========================================================
            ! Get current positive residual indices in original data
            ! R: idx_Hbar_pos2 <- idx_not_jl_or_jh[idpos]
            do i = 1, n_idpos
                idx_Hbar_pos2(i) = idx_not_jl_or_jh(idpos(i))
            end do

            ! Check which ones match previous round
            n_matched_pos = 0
            n_unmatched_pos = 0
            do i = 1, n_idpos
                ! R: matched_rows_pos <- idx_Hbar_pos2 %in% idx_Hbar_pos
                matched_rows_pos(i) = .false.
                do j = 1, n_pos_prev
                    if (idx_Hbar_pos2(i) == idx_Hbar_pos(j)) then
                        matched_rows_pos(i) = .true.
                        row_map_pos(i) = j  ! Store which row in gammaxs_pos to reuse
                        n_matched_pos = n_matched_pos + 1
                        exit
                    end if
                end do
                if (.not. matched_rows_pos(i)) then
                    n_unmatched_pos = n_unmatched_pos + 1
                    unmatched_idx_pos(n_unmatched_pos) = i
                end if
            end do

            ! DEBUG: Print for rd=3
!DEBUG            if (rd == 3) then
!DEBUG                write(*,*) 'Positive Hbar matching:'
!DEBUG                write(*,*) '  n_matched_pos:', n_matched_pos
!DEBUG                write(*,*) '  n_unmatched_pos:', n_unmatched_pos
!DEBUG                if (n_idpos > 0) write(*,*) '  idx_Hbar_pos2 (first 5):', idx_Hbar_pos2(1:min(5, n_idpos))
!DEBUG                if (n_pos_prev > 0) write(*,*) '  idx_Hbar_pos (first 5):', idx_Hbar_pos(1:min(5, n_pos_prev))
!DEBUG            end if

            ! Copy matched rows from previous gammaxs_pos
            ! R lines 673-675: indices <- (nvar+2):(nvar+1+length(u.in.IBs))
            !                  gammaxs.temp[indices[rows.pos2], ] <- gammaxs.pos[rows.pos, ]
            do i = 1, n_idpos
                curr_idx = nvar + 1 + i  ! Position in gammaxs_temp
                if (matched_rows_pos(i)) then
                    ! Reuse stored row
                    do j = 1, nvar+1
                        gammaxs_temp(curr_idx, j) = gammaxs_pos(row_map_pos(i), j)
                    end do
                    bs_temp(curr_idx) = bs_pos(row_map_pos(i))
                else
                    ! Compute new row: -A[idx,:] %*% xhinv (R lines 679-681)
                    ! R: Pxhbarxhinv.pos <- A[rows.pos2.nm,] %*% xhinv
                    ! R: gammaxs.temp[indices[!matched_rows_pos],] <- - Pxhbarxhinv.pos
                    r_idx = idx_Hbar_pos2(i)
                    do j = 1, nvar+1
                        gammaxs_temp(curr_idx, j) = 0.0d0
                        do k = 1, nvar+1
                            gammaxs_temp(curr_idx, j) = gammaxs_temp(curr_idx, j) - A(r_idx, k) * xhinv(k, j)
                        end do
                    end do
                    ! R: bs.temp[indices[!matched_rows_pos]] <- - Pxhbarxhinv.pos %*% y[idx_not_jl_or_jh[H]] + y[rows.pos2.nm]
                    ! Since gammaxs_temp = -Pxhbarxhinv.pos, we want gammaxs_temp %*% y[H] + y[rows]
                    bs_temp(curr_idx) = 0.0d0
                    do j = 1, nvar+1
                        bs_temp(curr_idx) = bs_temp(curr_idx) + gammaxs_temp(curr_idx, j) * y(idx_not_jl_or_jh(H_subsample(j)))
                    end do
                    bs_temp(curr_idx) = bs_temp(curr_idx) + y(r_idx)

                    ! DEBUG: Print first few new positive rows for rd=3
!DEBUG                    if (rd == 3 .and. i <= 3) then
!DEBUG                        write(*,'(A,I2,A,I5)') 'Step 2: New positive Hbar row ', i, ', orig_idx=', r_idx
!DEBUG                        write(*,'(A,2F15.8)') '  A[orig_idx,]:', A(r_idx, :)
!DEBUG                        write(*,'(A,2F15.8)') '  gammaxs_temp row:', gammaxs_temp(curr_idx, :)
!DEBUG                        write(*,'(A,F15.8)') '  bs_temp:', bs_temp(curr_idx)
!DEBUG                    end if
                end if
            end do

            ! ========================================================
            ! Incremental update of negative Hbar rows (R lines 684-698)
            ! ========================================================
            do i = 1, n_idneg
                idx_Hbar_neg2(i) = idx_not_jl_or_jh(idneg(i))
            end do

            n_matched_neg = 0
            n_unmatched_neg = 0
            do i = 1, n_idneg
                matched_rows_neg(i) = .false.
                do j = 1, n_neg_prev
                    if (idx_Hbar_neg2(i) == idx_Hbar_neg(j)) then
                        matched_rows_neg(i) = .true.
                        row_map_neg(i) = j
                        n_matched_neg = n_matched_neg + 1
                        exit
                    end if
                end do
                if (.not. matched_rows_neg(i)) then
                    n_unmatched_neg = n_unmatched_neg + 1
                    unmatched_idx_neg(n_unmatched_neg) = i
                end if
            end do

            ! Copy matched rows or compute new ones
            do i = 1, n_idneg
                curr_idx = nvar + 1 + n_idpos + i  ! Position in gammaxs_temp
                if (matched_rows_neg(i)) then
                    do j = 1, nvar+1
                        gammaxs_temp(curr_idx, j) = gammaxs_neg(row_map_neg(i), j)
                    end do
                    bs_temp(curr_idx) = bs_neg(row_map_neg(i))
                else
                    ! R: Pxhbarxhinv.neg <- - A[rows.neg2.nm,] %*% xhinv (line 695)
                    r_idx = idx_Hbar_neg2(i)
                    do j = 1, nvar+1
                        gammaxs_temp(curr_idx, j) = 0.0d0
                        do k = 1, nvar+1
                            gammaxs_temp(curr_idx, j) = gammaxs_temp(curr_idx, j) + A(r_idx, k) * xhinv(k, j)
                        end do
                    end do
                    ! R line 697: bs.temp[...] <- - Pxhbarxhinv.neg %*% y[...] - y[rows.neg2.nm]
                    ! Pxhbarxhinv.neg = -A*xhinv, gammaxs_temp = A*xhinv
                    ! R wants: -Pxhbarxhinv.neg*y[H] - y[rows] = A*xhinv*y[H] - y[rows] = gammaxs_temp*y[H] - y[rows]
                    bs_temp(curr_idx) = 0.0d0
                    do j = 1, nvar+1
                        bs_temp(curr_idx) = bs_temp(curr_idx) + gammaxs_temp(curr_idx, j) * y(idx_not_jl_or_jh(H_subsample(j)))
                    end do
                    bs_temp(curr_idx) = bs_temp(curr_idx) - y(r_idx)

                    ! DEBUG: Print first few new negative rows for rd=3
!DEBUG                    if (rd == 3 .and. i <= 3) then
!DEBUG                        write(*,'(A,I2,A,I5)') 'Step 3: New negative Hbar row ', i, ', orig_idx=', r_idx
!DEBUG                        write(*,'(A,2F15.8)') '  A[orig_idx,]:', A(r_idx, :)
!DEBUG                        write(*,'(A,2F15.8)') '  -A[orig_idx,] %*% xhinv:', gammaxs_temp(curr_idx, :)
!DEBUG                        write(*,'(A,F15.8)') '  bs_temp:', bs_temp(curr_idx)
!DEBUG                    end if
                end if
            end do

            ! ========================================================
            ! Update aggregate rows (R lines 700-724)
            ! ========================================================
            if (any(sl) .and. any(sh)) then
                ! R lines 701-704: Update both sl and sh aggregates
                ! gammaxs.temp[m+1,] <- gammaxs.temp[m+1,] %*% xhinv
                do j = 1, nvar+1
                    temp_vec(j) = 0.0d0
                    do k = 1, nvar+1
                        temp_vec(j) = temp_vec(j) + gammaxs_temp(m+1, k) * xhinv(k, j)
                    end do
                end do
                do j = 1, nvar+1
                    gammaxs_temp(m+1, j) = temp_vec(j)
                end do
                ! bs.temp[m+1] <- gammaxs.temp[m+1,] %*% y[idx_not_jl_or_jh[H]] - bs.temp[m+1]
                min_k = 0.0d0
                do j = 1, nvar+1
                    min_k = min_k + gammaxs_temp(m+1, j) * y(idx_not_jl_or_jh(H_subsample(j)))
                end do
                bs_temp(m+1) = min_k - bs_temp(m+1)

                ! gammaxs.temp[m+2,] <- - gammaxs.temp[m+2,] %*% xhinv
                do j = 1, nvar+1
                    temp_vec(j) = 0.0d0
                    do k = 1, nvar+1
                        temp_vec(j) = temp_vec(j) - gammaxs_temp(m+2, k) * xhinv(k, j)
                    end do
                end do
                do j = 1, nvar+1
                    gammaxs_temp(m+2, j) = temp_vec(j)
                end do
                ! bs.temp[m+2] <- gammaxs.temp[m+2,] %*% y[...] + bs.temp[m+2]
                min_k = 0.0d0
                do j = 1, nvar+1
                    min_k = min_k + gammaxs_temp(m+2, j) * y(idx_not_jl_or_jh(H_subsample(j)))
                end do
                bs_temp(m+2) = min_k + bs_temp(m+2)

                ! DEBUG: Print aggregate transformations for rd=3
!DEBUG                if (rd == 3) then
!DEBUG                    write(*,*) 'Step 4: Aggregate transformations (both sl and sh)'
!DEBUG                    write(*,'(A,2F15.8)') '  gammaxs_temp[m+1,]:', gammaxs_temp(m+1, :)
!DEBUG                    write(*,'(A,F15.8)') '  bs_temp[m+1]:', bs_temp(m+1)
!DEBUG                    write(*,'(A,2F15.8)') '  gammaxs_temp[m+2,]:', gammaxs_temp(m+2, :)
!DEBUG                    write(*,'(A,F15.8)') '  bs_temp[m+2]:', bs_temp(m+2)
!DEBUG                end if
            else if (any(sl)) then
                ! Only sl aggregate (R lines 709-710)
                do j = 1, nvar+1
                    temp_vec(j) = 0.0d0
                    do k = 1, nvar+1
                        temp_vec(j) = temp_vec(j) + gammaxs_temp(m+1, k) * xhinv(k, j)
                    end do
                end do
                do j = 1, nvar+1
                    gammaxs_temp(m+1, j) = temp_vec(j)
                end do
                min_k = 0.0d0
                do j = 1, nvar+1
                    min_k = min_k + gammaxs_temp(m+1, j) * y(idx_not_jl_or_jh(H_subsample(j)))
                end do
                bs_temp(m+1) = min_k - bs_temp(m+1)

                ! DEBUG: Print aggregate transformations for rd=3
!DEBUG                if (rd == 3) then
!DEBUG                    write(*,*) 'Step 4: Aggregate transformations (sl only)'
!DEBUG                    write(*,'(A,2F15.8)') '  gammaxs_temp[m+1,]:', gammaxs_temp(m+1, :)
!DEBUG                    write(*,'(A,F15.8)') '  bs_temp[m+1]:', bs_temp(m+1)
!DEBUG                end if
            else if (any(sh)) then
                ! Only sh aggregate (R lines 715-716)
                do j = 1, nvar+1
                    temp_vec(j) = 0.0d0
                    do k = 1, nvar+1
                        temp_vec(j) = temp_vec(j) - gammaxs_temp(m+2, k) * xhinv(k, j)
                    end do
                end do
                do j = 1, nvar+1
                    gammaxs_temp(m+2, j) = temp_vec(j)
                end do
                min_k = 0.0d0
                do j = 1, nvar+1
                    min_k = min_k + gammaxs_temp(m+2, j) * y(idx_not_jl_or_jh(H_subsample(j)))
                end do
                bs_temp(m+2) = min_k + bs_temp(m+2)

                ! DEBUG: Print aggregate transformations for rd=3
!DEBUG                if (rd == 3) then
!DEBUG                    write(*,*) 'Step 4: Aggregate transformations (sh only)'
!DEBUG                    write(*,'(A,2F15.8)') '  gammaxs_temp[m+2,]:', gammaxs_temp(m+2, :)
!DEBUG                    write(*,'(A,F15.8)') '  bs_temp[m+2]:', bs_temp(m+2)
!DEBUG                end if
            end if

            ! ========================================================
            ! Build final gammaxs and bs from gammaxs_temp (R lines 705-724)
            ! ========================================================
            ! For rd>2, gammaxs_temp contains:
            !   Rows 1:(nvar+1): xhinv
            !   Rows (nvar+2):(nvar+1+n_idpos): positive Hbar rows
            !   Rows (nvar+2+n_idpos):(nvar+1+n_idpos+n_idneg): negative Hbar rows
            !   Rows m+1, m+2: aggregates (if they exist)
            ! R: ms <- nvar + 1 + length(u.in.IBs) + length(v.in.IBs) + (1 if sl) + (1 if sh)
            ! R: gammaxs <- gammaxs.temp[c(1:(ms-2),m+1,m+2),]
            n_xhinv_hbar = nvar + 1 + n_idpos + n_idneg  ! This is ms-2 when both aggregates exist

            ! DEBUG: Print for rd=3
!DEBUG            if (rd == 3) then
!DEBUG                write(*,*) 'Final tableau assembly:'
!DEBUG                write(*,*) '  n_xhinv_hbar:', n_xhinv_hbar
!DEBUG                write(*,*) '  any(sl):', any(sl)
!DEBUG                write(*,*) '  any(sh):', any(sh)
!DEBUG                write(*,*) '  gammaxs_temp[1,]:', gammaxs_temp(1, 1), gammaxs_temp(1, 2)
!DEBUG                write(*,*) '  gammaxs_temp[2,]:', gammaxs_temp(2, 1), gammaxs_temp(2, 2)
!DEBUG                if (n_idpos > 0) then
!DEBUG                    write(*,*) '  gammaxs_temp[3,] (first pos Hbar):', gammaxs_temp(3, 1), gammaxs_temp(3, 2)
!DEBUG                end if
!DEBUG                if (any(sl)) write(*,*) '  gammaxs_temp[m+1,]:', gammaxs_temp(m+1, 1), gammaxs_temp(m+1, 2)
!DEBUG                if (any(sh)) write(*,*) '  gammaxs_temp[m+2,]:', gammaxs_temp(m+2, 1), gammaxs_temp(m+2, 2)
!DEBUG            end if

            if (any(sl) .and. any(sh)) then
                ! Copy xhinv and Hbar rows 1:(ms-2), then aggregates m+1, m+2
                do i = 1, n_xhinv_hbar
                    do j = 1, nvar+1
                        gammaxs(i, j) = gammaxs_temp(i, j)
                    end do
                    bs(i) = bs_temp(i)
                end do
                do j = 1, nvar+1
                    gammaxs(n_xhinv_hbar + 1, j) = gammaxs_temp(m+1, j)
                    gammaxs(n_xhinv_hbar + 2, j) = gammaxs_temp(m+2, j)
                end do
                bs(n_xhinv_hbar + 1) = bs_temp(m+1)
                bs(n_xhinv_hbar + 2) = bs_temp(m+2)
                ! Build lambda (R line 707)
                do i = 1, n_idpos
                    lambda(i) = tau * ws(idpos(i))
                end do
                do i = 1, n_idneg
                    lambda(n_idpos + i) = (1.0d0 - tau) * ws(idneg(i))
                end do
                lambda(n_idpos + n_idneg + 1) = 1.0d0 - tau
                lambda(n_idpos + n_idneg + 2) = tau
            else if (any(sl)) then
                ! Copy xhinv and Hbar rows 1:(ms-1), then aggregate m+1
                do i = 1, n_xhinv_hbar
                    do j = 1, nvar+1
                        gammaxs(i, j) = gammaxs_temp(i, j)
                    end do
                    bs(i) = bs_temp(i)
                end do
                do j = 1, nvar+1
                    gammaxs(n_xhinv_hbar + 1, j) = gammaxs_temp(m+1, j)
                end do
                bs(n_xhinv_hbar + 1) = bs_temp(m+1)
                do i = 1, n_idpos
                    lambda(i) = tau * ws(idpos(i))
                end do
                do i = 1, n_idneg
                    lambda(n_idpos + i) = (1.0d0 - tau) * ws(idneg(i))
                end do
                lambda(n_idpos + n_idneg + 1) = 1.0d0 - tau
            else if (any(sh)) then
                ! Copy xhinv and Hbar rows 1:(ms-1), then aggregate m+2
                do i = 1, n_xhinv_hbar
                    do j = 1, nvar+1
                        gammaxs(i, j) = gammaxs_temp(i, j)
                    end do
                    bs(i) = bs_temp(i)
                end do
                do j = 1, nvar+1
                    gammaxs(n_xhinv_hbar + 1, j) = gammaxs_temp(m+2, j)
                end do
                bs(n_xhinv_hbar + 1) = bs_temp(m+2)
                do i = 1, n_idpos
                    lambda(i) = tau * ws(idpos(i))
                end do
                do i = 1, n_idneg
                    lambda(n_idpos + i) = (1.0d0 - tau) * ws(idneg(i))
                end do
                lambda(n_idpos + n_idneg + 1) = tau
            else
                ! No aggregates - copy all xhinv and Hbar rows (1:ms)
                do i = 1, n_xhinv_hbar
                    do j = 1, nvar+1
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

            ! Compute Pxhbarxhinv and objective row (R lines 726-731)
            ! R: Pxhbarxhinv <- - gammaxs[(nvar+2):ms,]
            ! This extracts the Hbar rows (which are already stored in gammaxs after xhinv rows)
            ! R: gammaxs <- rbind(gammaxs, tau * ws[H] + t(lambda) %*% Pxhbarxhinv)
            ! Compute current ms (number of rows in gammaxs after assembly)
            if (any(sl) .and. any(sh)) then
                ms_current = n_xhinv_hbar + 2
                n_lambda = n_idpos + n_idneg + 2
            else if (any(sl) .or. any(sh)) then
                ms_current = n_xhinv_hbar + 1
                n_lambda = n_idpos + n_idneg + 1
            else
                ms_current = n_xhinv_hbar
                n_lambda = n_idpos + n_idneg
            end if

            ! Compute objective row
            do j = 1, nvar+1
                obj_row(j) = tau * ws(H_subsample(j))
                ! Pxhbarxhinv are rows (nvar+2):ms_current
                do i = 1, n_lambda
                    ! Pxhbarxhinv is -gammaxs[nvar+1+i,]
                    obj_row(j) = obj_row(j) - lambda(i) * gammaxs(nvar + 1 + i, j)
                end do
            end do
            do j = 1, nvar+1
                gammaxs(ms_current + 1, j) = obj_row(j)
            end do
            bs(ms_current + 1) = 0.0d0

            ! DEBUG: Print objective row for rd=3
!DEBUG            if (rd == 3) then
!DEBUG                write(*,*) 'Step 5: Objective row computation'
!DEBUG                write(*,*) '  ms_current:', ms_current
!DEBUG                write(*,*) '  n_lambda:', n_lambda
!DEBUG                write(*,'(A,2F15.8)') '  ws[H]:', ws(H_subsample(1)), ws(H_subsample(2))
!DEBUG                write(*,'(A,2F15.8)') '  lambda[1:2]:', lambda(1), lambda(2)
!DEBUG                write(*,'(A,2F15.8)') '  obj_row:', obj_row(1:2)
!DEBUG                write(*,'(A,F15.8)') '  bs[ms+1]:', bs(ms_current + 1)
!DEBUG            end if

            ! Initialize IBs, freevarrows, r1s, r2s (same as rd==2)
            do i = 1, nvar+1
                IBs(i) = i
            end do
            do i = 1, n_idpos
                u_in_IBs(i) = idpos(i) + nvar + 1
                IBs(nvar + 1 + i) = u_in_IBs(i)
            end do
            do i = 1, n_idneg
                v_in_IBs(i) = idneg(i) + nvar + 1 + ms
                IBs(nvar + 1 + n_idpos + i) = v_in_IBs(i)
            end do
            if (any(sl) .and. any(sh)) then
                IBs(nvar + 1 + n_idpos + n_idneg + 1) = nvar + 1 + 2*ms - 1
                IBs(nvar + 1 + n_idpos + n_idneg + 2) = nvar + 1 + ms
            else if (any(sl)) then
                IBs(nvar + 1 + n_idpos + n_idneg + 1) = nvar + 1 + 2*ms
            else if (any(sh)) then
                IBs(nvar + 1 + n_idpos + n_idneg + 1) = nvar + 1 + ms
            end if

            do i = 1, nvar+1
                freevarrows(i) = .true.
            end do
            do i = nvar+2, nvar+1+n_hbar
                freevarrows(i) = .false.
            end do
            if (any(sl) .and. any(sh)) then
                freevarrows(nvar + 1 + n_hbar + 1) = .true.
                freevarrows(nvar + 1 + n_hbar + 2) = .true.
                freevarrows(ms + 1) = .true.
            else if (any(sl)) then
                freevarrows(nvar + 1 + n_hbar + 1) = .true.
                freevarrows(ms + 1) = .true.
            else if (any(sh)) then
                freevarrows(nvar + 1 + n_hbar + 1) = .true.
                freevarrows(ms + 1) = .true.
            else
                freevarrows(ms + 1) = .true.
            end if

            do i = 1, nvar+1
                r1s(i) = H_subsample(i) + nvar + 1
                r2s(i) = r1s(i) + ms
            end do

            ! Copy to simplex arrays (same as rd==2)
            do i = 1, ms+1
                do j = 1, nvar+1
                    gammaxs_simplex(i, j) = gammaxs(i, j)
                end do
                bs_simplex(i) = bs(i)
            end do

            ! Run simplex (same as rd==2)
            call run_simplex_ppro(gammaxs_simplex, bs_simplex, IBs, freevarrows, r1s, r2s, rr, ws, &
                                  m, ms, nvar, tau, tol, maxit, bland, iter)

            ! Extract solution (same as rd==2)
            call extract_solution(gammaxs_simplex, bs_simplex, IBs, ms, nvar, estimate, &
                                  u_subsample, v_subsample, r1s)

            ! Store results (same as rd==2)
            ll_est_sorted(rd) = estimate(1) + estimate(2) * z_sorted(rd)
            d_ll_est_sorted(rd) = estimate(2)
            it_num_sorted(rd) = iter

            ! Compute full residuals (same as rd==2)
            do i = 1, m
                r(i) = y(i) - (A(i, 1) * estimate(1) + A(i, 2) * estimate(2))
            end do

            ! Map H back to full indices (exactly as rd==2 and R)
            ! R code: H <- r1 - 1 - nvar; H <- idx_not_jl_or_jh[H]
            ! NO fallback to H_prev - always use r1s unconditionally
            do i = 1, nvar+1
                H_mat_sorted(rd, i) = idx_not_jl_or_jh(r1s(i) - nvar - 1)
            end do

            ! Set H residuals to 0 (same as rd==2)
            do i = 1, nvar+1
                if (H_mat_sorted(rd, i) > 0 .and. H_mat_sorted(rd, i) <= m) then
                    r(H_mat_sorted(rd, i)) = 0.0d0
                end if
            end do

            ! Store full residuals
            residual_est_sorted(rd, :) = r

            ! Check for bad signs (same as rd==2)
            if (count(sl) > 0 .or. count(sh) > 0) then
                n_bad_signs = 0
                do i = 1, m
                    if ((sh(i) .and. r(i) < 0.0d0) .or. (sl(i) .and. r(i) > 0.0d0)) then
                        n_bad_signs = n_bad_signs + 1
                    end if
                end do

                if (n_bad_signs > 0) then
                    if (dble(n_bad_signs) > 0.1d0 * dble(ms)) then
                        mmm = mmm * 2.0d0
                        not_new_sl_sh = .true.
                    else
                        do i = 1, m
                            if (sh(i) .and. r(i) < 0.0d0) sh(i) = .false.
                            if (sl(i) .and. r(i) > 0.0d0) sl(i) = .false.
                        end do
                        not_new_sl_sh = .false.
                    end if
                else
                    ! No bad signs: success! Store for next round
                    do i = 1, nvar+1
                        do j = 1, nvar+1
                            xhinv_stored(i, j) = gammaxs_simplex(i, j)
                        end do
                        bs_stored(i) = bs_simplex(i)
                    end do

                    ms_org = ms
                    if (any(sl)) ms_org = ms_org - 1
                    if (any(sh)) ms_org = ms_org - 1

                    n_pos_prev = 0
                    n_neg_prev = 0
                    do i = nvar+2, ms_org
                        if (IBs(i) > nvar + 1 + ms) then
                            n_neg_prev = n_neg_prev + 1
                            curr_idx = IBs(i) - nvar - 1 - ms
                            idx_Hbar_neg(n_neg_prev) = idx_not_jl_or_jh(curr_idx)
                            do j = 1, nvar+1
                                gammaxs_neg(n_neg_prev, j) = gammaxs_simplex(i, j)
                            end do
                            bs_neg(n_neg_prev) = bs_simplex(i)
                        else
                            n_pos_prev = n_pos_prev + 1
                            curr_idx = IBs(i) - nvar - 1
                            idx_Hbar_pos(n_pos_prev) = idx_not_jl_or_jh(curr_idx)
                            do j = 1, nvar+1
                                gammaxs_pos(n_pos_prev, j) = gammaxs_simplex(i, j)
                            end do
                            bs_pos(n_pos_prev) = bs_simplex(i)
                        end if
                    end do

                    not_optimal = .false.
                end if
            else
                ! No sl or sh: automatically good, store for next round
                do i = 1, nvar+1
                    do j = 1, nvar+1
                        xhinv_stored(i, j) = gammaxs_simplex(i, j)
                    end do
                    bs_stored(i) = bs_simplex(i)
                end do

                ms_org = ms
                n_pos_prev = 0
                n_neg_prev = 0
                do i = nvar+2, ms_org
                    if (IBs(i) > nvar + 1 + ms) then
                        n_neg_prev = n_neg_prev + 1
                        curr_idx = IBs(i) - nvar - 1 - ms
                        idx_Hbar_neg(n_neg_prev) = idx_not_jl_or_jh(curr_idx)
                        do j = 1, nvar+1
                            gammaxs_neg(n_neg_prev, j) = gammaxs_simplex(i, j)
                        end do
                        bs_neg(n_neg_prev) = bs_simplex(i)
                    else
                        n_pos_prev = n_pos_prev + 1
                        curr_idx = IBs(i) - nvar - 1
                        idx_Hbar_pos(n_pos_prev) = idx_not_jl_or_jh(curr_idx)
                        do j = 1, nvar+1
                            gammaxs_pos(n_pos_prev, j) = gammaxs_simplex(i, j)
                        end do
                        bs_pos(n_pos_prev) = bs_simplex(i)
                    end if
                end do

                not_optimal = .false.
            end if

        end if  ! rd == 2 or rd > 2

        end do  ! while (not_optimal)

    end do  ! rd = 2, rounds

    ! ============================================================
    ! Unsort results back to original order (matching R track_order=TRUE)
    ! R code: if (track_order) ll_est <- ll_est[order(original_order)]
    ! ============================================================
    do rd = 1, rounds
        k = z_order(rd)  ! Original position of this sorted element
        ll_est(k) = ll_est_sorted(rd)
        d_ll_est(k) = d_ll_est_sorted(rd)
        it_num(k) = it_num_sorted(rd)
        do i = 1, m
            residual_est(k, i) = residual_est_sorted(rd, i)
        end do
        do i = 1, nvar+1
            H_mat(k, i) = H_mat_sorted(rd, i)
        end do
    end do

contains

    ! PPRO-specific simplex algorithm (matches R's llqr_tau_seq_ppro lines 731-831)
    subroutine run_simplex_ppro(gx, bv, IBv, fvr, r1v, r2v, rrv, wv, ldgx, mv, nvr, &
                                tv, tl, mxit, bld, iters)
        implicit none
        integer, intent(in) :: ldgx, mv, nvr, mxit
        double precision, intent(inout) :: gx(ldgx, nvr+1), bv(mv+1)
        integer, intent(inout) :: IBv(mv+1), r1v(nvr+1), r2v(nvr+1)
        logical, intent(inout) :: fvr(mv+1)
        double precision, intent(inout) :: rrv(2, nvr+1)
        double precision, intent(in) :: wv(mv), tv, tl
        logical, intent(in) :: bld
        integer, intent(out) :: iters

        ! Local variables
        double precision :: yyv(mv+1), eev(mv+1), k_valsv(mv+1)
        integer :: ii, jj, kk, t_rrv, tsepv, tv_val
        integer :: idx_offset
        double precision :: rrlv, min_kv, pivot_val

        iters = 0


        do while (iters < mxit)
            ! Step 2: Compute reduced costs (PPRO version)
            ! R code line 735-736:
            ! rr[1, ] <- gammaxs[ms+1, ]
            ! rr[2, ] <- (ws[r1 - 1 - nvar] - rr[1, ])


            do ii = 1, nvr+1
                rrv(1, ii) = gx(mv+1, ii)

                ! PPRO formula: ws[r1 - 1 - nvar] - rr[1, ]
                idx_offset = r1v(ii) - 1 - nvr
                if (idx_offset >= 1 .and. idx_offset <= mv) then
                    rrv(2, ii) = wv(idx_offset) - rrv(1, ii)
                else
                    rrv(2, ii) = -rrv(1, ii)  ! If out of range, just negate
                end if

                ! R code: rr[1, r2==0] <- -abs(rr[1, r2==0])
                ! When r2[i] == 0, negate rr[1,i] to make it negative
                if (r2v(ii) == 0) then
                    rrv(1, ii) = -abs(rrv(1, ii))
                end if
            end do


            ! Check optimality
            rrlv = minval(rrv)
            if (rrlv >= -tl) then
                exit
            end if

            ! Step 3: Choose entering variable
            if (bld) then
                ! Bland's rule implementation
                if (any(rrv(1,:) < -tl)) then
                    do ii = 1, nvr+1
                        if (rrv(1,ii) < -tl) then
                            tv_val = r1v(ii)
                            t_rrv = ii
                            tsepv = 1
                            exit
                        end if
                    end do
                else
                    do ii = 1, nvr+1
                        if (rrv(2,ii) < -tl) then
                            tv_val = r2v(ii)
                            t_rrv = ii
                            tsepv = 2
                            exit
                        end if
                    end do
                end if
            else
                ! Standard rule: most negative
                do ii = 1, 2
                    do jj = 1, nvr+1
                        if (abs(rrv(ii,jj) - rrlv) < tl) then
                            t_rrv = jj
                            tsepv = ii
                            if (tsepv == 1) then
                                tv_val = r1v(t_rrv)
                            else
                                tv_val = r2v(t_rrv)
                            end if
                            goto 100
                        end if
                    end do
                end do
100             continue
            end if


            ! Step 4 & 5: Choose leaving variable (PPRO version)
            ! R code line 776-780
            ! IMPORTANT: Only copy rows 1:mv+1, not all ldgx rows!
            if (tsepv == 1) then
                do ii = 1, mv+1
                    yyv(ii) = gx(ii, t_rrv)
                end do
            else
                do ii = 1, mv+1
                    yyv(ii) = -gx(ii, t_rrv)
                end do
            end if


            ! Ratio test (R code line 783-800)
            ! k <- bs / yy
            ! k <- which((min(k[yy > 0 & !freevarrow]) == k) & (yy > 0))
            min_kv = huge(1.0d0)
            kk = 0

            do ii = 1, mv+1
                if (yyv(ii) > tl .and. .not. fvr(ii)) then
                    k_valsv(ii) = bv(ii) / yyv(ii)
                    if (k_valsv(ii) < min_kv - tl) then
                        min_kv = k_valsv(ii)
                        kk = ii
                    else if (abs(k_valsv(ii) - min_kv) < tl .and. bld) then
                        ! Bland's rule tie-breaking: choose smallest index
                        if (IBv(ii) < IBv(kk)) then
                            kk = ii
                        end if
                    end if
                end if
            end do

            if (kk == 0) then
                exit
            end if


            ! Step 6': Adjust yy if entering variable is v_i (R code line 805-807)
            ! if (tsep != 1){ yy[ms + 1] <- yy[ms + 1] + ws[r1[t_rr] - 1 - nvar] }
            if (tsepv /= 1) then
                idx_offset = r1v(t_rrv) - 1 - nvr
                if (idx_offset >= 1 .and. idx_offset <= mv) then
                    yyv(mv+1) = yyv(mv+1) + wv(idx_offset)
                end if
            end if

            ! Step 6: Pivoting (R code line 809-828)
            ! ee <- yy / yy[k]; ee[k] <- 1 - 1 / yy[k]
            do ii = 1, mv+1
                if (ii == kk) then
                    eev(ii) = 1.0d0 - 1.0d0 / yyv(kk)
                else
                    eev(ii) = yyv(ii) / yyv(kk)
                end if
            end do

            ! Update pivot column and r1, r2 (R code line 812-823)
            if (IBv(kk) <= (mv + nvr + 1)) then
                ! R: if (IBs[k] <= (ms + nvar + 1))
                ! IMPORTANT: Only zero out rows 1:mv+1, not all ldgx rows!
                do ii = 1, mv+1
                    gx(ii, t_rrv) = 0.0d0
                end do
                gx(kk, t_rrv) = 1.0d0
                r1v(t_rrv) = IBv(kk)
                r2v(t_rrv) = IBv(kk) + mv
            else
                ! R: else { gammaxs[k, t_rr] <- -1; gammaxs[(ms + 1), t_rr] <- ws[IBs[k] - ms - nvar - 1] }
                ! IMPORTANT: Only zero out rows 1:mv+1, not all ldgx rows!
                do ii = 1, mv+1
                    gx(ii, t_rrv) = 0.0d0
                end do
                gx(kk, t_rrv) = -1.0d0
                idx_offset = IBv(kk) - mv - nvr - 1
                if (idx_offset >= 1 .and. idx_offset <= mv) then
                    gx(mv+1, t_rrv) = wv(idx_offset)
                end if
                r1v(t_rrv) = IBv(kk) - mv
                r2v(t_rrv) = IBv(kk)
            end if

            ! Update all columns (tcrossprod)
            do jj = 1, nvr+1
                pivot_val = gx(kk, jj)
                do ii = 1, mv+1
                    gx(ii, jj) = gx(ii, jj) - eev(ii) * pivot_val
                end do
            end do

            ! Update b
            pivot_val = bv(kk)
            do ii = 1, mv+1
                bv(ii) = bv(ii) - eev(ii) * pivot_val
            end do

            IBv(kk) = tv_val

            ! R code: if (t <= nvar+1) { freevarrow[k] <- TRUE }
            ! Mark row as free variable if a coefficient entered
            if (tv_val <= nvr + 1) then
                fvr(kk) = .true.
            end if

            iters = iters + 1

        end do

    end subroutine run_simplex_ppro

    ! Extract solution from basis
    subroutine extract_solution(gx, bv, IBv, mv, nvr, est, uv, vv, r1v)
        implicit none
        integer, intent(in) :: mv, nvr
        double precision, intent(in) :: gx(mv+1, nvr+1), bv(mv+1)
        integer, intent(in) :: IBv(mv+1), r1v(nvr+1)
        double precision, intent(out) :: est(nvr+1), uv(mv), vv(mv)
        integer :: ii, jj, u_idx, v_idx
        logical :: u_tmp(mv), v_tmp(mv)

        ! Initialize
        est = 0.0d0
        uv = 0.0d0
        vv = 0.0d0

        ! Extract from basis (only loop through m, not m+1)
        do ii = 1, mv
            if (IBv(ii) > nvr + 1 .and. IBv(ii) <= nvr + 1 + mv) then
                ! u variable is basic
                u_idx = IBv(ii) - nvr - 1
                if (u_idx >= 1 .and. u_idx <= mv) then
                    uv(u_idx) = bv(ii)
                end if
            else if (IBv(ii) > nvr + 1 + mv) then
                ! v variable is basic
                v_idx = IBv(ii) - nvr - 1 - mv
                if (v_idx >= 1 .and. v_idx <= mv) then
                    vv(v_idx) = bv(ii)
                end if
            else if (IBv(ii) >= 1 .and. IBv(ii) <= nvr + 1) then
                ! Coefficient variable is basic
                est(IBv(ii)) = bv(ii)
            end if
        end do
    end subroutine extract_solution

end subroutine llqr_ppro_fortran
