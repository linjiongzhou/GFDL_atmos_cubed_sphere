!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the FV3 dynamical core.
!*
!* The FV3 dynamical core is free software: you can redistribute it
!* and/or modify it under the terms of the
!* GNU Lesser General Public License as published by the
!* Free Software Foundation, either version 3 of the License, or
!* (at your option) any later version.
!*
!* The FV3 dynamical core is distributed in the hope that it will be
!* useful, but WITHOUT ANYWARRANTY; without even the implied warranty
!* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!* See the GNU General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with the FV3 dynamical core.
!* If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

! =======================================================================
! Fast Physics Interface
! Developer: Linjiong Zhou
! Initial Development: 3/5/2021
! =======================================================================

module fast_phys_mod

    use constants_mod, only: rdgas, grav
    use fv_grid_utils_mod, only: cubed_to_latlon, update_dwinds_phys
    use fv_arrays_mod, only: fv_grid_type, fv_grid_bounds_type, inline_mp_type
    use mpp_domains_mod, only: domain2d, mpp_update_domains
    use fv_timing_mod, only: timing_on, timing_off
    use tracer_manager_mod, only: get_tracer_index
    use field_manager_mod, only: model_atmos
    use gfdl_mp_mod, only: gfdl_mp_driver, fast_sat_adj, iqs, c_liq, c_ice, cv_air, cv_vap
    use mpp_mod, only: stdout, mpp_chksum
    use module_mp_fast_sbm, only: fast_sbm
    
    implicit none
    
    private

    real, parameter :: consv_min = 0.001

    public :: fast_phys

contains

subroutine fast_phys (is, ie, js, je, isd, ied, jsd, jed, km, npx, npy, &
               c2l_ord, mdt, consv, akap, pfull, hs, te0_2d, ua, va, u, &
               v, w, pt, delp, delz, q_con, cappa, q, pkz, te, peln, pe, pk, ps, &
               inline_mp, gridstruct, domain, bd, hydrostatic, do_adiabatic_init, &
               do_inline_mp, do_sat_adj, last_step, omga, r_vir, &
               do_fsbm, a_step, fsbm_bin, fsbm_dx, fsbm_dy, pt_old, q_old, warm_start)
    
    implicit none
    
    ! -----------------------------------------------------------------------
    ! input / output arguments
    ! -----------------------------------------------------------------------

    integer, intent (in) :: is, ie, js, je, isd, ied, jsd, jed, km, npx, npy, c2l_ord
    integer, intent (in) :: a_step, fsbm_bin

    logical, intent (in) :: hydrostatic, do_adiabatic_init, do_inline_mp, do_sat_adj, last_step
    logical, intent (in) :: do_fsbm, warm_start

    real, intent (in) :: consv, mdt, akap, fsbm_dx, fsbm_dy, r_vir

    real, intent (in), dimension (km) :: pfull

    real, intent (in), dimension (isd:ied, jsd:jed) :: hs

    real, intent (inout), dimension (isd:ied, jsd:jed) :: ps

    real, intent (inout), dimension (is:ie, js:je) :: te0_2d

    real, intent (inout), dimension (is:ie, js:je, km+1) :: pk

    real, intent (inout), dimension (is:ie, km+1, js:je) :: peln

    real, intent (inout), dimension (is:, js:, 1:) :: delz
    
    real, intent (inout), dimension (isd:, jsd:, 1:) :: q_con, cappa, w
    
    real, intent (inout), dimension (isd:ied, jsd:jed, km) :: pt, ua, va, delp, pt_old, q_old, omga

    real, intent (inout), dimension (isd:ied, jsd:jed, km, *) :: q

    real, intent (inout), dimension (isd:ied, jsd:jed+1, km) :: u

    real, intent (inout), dimension (isd:ied+1, jsd:jed, km) :: v

    real, intent (inout), dimension (is-1:ie+1, km+1, js-1:je+1) :: pe

    real, intent (out), dimension (is:ie, js:je, km) :: pkz

    real, intent (out), dimension (isd:ied, jsd:jed, km) :: te

    type (fv_grid_type), intent (in), target :: gridstruct

    type (fv_grid_bounds_type), intent (in) :: bd

    type (domain2d), intent (inout) :: domain

    type (inline_mp_type), intent (inout) :: inline_mp

    ! -----------------------------------------------------------------------
    ! local variables
    ! -----------------------------------------------------------------------

    integer :: i, j, k, kmp, n_chem, num_sbmradar, itimestep, unit, n
    integer :: sphum, liq_wat, ice_wat, rainwat, snowwat, graupel, cld_amt, ccn_cm3, cin_cm3
    integer :: ql_num, qr_num, qi_num, qs_num, qg_num, qa_num, qn_num

    integer, dimension (fsbm_bin) :: qlr_ind, qis_ind, qg_ind, qa_ind, qn_ind

    logical :: diagflag = .false.

    real :: rrg, dqv, dql, dqr, dqi, dqs, dqg, ps_dt, cvm
    real :: qliq, qsol, f_sum, mu, sigma, alpha, beta, qsat, rh

    real, parameter :: xr_a = 0.25 ! p value in Xu and Randall (1996)
    real, parameter :: xr_b = 100. ! alpha_0 value in Xu and Randall (1996)
    real, parameter :: xr_c = 0.49 ! gamma value in Xu and Randall (1996)

    real, dimension (fsbm_bin) :: f

    real, dimension (is:ie) :: gsize

    real, dimension (is:ie, km) :: q2, q3

    real, dimension (is:ie, km+1) :: phis

    real, allocatable, dimension (:,:) :: dz, wa

    real, allocatable, dimension (:,:,:) :: u_dt, v_dt, dp0, u0, v0
    
    real, allocatable, dimension (:,:) :: xland, rainnc, rainncv, snownc, snowncv, graupelnc, graupelncv
    real, allocatable, dimension (:,:,:) :: ur, vr, wr, dz8w, p_phy, pi_phy, rho_phy, th_phy, sbqna, sbqnn
    real, allocatable, dimension (:,:,:) :: sbqv, sbqc, sbqr, sbqi, sbqs, sbqg, sbqnc, sbqnr, sbqni, sbqns, sbqng
    real, allocatable, dimension (:,:,:) :: ma, lh_rate, ce_rate, ds_rate, melt_rate, frz_rate, th_old, qv_old
    real, allocatable, dimension (:,:,:) :: cldnucl_rate, icenucl_rate, n_reg_ccn, pkz0, delz0, dlnp
    real, allocatable, dimension (:,:,:,:) :: chem_new
    real, allocatable, dimension (:,:,:,:) :: sbmradar

    character (len = 4) :: ind

    sphum = get_tracer_index (model_atmos, 'sphum')
    liq_wat = get_tracer_index (model_atmos, 'liq_wat')
    ice_wat = get_tracer_index (model_atmos, 'ice_wat')
    rainwat = get_tracer_index (model_atmos, 'rainwat')
    snowwat = get_tracer_index (model_atmos, 'snowwat')
    graupel = get_tracer_index (model_atmos, 'graupel')
    cld_amt = get_tracer_index (model_atmos, 'cld_amt')
    ccn_cm3 = get_tracer_index (model_atmos, 'ccn_cm3')
    cin_cm3 = get_tracer_index (model_atmos, 'cin_cm3')
    ql_num = get_tracer_index (MODEL_ATMOS, 'ql_num')
    qr_num = get_tracer_index (MODEL_ATMOS, 'qr_num')
    qi_num = get_tracer_index (MODEL_ATMOS, 'qi_num')
    qs_num = get_tracer_index (MODEL_ATMOS, 'qs_num')
    qg_num = get_tracer_index (MODEL_ATMOS, 'qg_num')
    qa_num = get_tracer_index (MODEL_ATMOS, 'qa_num')
    qn_num = get_tracer_index (MODEL_ATMOS, 'qn_num')

    rrg = - rdgas / grav

    ! time saving trick
    if (last_step) then
        kmp = 1
    else
        do k = 1, km
            kmp = k
            if (pfull (k) .gt. 50.E2) exit
        enddo
    endif

    !-----------------------------------------------------------------------
    ! Fast Saturation Adjustment >>>
    !-----------------------------------------------------------------------

    ! Note: pt at this stage is T_v
    if (do_adiabatic_init .or. do_sat_adj) then

        call timing_on ('fast_sat_adj')

        allocate (dz (is:ie, kmp:km))

!$OMP parallel do default (none) shared (is, ie, js, je, isd, jsd, kmp, km, te, &
!$OMP                                    delp, hydrostatic, hs, pt, peln, delz, rainwat, &
!$OMP                                    liq_wat, ice_wat, snowwat, graupel, q_con, &
!$OMP                                    sphum, pkz, last_step, consv, te0_2d, gridstruct, &
!$OMP                                    q, mdt, cld_amt, cappa, rrg, akap, ccn_cm3, &
!$OMP                                    cin_cm3, inline_mp) &
!$OMP                           private (q2, q3, gsize, dz)

        do j = js, je

            gsize (is:ie) = sqrt (gridstruct%area_64 (is:ie, j))

            if (ccn_cm3 .gt. 0) then
                q2 (is:ie, kmp:km) = q (is:ie, j, kmp:km, ccn_cm3)
            else
                q2 (is:ie, kmp:km) = 0.0
            endif
            if (cin_cm3 .gt. 0) then
                q3 (is:ie, kmp:km) = q (is:ie, j, kmp:km, cin_cm3)
            else
                q3 (is:ie, kmp:km) = 0.0
            endif
 
            if (.not. hydrostatic) then
                dz (is:ie, kmp:km) = delz (is:ie, j, kmp:km)
            else
                dz (is:ie, kmp:km) = (peln (is:je, kmp:km, j) - peln (is:ie, kmp+1:km+1, j)) * &
                    rdgas * pt (is:ie, j, kmp:km) / grav
            endif

            call fast_sat_adj (abs (mdt), is, ie, kmp, km, hydrostatic, consv .gt. consv_min, &
                     te (is:ie, j, kmp:km), q (is:ie, j, kmp:km, sphum), q (is:ie, j, kmp:km, liq_wat), &
                     q (is:ie, j, kmp:km, rainwat), q (is:ie, j, kmp:km, ice_wat), &
                     q (is:ie, j, kmp:km, snowwat), q (is:ie, j, kmp:km, graupel), &
                     q (is:ie, j, kmp:km, cld_amt), q2 (is:ie, kmp:km), q3 (is:ie, kmp:km), &
                     hs (is:ie, j), dz (is:ie, kmp:km), pt (is:ie, j, kmp:km), &
                     delp (is:ie, j, kmp:km), &
#ifdef USE_COND
                     q_con (is:ie, j, kmp:km), &
#else
                     q_con (isd:, jsd, 1:), &
#endif
#ifdef MOIST_CAPPA
                     cappa (is:ie, j, kmp:km), &
#else
                     cappa (isd:, jsd, 1:), &
#endif
                     gsize, last_step, inline_mp%cond (is:ie, j), inline_mp%reevap (is:ie, j), &
                     inline_mp%dep (is:ie, j), inline_mp%sub (is:ie, j))

            ! update pkz
            if (.not. hydrostatic) then
#ifdef MOIST_CAPPA
                pkz (is:ie, j, kmp:km) = exp (cappa (is:ie, j, kmp:km) * &
                    log (rrg * delp (is:ie, j, kmp:km) / &
                    delz (is:ie, j, kmp:km) * pt (is:ie, j, kmp:km)))
#else
                pkz (is:ie, j, kmp:km) = exp (akap * log (rrg * delp (is:ie, j, kmp:km) / &
                    delz (is:ie, j, kmp:km) * pt (is:ie, j, kmp:km)))
#endif
            endif
 
            if (consv .gt. consv_min) then
                do i = is, ie
                    do k = kmp, km
                        te0_2d (i, j) = te0_2d (i, j) + te (i, j, k)
                    enddo
                enddo
            endif

        enddo

        deallocate (dz)

        call timing_off ('fast_sat_adj')

    endif

    !-----------------------------------------------------------------------
    ! <<< Fast Saturation Adjustment
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! Inline GFDL MP >>>
    !-----------------------------------------------------------------------

    if ((.not. do_adiabatic_init) .and. do_inline_mp) then

        call timing_on ('gfdl_mp')

        allocate (u_dt (isd:ied, jsd:jed, km))
        allocate (v_dt (isd:ied, jsd:jed, km))

        do k = 1, km
            do j = jsd, jed
                do i = isd, ied
                    u_dt (i, j, k) = 0.
                    v_dt (i, j, k) = 0.
                enddo
            enddo
        enddo

        ! save D grid u and v
        if (consv .gt. consv_min) then
            allocate (u0 (isd:ied, jsd:jed+1, km))
            allocate (v0 (isd:ied+1, jsd:jed, km))
            u0 = u
            v0 = v
        endif

        ! D grid wind to A grid wind remap
        call cubed_to_latlon (u, v, ua, va, gridstruct, npx, npy, km, 1, gridstruct%grid_type, &
                 domain, gridstruct%bounded_domain, c2l_ord, bd)

        ! save delp
        if (consv .gt. consv_min) then
            allocate (dp0 (isd:ied, jsd:jed, km))
            dp0 = delp
        endif

        allocate (dz (is:ie, kmp:km))
        allocate (wa (is:ie, kmp:km))

!$OMP parallel do default (none) shared (is, ie, js, je, isd, jsd, kmp, km, pe, ua, va, &
!$OMP                                    te, delp, hydrostatic, hs, pt, peln, delz, &
!$OMP                                    rainwat, liq_wat, ice_wat, snowwat, graupel, q_con, &
!$OMP                                    sphum, w, pk, pkz, last_step, consv, te0_2d, &
!$OMP                                    gridstruct, q, mdt, cld_amt, cappa, rrg, akap, &
!$OMP                                    ccn_cm3, cin_cm3, inline_mp, do_inline_mp, ps) &
!$OMP                           private (u_dt, v_dt, q2, q3, gsize, dz, wa)

        do j = js, je

            gsize (is:ie) = sqrt (gridstruct%area_64 (is:ie, j))

            if (ccn_cm3 .gt. 0) then
                q2 (is:ie, kmp:km) = q (is:ie, j, kmp:km, ccn_cm3)
            else
                q2 (is:ie, kmp:km) = 0.0
            endif
            if (cin_cm3 .gt. 0) then
                q3 (is:ie, kmp:km) = q (is:ie, j, kmp:km, cin_cm3)
            else
                q3 (is:ie, kmp:km) = 0.0
            endif
 
            ! note: ua and va are A-grid variables
            ! note: pt is virtual temperature at this point
            ! note: w is vertical velocity (m/s)
            ! note: delz is negative, delp is positive, delz doesn't change in constant volume situation
            ! note: hs is geopotential height (m^2/s^2)
            ! note: the unit of q2 or q3 is #/cm^3
            ! note: the unit of area is m^2
            ! note: the unit of prew, prer, prei, pres, preg is mm/day
            ! note: the unit of cond, dep, reevap, sub is mm/day

            ! save ua, va for wind tendency calculation
            u_dt (is:ie, j, kmp:km) = ua (is:ie, j, kmp:km)
            v_dt (is:ie, j, kmp:km) = va (is:ie, j, kmp:km)

            if (allocated (inline_mp%liq_wat_dt)) inline_mp%liq_wat_dt (is:ie, j, kmp:km) = &
                inline_mp%liq_wat_dt (is:ie, j, kmp:km) - q (is:ie, j, kmp:km, liq_wat)
            if (allocated (inline_mp%ice_wat_dt)) inline_mp%ice_wat_dt (is:ie, j, kmp:km) = &
                inline_mp%ice_wat_dt (is:ie, j, kmp:km) - q (is:ie, j, kmp:km, ice_wat)
            if (allocated (inline_mp%qv_dt)) inline_mp%qv_dt (is:ie, j, kmp:km) = &
                inline_mp%qv_dt (is:ie, j, kmp:km) - q (is:ie, j, kmp:km, sphum)
            if (allocated (inline_mp%ql_dt)) inline_mp%ql_dt (is:ie, j, kmp:km) = &
                inline_mp%ql_dt (is:ie, j, kmp:km) - (q (is:ie, j, kmp:km, liq_wat) + &
                q (is:ie, j, kmp:km, rainwat))
            if (allocated (inline_mp%qi_dt)) inline_mp%qi_dt (is:ie, j, kmp:km) = &
                inline_mp%qi_dt (is:ie, j, kmp:km) - (q (is:ie, j, kmp:km, ice_wat) + &
                q (is:ie, j, kmp:km, snowwat) + q (is:ie, j, kmp:km, graupel))
            if (allocated (inline_mp%qr_dt)) inline_mp%qr_dt (is:ie, j, kmp:km) = &
                inline_mp%qr_dt (is:ie, j, kmp:km) - q (is:ie, j, kmp:km, rainwat)
            if (allocated (inline_mp%qs_dt)) inline_mp%qs_dt (is:ie, j, kmp:km) = &
                inline_mp%qs_dt (is:ie, j, kmp:km) - q (is:ie, j, kmp:km, snowwat)
            if (allocated (inline_mp%qg_dt)) inline_mp%qg_dt (is:ie, j, kmp:km) = &
                inline_mp%qg_dt (is:ie, j, kmp:km) - q (is:ie, j, kmp:km, graupel)
            if (allocated (inline_mp%t_dt)) inline_mp%t_dt (is:ie, j, kmp:km) = &
                inline_mp%t_dt (is:ie, j, kmp:km) - pt (is:ie, j, kmp:km)
            if (allocated (inline_mp%u_dt)) inline_mp%u_dt (is:ie, j, kmp:km) = &
                inline_mp%u_dt (is:ie, j, kmp:km) - ua (is:ie, j, kmp:km)
            if (allocated (inline_mp%v_dt)) inline_mp%v_dt (is:ie, j, kmp:km) = &
                inline_mp%v_dt (is:ie, j, kmp:km) - va (is:ie, j, kmp:km)

            if (.not. hydrostatic) then
                wa (is:ie, kmp:km) = w (is:ie, j, kmp:km)
                dz (is:ie, kmp:km) = delz (is:ie, j, kmp:km)
            else
                dz (is:ie, kmp:km) = (peln (is:je, kmp:km, j) - peln (is:ie, kmp+1:km+1, j)) * &
                    rdgas * pt (is:ie, j, kmp:km) / grav
            endif

            call gfdl_mp_driver (q (is:ie, j, kmp:km, sphum), q (is:ie, j, kmp:km, liq_wat), &
                     q (is:ie, j, kmp:km, rainwat), q (is:ie, j, kmp:km, ice_wat), &
                     q (is:ie, j, kmp:km, snowwat), q (is:ie, j, kmp:km, graupel), &
                     q (is:ie, j, kmp:km, cld_amt), q2 (is:ie, kmp:km), &
                     q3 (is:ie, kmp:km), pt (is:ie, j, kmp:km), wa (is:ie, kmp:km), &
                     ua (is:ie, j, kmp:km), va (is:ie, j, kmp:km), dz (is:ie, kmp:km), &
                     delp (is:ie, j, kmp:km), gsize, abs (mdt), hs (is:ie, j), &
                     inline_mp%prew (is:ie, j), inline_mp%prer (is:ie, j), &
                     inline_mp%prei (is:ie, j), inline_mp%pres (is:ie, j), &
                     inline_mp%preg (is:ie, j), hydrostatic, is, ie, kmp, km, &
#ifdef USE_COND
                     q_con (is:ie, j, kmp:km), &
#else
                     q_con (isd:, jsd, 1:), &
#endif
#ifdef MOIST_CAPPA
                     cappa (is:ie, j, kmp:km), &
#else
                     cappa (isd:, jsd, 1:), &
#endif
                     consv .gt. consv_min, te (is:ie, j, kmp:km), inline_mp%cond (is:ie, j), &
                     inline_mp%dep (is:ie, j), inline_mp%reevap (is:ie, j), inline_mp%sub (is:ie, j), &
                     last_step, do_inline_mp)

            if (.not. hydrostatic) then
                w (is:ie, j, kmp:km) = wa (is:ie, kmp:km)
            endif

            ! compute wind tendency at A grid fori D grid wind update
            u_dt (is:ie, j, kmp:km) = (ua (is:ie, j, kmp:km) - u_dt (is:ie, j, kmp:km)) / abs (mdt)
            v_dt (is:ie, j, kmp:km) = (va (is:ie, j, kmp:km) - v_dt (is:ie, j, kmp:km)) / abs (mdt)

            if (allocated (inline_mp%liq_wat_dt)) inline_mp%liq_wat_dt (is:ie, j, kmp:km) = &
                inline_mp%liq_wat_dt (is:ie, j, kmp:km) + q (is:ie, j, kmp:km, liq_wat)
            if (allocated (inline_mp%ice_wat_dt)) inline_mp%ice_wat_dt (is:ie, j, kmp:km) = &
                inline_mp%ice_wat_dt (is:ie, j, kmp:km) + q (is:ie, j, kmp:km, ice_wat)
            if (allocated (inline_mp%qv_dt)) inline_mp%qv_dt (is:ie, j, kmp:km) = &
                inline_mp%qv_dt (is:ie, j, kmp:km) + q (is:ie, j, kmp:km, sphum)
            if (allocated (inline_mp%ql_dt)) inline_mp%ql_dt (is:ie, j, kmp:km) = &
                inline_mp%ql_dt (is:ie, j, kmp:km) + (q (is:ie, j, kmp:km, liq_wat) + &
                q (is:ie, j, kmp:km, rainwat))
            if (allocated (inline_mp%qi_dt)) inline_mp%qi_dt (is:ie, j, kmp:km) = &
                inline_mp%qi_dt (is:ie, j, kmp:km) + (q (is:ie, j, kmp:km, ice_wat) + &
                q (is:ie, j, kmp:km, snowwat) + q (is:ie, j, kmp:km, graupel))
            if (allocated (inline_mp%qr_dt)) inline_mp%qr_dt (is:ie, j, kmp:km) = &
                inline_mp%qr_dt (is:ie, j, kmp:km) + q (is:ie, j, kmp:km, rainwat)
            if (allocated (inline_mp%qs_dt)) inline_mp%qs_dt (is:ie, j, kmp:km) = &
                inline_mp%qs_dt (is:ie, j, kmp:km) + q (is:ie, j, kmp:km, snowwat)
            if (allocated (inline_mp%qg_dt)) inline_mp%qg_dt (is:ie, j, kmp:km) = &
                inline_mp%qg_dt (is:ie, j, kmp:km) + q (is:ie, j, kmp:km, graupel)
            if (allocated (inline_mp%t_dt)) inline_mp%t_dt (is:ie, j, kmp:km) = &
                inline_mp%t_dt (is:ie, j, kmp:km) + pt (is:ie, j, kmp:km)
            if (allocated (inline_mp%u_dt)) inline_mp%u_dt (is:ie, j, kmp:km) = &
                inline_mp%u_dt (is:ie, j, kmp:km) + ua (is:ie, j, kmp:km)
            if (allocated (inline_mp%v_dt)) inline_mp%v_dt (is:ie, j, kmp:km) = &
                inline_mp%v_dt (is:ie, j, kmp:km) + va (is:ie, j, kmp:km)

            ! update pe, peln, pk, ps
            do k = kmp + 1, km + 1
                pe (is:ie, k, j) = pe (is:ie, k-1, j) + delp (is:ie, j, k-1)
                peln (is:ie, k, j) = log (pe (is:ie, k, j))
                pk (is:ie, j, k) = exp (akap * peln (is:ie, k, j))
            enddo

            ps (is:ie, j) = pe (is:ie, km+1, j)

            ! update pkz
            if (.not. hydrostatic) then
#ifdef MOIST_CAPPA
                pkz (is:ie, j, kmp:km) = exp (cappa (is:ie, j, kmp:km) * &
                    log (rrg * delp (is:ie, j, kmp:km) / &
                    delz (is:ie, j, kmp:km) * pt (is:ie, j, kmp:km)))
#else
                pkz (is:ie, j, kmp:km) = exp (akap * log (rrg * delp (is:ie, j, kmp:km) / &
                    delz (is:ie, j, kmp:km) * pt (is:ie, j, kmp:km)))
#endif
            endif
 
            if (consv .gt. consv_min) then
                do i = is, ie
                    do k = kmp, km
                        te0_2d (i, j) = te0_2d (i, j) + te (i, j, k)
                    enddo
                enddo
            endif

        enddo

        deallocate (dz)
        deallocate (wa)

        ! Note: (ua, va) are *lat-lon* wind tendenies on cell centers
        if ( gridstruct%square_domain ) then
            call mpp_update_domains (u_dt, domain, whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.false.)
            call mpp_update_domains (v_dt, domain, whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.true.)
        else
            call mpp_update_domains (u_dt, domain, complete=.false.)
            call mpp_update_domains (v_dt, domain, complete=.true.)
        endif
        ! update u_dt and v_dt in halo
        call mpp_update_domains (u_dt, v_dt, domain)

        ! update D grid wind
        call update_dwinds_phys (is, ie, js, je, isd, ied, jsd, jed, abs (mdt), u_dt, v_dt, u, v, &
                 gridstruct, npx, npy, km, domain)

        ! update dry total energy
        if (consv .gt. consv_min) then
!$OMP parallel do default (none) shared (is, ie, js, je, km, te0_2d, hydrostatic, delp, &
!$OMP                                    gridstruct, u, v, dp0, u0, v0, hs, delz, w) &
!$OMP                           private (phis)
            do j = js, je
                if (hydrostatic) then
                    do k = 1, km
                        do i = is, ie
                            te0_2d (i, j) = te0_2d (i, j) + delp (i, j, k) * &
                                (0.25 * gridstruct%rsin2 (i, j) * (u (i, j, k) ** 2 + &
                                u (i, j+1, k) ** 2 + v (i, j, k) ** 2 + v (i+1, j, k) ** 2 - &
                                (u (i, j, k) + u (i, j+1, k)) * (v (i, j, k) + v (i+1, j, k)) * &
                                gridstruct%cosa_s (i, j))) - dp0 (i, j, k) * &
                                (0.25 * gridstruct%rsin2 (i, j) * (u0 (i, j, k) ** 2 + &
                                u0 (i, j+1, k) ** 2 + v0 (i, j, k) ** 2 + v0 (i+1, j, k) ** 2 - &
                                (u0 (i, j, k) + u0 (i, j+1, k)) * (v0 (i, j, k) + v0 (i+1, j, k)) * &
                                gridstruct%cosa_s (i, j)))
                        enddo
                    enddo
                else
                    do i = is, ie
                        phis (i, km+1) = hs (i, j)
                    enddo
                    do k = km, 1, -1
                        do i = is, ie
                            phis (i, k) = phis (i, k+1) - grav * delz (i, j, k)
                        enddo
                    enddo
                    do k = 1, km
                        do i = is, ie
                            te0_2d (i, j) = te0_2d (i, j) + delp (i, j, k) * &
                                (0.5 * (phis (i, k) + phis (i, k+1) + w (i, j, k) ** 2 + 0.5 * &
                                gridstruct%rsin2 (i, j) * (u (i, j, k) ** 2 + u (i, j+1, k) ** 2 + &
                                v (i, j, k) ** 2 + v (i+1, j, k) ** 2 - (u (i, j, k) + &
                                u (i, j+1, k)) * (v (i, j, k) + v (i+1, j, k)) * &
                                gridstruct%cosa_s (i, j)))) - dp0 (i, j, k) * &
                                (0.5 * (phis (i, k) + phis (i, k+1) + w (i, j, k) ** 2 + &
                                0.5 * gridstruct%rsin2 (i, j) * (u0 (i, j, k) ** 2 + &
                                u0 (i, j+1, k) ** 2 + v0 (i, j, k) ** 2 + v0 (i+1, j, k) ** 2 - &
                                (u0 (i, j, k) + u0 (i, j+1, k)) * (v0 (i, j, k) + v0 (i+1, j, k)) * &
                                gridstruct%cosa_s (i, j))))
                        enddo
                    enddo
                endif
            enddo
        end if

        deallocate (u_dt)
        deallocate (v_dt)
        if (consv .gt. consv_min) then
            deallocate (u0)
            deallocate (v0)
            deallocate (dp0)
        endif

        call timing_off ('gfdl_mp')

    endif

    !-----------------------------------------------------------------------
    ! <<< Inline GFDL MP
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! Fast Spectral-Bin Microphysics >>>
    !-----------------------------------------------------------------------

    if ((.not. do_adiabatic_init) .and. do_inline_mp .and. do_fsbm) then

        call timing_on ('fsbm')

        f_sum = 0
        do n = 1, fsbm_bin

            ! normal distribution
            !mu = (1 + fsbm_bin) / 2.
            !sigma = 10.
            !f (n) = 1. / (sigma * sqrt (2. * pi)) * exp (- (n - mu) ** 2. / (2. * sigma ** 2.))

            ! gamma distribution
            alpha = 3.
            beta = 3.
            f (n) = 1. / (gamma (alpha) * beta ** alpha) * n ** (alpha - 1.) * exp (- n / beta)

            ! sum up
            f_sum = f_sum + f (n)

            ! get tracer index
            if (n .lt. 10) then
                write (ind, '(I1)') n
            else
                write (ind, '(I2)') n
            endif
            qlr_ind (n) = get_tracer_index (MODEL_ATMOS, 'qlr_' // trim (ind))
            qis_ind (n) = get_tracer_index (MODEL_ATMOS, 'qis_' // trim (ind))
            qg_ind (n) = get_tracer_index (MODEL_ATMOS, 'qg_' // trim (ind))
            qa_ind (n) = get_tracer_index (MODEL_ATMOS, 'qa_' // trim (ind))
            qn_ind (n) = get_tracer_index (MODEL_ATMOS, 'qn_' // trim (ind))

        enddo
        f = f / f_sum

        n_chem = fsbm_bin * 5
        num_sbmradar = 55

        allocate (xland (is:ie, js:je), rainnc (is:ie, js:je), rainncv (is:ie, js:je), snownc (is:ie, js:je))
        allocate (snowncv (is:ie, js:je), graupelnc (is:ie, js:je), graupelncv (is:ie, js:je))
        allocate (ur (is-1:ie+1, km, js-1:je+1), vr (is-1:ie+1, km, js-1:je+1), wr (is-1:ie+1, km, js-1:je+1))
        allocate (dz8w (is-1:ie+1, km, js-1:je+1), n_reg_ccn (is:ie, km, js:je))
        allocate (p_phy (is-1:ie+1, km, js-1:je+1), pi_phy (is-1:ie+1, km, js-1:je+1), sbqnn (is:ie, km, js:je))
        allocate (th_phy (is-1:ie+1, km, js-1:je+1), sbqv (is-1:ie+1, km, js-1:je+1), sbqc (is:ie, km, js:je))
        allocate (sbqr (is:ie, km, js:je), sbqi (is:ie, km, js:je), sbqs (is:ie, km, js:je))
        allocate (sbqg (is:ie, km, js:je), sbqnc (is:ie, km, js:je), sbqnr (is:ie, km, js:je))
        allocate (sbqni (is:ie, km, js:je), sbqns (is:ie, km, js:je), sbqng (is:ie, km, js:je))
        allocate (sbqna (is:ie, km, js:je), ma (is:ie, km, js:je), lh_rate (is:ie, km, js:je))
        allocate (ce_rate (is:ie, km, js:je), ds_rate (is:ie, km, js:je), melt_rate (is:ie, km, js:je))
        allocate (frz_rate (is:ie, km, js:je), cldnucl_rate (is:ie, km, js:je), icenucl_rate (is:ie, km, js:je))
        allocate (th_old (is-1:ie+1, km, js-1:je+1), qv_old (is-1:ie+1, km, js-1:je+1))
        allocate (pkz0 (isd:ied, jsd:jed, km), delz0 (isd:ied, jsd:jed, km))
        allocate (sbmradar (is:ie, km, js:je, num_sbmradar), rho_phy (is-1:ie+1, km, js-1:je+1))
        allocate (chem_new (is:ie, km, js:je, n_chem), dlnp (isd:ied, jsd:jed, km))

!$OMP parallel do default (none) shared (is, ie, js, je, km, pkz, pkz0, delz, delz0, dlnp, peln, hydrostatic)

        do k = 1, km
            do j = js, je
                do i = is, ie
                    pkz0 (i, j, k) = pkz (i, j, k)
                    if (hydrostatic) then
                        dlnp (i, j, k) = peln (i, k+1, j) - peln (i, k, j)
                    else
                        delz0 (i, j, k) = delz (i, j, k)
                    endif
                enddo
            enddo
        enddo

        call mpp_update_domains (pt, domain)
        call mpp_update_domains (ua, domain)
        call mpp_update_domains (va, domain)
        call mpp_update_domains (q (:,:,:, sphum), domain)
        call mpp_update_domains (delp, domain)
        if (hydrostatic) then
            call mpp_update_domains (dlnp, domain)
        else
            call mpp_update_domains (w, domain)
            call mpp_update_domains (delz0, domain)
        endif
        call mpp_update_domains (pkz0, domain)
        call mpp_update_domains (pt_old, domain)
        call mpp_update_domains (q_old, domain)

!$OMP parallel do default (none) shared (is, ie, js, je, km, ua, va, w, ur, vr, wr, dz8w, delz0, rho_phy, delp, &
!$OMP                                    p_phy, pt, pi_phy, pkz0, th_phy, pt_old, th_old, q_old, qv_old, &
!$OMP                                    sbqv, q, sphum, dlnp, hydrostatic, omga)

        do k = 1, km
            do j = js-1, je+1
                do i = is-1, ie+1

                    ur (i, k, j) = ua (i, j, km+1-k)
                    vr (i, k, j) = va (i, j, km+1-k)

                    if (hydrostatic) then
                        p_phy (i, k, j) = delp (i, j, km+1-k) / dlnp (i, j, km+1-k)
                        rho_phy (i, k, j) = p_phy (i, k, j) / (rdgas * pt (i, j, km+1-k))
                        dz8w (i, k, j) = delp (i, j, km+1-k) / (rho_phy (i, k, j) * grav)
                        wr (i, k, j) = - omga (i, j, km+1-k) * dz8w (i, k, j) / delp (i, j, km+1-k)
                    else
                        dz8w (i, k, j) = - delz0 (i, j, km+1-k)
                        rho_phy (i, k, j) = - delp (i, j, km+1-k) / delz0 (i, j, km+1-k) / grav
                        p_phy (i, k, j) = rho_phy (i, k, j) * rdgas * pt (i, j, km+1-k)
                        wr (i, k, j) = w (i, j, km+1-k)
                    endif

                    pi_phy (i, k, j) = pkz0 (i, j, km+1-k)
                    th_phy (i, k, j) = pt (i, j, km+1-k) / pi_phy (i, k, j)
                    th_old (i, k, j) = pt_old (i, j, km+1-k) / pi_phy (i, k, j)
                    qv_old (i, k, j) = q_old (i, j, km+1-k)
                    sbqv (i, k, j) = q (i, j, km+1-k, sphum)

                enddo
            enddo
        enddo

        do j = js, je
            do i = is, ie
                if (hs (i, j) .gt. 0) then
                    xland (i, j) = 1
                else
                    xland (i, j) = 0
                endif
            enddo
        enddo

!$OMP parallel do default (none) shared (is, ie, js, je, km, hs, rho_phy, delp, pt, q, &
!$OMP                                    chem_new, te, sphum, liq_wat, ice_wat, rainwat, &
!$OMP                                    snowwat, graupel, ma, lh_rate, ce_rate, ds_rate, melt_rate, &
!$OMP                                    frz_rate, consv, f, qlr_ind, qis_ind, qg_ind, qa_ind, qn_ind, a_step, &
!$OMP                                    fsbm_bin, r_vir, warm_start, itimestep) &
!$OMP                           private (qliq, qsol, cvm)

        do k = 1, km
            do j = js, je
                do i = is, ie

                    do n = 1, fsbm_bin
                        if (.not. warm_start .and. a_step .eq. 1) then
                            itimestep = 1
                            chem_new (i, k, j, fsbm_bin*0+n) = (q (i, j, km+1-k, liq_wat) + &
                                q (i, j, km+1-k, rainwat)) * f (n)
                            chem_new (i, k, j, fsbm_bin*1+n) = (q (i, j, km+1-k, ice_wat) + &
                                q (i, j, km+1-k, snowwat)) * f (n)
                            chem_new (i, k, j, fsbm_bin*2+n) = q (i, j, km+1-k, graupel) * f (n)
                            chem_new (i, k, j, fsbm_bin*3+n) = 1.e8
                            chem_new (i, k, j, fsbm_bin*4+n) = 0.0
                        else
                            itimestep = 2
                            chem_new (i, k, j, fsbm_bin*0+n) = q (i, j, km+1-k, qlr_ind (n))
                            chem_new (i, k, j, fsbm_bin*1+n) = q (i, j, km+1-k, qis_ind (n))
                            chem_new (i, k, j, fsbm_bin*2+n) = q (i, j, km+1-k, qg_ind (n))
                            chem_new (i, k, j, fsbm_bin*3+n) = q (i, j, km+1-k, qa_ind (n)) 
                            chem_new (i, k, j, fsbm_bin*4+n) = q (i, j, km+1-k, qn_ind (n)) 
                        endif
                    enddo

                    ma (i, k, j) = 0.0
                    lh_rate (i, k, j) = 0.0
                    ce_rate (i, k, j) = 0.0
                    ds_rate (i, k, j) = 0.0
                    melt_rate (i, k, j) = 0.0
                    frz_rate (i, k, j) = 0.0

                    if (consv .gt. consv_min) then
                        qliq = q (i, j, k, liq_wat) + q (i, j, k, rainwat)
                        qsol = q (i, j, k, ice_wat) + q (i, j, k, snowwat) + q (i, j, k, graupel)
                        cvm = (1 - (q (i, j, k, sphum) + qliq + qsol)) * cv_air + &
                            q (i, j, k, sphum) * cv_vap + qliq * c_liq + qsol * c_ice
                        te (i, j, k) = - cvm * pt (i, j, k) / ((1. + r_vir * q (i, j, k, sphum)) * &
                            (1. - (qliq + qsol))) * delp (i, j, k)
                    endif

                enddo
            enddo
        enddo

        !unit = stdout ()
        !write (unit,*) 'fsbm chksum before: temp', mpp_chksum (pt (is:ie, js:je,:))
        !write (unit,*) 'fsbm chksum before: qv', mpp_chksum (q (is:ie, js:je,:, sphum))
        !write (unit,*) 'fsbm chksum before: qml', mpp_chksum (q (is:ie, js:je,:, liq_wat))
        !write (unit,*) 'fsbm chksum before: qmr', mpp_chksum (q (is:ie, js:je,:, rainwat))
        !write (unit,*) 'fsbm chksum before: qmi', mpp_chksum (q (is:ie, js:je,:, ice_wat))
        !write (unit,*) 'fsbm chksum before: qms', mpp_chksum (q (is:ie, js:je,:, snowwat))
        !write (unit,*) 'fsbm chksum before: qmg', mpp_chksum (q (is:ie, js:je,:, graupel))
        !write (unit,*) 'fsbm chksum before: qnl', mpp_chksum (q (is:ie, js:je,:, ql_num))
        !write (unit,*) 'fsbm chksum before: qnr', mpp_chksum (q (is:ie, js:je,:, qr_num))
        !write (unit,*) 'fsbm chksum before: qni', mpp_chksum (q (is:ie, js:je,:, qi_num))
        !write (unit,*) 'fsbm chksum before: qns', mpp_chksum (q (is:ie, js:je,:, qs_num))
        !write (unit,*) 'fsbm chksum before: qng', mpp_chksum (q (is:ie, js:je,:, qg_num))
        !write (unit,*) 'fsbm chksum before: qlr', mpp_chksum (q (is:ie, js:je,:, 7+1:7+33))
        !write (unit,*) 'fsbm chksum before: qis', mpp_chksum (q (is:ie, js:je,:, 7+34:7+66))
        !write (unit,*) 'fsbm chksum before: qgg', mpp_chksum (q (is:ie, js:je,:, 7+67:7+99))
        !write (unit,*) 'fsbm chksum before: ccn', mpp_chksum (q (is:ie, js:je,:, 7+100:7+132))
  
!$OMP parallel do default (none) shared (is, ie, js, je, km, inline_mp, q, sphum, liq_wat, &
!$OMP                                    rainwat, ice_wat, snowwat, graupel, pt)

        do k = 1, km
            do j = js, je
                do i = is, ie
                    if (allocated (inline_mp%liq_wat_dt)) &
                        inline_mp%liq_wat_dt (i, j, k) = inline_mp%liq_wat_dt (i, j, k) - q (i, j, k, liq_wat)
                    if (allocated (inline_mp%ice_wat_dt)) &
                        inline_mp%ice_wat_dt (i, j, k) = inline_mp%ice_wat_dt (i, j, k) - q (i, j, k, ice_wat)
                    if (allocated (inline_mp%qv_dt)) &
                        inline_mp%qv_dt (i, j, k) = inline_mp%qv_dt (i, j, k) - q (i, j, k, sphum)
                    if (allocated (inline_mp%ql_dt)) &
                        inline_mp%ql_dt (i, j, k) = inline_mp%ql_dt (i, j, k) - &
                        (q (i, j, k, liq_wat) + q (i, j, k, rainwat))
                    if (allocated (inline_mp%qi_dt)) &
                        inline_mp%qi_dt (i, j, k) = inline_mp%qi_dt (i, j, k) - &
                        (q (i, j, k, ice_wat) + q (i, j, k, snowwat) + q (i, j, k, graupel))
                    if (allocated (inline_mp%qr_dt)) &
                        inline_mp%qr_dt (i, j, k) = inline_mp%qr_dt (i, j, k) - q (i, j, k, rainwat)
                    if (allocated (inline_mp%qs_dt)) &
                        inline_mp%qs_dt (i, j, k) = inline_mp%qs_dt (i, j, k) - q (i, j, k, snowwat)
                    if (allocated (inline_mp%qg_dt)) &
                        inline_mp%qg_dt (i, j, k) = inline_mp%qg_dt (i, j, k) - q (i, j, k, graupel)
                    if (allocated (inline_mp%t_dt)) &
                        inline_mp%t_dt (i, j, k) = inline_mp%t_dt (i, j, k) - pt (i, j, k)
                enddo
            enddo
        enddo

        call fast_sbm (wr, ur, vr, th_old, chem_new, n_chem, itimestep, abs (mdt), fsbm_dx, &
            fsbm_dy, dz8w, rho_phy, p_phy, pi_phy, th_phy, xland, sbqv, sbqc, sbqr, sbqi, &
            sbqs, sbqg, qv_old, sbqnc, sbqnr, sbqni, sbqns, sbqng, sbqna, sbqnn, 1, npx, &
            1, npy, 1, km, is, ie, js, je, 1, km, is, ie, js, je, 1, km, diagflag, sbmradar, &
            num_sbmradar, rainnc, rainncv, snownc, snowncv, graupelnc, graupelncv, ma, &
            lh_rate, ce_rate, ds_rate, melt_rate, frz_rate, cldnucl_rate, icenucl_rate, &
            n_reg_ccn)

        do j = js, je
            do i = is, ie
                if (do_inline_mp) then
                    inline_mp%prer (i, j) = inline_mp%prer (i, j) + rainncv (i, j) / abs (mdt) * 86400
                    inline_mp%pres (i, j) = inline_mp%pres (i, j) + snowncv (i, j) / abs (mdt) * 86400
                    inline_mp%preg (i, j) = inline_mp%preg (i, j) + graupelncv (i, j) / abs (mdt) * 86400
                endif
            enddo
        enddo

!$OMP parallel do default (none) shared (is, ie, js, je, km, sbqv, sbqc, sbqr, sbqi, sbqs, sbqg, q, th_phy, &
!$OMP                                    pi_phy, th_old, pt, pt_old, q_old, qv_old, q_con, cappa, r_vir, &
!$OMP                                    te, delp, sphum, liq_wat, ice_wat, rainwat, snowwat, graupel, &
!$OMP                                    consv, cld_amt, rho_phy, qlr_ind, qis_ind, qg_ind, qa_ind, qn_ind,&
!$OMP                                    chem_new, fsbm_bin, te0_2d, ql_num, qr_num, qi_num, qs_num, &
!$OMP                                    qg_num, qa_num, qn_num, sbqnc, sbqnr, sbqni, sbqns, sbqng,&
!$OMP                                    sbqna, sbqnn) &
!$OMP                           private (qliq, qsol, cvm, dqv, dql, dqr, dqi, dqs, dqg, rh, qsat, ps_dt)

        do k = 1, km
            do j = js, je
                do i = is, ie

                    dqv = sbqv (i, km+1-k, j) - q (i, j, k, sphum)
                    dql = sbqc (i, km+1-k, j) - q (i, j, k, liq_wat)
                    dqr = sbqr (i, km+1-k, j) - q (i, j, k, rainwat)
                    dqi = sbqi (i, km+1-k, j) - q (i, j, k, ice_wat)
                    dqs = sbqs (i, km+1-k, j) - q (i, j, k, snowwat)
                    dqg = sbqg (i, km+1-k, j) - q (i, j, k, graupel)
                    ps_dt = 1 + dqv + dql + dqr + dqi + dqs + dqg
                    q (i, j, k, sphum) = sbqv (i, km+1-k, j) / ps_dt
                    q (i, j, k, liq_wat) = sbqc (i, km+1-k, j) / ps_dt
                    q (i, j, k, rainwat) = sbqr (i, km+1-k, j) / ps_dt
                    q (i, j, k, ice_wat) = sbqi (i, km+1-k, j) / ps_dt
                    q (i, j, k, snowwat) = sbqs (i, km+1-k, j) / ps_dt
                    q (i, j, k, graupel) = sbqg (i, km+1-k, j) / ps_dt
                    q (i, j, k, ql_num) = sbqnc (i, km+1-k, j) / ps_dt
                    q (i, j, k, qr_num) = sbqnr (i, km+1-k, j) / ps_dt
                    q (i, j, k, qi_num) = sbqni (i, km+1-k, j) / ps_dt
                    q (i, j, k, qs_num) = sbqns (i, km+1-k, j) / ps_dt
                    q (i, j, k, qg_num) = sbqng (i, km+1-k, j) / ps_dt
                    q (i, j, k, qa_num) = sbqna (i, km+1-k, j) / ps_dt
                    q (i, j, k, qn_num) = sbqnn (i, km+1-k, j) / ps_dt
                    pt (i, j, k) = th_phy (i, km+1-k, j) * pi_phy (i, km+1-k, j)
                    pt_old (i, j, k) = th_old (i, km+1-k, j) * pi_phy (i, km+1-k, j)
                    q_old (i, j, k) = qv_old (i, km+1-k, j)
                    do n = 1, fsbm_bin
                        q (i, j, k, qlr_ind (n)) = chem_new (i, km+1-k, j, fsbm_bin*0+n)
                        q (i, j, k, qis_ind (n)) = chem_new (i, km+1-k, j, fsbm_bin*1+n)
                        q (i, j, k, qg_ind (n)) = chem_new (i, km+1-k, j, fsbm_bin*2+n)
                        q (i, j, k, qa_ind (n)) = chem_new (i, km+1-k, j, fsbm_bin*3+n)
                        q (i, j, k, qn_ind (n)) = chem_new (i, km+1-k, j, fsbm_bin*4+n)
                    enddo

                    qliq = q (i, j, k, liq_wat) + q (i, j, k, rainwat)
                    qsol = q (i, j, k, ice_wat) + q (i, j, k, snowwat) + q (i, j, k, graupel)
                    qsat = iqs (pt (i, j, k), rho_phy (i, km+1-k, j), qsat)
                    rh = q (i, j, k, sphum) / qsat

                    if (rh >= 1.0) then
                        q (i, j, k, cld_amt) = 1.0
                    elseif (rh > 0.75 .and. qliq + qsol > 1.e-6) then
                        q (i, j, k, cld_amt) = rh ** xr_a * (1.0 - exp (- xr_b * max (0.0, qliq + qsol) / &
                            max (1.e-5, (max (1.e-10, 1.0 - rh) * qsat) ** xr_c)))
                        q (i, j, k, cld_amt) = max (0.0, min (1., q (i, j, k, cld_amt)))
                    else
                        q (i, j, k, cld_amt) = 0.0
                    endif

                    cvm = (1 - (q (i, j, k, sphum) + qliq + qsol)) * cv_air + &
                        q (i, j, k, sphum) * cv_vap + qliq * c_liq + qsol * c_ice
#ifdef USE_COND
                    q_con (i, j, k) = qliq + qsol
#endif
#ifdef MOIST_CAPPA
                    cappa (i, j, k) = rdgas / (rdgas + cvm / (1. + r_vir * q (i, j, k, sphum)))
#endif

                    delp (i, j, k) = delp (i, j, k) * ps_dt
                    if (consv .gt. consv_min) then
                        te (i, j, k) = te (i, j, k) + &
                            cvm * pt (i, j, k) / ((1. + r_vir * q (i, j, k, sphum)) * &
                            (1. - (qliq + qsol))) * delp (i, j, k)
                        te0_2d (i, j) = te0_2d (i, j) + te (i, j, k)
                    endif

                enddo
            enddo
        enddo

!$OMP parallel do default (none) shared (is, ie, js, je, km, inline_mp, q, sphum, liq_wat, &
!$OMP                                    rainwat, ice_wat, snowwat, graupel, pt, mdt)

        do k = 1, km
            do j = js, je
                do i = is, ie
                    if (allocated (inline_mp%liq_wat_dt)) &
                        inline_mp%liq_wat_dt (i, j, k) = inline_mp%liq_wat_dt (i, j, k) + q (i, j, k, liq_wat)
                    if (allocated (inline_mp%ice_wat_dt)) &
                        inline_mp%ice_wat_dt (i, j, k) = inline_mp%ice_wat_dt (i, j, k) + q (i, j, k, ice_wat)
                    if (allocated (inline_mp%qv_dt)) &
                        inline_mp%qv_dt (i, j, k) = inline_mp%qv_dt (i, j, k) + q (i, j, k, sphum)
                    if (allocated (inline_mp%ql_dt)) &
                        inline_mp%ql_dt (i, j, k) = inline_mp%ql_dt (i, j, k) + &
                        (q (i, j, k, liq_wat) + q (i, j, k, rainwat))
                    if (allocated (inline_mp%qi_dt)) &
                        inline_mp%qi_dt (i, j, k) = inline_mp%qi_dt (i, j, k) + &
                        (q (i, j, k, ice_wat) + q (i, j, k, snowwat) + q (i, j, k, graupel))
                    if (allocated (inline_mp%qr_dt)) &
                        inline_mp%qr_dt (i, j, k) = inline_mp%qr_dt (i, j, k) + q (i, j, k, rainwat)
                    if (allocated (inline_mp%qs_dt)) &
                        inline_mp%qs_dt (i, j, k) = inline_mp%qs_dt (i, j, k) + q (i, j, k, snowwat)
                    if (allocated (inline_mp%qg_dt)) &
                        inline_mp%qg_dt (i, j, k) = inline_mp%qg_dt (i, j, k) + q (i, j, k, graupel)
                    if (allocated (inline_mp%t_dt)) &
                        inline_mp%t_dt (i, j, k) = inline_mp%t_dt (i, j, k) + pt (i, j, k)
                enddo
            enddo
        enddo

        !unit = stdout ()
        !write (unit,*) 'fsbm chksum after: temp', mpp_chksum (pt (is:ie, js:je,:))
        !write (unit,*) 'fsbm chksum after: qv', mpp_chksum (q (is:ie, js:je,:, sphum))
        !write (unit,*) 'fsbm chksum after: qml', mpp_chksum (q (is:ie, js:je,:, liq_wat))
        !write (unit,*) 'fsbm chksum after: qmr', mpp_chksum (q (is:ie, js:je,:, rainwat))
        !write (unit,*) 'fsbm chksum after: qmi', mpp_chksum (q (is:ie, js:je,:, ice_wat))
        !write (unit,*) 'fsbm chksum after: qms', mpp_chksum (q (is:ie, js:je,:, snowwat))
        !write (unit,*) 'fsbm chksum after: qmg', mpp_chksum (q (is:ie, js:je,:, graupel))
        !write (unit,*) 'fsbm chksum after: qnl', mpp_chksum (q (is:ie, js:je,:, ql_num))
        !write (unit,*) 'fsbm chksum after: qnr', mpp_chksum (q (is:ie, js:je,:, qr_num))
        !write (unit,*) 'fsbm chksum after: qni', mpp_chksum (q (is:ie, js:je,:, qi_num))
        !write (unit,*) 'fsbm chksum after: qns', mpp_chksum (q (is:ie, js:je,:, qs_num))
        !write (unit,*) 'fsbm chksum after: qng', mpp_chksum (q (is:ie, js:je,:, qg_num))
        !write (unit,*) 'fsbm chksum after: qlr', mpp_chksum (q (is:ie, js:je,:, 7+1:7+33))
        !write (unit,*) 'fsbm chksum after: qis', mpp_chksum (q (is:ie, js:je,:, 7+34:7+66))
        !write (unit,*) 'fsbm chksum after: qgg', mpp_chksum (q (is:ie, js:je,:, 7+67:7+99))
        !write (unit,*) 'fsbm chksum after: ccn', mpp_chksum (q (is:ie, js:je,:, 7+100:7+132))
  
        deallocate (xland, rainnc, rainncv, snownc, snowncv, graupelnc, graupelncv, sbmradar, chem_new)
        deallocate (ur, vr, wr, dz8w, p_phy, pi_phy, rho_phy, th_phy, pkz0, delz0, dlnp, th_old, qv_old)
        deallocate (sbqv, sbqc, sbqr, sbqi, sbqs, sbqg, sbqnc, sbqnr, sbqni, sbqns, sbqng, sbqna, sbqnn)
        deallocate (ma, lh_rate, ce_rate, ds_rate, melt_rate, frz_rate, cldnucl_rate, icenucl_rate, n_reg_ccn)

        call timing_off ('fsbm')

    endif

    !-----------------------------------------------------------------------
    ! <<< Fast Spectral-Bin Microphysics
    !-----------------------------------------------------------------------

end subroutine fast_phys

end module fast_phys_mod
