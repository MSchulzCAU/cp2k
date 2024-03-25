!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
!> \brief TDA-Bethe Salpeter Equation(BSE).
! **************************************************************************************************
MODULE qs_tddfpt2_bse_utils

   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_control_types,                ONLY: stda_control_type,&
                                              tddfpt2_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              cp_dbcsr_plus_fm_fm_t,&
                                              cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_allocate_matrix_set
   USE cp_fm_basic_linalg,              ONLY: cp_fm_row_scale,&
                                              cp_fm_schur_product
   USE cp_fm_diag,                      ONLY: choose_eigv_solver,&
                                              cp_fm_power
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_set_submatrix,&
                                              cp_fm_to_fm,&
                                              cp_fm_type,&
                                              cp_fm_vectorssum
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add_on_diag, dbcsr_create, dbcsr_distribution_type, dbcsr_filter, dbcsr_finalize, &
        dbcsr_get_block_p, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_p_type, &
        dbcsr_release, dbcsr_set, dbcsr_type, dbcsr_type_antisymmetric, dbcsr_type_no_symmetry, &
        dbcsr_type_symmetric
   USE ewald_environment_types,         ONLY: ewald_env_create,&
                                              ewald_env_get,&
                                              ewald_env_set,&
                                              ewald_environment_type,&
                                              read_ewald_section_tb
   USE ewald_methods_tb,                ONLY: tb_ewald_overlap,&
                                              tb_spme_evaluate
   USE ewald_pw_types,                  ONLY: ewald_pw_create,&
                                              ewald_pw_type
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE iterate_matrix,                  ONLY: matrix_sqrt_Newton_Schulz
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: oorootpi
   USE message_passing,                 ONLY: mp_para_env_type
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind_set,&
                                              qs_kind_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
   USE qs_tddfpt2_bse_types,            ONLY: bse_env_type
   USE qs_tddfpt2_subgroups,            ONLY: tddfpt_subgroup_env_type
   USE qs_tddfpt2_types,                ONLY: tddfpt_work_matrices
   USE scf_control_types,               ONLY: scf_control_type
   USE util,                            ONLY: get_limit
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddfpt2_bse_utils'

   PUBLIC:: bse_init_matrices!, stda_calculate_kernel, get_lowdin_x, setup_gamma

CONTAINS

! **************************************************************************************************
!> \brief Calculate sTDA matrices
!> \param qs_env ...
!> \param bse_kernel ...
!> \param sub_env ...
!> \param work ...
!> \param tddfpt_control ...
! **************************************************************************************************
   SUBROUTINE bse_init_matrices(qs_env, bse_kernel, sub_env, work, tddfpt_control)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(bse_env_type)                                 :: bse_kernel
      TYPE(tddfpt_subgroup_env_type), INTENT(in)         :: sub_env
      TYPE(tddfpt_work_matrices)                         :: work
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control

      CHARACTER(len=*), PARAMETER :: routineN = 'bse_init_matrices'

      INTEGER                                            :: handle
      LOGICAL                                            :: do_coulomb
      TYPE(cell_type), POINTER                           :: cell, cell_ref
      TYPE(ewald_environment_type), POINTER              :: ewald_env
      TYPE(ewald_pw_type), POINTER                       :: ewald_pw
      TYPE(section_vals_type), POINTER                   :: ewald_section, poisson_section, &
                                                            print_section
      TYPE(dbcsr_type)                                   :: shalf 

      CALL timeset(routineN, handle)

      do_coulomb = .NOT. tddfpt_control%rks_triplets
      IF (do_coulomb) THEN
         ! calculate exchange gamma matrix
         CALL setup_gamma(qs_env, bse_kernel, sub_env, work%gamma_exchange)
      END IF

      ! calculate S_half and Lowdin MO coefficients
      CALL get_lowdin_mo_coefficients(qs_env, sub_env, work)

      shalf = work%shalf
      CALL timestop(handle)

   END SUBROUTINE bse_init_matrices

! **************************************************************************************************
!> \brief Calculate sTDA exchange-type contributions
!> \param qs_env ...
!> \param stda_env ...
!> \param sub_env ...
!> \param gamma_matrix sTDA exchange-type contributions
!> \param ndim ...
!> \note  Note the specific sTDA notation exchange-type integrals (ia|jb) refer to Coulomb interaction
! **************************************************************************************************
   SUBROUTINE setup_gamma(qs_env, stda_env, sub_env, gamma_matrix, ndim)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(bse_env_type)                                :: stda_env
      TYPE(tddfpt_subgroup_env_type), INTENT(in)         :: sub_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: gamma_matrix
      INTEGER, INTENT(IN), OPTIONAL                      :: ndim

      CHARACTER(len=*), PARAMETER                        :: routineN = 'setup_gamma'
      REAL(KIND=dp), PARAMETER                           :: rsmooth = 1.0_dp

      INTEGER                                            :: handle, i, iatom, icol, ikind, imat, &
                                                            irow, jatom, jkind, natom, nmat
      INTEGER, DIMENSION(:), POINTER                     :: row_blk_sizes
      LOGICAL                                            :: found
      REAL(KIND=dp)                                      :: dfcut, dgb, dr, eta, fcut, r, rcut, &
                                                            rcuta, rcutb, x
      REAL(KIND=dp), DIMENSION(3)                        :: rij
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: dgblock, gblock
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: n_list

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, natom=natom)
      dbcsr_dist => sub_env%dbcsr_dist
      ! Using the overlap list here can have a considerable effect on the number of
      ! terms calculated. This makes gamma also dependent on EPS_DEFAULT -> Overlap
      n_list => sub_env%sab_orb

      IF (PRESENT(ndim)) THEN
         nmat = ndim
      ELSE
         nmat = 1
      END IF
      CPASSERT(nmat == 1 .OR. nmat == 4)
      CPASSERT(.NOT. ASSOCIATED(gamma_matrix))
      CALL dbcsr_allocate_matrix_set(gamma_matrix, nmat)

      ALLOCATE (row_blk_sizes(natom))
      row_blk_sizes(1:natom) = 1
      DO imat = 1, nmat
         ALLOCATE (gamma_matrix(imat)%matrix)
      END DO

      CALL dbcsr_create(gamma_matrix(1)%matrix, name="gamma", dist=dbcsr_dist, &
                        matrix_type=dbcsr_type_symmetric, row_blk_size=row_blk_sizes, &
                        col_blk_size=row_blk_sizes, nze=0)
      DO imat = 2, nmat
         CALL dbcsr_create(gamma_matrix(imat)%matrix, name="dgamma", dist=dbcsr_dist, &
                           matrix_type=dbcsr_type_antisymmetric, row_blk_size=row_blk_sizes, &
                           col_blk_size=row_blk_sizes, nze=0)
      END DO

      DEALLOCATE (row_blk_sizes)

      ! setup the matrices using the neighbor list
      DO imat = 1, nmat
         CALL cp_dbcsr_alloc_block_from_nbl(gamma_matrix(imat)%matrix, n_list)
         CALL dbcsr_set(gamma_matrix(imat)%matrix, 0.0_dp)
      END DO

      NULLIFY (nl_iterator)
      CALL neighbor_list_iterator_create(nl_iterator, n_list)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                iatom=iatom, jatom=jatom, r=rij)

         dr = SQRT(SUM(rij(:)**2)) ! interatomic distance

         eta = (stda_env%kind_param_set(ikind)%kind_param%hardness_param + &
                stda_env%kind_param_set(jkind)%kind_param%hardness_param)/2.0_dp

         icol = MAX(iatom, jatom)
         irow = MIN(iatom, jatom)

         NULLIFY (gblock)
         CALL dbcsr_get_block_p(matrix=gamma_matrix(1)%matrix, &
                                row=irow, col=icol, BLOCK=gblock, found=found)
         CPASSERT(found)

         ! get rcuta and rcutb
         rcuta = stda_env%kind_param_set(ikind)%kind_param%rcut
         rcutb = stda_env%kind_param_set(jkind)%kind_param%rcut
         rcut = rcuta + rcutb

         !>   Computes the short-range gamma parameter from
         !>   Nataga-Mishimoto-Ohno-Klopman formula equivalently as it is done for xTB
         IF (dr < 1.e-6) THEN
            ! on site terms
            gblock(:, :) = gblock(:, :) + eta
         ELSEIF (dr > rcut) THEN
            ! do nothing
         ELSE
            IF (dr < rcut - rsmooth) THEN
               fcut = 1.0_dp
            ELSE
               r = dr - (rcut - rsmooth)
               x = r/rsmooth
               fcut = -6._dp*x**5 + 15._dp*x**4 - 10._dp*x**3 + 1._dp
            END IF
            gblock(:, :) = gblock(:, :) + &
                           fcut*(1._dp/(dr**(stda_env%alpha_param) + eta**(-stda_env%alpha_param))) &
                           **(1._dp/stda_env%alpha_param) - fcut/dr
         END IF

         IF (nmat > 1) THEN
            !>   Computes the short-range gamma parameter from
            !>   Nataga-Mishimoto-Ohno-Klopman formula equivalently as it is done for xTB
            !>   Derivatives
            IF (dr < 1.e-6 .OR. dr > rcut) THEN
               ! on site terms or beyond cutoff
               dgb = 0.0_dp
            ELSE
               IF (dr < rcut - rsmooth) THEN
                  fcut = 1.0_dp
                  dfcut = 0.0_dp
               ELSE
                  r = dr - (rcut - rsmooth)
                  x = r/rsmooth
                  fcut = -6._dp*x**5 + 15._dp*x**4 - 10._dp*x**3 + 1._dp
                  dfcut = -30._dp*x**4 + 60._dp*x**3 - 30._dp*x**2
                  dfcut = dfcut/rsmooth
               END IF
               dgb = dfcut*(1._dp/(dr**(stda_env%alpha_param) + eta**(-stda_env%alpha_param))) &
                     **(1._dp/stda_env%alpha_param)
               dgb = dgb - dfcut/dr + fcut/dr**2
               dgb = dgb - fcut*(1._dp/(dr**(stda_env%alpha_param) + eta**(-stda_env%alpha_param))) &
                     **(1._dp/stda_env%alpha_param + 1._dp)*dr**(stda_env%alpha_param - 1._dp)
            END IF
            DO imat = 2, nmat
               NULLIFY (dgblock)
               CALL dbcsr_get_block_p(matrix=gamma_matrix(imat)%matrix, &
                                      row=irow, col=icol, BLOCK=dgblock, found=found)
               IF (found) THEN
                  IF (dr > 1.e-6) THEN
                     i = imat - 1
                     IF (irow == iatom) THEN
                        dgblock(:, :) = dgblock(:, :) + dgb*rij(i)/dr
                     ELSE
                        dgblock(:, :) = dgblock(:, :) - dgb*rij(i)/dr
                     END IF
                  END IF
               END IF
            END DO
         END IF

      END DO

      CALL neighbor_list_iterator_release(nl_iterator)

      DO imat = 1, nmat
         CALL dbcsr_finalize(gamma_matrix(imat)%matrix)
      END DO

      CALL timestop(handle)

   END SUBROUTINE setup_gamma

! **************************************************************************************************
!> \brief Calculate Lowdin MO coefficients
!> \param qs_env ...
!> \param sub_env ...
!> \param work ...
! **************************************************************************************************
   SUBROUTINE get_lowdin_mo_coefficients(qs_env, sub_env, work)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddfpt_subgroup_env_type), INTENT(in)         :: sub_env
      TYPE(tddfpt_work_matrices)                         :: work

      CHARACTER(len=*), PARAMETER :: routineN = 'get_lowdin_mo_coefficients'

      INTEGER                                            :: handle, i, iounit, ispin, j, &
                                                            max_iter_lanczos, nactive, ndep, nsgf, &
                                                            nspins, order_lanczos
      LOGICAL                                            :: converged
      REAL(KIND=dp)                                      :: eps_lanczos, sij, threshold
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: slam
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         POINTER                                         :: local_data
      TYPE(cp_fm_struct_type), POINTER                   :: fmstruct
      TYPE(cp_fm_type)                                   :: fm_s_half, fm_work1
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrixkp_s
      TYPE(dbcsr_type)                                   :: sm_hinv
      TYPE(dbcsr_type), POINTER                          :: sm_h, sm_s
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(scf_control_type), POINTER                    :: scf_control

      CALL timeset(routineN, handle)

      NULLIFY (logger) !get output_unit
      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)

      ! Calculate S^1/2 matrix
      IF (iounit > 0) THEN
         WRITE (iounit, "(1X,A)") "", &
            "-------------------------------------------------------------------------------", &
            "-                      sTDA Kernel: Create Matrix SQRT(S)                     -", &
            "-------------------------------------------------------------------------------"
      END IF

      IF (sub_env%is_split) THEN
         CPABORT('SPLIT')
      ELSE
         CALL get_qs_env(qs_env=qs_env, matrix_s_kp=matrixkp_s)
         CPASSERT(ASSOCIATED(matrixkp_s))
         IF (SIZE(matrixkp_s, 2) > 1) CPWARN("not implemented for k-points.")
         sm_s => matrixkp_s(1, 1)%matrix
      END IF
      sm_h => work%shalf

      CALL dbcsr_create(sm_hinv, template=sm_s)
      CALL dbcsr_add_on_diag(sm_h, 1.0_dp)
      threshold = 1.0e-8_dp
      order_lanczos = 3
      eps_lanczos = 1.0e-4_dp
      max_iter_lanczos = 40
      CALL matrix_sqrt_Newton_Schulz(sm_h, sm_hinv, sm_s, &
                                     threshold, order_lanczos, eps_lanczos, max_iter_lanczos, &
                                     converged=converged)
      CALL dbcsr_release(sm_hinv)
      !
      NULLIFY (qs_kind_set)
      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set)
      ! Get the total number of contracted spherical Gaussian basis functions
      CALL get_qs_kind_set(qs_kind_set, nsgf=nsgf)
      !
      IF (.NOT. converged) THEN
         IF (iounit > 0) THEN
            WRITE (iounit, "(T3,A)") "STDA| Newton-Schulz iteration did not converge"
            WRITE (iounit, "(T3,A)") "STDA| Calculate SQRT(S) from diagonalization"
         END IF
         CALL get_qs_env(qs_env=qs_env, scf_control=scf_control)
         ! Provide full size work matrices
         CALL cp_fm_struct_create(fmstruct=fmstruct, &
                                  para_env=sub_env%para_env, &
                                  context=sub_env%blacs_env, &
                                  nrow_global=nsgf, &
                                  ncol_global=nsgf)
         CALL cp_fm_create(matrix=fm_s_half, matrix_struct=fmstruct, name="S^(1/2) MATRIX")
         CALL cp_fm_create(matrix=fm_work1, matrix_struct=fmstruct, name="TMP MATRIX")
         CALL cp_fm_struct_release(fmstruct=fmstruct)
         CALL copy_dbcsr_to_fm(sm_s, fm_s_half)
         CALL cp_fm_power(fm_s_half, fm_work1, 0.5_dp, scf_control%eps_eigval, ndep)
         IF (ndep /= 0) &
            CALL cp_warn(__LOCATION__, &
                         "Overlap matrix exhibits linear dependencies. At least some "// &
                         "eigenvalues have been quenched.")
         CALL copy_fm_to_dbcsr(fm_s_half, sm_h)
         CALL cp_fm_release(fm_s_half)
         CALL cp_fm_release(fm_work1)
         IF (iounit > 0) WRITE (iounit, *)
      END IF

      nspins = SIZE(sub_env%mos_occ)

      DO ispin = 1, nspins
         CALL cp_fm_get_info(work%ctransformed(ispin), ncol_global=nactive)
         CALL cp_dbcsr_sm_fm_multiply(work%shalf, sub_env%mos_occ(ispin), &
                                      work%ctransformed(ispin), nactive, alpha=1.0_dp, beta=0.0_dp)
      END DO

      ! for Lowdin forces
      CALL cp_fm_create(matrix=fm_work1, matrix_struct=work%S_eigenvectors%matrix_struct, name="TMP MATRIX")
      CALL copy_dbcsr_to_fm(sm_s, fm_work1)
      CALL choose_eigv_solver(fm_work1, work%S_eigenvectors, work%S_eigenvalues)
      CALL cp_fm_release(fm_work1)
      !
      ALLOCATE (slam(nsgf, 1))
      DO i = 1, nsgf
         IF (work%S_eigenvalues(i) > 0._dp) THEN
            slam(i, 1) = SQRT(work%S_eigenvalues(i))
         ELSE
            CPABORT("S matrix not positive definit")
         END IF
      END DO
      DO i = 1, nsgf
         CALL cp_fm_set_submatrix(work%slambda, slam, 1, i, nsgf, 1, 1.0_dp, 0.0_dp)
      END DO
      DO i = 1, nsgf
         CALL cp_fm_set_submatrix(work%slambda, slam, i, 1, 1, nsgf, 1.0_dp, 1.0_dp, .TRUE.)
      END DO
      CALL cp_fm_get_info(work%slambda, local_data=local_data)
      DO i = 1, SIZE(local_data, 2)
         DO j = 1, SIZE(local_data, 1)
            sij = local_data(j, i)
            IF (sij > 0.0_dp) sij = 1.0_dp/sij
            local_data(j, i) = sij
         END DO
      END DO
      DEALLOCATE (slam)

      CALL timestop(handle)

   END SUBROUTINE get_lowdin_mo_coefficients


END MODULE qs_tddfpt2_bse_utils
