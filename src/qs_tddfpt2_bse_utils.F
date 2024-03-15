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
!> \param stda_kernel ...
!> \param sub_env ...
!> \param work ...
!> \param tddfpt_control ...
! **************************************************************************************************
   SUBROUTINE bse_init_matrices(qs_env, stda_kernel, sub_env, work, tddfpt_control)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(bse_env_type)                                :: stda_kernel
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

      CALL timeset(routineN, handle)

      !do_coulomb = .NOT. tddfpt_control%rks_triplets
      !IF (do_coulomb) THEN
      !   ! calculate exchange gamma matrix
      !   CALL setup_gamma(qs_env, stda_kernel, sub_env, work%gamma_exchange)
      !END IF

      ! calculate S_half and Lowdin MO coefficients
      !CALL get_lowdin_mo_coefficients(qs_env, sub_env, work)

      CALL timestop(handle)

   END SUBROUTINE bse_init_matrices


END MODULE qs_tddfpt2_bse_utils
