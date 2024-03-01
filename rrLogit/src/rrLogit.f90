!#####################################################################
!# wrapper functions for rrLogit shared object
!#####################################################################
subroutine rrlogit( dim_vec, &
     model_matrix, fitweight_row_input_data, fitweight_cov_patt, &
     fitweight_data_patt, freq_for_data_patt, cov_patt, data_patt, &
     cov_patt_for_data_patt, response_for_data_patt, &
     survey_mode_int, baseline_int, pert_mat, pert_mat_inv, &
     prior_int, prior_freq_tot, prior_alloc_supplied_int, prior_alloc, &
     saturated_int, method_int, ctrl_int, ctrl_real, special_codes, &
     dim_vec_mcmc, vhat_coef_rwm, &
     micro_data_int, freq_row_input_data_int, &
     f_mat_int, response_row_input_data, freq_for_data_patt_int, &
     start_val_source_int, &
     outputs_int, outputs_real, score, hess, &
     coefficients, coef_vec, vhat_coef, fitted_pi, fitted_pistar, &
     marg_proportions_ystar, marg_proportions_y, fhat_mat, hessA, &
     fstar_mat, fitweight_mat, &
     coef_mcmc, coef_vec_mcmc, vhat_coef_mcmc, &
     fitted_pi_mcmc, fitted_pistar_mcmc, &
     n_actual, mh_accept_rate, start_logP, &
     logP_series, coef_vec_series, imp_mat_series, imp_vec_series, &
     pi_marg_mcmc, pi_marg_series, approx_bayes_log_imp_ratios, &
     status, msg_len_max, msg_codes, msg_len_actual )
   !#############################################################
   ! This is a wrapper function for run_rrlogit
   ! status = 0 means everything ran OK, or run was aborted for
   !    some reason (see msg)
   ! Other value of status indicates a fatal error.
   !#############################################################
   use error_handler
   use program_constants
   use dynalloc
   use math_R
   use rrlogit_engine
   implicit none
   ! declare input arguments
   integer(kind=our_int), intent(in) :: dim_vec(8)
   real(kind=our_dble), intent(in) :: model_matrix( dim_vec(3), dim_vec(5) )
   real(kind=our_dble), intent(in) :: fitweight_row_input_data( dim_vec(1) )
   real(kind=our_dble), intent(in) :: fitweight_cov_patt( dim_vec(3) )
   real(kind=our_dble), intent(in) :: fitweight_data_patt( dim_vec(4) )
   real(kind=our_dble), intent(in) :: freq_for_data_patt( dim_vec(4) )
   integer(kind=our_int), intent(in) :: cov_patt( dim_vec(1) )
   integer(kind=our_int), intent(in) :: &
        data_patt( dim_vec(1) * ( 1 - dim_vec(8) ) )
   integer(kind=our_int), intent(in) :: cov_patt_for_data_patt( dim_vec(4) )
   integer(kind=our_int), intent(in) :: response_for_data_patt( dim_vec(4) )
   integer(kind=our_int), intent(in) :: survey_mode_int
   integer(kind=our_int), intent(in) :: baseline_int
   real(kind=our_dble), intent(in) :: pert_mat( dim_vec(2), dim_vec(2) )
   real(kind=our_dble), intent(in) :: pert_mat_inv( dim_vec(2), dim_vec(2) )
   integer(kind=our_int), intent(in) :: prior_int
   real(kind=our_dble), intent(in) :: prior_freq_tot
   integer(kind=our_int), intent(in) :: prior_alloc_supplied_int
   real(kind=our_dble), intent(inout) :: prior_alloc( dim_vec(2) )
   integer(kind=our_int), intent(in) :: saturated_int
   integer(kind=our_int), intent(in) :: method_int
   integer(kind=our_int), intent(in) :: ctrl_int(12)
   real(kind=our_dble), intent(in) :: ctrl_real(8)
   real(kind=our_dble), intent(in) :: special_codes(4)
   ! mcmc inputs
   integer(kind=our_int), intent(in) :: dim_vec_mcmc(2)
   real(kind=our_dble), intent(out) :: vhat_coef_rwm( &
        (1-saturated_int) * dim_vec(6), (1-saturated_int) * dim_vec(6) )
   integer(kind=our_int), intent(in) :: micro_data_int
   integer(kind=our_int), intent(in) :: freq_row_input_data_int( dim_vec(1) )
   integer(kind=our_int), intent(in) :: f_mat_int( dim_vec(1), dim_vec(2) )
   integer(kind=our_int), intent(in) :: &
        response_row_input_data( ( 1 - dim_vec(8) ) * dim_vec(1) )
   integer(kind=our_int), intent(in) :: freq_for_data_patt_int( dim_vec(4) )
   ! starting values indicator
   integer(kind=our_int), intent(in) :: start_val_source_int
   ! outputs
   integer(kind=our_int), intent(out) :: outputs_int(6)
   real(kind=our_dble), intent(out) :: outputs_real(2)
   real(kind=our_dble), intent(out) :: score( dim_vec(6) )
   real(kind=our_dble), intent(out) :: hess( dim_vec(6), dim_vec(6) )
   real(kind=our_dble), intent(inout) :: coefficients( dim_vec(5), dim_vec(2) ) 
   real(kind=our_dble), intent(out) :: coef_vec( dim_vec(6) )
   real(kind=our_dble), intent(out) :: vhat_coef( dim_vec(6), dim_vec(6) )
   real(kind=our_dble), intent(inout) :: fitted_pi( dim_vec(3), dim_vec(2) )
   real(kind=our_dble), intent(out) :: fitted_pistar( dim_vec(3), dim_vec(2) )
   real(kind=our_dble), intent(out) :: marg_proportions_ystar( dim_vec(2) )
   real(kind=our_dble), intent(out) :: marg_proportions_y( dim_vec(2) )
   real(kind=our_dble), intent(out) :: fhat_mat( dim_vec(3), dim_vec(2) )
   real(kind=our_dble), intent(out) :: hessA( dim_vec(6), dim_vec(6) )
   real(kind=our_dble), intent(out) :: fstar_mat( dim_vec(3), dim_vec(2) )
   real(kind=our_dble), intent(out) :: fitweight_mat( dim_vec(3), dim_vec(2) )
   ! mcmc outputs
   real(kind=our_dble), intent(out) :: coef_mcmc( &
        (1-saturated_int) * dim_vec(5), (1-saturated_int) * dim_vec(2) )
   real(kind=our_dble), intent(out) :: coef_vec_mcmc( &
        (1-saturated_int) * dim_vec(6) )
   real(kind=our_dble), intent(out) :: vhat_coef_mcmc( &
        (1-saturated_int) * dim_vec(6), (1-saturated_int) * dim_vec(6) )
   real(kind=our_dble), intent(out) :: fitted_pi_mcmc( &
        dim_vec(3), dim_vec(2) )
   real(kind=our_dble), intent(out) :: fitted_pistar_mcmc( &
        dim_vec(3), dim_vec(2) )
   integer(kind=our_int), intent(out) :: n_actual(3)
   real(kind=our_dble), intent(out) :: mh_accept_rate
   real(kind=our_dble), intent(out) :: start_logP
   real(kind=our_dble), intent(out) :: logP_series( dim_vec_mcmc(1) )
   real(kind=our_dble), intent(out) :: coef_vec_series( dim_vec_mcmc(1), &
        (1-saturated_int) * dim_vec(6) )
   integer(kind=our_int), intent(out) :: &
        imp_mat_series( dim_vec(1), dim_vec(2), dim_vec_mcmc(2) )
   integer(kind=our_int), intent(out) :: &
        imp_vec_series( dim_vec(1), dim_vec_mcmc(2) )
   real(kind=our_dble), intent(out) :: pi_marg_mcmc( dim_vec(2) )
   real(kind=our_dble), intent(out) :: &
        pi_marg_series( dim_vec_mcmc(1), dim_vec(2) )
   real(kind=our_dble), intent(out) :: &
        approx_bayes_log_imp_ratios( dim_vec_mcmc(1) )
   ! messaging outputs
   integer(kind=our_int), intent(out) :: status
   integer(kind=our_int), intent(in) :: msg_len_max
   integer(kind=our_int), intent(out) :: msg_codes( msg_len_max, 17 )
   integer(kind=our_int), intent(out) :: msg_len_actual
   ! declare locals
   integer(kind=our_int) :: ijunk
   type(workspace_type_rrlogit) :: work
   type(error_type) :: err
   !
   status = 1
   call err_reset(err)
   if( get_randgen_state_R(err) == RETURN_FAIL ) goto 800
   if( run_rrlogit(  &
        dim_vec(1), &   ! nrow_input_data
        dim_vec(2), &   ! n_levels
        dim_vec(3), &   ! n_cov_patt
        dim_vec(4), &   ! n_data_patt
        dim_vec(5), &   ! ncol_model_matrix
        dim_vec(6), &   ! nparam_this_model
        dim_vec(7), &   ! nparam_sat_model
        dim_vec(8), &   ! wide_format_int
        model_matrix, fitweight_row_input_data, fitweight_cov_patt, &
        fitweight_data_patt, freq_for_data_patt, cov_patt, data_patt, &
        cov_patt_for_data_patt, response_for_data_patt, &
        survey_mode_int, baseline_int, pert_mat, pert_mat_inv, &
        prior_int, prior_freq_tot, prior_alloc_supplied_int, prior_alloc, &
        saturated_int, method_int, &
        ctrl_int(1), &  ! iter_max_nr
        ctrl_int(2), &  ! iter_max_fs
        ctrl_int(3), &  ! iter_max_em
        ctrl_int(4), &  ! iter_max_mstep
        ctrl_int(5), &  ! iter_approx_bayes
        ctrl_int(6), &  ! impute_approx_bayes_int
        ctrl_int(7), &  ! iter_mcmc
        ctrl_int(8), &  ! burn_mcmc
        ctrl_int(9), &  ! thin_mcmc
        ctrl_int(10), & ! impute_every
        ctrl_int(11), & ! type_mcmc_int
        ctrl_int(12), & ! stuck_limit
        ctrl_real(1), & ! crit_converge
        ctrl_real(2), & ! crit_boundary
        ctrl_real(3), & ! df_da
        ctrl_real(4), & ! step_size_da
        ctrl_real(5), & ! scale_fac_da
        ctrl_real(6), & ! df_rwm
        ctrl_real(7), & ! scale_fac_rwm
        ctrl_real(8), & ! start_val_jitter
        special_codes(1), & ! mvcode
        special_codes(2), & ! nancode
        special_codes(3), & ! infcode
        special_codes(4), & ! neginfcode
        dim_vec_mcmc(1), &  ! series_length
        dim_vec_mcmc(2), &  ! n_impute
        vhat_coef_rwm, &
        micro_data_int, freq_row_input_data_int, &
        f_mat_int, response_row_input_data, freq_for_data_patt_int, &
        start_val_source_int, &
        outputs_int(1), & ! iter
        outputs_int(2), & ! converged_int
        outputs_int(3), & ! boundary_int
        outputs_int(4), & ! aborted_int
        outputs_int(5), & ! vhat_failed_int
        outputs_int(6), & ! dap_failed_int
        outputs_real(1), & ! loglik
        outputs_real(2), & ! logprior
        score, hess, &
        coefficients, coef_vec, vhat_coef, fitted_pi, fitted_pistar, &
        marg_proportions_ystar, marg_proportions_y, fhat_mat, hessA, &
        fstar_mat, fitweight_mat, &
        coef_mcmc, coef_vec_mcmc, vhat_coef_mcmc, fitted_pi_mcmc, &
        fitted_pistar_mcmc, &
        n_actual(1), & ! n_iter_actual
        n_actual(2), & ! n_sample_actual
        n_actual(3), & ! n_imp_actual
        mh_accept_rate, start_logP, logP_series, coef_vec_series, &
        imp_mat_series, imp_vec_series, pi_marg_mcmc, pi_marg_series, &
        approx_bayes_log_imp_ratios, &
        work, err ) == RETURN_FAIL ) goto 800
   ! normal exit
   status = 0
800 continue
   ! report message if present
   msg_codes(:,:) = 0
   msg_len_actual = 0
   if( err_msg_present(err) ) call err_get_codes(err, &
        msg_codes, msg_len_actual)
   ! cleanup
   call err_reset(err)
   ijunk = nullify_workspace_type_rrlogit(work, err)
   ijunk = put_randgen_state_R(err)
end subroutine rrlogit
!#####################################################################
subroutine rrlogit_predict( dim_vec, &
     model_matrix, fitweight_row_input_data, fitweight_cov_patt, &
     fitweight_data_patt, freq_for_data_patt, cov_patt, data_patt, &
     cov_patt_for_data_patt, response_for_data_patt, &
     survey_mode_int, baseline_int, pert_mat, pert_mat_inv, &
     prior_int, prior_freq_tot, prior_alloc_supplied_int, prior_alloc, &
     saturated_int, method_int, ctrl_int, ctrl_real, special_codes, &
     type_int, noisy_int, se_fit_int, freq_row_input_data, &
     coefficients, vhat_coef, &
     ! args on next line are inout
     fitted_pi, fitted_pistar, &
     ! outputs
     fitted, se_mat, vhat_fitted_array, &
     status, msg_len_max, msg_codes, msg_len_actual )
   !#############################################################
   ! This is a wrapper function for run_rrlogit_predict
   ! status = 0 means everything ran OK, or run was aborted for
   !    some reason (see msg)
   ! Other value of status indicates a fatal error.
   !#############################################################
   use error_handler
   use program_constants
   use dynalloc
   use math_R
   use rrlogit_engine
   implicit none
   ! declare input arguments
   integer(kind=our_int), intent(in) :: dim_vec(8)
   real(kind=our_dble), intent(in) :: model_matrix( dim_vec(3), dim_vec(5) )
   real(kind=our_dble), intent(in) :: fitweight_row_input_data( dim_vec(1) )
   real(kind=our_dble), intent(in) :: fitweight_cov_patt( dim_vec(3) )
   real(kind=our_dble), intent(in) :: fitweight_data_patt( dim_vec(4) )
   real(kind=our_dble), intent(in) :: freq_for_data_patt( dim_vec(4) )
   integer(kind=our_int), intent(in) :: cov_patt( dim_vec(1) )
   integer(kind=our_int), intent(in) :: &
        data_patt( dim_vec(1) * ( 1 - dim_vec(8) ) )
   integer(kind=our_int), intent(in) :: cov_patt_for_data_patt( dim_vec(4) )
   integer(kind=our_int), intent(in) :: response_for_data_patt( dim_vec(4) )
   integer(kind=our_int), intent(in) :: survey_mode_int
   integer(kind=our_int), intent(in) :: baseline_int
   real(kind=our_dble), intent(in) :: pert_mat( dim_vec(2), dim_vec(2) )
   real(kind=our_dble), intent(in) :: pert_mat_inv( dim_vec(2), dim_vec(2) )
   integer(kind=our_int), intent(in) :: prior_int
   real(kind=our_dble), intent(in) :: prior_freq_tot
   integer(kind=our_int), intent(in) :: prior_alloc_supplied_int
   real(kind=our_dble), intent(inout) :: prior_alloc( dim_vec(2) )
   integer(kind=our_int), intent(in) :: saturated_int
   integer(kind=our_int), intent(in) :: method_int
   integer(kind=our_int), intent(in) :: ctrl_int(12)
   real(kind=our_dble), intent(in) :: ctrl_real(8)
   real(kind=our_dble), intent(in) :: special_codes(4)
   integer(kind=our_int), intent(in) :: type_int
   integer(kind=our_int), intent(in) :: noisy_int
   integer(kind=our_int), intent(in) :: se_fit_int
   real(kind=our_dble), intent(in) :: freq_row_input_data( dim_vec(1) )
   real(kind=our_dble), intent(in) :: coefficients( dim_vec(5), dim_vec(2) ) 
   real(kind=our_dble), intent(in) :: vhat_coef( dim_vec(6), dim_vec(6) )
   ! inouts
   real(kind=our_dble), intent(inout) :: fitted_pi( dim_vec(3), dim_vec(2) )
   real(kind=our_dble), intent(inout) :: fitted_pistar( dim_vec(3), dim_vec(2) )
   ! outputs
   real(kind=our_dble), intent(out) :: fitted( dim_vec(1), dim_vec(2) )
   real(kind=our_dble), intent(out) :: se_mat( dim_vec(1) * se_fit_int, &
        dim_vec(2) * se_fit_int )
   real(kind=our_dble), intent(out) :: &
        vhat_fitted_array( dim_vec(1) * se_fit_int, &
        dim_vec(2) * se_fit_int, dim_vec(2) * se_fit_int )
   ! messaging outputs
   integer(kind=our_int), intent(out) :: status
   integer(kind=our_int), intent(in) :: msg_len_max
   integer(kind=our_int), intent(out) :: msg_codes( msg_len_max, 17 )
   integer(kind=our_int), intent(out) :: msg_len_actual
   ! declare locals
   integer(kind=our_int) :: ijunk
   type(workspace_type_rrlogit) :: work
   type(error_type) :: err
   ! begin
   status = 1
   call err_reset(err)
   if( get_randgen_state_R(err) == RETURN_FAIL ) goto 800
   if( run_rrlogit_predict(  &
        dim_vec(1), &   ! nrow_input_data
        dim_vec(2), &   ! n_levels
        dim_vec(3), &   ! n_cov_patt
        dim_vec(4), &   ! n_data_patt
        dim_vec(5), &   ! ncol_model_matrix
        dim_vec(6), &   ! nparam_this_model
        dim_vec(7), &   ! nparam_sat_model
        dim_vec(8), &   ! wide_format_int
        model_matrix, fitweight_row_input_data, fitweight_cov_patt, &
        fitweight_data_patt, freq_for_data_patt, cov_patt, data_patt, &
        cov_patt_for_data_patt, response_for_data_patt, &
        survey_mode_int, baseline_int, pert_mat, pert_mat_inv, &
        prior_int, prior_freq_tot, prior_alloc_supplied_int, prior_alloc, &
        saturated_int, method_int, &
        ctrl_int(1), &  ! iter_max_nr
        ctrl_int(2), &  ! iter_max_fs
        ctrl_int(3), &  ! iter_max_em
        ctrl_int(4), &  ! iter_max_mstep
        ctrl_real(1), & ! crit_converge
        ctrl_real(2), & ! crit_boundary
        special_codes(1), & ! mvcode
        special_codes(2), & ! nancode
        special_codes(3), & ! infcode
        special_codes(4), & ! neginfcode
        type_int, noisy_int, se_fit_int, freq_row_input_data, &
        coefficients, vhat_coef, &
        fitted_pi, fitted_pistar, &
        fitted, se_mat, vhat_fitted_array, &
        work, err ) == RETURN_FAIL ) goto 800
   ! normal exit
   status = 0
800 continue
   ! report message if present
   msg_codes(:,:) = 0
   msg_len_actual = 0
   if( err_msg_present(err) ) call err_get_codes(err, &
        msg_codes, msg_len_actual)
   ! cleanup
   call err_reset(err)
   ijunk = nullify_workspace_type_rrlogit(work, err)
   ijunk = put_randgen_state_R(err)
 end subroutine rrlogit_predict
!#####################################################################
subroutine rrlogit_impute( dim_vec, &
     model_matrix, fitweight_row_input_data, fitweight_cov_patt, &
     fitweight_data_patt, freq_for_data_patt, cov_patt, data_patt, &
     cov_patt_for_data_patt, response_for_data_patt, &
     survey_mode_int, baseline_int, pert_mat, pert_mat_inv, &
     prior_int, prior_freq_tot, prior_alloc_supplied_int, prior_alloc, &
     saturated_int, method_int, ctrl_int, ctrl_real, special_codes, &
     type_int, cond_mean_int, micro_data_int, &
     freq_row_input_data, freq_row_input_data_int, f_mat, f_mat_int, &
     response_row_input_data, &
     coefficients, &
     ! args on next line are inout
     fitted_pi, fitted_pistar, &
     ! outputs
     cond_means, imp_mat, imp_vec, &
     status, msg_len_max, msg_codes, msg_len_actual )
   !#############################################################
   ! This is a wrapper function for run_rrlogit_predict
   ! status = 0 means everything ran OK, or run was aborted for
   !    some reason (see msg)
   ! Other value of status indicates a fatal error.
   !#############################################################
   use error_handler
   use program_constants
   use dynalloc
   use math_R
   use rrlogit_engine
   implicit none
   ! declare input arguments
   integer(kind=our_int), intent(in) :: dim_vec(8)
   real(kind=our_dble), intent(in) :: model_matrix( dim_vec(3), dim_vec(5) )
   real(kind=our_dble), intent(in) :: fitweight_row_input_data( dim_vec(1) )
   real(kind=our_dble), intent(in) :: fitweight_cov_patt( dim_vec(3) )
   real(kind=our_dble), intent(in) :: fitweight_data_patt( dim_vec(4) )
   real(kind=our_dble), intent(in) :: freq_for_data_patt( dim_vec(4) )
   integer(kind=our_int), intent(in) :: cov_patt( dim_vec(1) )
   integer(kind=our_int), intent(in) :: &
        data_patt( dim_vec(1) * ( 1 - dim_vec(8) ) )
   integer(kind=our_int), intent(in) :: cov_patt_for_data_patt( dim_vec(4) )
   integer(kind=our_int), intent(in) :: response_for_data_patt( dim_vec(4) )
   integer(kind=our_int), intent(in) :: survey_mode_int
   integer(kind=our_int), intent(in) :: baseline_int
   real(kind=our_dble), intent(in) :: pert_mat( dim_vec(2), dim_vec(2) )
   real(kind=our_dble), intent(in) :: pert_mat_inv( dim_vec(2), dim_vec(2) )
   integer(kind=our_int), intent(in) :: prior_int
   real(kind=our_dble), intent(in) :: prior_freq_tot
   integer(kind=our_int), intent(in) :: prior_alloc_supplied_int
   real(kind=our_dble), intent(inout) :: prior_alloc( dim_vec(2) )
   integer(kind=our_int), intent(in) :: saturated_int
   integer(kind=our_int), intent(in) :: method_int
   integer(kind=our_int), intent(in) :: ctrl_int(12)
   real(kind=our_dble), intent(in) :: ctrl_real(8)
   real(kind=our_dble), intent(in) :: special_codes(4)
   integer(kind=our_int), intent(in) :: type_int
   integer(kind=our_int), intent(in) :: cond_mean_int ! for dimensioning
   integer(kind=our_int), intent(in) :: micro_data_int
   real(kind=our_dble), intent(in) :: &
        freq_row_input_data( cond_mean_int * dim_vec(1) )
   integer(kind=our_int), intent(in) :: &
        freq_row_input_data_int( ( 1 - cond_mean_int ) * dim_vec(1) )
   real(kind=our_dble), intent(in) :: &
        f_mat( cond_mean_int * dim_vec(1), cond_mean_int * dim_vec(2) )
   integer(kind=our_int), intent(in) :: &
        f_mat_int( ( 1 - cond_mean_int ) * dim_vec(1), &
        ( 1 - cond_mean_int ) * dim_vec(2) )
   integer(kind=our_int), intent(in) :: &
        response_row_input_data( ( 1 - dim_vec(8) ) * dim_vec(1) )
   real(kind=our_dble), intent(in) :: coefficients( dim_vec(5), dim_vec(2) ) 
   ! inouts
   real(kind=our_dble), intent(inout) :: fitted_pi( dim_vec(3), dim_vec(2) )
   real(kind=our_dble), intent(inout) :: fitted_pistar( dim_vec(3), dim_vec(2) )
   ! outputs
   real(kind=our_dble), intent(out) :: &
        cond_means( cond_mean_int * dim_vec(1), cond_mean_int * dim_vec(2) )
   integer(kind=our_int), intent(out) :: &
        imp_mat( ( 1 - cond_mean_int ) * dim_vec(1), &
        ( 1 - cond_mean_int ) * dim_vec(2) )
   integer(kind=our_int), intent(out) :: &
        imp_vec( ( 1 - cond_mean_int ) * dim_vec(1) )
   ! messaging outputs
   integer(kind=our_int), intent(out) :: status
   integer(kind=our_int), intent(in) :: msg_len_max
   integer(kind=our_int), intent(out) :: msg_codes( msg_len_max, 17 )
   integer(kind=our_int), intent(out) :: msg_len_actual
   ! declare locals
   integer(kind=our_int) :: ijunk
   type(workspace_type_rrlogit) :: work
   type(error_type) :: err
   ! begin
   status = 1
   call err_reset(err)
   if( get_randgen_state_R(err) == RETURN_FAIL ) goto 800
   if( run_rrlogit_impute(  &
        dim_vec(1), &   ! nrow_input_data
        dim_vec(2), &   ! n_levels
        dim_vec(3), &   ! n_cov_patt
        dim_vec(4), &   ! n_data_patt
        dim_vec(5), &   ! ncol_model_matrix
        dim_vec(6), &   ! nparam_this_model
        dim_vec(7), &   ! nparam_sat_model
        dim_vec(8), &   ! wide_format_int
        model_matrix, fitweight_row_input_data, fitweight_cov_patt, &
        fitweight_data_patt, freq_for_data_patt, cov_patt, data_patt, &
        cov_patt_for_data_patt, response_for_data_patt, &
        survey_mode_int, baseline_int, pert_mat, pert_mat_inv, &
        prior_int, prior_freq_tot, prior_alloc_supplied_int, prior_alloc, &
        saturated_int, method_int, &
        ctrl_int(1), &  ! iter_max_nr
        ctrl_int(2), &  ! iter_max_fs
        ctrl_int(3), &  ! iter_max_em
        ctrl_int(4), &  ! iter_max_mstep
        ctrl_real(1), & ! crit_converge
        ctrl_real(2), & ! crit_boundary
        special_codes(1), & ! mvcode
        special_codes(2), & ! nancode
        special_codes(3), & ! infcode
        special_codes(4), & ! neginfcode
        type_int, micro_data_int, &
        freq_row_input_data, freq_row_input_data_int, f_mat, f_mat_int, &
        response_row_input_data, &
        coefficients, &
        fitted_pi, fitted_pistar, &
        cond_means, imp_mat, imp_vec, &
        work, err ) == RETURN_FAIL ) goto 800
   ! normal exit
   status = 0
800 continue
   ! report message if present
   msg_codes(:,:) = 0
   msg_len_actual = 0
   if( err_msg_present(err) ) call err_get_codes(err, &
        msg_codes, msg_len_actual)
   ! cleanup
   call err_reset(err)
   ijunk = nullify_workspace_type_rrlogit(work, err)
   ijunk = put_randgen_state_R(err)
 end subroutine rrlogit_impute
!#####################################################################
