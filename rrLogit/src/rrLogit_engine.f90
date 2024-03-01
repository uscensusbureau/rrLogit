!#####################################################################
module rrlogit_engine
   ! data structures and computational routines for rrLogit modeling
   use program_constants
   use error_handler
   use dynalloc
   use quick_sort
   use tabulate
   use matrix_methods
   use math_R
   use math_funcs
   implicit none
   private ! by default
   ! declare public types (contents are still private)
   public :: workspace_type_rrlogit
   ! declare public parameters
   public :: model_type_str_len
   ! declare public procedures
   public :: nullify_workspace_type_rrlogit, run_rrlogit, &
        run_rrlogit_predict, run_rrlogit_impute
   ! parameters private to this module
   real(kind=our_dble), parameter :: &
        log_huge = log( huge( real(0,kind=our_dble) ) ), &
        log_tiny = log( tiny( real(0,kind=our_dble) ) )
   character(len=*), parameter :: modname = "rrlogit_engine"
   integer(kind=our_int), parameter :: model_type_str_len = 30, &
        prior_type_str_len = 30, method_str_len = 30
   ! private data types
   !###################################################################
   type :: workspace_type_rrlogit
      ! for holding data and fitting a model
      integer(kind=our_int) :: nrow_input_data = 0
      integer(kind=our_int) :: n_levels = 0
      integer(kind=our_int) :: n_cov_patt = 0
      integer(kind=our_int) :: n_data_patt = 0
      integer(kind=our_int) :: ncol_model_matrix = 0
      integer(kind=our_int) :: nparam_this_model = 0
      integer(kind=our_int) :: nparam_sat_model = 0
      integer(kind=our_int) :: baseline = 0
      logical :: wide_format = .false.
      logical :: survey_mode = .false.
      logical :: saturated = .false.
      character(len=method_str_len) :: method = ""
      integer(kind=our_int) :: iter_max = 0
      integer(kind=our_int) :: iter_max_nr = 0
      integer(kind=our_int) :: iter_max_fs = 0
      integer(kind=our_int) :: iter_max_em = 0
      integer(kind=our_int) :: iter_max_mstep = 0
      real(kind=our_dble) :: crit_converge = 0.D0
      real(kind=our_dble) :: crit_boundary = 0.D0
      real(kind=our_dble) :: mvcode = 0.D0
      real(kind=our_dble) :: nancode = 0.D0
      real(kind=our_dble) :: infcode = 0.D0
      real(kind=our_dble) :: neginfcode = 0.D0
      character(len=method_str_len) :: start_val_source = ""
      real(kind=our_dble) :: start_val_jitter = 0.D0
      ! "fitweight" refers to either a frequency or survey weight,
      ! whichever plays the role of frequency in model fitting
      real(kind=our_dble), allocatable :: model_matrix(:,:)
      real(kind=our_dble), allocatable :: fitweight_row_input_data(:)
      real(kind=our_dble), allocatable :: fitweight_cov_patt(:)
      real(kind=our_dble), allocatable :: fitweight_data_patt(:)
      integer(kind=our_int), allocatable :: cov_patt(:)
      integer(kind=our_int), allocatable :: data_patt(:)
      integer(kind=our_int), allocatable :: cov_patt_for_data_patt(:)
      integer(kind=our_int), allocatable :: response_for_data_patt(:)
      real(kind=our_dble), allocatable :: pert_mat(:,:)
      real(kind=our_dble), allocatable :: pert_mat_inv(:,:)
      logical :: pert_mat_identity = .false.
      character(len=prior_type_str_len) :: prior = ""
      real(kind=our_dble) :: prior_freq_tot = 0.D0
      logical :: prior_alloc_supplied = .false.
      real(kind=our_dble), allocatable :: prior_alloc(:)
      logical :: neg_prob = .false.
      integer(kind=our_int), allocatable :: data_patt_st(:)
      integer(kind=our_int), allocatable :: data_patt_fin(:)
      logical, allocatable :: active_ystar(:,:)
      logical, allocatable :: empty_cov_patt(:)
      integer(kind=our_int) :: n_cov_patt_empty = 0
      integer(kind=our_int), allocatable :: cov_patt_order(:)
      integer(kind=our_int), allocatable :: cov_patt_st(:)
      integer(kind=our_int), allocatable :: cov_patt_fin(:)
      ! for convenience, we use shorthand symbols for key dimensions
      integer(kind=our_int) :: n = 0  ! same as nrow_input_data
      integer(kind=our_int) :: r = 0  ! same as n_levels
      integer(kind=our_int) :: p = 0  ! same as ncol_model_matrix
      integer(kind=our_int) :: d = 0  ! same as nparam_this_model
      real(kind=our_dble), allocatable :: beta(:,:), &
           beta_vec(:), beta_vec_old(:), beta_vec_new(:), &
           beta_mstep(:), beta_mstep_old(:), beta_mstep_new(:), &
           pi_mat(:,:), pistar_mat(:,:), score(:), hess(:,:), &
           pi_mat_old(:,:), pistar_mat_old(:,:), &
           pi_mat_mstep(:,:), pi_mat_mstep_old(:,:), &
           phi_mat(:,:), fhat_mat(:,:), &
           fstar_mat(:,:), fitweight_mat(:,:), &
           scoreA(:), hessA(:,:), &
           prior_freq(:,:), prior_freq_tot_cov_patt(:), &
           marg_proportions_ystar(:), marg_proportions_y(:), &
           mini_em_prior_freq(:), mini_em_result(:), pi_marg(:), &
           vhat_coef(:,:), pi_vec(:), pi_vec_old(:), pistar_vec(:), &
           fhat_vec(:), se_fit_mat(:,:), vhat_fit_array(:,:,:)
      real(kind=our_dble) :: mini_em_prior_freq_tot = 0.D0
      integer(kind=our_int), allocatable :: fimp_vec(:), fimp_mat(:,:)
      real(kind=our_dble) :: loglik = 0.D0, logP = 0.D0, &
           loglikA = 0.D0, logprior = 0.D0, logpriorA = 0.D0, &
           logprior_new = 0.D0, loglik_new = 0.D0, logpriorA_new = 0.D0, &
           loglikA_new = 0.D0
      real(kind=our_dble), allocatable :: loglik_vec(:), logP_vec(:)
      integer(kind=our_int) :: iter = 0, iter_mstep = 0, iter_em = 0, &
           step_halving = 0
      integer(kind=our_int), allocatable :: iter_em_saturated(:)
      logical, allocatable :: converged_em_saturated(:), &
           boundary_em_saturated(:)
      logical :: converged = .false., aborted = .false., &
           aborted_estep = .false., aborted_mstep = .false., &
           aborted_em = .false., &
           boundary = .false., converged_mstep = .false., &
           vhat_failed = .false., dap_failed = .false., &
           any_zero_pi = .false., any_zero_pi_mstep = .false.
      real(kind=our_dble), allocatable :: dldpistar(:), dldpi(:), &
           d2ldpi2(:,:), &
           Edldpistar(:), Edldpi(:), &
           dpidbeta(:,:), d2pidbeta2(:,:,:)
      real(kind=our_dble), allocatable :: wkrA(:), wkrB(:), &
           wkrdA(:,:), wkrdB(:,:), wkrdC(:,:), wkddA(:,:), wkddB(:,:), &
           wkdA(:), wkdB(:), wkdC(:), &
           wkppA(:,:), wkppB(:,:), wkpA(:), wkpB(:), &
           wkprA(:,:), wkprB(:,:)
      ! objects below pertain to MCMC and approxBayes
      integer(kind=our_int) :: iter_mcmc_nominal = 0
      integer(kind=our_int) :: iter_mcmc = 0
      integer(kind=our_int) :: burn_mcmc = 0
      integer(kind=our_int) :: thin_mcmc = 0
      integer(kind=our_int) :: impute_every = 0
      character(len=method_str_len) :: type_mcmc = ""
      integer(kind=our_int) :: stuck_limit = 0
      integer(kind=our_int) :: iter_approx_bayes = 0
      logical :: impute_approx_bayes = .false.
      real(kind=our_dble) :: df_da = 0.D0
      real(kind=our_dble) :: step_size_da = 0.D0
      real(kind=our_dble) :: scale_fac_da = 0.D0
      real(kind=our_dble) :: df_rwm = 0.D0
      real(kind=our_dble) :: scale_fac_rwm = 0.D0
      integer(kind=our_int) :: series_length = 0
      integer(kind=our_int) :: n_impute = 0
      integer(kind=our_int) :: beta_accept_count = 0
      integer(kind=our_int) :: beta_current_reject_run = 0
      real(kind=our_dble) :: beta_accept_rate = 0.D0
      real(kind=our_dble) :: beta_vec_log_dens = 0.D0
      real(kind=our_dble), allocatable :: mh_ratios_beta(:)
      logical, allocatable :: mh_accept_beta(:)
      real(kind=our_dble), allocatable :: beta_can(:)  ! candidate
      real(kind=our_dble), allocatable :: beta_center(:)
      real(kind=our_dble), allocatable :: vhat_beta_rwm(:,:)
      real(kind=our_dble), allocatable :: beta_scale(:,:)
      real(kind=our_dble), allocatable :: beta_scale_inv(:,:)
      real(kind=our_dble), allocatable :: beta_scale_inv_sqrt(:,:)
      real(kind=our_dble), allocatable :: beta_scale_sqrt(:,:)
      real(kind=our_dble), allocatable :: beta_mean(:)
      real(kind=our_dble), allocatable :: beta_cov_mat(:,:)
      real(kind=our_dble), allocatable :: pi_mat_mean(:,:)
      real(kind=our_dble), allocatable :: pistar_mat_mean(:,:)
      real(kind=our_dble), allocatable :: pi_marg_mean(:)
      integer(kind=our_int) :: iter_past_burn_in = 0
      integer(kind=our_int) :: store_count = 0
      integer(kind=our_int) :: imp_count = 0
      logical :: store_this_iter = .false.
      logical :: imp_this_iter = .false.
      logical :: micro_data = .false.
      integer(kind=our_int), allocatable :: freq_row_input_data_int(:)
      integer(kind=our_int), allocatable :: f_mat_int(:,:)
      integer(kind=our_int), allocatable :: response_row_input_data(:)
      integer(kind=our_int), allocatable :: freq_for_data_patt_int(:)
   end type workspace_type_rrlogit
   !###################################################################
   contains
   !###################################################################
   integer(kind=our_int) function nullify_workspace_type_rrlogit( &
        work, err ) result( answer )
      implicit none
      ! args
      type(workspace_type_rrlogit), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! locals
      character(len=*), parameter :: &
           subname = "nullify_workspace_type_rrlogit"
      integer(kind=our_int) :: status
      ! begin
      answer = RETURN_FAIL
      work%nrow_input_data = 0
      work%n_levels = 0
      work%n_cov_patt = 0
      work%n_data_patt = 0
      work%ncol_model_matrix = 0
      work%nparam_this_model = 0
      work%nparam_sat_model = 0
      work%baseline = 0
      work%wide_format = .false.
      work%survey_mode = .false.
      work%saturated = .false.
      work%method = ""
      work%iter_max = 0
      work%iter_max_nr = 0
      work%iter_max_fs = 0
      work%iter_max_em = 0
      work%iter_max_mstep= 0
      work%crit_converge = 0.D0
      work%crit_boundary = 0.D0
      work%mvcode = 0.D0
      work%nancode = 0.D0
      work%infcode = 0.D0
      work%neginfcode = 0.D0
      work%start_val_source = ""
      work%start_val_jitter = 0.D0
      if( allocated( work%model_matrix ) ) then
         deallocate( work%model_matrix, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%fitweight_row_input_data ) ) then
         deallocate( work%fitweight_row_input_data, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%fitweight_cov_patt ) ) then
         deallocate( work%fitweight_cov_patt, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%fitweight_data_patt ) ) then
         deallocate( work%fitweight_data_patt, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%cov_patt ) ) then
         deallocate( work%cov_patt, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%data_patt ) ) then
         deallocate( work%data_patt, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%cov_patt_for_data_patt ) ) then
         deallocate( work%cov_patt_for_data_patt, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%response_for_data_patt ) ) then
         deallocate( work%response_for_data_patt, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%pert_mat ) ) then
         deallocate( work%pert_mat, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%pert_mat_inv ) ) then
         deallocate( work%pert_mat_inv, stat=status )
         if( status /= 0 ) goto 800
      end if
      work%pert_mat_identity = .false.
      work%prior = ""
      work%prior_freq_tot = 0.D0
      work%prior_alloc_supplied = .false.
      if( allocated( work%prior_alloc ) ) then
         deallocate( work%prior_alloc, stat=status )
         if( status /= 0 ) goto 800
      end if
      work%neg_prob = .false.
      if( allocated( work%data_patt_st ) ) then
         deallocate( work%data_patt_st, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%data_patt_fin ) ) then
         deallocate( work%data_patt_fin, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%active_ystar ) ) then
         deallocate( work%active_ystar, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%empty_cov_patt ) ) then
         deallocate( work%empty_cov_patt, stat=status )
         if( status /= 0 ) goto 800
      end if
      work%n_cov_patt_empty = 0
      if( allocated( work%cov_patt_order ) ) then
         deallocate( work%cov_patt_order, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%cov_patt_st ) ) then
         deallocate( work%cov_patt_st, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%cov_patt_fin ) ) then
         deallocate( work%cov_patt_fin, stat=status )
         if( status /= 0 ) goto 800
      end if
      !
      work%n = 0
      work%r = 0
      work%p = 0
      work%d = 0
      if( allocated( work%beta ) ) then
         deallocate( work%beta, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%beta_vec ) ) then
         deallocate( work%beta_vec, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%beta_vec_old ) ) then
         deallocate( work%beta_vec_old, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%beta_vec_new ) ) then
         deallocate( work%beta_vec_new, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%beta_mstep ) ) then
         deallocate( work%beta_mstep, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%beta_mstep_old ) ) then
         deallocate( work%beta_mstep_old, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%beta_mstep_new ) ) then
         deallocate( work%beta_mstep_new, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%pi_mat ) ) then
         deallocate( work%pi_mat, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%pistar_mat ) ) then
         deallocate( work%pistar_mat, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%score ) ) then
         deallocate( work%score, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%hess ) ) then
         deallocate( work%hess, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%pi_mat_old ) ) then
         deallocate( work%pi_mat_old, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%pistar_mat_old ) ) then
         deallocate( work%pistar_mat_old, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%pi_mat_mstep ) ) then
         deallocate( work%pi_mat_mstep, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%pi_mat_mstep_old ) ) then
         deallocate( work%pi_mat_mstep_old, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%phi_mat ) ) then
         deallocate( work%phi_mat, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%fhat_mat ) ) then
         deallocate( work%fhat_mat, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%fstar_mat ) ) then
         deallocate( work%fstar_mat, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%fitweight_mat ) ) then
         deallocate( work%fitweight_mat, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%scoreA ) ) then
         deallocate( work%scoreA, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%hessA ) ) then
         deallocate( work%hessA, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%prior_freq ) ) then
         deallocate( work%prior_freq, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%prior_freq_tot_cov_patt ) ) then
         deallocate( work%prior_freq_tot_cov_patt, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%marg_proportions_ystar ) ) then
         deallocate( work%marg_proportions_ystar, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%marg_proportions_y ) ) then
         deallocate( work%marg_proportions_y, stat=status )
         if( status /= 0 ) goto 800
      end if
      work%mini_em_prior_freq_tot = 0.D0
      if( allocated( work%mini_em_prior_freq ) ) then
         deallocate( work%mini_em_prior_freq, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%mini_em_result ) ) then
         deallocate( work%mini_em_result, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%pi_marg ) ) then
         deallocate( work%pi_marg, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%vhat_coef ) ) then
         deallocate( work%vhat_coef, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%pi_vec ) ) then
         deallocate( work%pi_vec, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%pi_vec_old ) ) then
         deallocate( work%pi_vec_old, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%pistar_vec ) ) then
         deallocate( work%pistar_vec, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%fhat_vec ) ) then
         deallocate( work%fhat_vec, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%se_fit_mat ) ) then
         deallocate( work%se_fit_mat, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%vhat_fit_array ) ) then
         deallocate( work%vhat_fit_array, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%fimp_vec ) ) then
         deallocate( work%fimp_vec, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%fimp_mat ) ) then
         deallocate( work%fimp_mat, stat=status )
         if( status /= 0 ) goto 800
      end if
      work%loglik = 0.D0
      work%logP = 0.D0
      work%loglikA = 0.D0
      work%logprior = 0.D0
      work%logpriorA = 0.D0
      work%logprior_new = 0.D0
      work%loglik_new = 0.D0
      work%logpriorA_new = 0.D0
      work%loglikA_new = 0.D0
      if( allocated( work%loglik_vec ) ) then
         deallocate( work%loglik_vec, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%logP_vec ) ) then
         deallocate( work%logP_vec, stat=status )
         if( status /= 0 ) goto 800
      end if
      work%iter = 0
      work%step_halving = 0
      work%iter_mstep = 0
      work%iter_em = 0
      if( allocated( work%iter_em_saturated ) ) then
         deallocate( work%iter_em_saturated, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%converged_em_saturated ) ) then
         deallocate( work%converged_em_saturated, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%boundary_em_saturated ) ) then
         deallocate( work%boundary_em_saturated, stat=status )
         if( status /= 0 ) goto 800
      end if
      work%converged = .false.
      work%converged_mstep = .false.
      work%aborted = .false.
      work%aborted_estep = .false.
      work%aborted_mstep = .false.
      work%aborted_em = .false.
      work%boundary = .false.
      work%vhat_failed = .false.
      work%dap_failed = .false.
      work%any_zero_pi = .false.
      work%any_zero_pi_mstep = .false.
      if( allocated( work%dldpistar ) ) then
         deallocate( work%dldpistar, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%Edldpistar ) ) then
         deallocate( work%Edldpistar, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%dldpi ) ) then
         deallocate( work%dldpi, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%d2ldpi2 ) ) then
         deallocate( work%d2ldpi2, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%Edldpi ) ) then
         deallocate( work%Edldpi, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%dpidbeta ) ) then
         deallocate( work%dpidbeta, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%d2pidbeta2 ) ) then
         deallocate( work%d2pidbeta2, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%wkrA ) ) then
         deallocate( work%wkrA, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%wkrB ) ) then
         deallocate( work%wkrB, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%wkrdA ) ) then
         deallocate( work%wkrdA, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%wkrdB ) ) then
         deallocate( work%wkrdB, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%wkrdC ) ) then
         deallocate( work%wkrdC, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%wkddA ) ) then
         deallocate( work%wkddA, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%wkddB ) ) then
         deallocate( work%wkddB, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%wkdA ) ) then
         deallocate( work%wkdA, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%wkdB ) ) then
         deallocate( work%wkdB, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%wkdC ) ) then
         deallocate( work%wkdC, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%wkppA ) ) then
         deallocate( work%wkppA, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%wkppB ) ) then
         deallocate( work%wkppB, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%wkpA ) ) then
         deallocate( work%wkpA, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%wkpB ) ) then
         deallocate( work%wkpB, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%wkprA ) ) then
         deallocate( work%wkprA, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%wkprB ) ) then
         deallocate( work%wkprB, stat=status )
         if( status /= 0 ) goto 800
      end if
      !
      work%iter_mcmc_nominal = 0
      work%iter_mcmc = 0
      work%burn_mcmc = 0
      work%thin_mcmc = 0
      work%impute_every = 0
      work%type_mcmc = ""
      work%stuck_limit = 0
      work%iter_approx_bayes = 0
      work%impute_approx_bayes = .false.
      work%df_da = 0.D0
      work%step_size_da = 0.D0
      work%scale_fac_da = 0.D0
      work%df_rwm = 0.D0
      work%scale_fac_rwm = 0.D0
      !
      work%series_length = 0
      work%n_impute = 0
      work%beta_accept_count = 0
      work%beta_current_reject_run = 0
      work%beta_accept_rate = 0.D0
      work%beta_vec_log_dens = 0.D0
      if( allocated( work%mh_ratios_beta ) ) then
         deallocate( work%mh_ratios_beta, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%mh_accept_beta ) ) then
         deallocate( work%mh_accept_beta, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%beta_can ) ) then
         deallocate( work%beta_can, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%beta_center ) ) then
         deallocate( work%beta_center, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%vhat_beta_rwm ) ) then
         deallocate( work%vhat_beta_rwm, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%beta_scale ) ) then
         deallocate( work%beta_scale, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%beta_scale_inv ) ) then
         deallocate( work%beta_scale_inv, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%beta_scale_inv_sqrt ) ) then
         deallocate( work%beta_scale_inv_sqrt, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%beta_scale_sqrt ) ) then
         deallocate( work%beta_scale_sqrt, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%beta_mean ) ) then
         deallocate( work%beta_mean, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%beta_cov_mat ) ) then
         deallocate( work%beta_cov_mat, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%pi_mat_mean ) ) then
         deallocate( work%pi_mat_mean, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%pistar_mat_mean ) ) then
         deallocate( work%pistar_mat_mean, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%pi_marg_mean ) ) then
         deallocate( work%pi_marg_mean, stat=status )
         if( status /= 0 ) goto 800
      end if
      work%iter_past_burn_in = 0
      work%store_count = 0
      work%imp_count = 0
      work%store_this_iter = .false.
      work%imp_this_iter = .false.
      work%micro_data = .false.
      if( allocated( work%freq_row_input_data_int ) ) then
         deallocate( work%freq_row_input_data_int, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%f_mat_int ) ) then
         deallocate( work%f_mat_int, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%response_row_input_data ) ) then
         deallocate( work%response_row_input_data, stat=status )
         if( status /= 0 ) goto 800
      end if
      if( allocated( work%freq_for_data_patt_int ) ) then
         deallocate( work%freq_for_data_patt_int, stat=status )
         if( status /= 0 ) goto 800
      end if
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
   end function nullify_workspace_type_rrlogit
   !###################################################################
   integer(kind=our_int) function run_rrlogit( &
        nrow_input_data, n_levels, n_cov_patt, n_data_patt, &
        ncol_model_matrix, nparam_this_model, nparam_sat_model, &
        wide_format_int, &
        model_matrix, fitweight_row_input_data, fitweight_cov_patt, &
        fitweight_data_patt, freq_for_data_patt, cov_patt, data_patt, &
        cov_patt_for_data_patt, response_for_data_patt, &
        survey_mode_int, baseline_int, pert_mat, pert_mat_inv, &
        prior_int, prior_freq_tot, prior_alloc_supplied_int, prior_alloc, &
        saturated_int, method_int, &
        iter_max_nr, iter_max_fs, iter_max_em, iter_max_mstep, &
        iter_approx_bayes, impute_approx_bayes_int, &
        iter_mcmc, burn_mcmc, thin_mcmc, &
        impute_every, type_mcmc_int, stuck_limit, &
        crit_converge, crit_boundary, &
        df_da, step_size_da, scale_fac_da, df_rwm, scale_fac_rwm, &
        start_val_jitter, &
        mvcode, nancode, infcode, neginfcode, &
        series_length, n_impute, vhat_coef_rwm, &
        micro_data_int, freq_row_input_data_int, &
        f_mat_int, response_row_input_data, freq_for_data_patt_int, &
        start_val_source_int, &
        iter, converged_int, boundary_int, aborted_int, &
        vhat_failed_int, dap_failed_int, &
        loglik, logprior, score, hess, &
        coefficients, coef_vec, vhat_coef, fitted_pi, fitted_pistar, &
        marg_proportions_ystar, marg_proportions_y, fhat_mat, hessA, &
        fstar_mat, fitweight_mat, &
        coef_mcmc, coef_vec_mcmc, vhat_coef_mcmc, fitted_pi_mcmc, &
        fitted_pistar_mcmc, n_iter_actual, n_sample_actual, n_imp_actual, &
        mh_accept_rate, start_logP, logP_series, coef_vec_series, &
        imp_mat_series, imp_vec_series, pi_marg_mcmc, pi_marg_series, &
        approx_bayes_log_imp_ratios, &
        work, err ) result(answer)
      ! fits randomized response logistic model using Newton-Raphson,
      ! Fisher scoring or EM
      implicit none
      ! inputs
      integer(kind=our_int), intent(in) :: nrow_input_data
      integer(kind=our_int), intent(in) :: n_levels
      integer(kind=our_int), intent(in) :: n_cov_patt
      integer(kind=our_int), intent(in) :: n_data_patt
      integer(kind=our_int), intent(in) :: ncol_model_matrix
      integer(kind=our_int), intent(in) :: nparam_this_model
      integer(kind=our_int), intent(in) :: nparam_sat_model
      integer(kind=our_int), intent(in) :: wide_format_int
      real(kind=our_dble), intent(in) :: model_matrix(:,:)
      real(kind=our_dble), intent(in) :: fitweight_row_input_data(:)
      real(kind=our_dble), intent(in) :: fitweight_cov_patt(:)
      real(kind=our_dble), intent(in) :: fitweight_data_patt(:)
      real(kind=our_dble), intent(in) :: freq_for_data_patt(:)
      integer(kind=our_int), intent(in) :: cov_patt(:)
      integer(kind=our_int), intent(in) :: data_patt(:)
      integer(kind=our_int), intent(in) :: cov_patt_for_data_patt(:)
      integer(kind=our_int), intent(in) :: response_for_data_patt(:)
      integer(kind=our_int), intent(in) :: survey_mode_int
      integer(kind=our_int), intent(in) :: baseline_int
      real(kind=our_dble), intent(in) :: pert_mat(:,:)
      real(kind=our_dble), intent(in) :: pert_mat_inv(:,:)
      integer(kind=our_int), intent(in) :: prior_int
      real(kind=our_dble), intent(in) :: prior_freq_tot
      integer(kind=our_int), intent(in) :: prior_alloc_supplied_int
      real(kind=our_dble), intent(inout) :: prior_alloc(:)
      integer(kind=our_int), intent(in) :: saturated_int
      integer(kind=our_int), intent(in) :: method_int
      integer(kind=our_int), intent(in) :: iter_max_nr
      integer(kind=our_int), intent(in) :: iter_max_fs
      integer(kind=our_int), intent(in) :: iter_max_em
      integer(kind=our_int), intent(in) :: iter_max_mstep
      integer(kind=our_int), intent(in) :: iter_approx_bayes
      integer(kind=our_int), intent(in) :: impute_approx_bayes_int
      integer(kind=our_int), intent(in) :: iter_mcmc
      integer(kind=our_int), intent(in) :: burn_mcmc
      integer(kind=our_int), intent(in) :: thin_mcmc
      integer(kind=our_int), intent(in) :: impute_every
      integer(kind=our_int), intent(in) :: type_mcmc_int
      integer(kind=our_int), intent(in) :: stuck_limit
      real(kind=our_dble), intent(in) :: crit_converge
      real(kind=our_dble), intent(in) :: crit_boundary
      real(kind=our_dble), intent(in) :: df_da
      real(kind=our_dble), intent(in) :: step_size_da
      real(kind=our_dble), intent(in) :: scale_fac_da
      real(kind=our_dble), intent(in) :: df_rwm
      real(kind=our_dble), intent(in) :: scale_fac_rwm
      real(kind=our_dble), intent(in) :: start_val_jitter
      real(kind=our_dble), intent(in) :: mvcode
      real(kind=our_dble), intent(in) :: nancode
      real(kind=our_dble), intent(in) :: infcode
      real(kind=our_dble), intent(in) :: neginfcode
      ! mcmc inputs
      integer(kind=our_int), intent(in) :: series_length
      integer(kind=our_int), intent(in) :: n_impute
      real(kind=our_dble), intent(in) :: vhat_coef_rwm(:,:)
      integer(kind=our_int), intent(in) :: micro_data_int
      integer(kind=our_int), intent(in) :: freq_row_input_data_int(:)
      integer(kind=our_int), intent(in) :: f_mat_int(:,:)
      integer(kind=our_int), intent(in) :: response_row_input_data(:)
      integer(kind=our_int), intent(in) :: freq_for_data_patt_int(:)
      ! starting values indicator
      integer(kind=our_int), intent(in) :: start_val_source_int
      ! outputs
      integer(kind=our_int), intent(out) :: iter
      integer(kind=our_int), intent(out) :: converged_int
      integer(kind=our_int), intent(out) :: boundary_int
      integer(kind=our_int), intent(out) :: aborted_int
      integer(kind=our_int), intent(out) :: vhat_failed_int
      integer(kind=our_int), intent(out) :: dap_failed_int
      real(kind=our_dble), intent(out) :: loglik
      real(kind=our_dble), intent(out) :: logprior
      real(kind=our_dble), intent(out) :: score(:)
      real(kind=our_dble), intent(out) :: hess(:,:)
      real(kind=our_dble), intent(inout) :: coefficients(:,:)
      real(kind=our_dble), intent(out) :: coef_vec(:)
      real(kind=our_dble), intent(out) :: vhat_coef(:,:)
      real(kind=our_dble), intent(inout) :: fitted_pi(:,:)
      real(kind=our_dble), intent(out) :: fitted_pistar(:,:)
      real(kind=our_dble), intent(out) :: marg_proportions_ystar(:)
      real(kind=our_dble), intent(out) :: marg_proportions_y(:)
      real(kind=our_dble), intent(out) :: fhat_mat(:,:)
      real(kind=our_dble), intent(out) :: hessA(:,:)
      real(kind=our_dble), intent(out) :: fstar_mat(:,:)
      real(kind=our_dble), intent(out) :: fitweight_mat(:,:)
      ! mcmc outputs
      real(kind=our_dble), intent(out) :: coef_mcmc(:,:)
      real(kind=our_dble), intent(out) :: coef_vec_mcmc(:)
      real(kind=our_dble), intent(out) :: vhat_coef_mcmc(:,:)
      real(kind=our_dble), intent(out) :: fitted_pi_mcmc(:,:)
      real(kind=our_dble), intent(out) :: fitted_pistar_mcmc(:,:)
      integer(kind=our_int), intent(out) :: n_iter_actual
      integer(kind=our_int), intent(out) :: n_sample_actual
      integer(kind=our_int), intent(out) :: n_imp_actual
      real(kind=our_dble), intent(out) :: mh_accept_rate
      real(kind=our_dble), intent(out) :: start_logP
      real(kind=our_dble), intent(out) :: logP_series(:)
      real(kind=our_dble), intent(out) :: coef_vec_series(:,:)
      integer(kind=our_int), intent(out) :: imp_mat_series(:,:,:)
      integer(kind=our_int), intent(out) :: imp_vec_series(:,:)
      real(kind=our_dble), intent(out) :: pi_marg_mcmc(:)
      real(kind=our_dble), intent(out) :: pi_marg_series(:,:)
      real(kind=our_dble), intent(out) :: approx_bayes_log_imp_ratios(:)
      ! workspaces
      type(workspace_type_rrlogit), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! locals
      integer(kind=our_int) :: ijunk
      character(len=*), parameter :: &
           subname = "run_rrlogit"
      ! begin
      answer = RETURN_FAIL
      !####
      if( put_data_into_workspace_rrlogit( &
           nrow_input_data, n_levels, n_cov_patt, n_data_patt, &
           ncol_model_matrix, nparam_this_model, nparam_sat_model, &
           wide_format_int, &
           model_matrix, fitweight_row_input_data, fitweight_cov_patt, &
           fitweight_data_patt, freq_for_data_patt, cov_patt, data_patt, &
           cov_patt_for_data_patt, response_for_data_patt, &
           survey_mode_int, baseline_int, pert_mat, pert_mat_inv, &
           prior_int, prior_freq_tot, prior_alloc_supplied_int, prior_alloc, & 
           saturated_int, method_int, &
           iter_max_nr, iter_max_fs, iter_max_em, iter_max_mstep, &
           crit_converge, crit_boundary, &
           mvcode, nancode, infcode, neginfcode, &
           work, err ) == RETURN_FAIL ) goto 800
      ! handle prior
      if( compute_prior_freq_rrlogit( work, err ) == RETURN_FAIL ) goto 800
      if( work%dap_failed ) goto 600
      ! workspace objects for MCMC or approxBayes
      if( ( work%method == "MCMC" ) .or. &
           ( work%method == "approxBayes" ) ) then
         if( create_mcmc_objects_rrlogit( &
              iter_approx_bayes, impute_approx_bayes_int, &
              iter_mcmc, burn_mcmc, thin_mcmc, &
              impute_every, type_mcmc_int, stuck_limit, &
              df_da, step_size_da, scale_fac_da, df_rwm, scale_fac_rwm, &
              series_length, n_impute, vhat_coef_rwm, &
              micro_data_int, freq_row_input_data_int, &
              f_mat_int, response_row_input_data, freq_for_data_patt_int,&
              coef_mcmc, coef_vec_mcmc, vhat_coef_mcmc, fitted_pi_mcmc, &
              fitted_pistar_mcmc, &
              n_iter_actual, n_sample_actual, n_imp_actual, &
              mh_accept_rate, start_logP, logP_series, coef_vec_series, &
              imp_mat_series, imp_vec_series, &
              pi_marg_mcmc, pi_marg_series, approx_bayes_log_imp_ratios, &
              work, err ) == RETURN_FAIL ) goto 800
      end if
      !
      if( compute_start_val_rrlogit( coefficients, fitted_pi, &
           start_val_source_int, start_val_jitter, &
           work, err ) == RETURN_FAIL ) goto 800
      !
      ! estimation procedure
      if( work%saturated ) then
         if( ( work%method == "EM" ) .or. ( work%method == "FS" ) &
              .or. ( work%method == "NR" ) ) then
            if( run_rrlogit_em_saturated( work, err ) == RETURN_FAIL ) goto 800
         else if( work%method == "approxBayes" ) then
            goto 60
         else if( work%method == "MCMC" ) then
            if( run_rrlogit_da_saturated( start_logP, logP_series, &
                 imp_mat_series, imp_vec_series, pi_marg_series, &
                 n_iter_actual, n_sample_actual, n_imp_actual, &
                 work, err ) == RETURN_FAIL ) goto 800
         else
            goto 50
         end if
      else
         if( work%method == "NR" ) then
            if( run_rrlogit_nr( work, err ) == RETURN_FAIL ) goto 800
         else if( work%method == "FS" ) then
            if( run_rrlogit_fs( work, err ) == RETURN_FAIL ) goto 800
         else if( work%method == "EM" ) then
            if( run_rrlogit_em( work, err ) == RETURN_FAIL ) goto 800
         else if( work%method == "MCMC" ) then
            if( work%type_mcmc == "DA" ) then
               if( run_rrlogit_da_nonsaturated( start_logP, logP_series, &
                    coef_vec_series, imp_mat_series, imp_vec_series, &
                    pi_marg_series, &
                    n_iter_actual, n_sample_actual, n_imp_actual, &
                    work, err ) == RETURN_FAIL ) goto 800
            else if( work%type_mcmc == "RWM" ) then
               if( run_rrlogit_rwm_nonsaturated( start_logP, logP_series, &
                    coef_vec_series, imp_mat_series, imp_vec_series, &
                    pi_marg_series, &
                    n_iter_actual, n_sample_actual, n_imp_actual, &
                    work, err ) == RETURN_FAIL ) goto 800
            else
               goto 50
            end if
         else if( work%method == "approxBayes" ) then
            if( run_rrlogit_approx_bayes( start_logP, logP_series, &
                 coef_vec_series, imp_mat_series, imp_vec_series, &
                 pi_marg_series, approx_bayes_log_imp_ratios, &
                 n_iter_actual, n_sample_actual, n_imp_actual, &
                 work, err ) == RETURN_FAIL ) goto 800
         else
            goto 50
         end if
      end if
600   continue
      ! return results
      iter = work%iter
      converged_int = 0
      if( work%converged ) converged_int = 1
      boundary_int = 0
      if( work%boundary ) boundary_int = 1
      aborted_int = 0
      if( work%aborted ) aborted_int = 1
      dap_failed_int = 0
      if( work%dap_failed ) dap_failed_int = 1
      loglik = work%loglik
      logprior = work%logprior
      if( size( score, kind=our_int ) /= work%d ) goto 105
      score(:) = work%score(:)
      if( size( hess, 1, kind=our_int ) /= work%d ) goto 106
      if( size( hess, 2, kind=our_int ) /= work%d ) goto 106
      hess(:,:) = work%hess(:,:)
      if( unstack_coef( work%beta_vec, coefficients, work ) &
           == RETURN_FAIL ) goto 800
      if( size( coef_vec, kind=our_int ) /= work%d ) goto 110
      coef_vec(:) = work%beta_vec(:)
      vhat_failed_int = 0
      if( work%vhat_failed ) vhat_failed_int = 1
      if( size( vhat_coef, 1, kind=our_int ) /= work%d ) goto 120
      if( size( vhat_coef, 2, kind=our_int ) /= work%d ) goto 120
      if( work%vhat_failed ) then
         vhat_coef(:,:) = 0.D0
      else
         vhat_coef(:,:) = work%vhat_coef(:,:)
      end if
      if( size( fitted_pi, 1, kind=our_int ) /= work%n_cov_patt ) goto 130
      if( size( fitted_pi, 2, kind=our_int ) /= work%r ) goto 130
      fitted_pi(:,:) = work%pi_mat(:,:)
      if( size( fitted_pistar, 1, kind=our_int ) /= work%n_cov_patt ) goto 140
      if( size( fitted_pistar, 2, kind=our_int ) /= work%r ) goto 140
      fitted_pistar(:,:) = work%pistar_mat(:,:)
      if( size( marg_proportions_ystar, kind=our_int ) /= work%r ) goto 145
      marg_proportions_ystar(:) = work%marg_proportions_ystar(:)
      if( size( marg_proportions_y, kind=our_int ) /= work%r ) goto 146
      marg_proportions_y(:) = work%marg_proportions_y(:)
      prior_alloc(:) = work%prior_alloc(:)
      if( size( fhat_mat, 1, kind=our_int ) /= work%n_cov_patt ) goto 150
      if( size( fhat_mat, 2, kind=our_int ) /= work%r ) goto 150
      if( estep_rrlogit( work, err, &
           skip_compute_pi_mat = work%saturated ) == RETURN_FAIL ) goto 800
      fhat_mat(:,:) = work%fhat_mat(:,:)
      if( size( hessA, 1, kind=our_int ) /= work%d ) goto 160
      if( size( hessA, 2, kind=our_int ) /= work%d ) goto 160
      hessA(:,:) = work%hessA(:,:)
      if( size( fstar_mat, 1, kind=our_int ) /= work%n_cov_patt ) goto 170
      if( size( fstar_mat, 2, kind=our_int ) /= work%r ) goto 170
      fstar_mat(:,:) = work%fstar_mat(:,:)
      if( size( fitweight_mat, 1, kind=our_int ) /= work%n_cov_patt ) goto 180
      if( size( fitweight_mat, 2, kind=our_int ) /= work%r ) goto 180
      fitweight_mat(:,:) = work%fitweight_mat(:,:)
      ! mcmc outputs not already transferred
      if( ( work%method == "MCMC" ) .or. ( work%method == "approxBayes") ) then
         if( .not. work%saturated ) then
            if( unstack_coef( work%beta_mean, coef_mcmc, work ) &
                 == RETURN_FAIL ) goto 800
            coef_vec_mcmc(:) = work%beta_mean(:)
            vhat_coef_mcmc(:,:) = work%beta_cov_mat(:,:)
            mh_accept_rate = work%beta_accept_rate
         end if
         fitted_pi_mcmc(:,:) = work%pi_mat_mean(:,:)
         fitted_pistar_mcmc(:,:) = work%pistar_mat_mean(:,:)
         pi_marg_mcmc(:) = work%pi_marg_mean(:)
      end if
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
50    call err_handle(err, 1, &
            comment = "Method not recognized" )
      goto 800
60    call err_handle(err, 1, &
            comment = "approxBayes not available for saturated model" )
      goto 800
105   call err_handle(err, 1, &
            comment = "Array score has incorrect size" )
      goto 800
106   call err_handle(err, 1, &
            comment = "Array hess has incorrect size" )
      goto 800
110   call err_handle(err, 1, &
            comment = "Array coef_vec has incorrect size" )
      goto 800
120   call err_handle(err, 1, &
            comment = "Array vhat_coef has incorrect size" )
      goto 800
130   call err_handle(err, 1, &
            comment = "Array fitted_pi has incorrect size" )
      goto 800
140   call err_handle(err, 1, &
            comment = "Array fitted_pistar has incorrect size" )
      goto 800
145   call err_handle(err, 1, &
            comment = "Array marg_proportions_ystar has incorrect size" )
      goto 800
146   call err_handle(err, 1, &
            comment = "Array marg_proportions_y has incorrect size" )
      goto 800
150   call err_handle(err, 1, &
            comment = "Array fhat_mat has incorrect size" )
      goto 800
160   call err_handle(err, 1, &
            comment = "Array hessA has incorrect size" )
      goto 800
170   call err_handle(err, 1, &
            comment = "Array fstar_mat has incorrect size" )
      goto 800
180   call err_handle(err, 1, &
            comment = "Array fitweight_mat has incorrect size" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
      ijunk = nullify_workspace_type_rrlogit( work, err )
   end function run_rrlogit
   !###################################################################
   integer(kind=our_int) function put_data_into_workspace_rrlogit( &
        nrow_input_data, n_levels, n_cov_patt, n_data_patt, &
        ncol_model_matrix, nparam_this_model, nparam_sat_model, &
        wide_format_int, &
        model_matrix, fitweight_row_input_data, fitweight_cov_patt, &
        fitweight_data_patt, freq_for_data_patt, cov_patt, data_patt, &
        cov_patt_for_data_patt, response_for_data_patt, &
        survey_mode_int, baseline_int, pert_mat, pert_mat_inv, &
        prior_int, prior_freq_tot, prior_alloc_supplied_int, prior_alloc, &
        saturated_int, method_int, &
        iter_max_nr, iter_max_fs, iter_max_em, iter_max_mstep, &
        crit_converge, crit_boundary, &
        mvcode, nancode, infcode, neginfcode, &
        work, err ) result(answer)
      implicit none
      ! inputs
      integer(kind=our_int), intent(in) :: nrow_input_data
      integer(kind=our_int), intent(in) :: n_levels
      integer(kind=our_int), intent(in) :: n_cov_patt
      integer(kind=our_int), intent(in) :: n_data_patt
      integer(kind=our_int), intent(in) :: ncol_model_matrix
      integer(kind=our_int), intent(in) :: nparam_this_model
      integer(kind=our_int), intent(in) :: nparam_sat_model
      integer(kind=our_int), intent(in) :: wide_format_int
      real(kind=our_dble), intent(in) :: model_matrix(:,:)
      real(kind=our_dble), intent(in) :: fitweight_row_input_data(:)
      real(kind=our_dble), intent(in) :: fitweight_cov_patt(:)
      real(kind=our_dble), intent(in) :: fitweight_data_patt(:)
      real(kind=our_dble), intent(in) :: freq_for_data_patt(:)
      integer(kind=our_int), intent(in) :: cov_patt(:)
      integer(kind=our_int), intent(in) :: data_patt(:)
      integer(kind=our_int), intent(in) :: cov_patt_for_data_patt(:)
      integer(kind=our_int), intent(in) :: response_for_data_patt(:)
      integer(kind=our_int), intent(in) :: survey_mode_int
      integer(kind=our_int), intent(in) :: baseline_int
      real(kind=our_dble), intent(in) :: pert_mat(:,:)
      real(kind=our_dble), intent(in) :: pert_mat_inv(:,:)
      integer(kind=our_int), intent(in) :: prior_int
      real(kind=our_dble), intent(in) :: prior_freq_tot
      integer(kind=our_int), intent(in) :: prior_alloc_supplied_int
      real(kind=our_dble), intent(in) :: prior_alloc(:)
      integer(kind=our_int), intent(in) :: saturated_int
      integer(kind=our_int), intent(in) :: method_int
      integer(kind=our_int), intent(in) :: iter_max_nr
      integer(kind=our_int), intent(in) :: iter_max_fs
      integer(kind=our_int), intent(in) :: iter_max_em
      integer(kind=our_int), intent(in) :: iter_max_mstep
      real(kind=our_dble), intent(in) :: crit_converge
      real(kind=our_dble), intent(in) :: crit_boundary
      real(kind=our_dble), intent(in) :: mvcode
      real(kind=our_dble), intent(in) :: nancode
      real(kind=our_dble), intent(in) :: infcode
      real(kind=our_dble), intent(in) :: neginfcode
      ! workspaces
      type(workspace_type_rrlogit), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! locals
      integer(kind=our_int) :: status
      integer(kind=our_int) :: i, j, k, cur_cov_patt, posn, cp, dp, ystar
      real(kind=our_dble) :: total_fitweight, sum, rtmp
      real(kind=our_dble), parameter :: tol = 1.D-08
      character(len=*), parameter :: &
           subname = "put_data_into_workspace_rrlogit"
      ! begin
      answer = RETURN_FAIL
      !####
      if( nrow_input_data < 0 ) goto 110
      work%nrow_input_data = nrow_input_data
      if( n_levels < 2 ) goto 120
      work%n_levels = n_levels
      if( n_cov_patt > work%nrow_input_data ) goto 130
      work%n_cov_patt = n_cov_patt
      if( n_data_patt < work%n_cov_patt ) goto 140
      work%n_data_patt = n_data_patt
      if( ncol_model_matrix < 1 ) goto 150
      work%ncol_model_matrix = ncol_model_matrix
      if( nparam_this_model < 1 ) goto 160
      work%nparam_this_model = nparam_this_model
      if( nparam_sat_model < 1 ) goto 170
      work%nparam_sat_model = nparam_sat_model
      if( work%nparam_this_model > work%nparam_sat_model ) goto 180
      if( baseline_int < 0 ) goto 190
      if( baseline_int > work%n_levels ) goto 190
      work%baseline = baseline_int
      if( survey_mode_int == 0 ) then
         work%survey_mode = .false.
      else if( survey_mode_int == 1 ) then
         work%survey_mode = .true.
      else
         goto 200
      end if
      if( wide_format_int == 0 ) then
         work%wide_format = .false.
      else if( wide_format_int == 1 ) then
         work%wide_format = .true.
      else
         goto 210
      end if
      !
      if( saturated_int == 0 ) then
         work%saturated = .false.
      else if( saturated_int == 1 ) then
         work%saturated = .true.
      else
         goto 205
      end if
      !
      if( method_int == 1 ) then
         work%method = "EM"
      else if( method_int == 2 ) then
         work%method = "NR"
      else if( method_int == 3 ) then
         work%method = "FS"
      else if( method_int == 4 ) then
         work%method = "MCMC"
      else if( method_int == 5 ) then
         work%method = "approxBayes"
      else
         goto 211
      end if
      if( iter_max_nr < 0 ) goto 212
      work%iter_max_nr = iter_max_nr
      if( iter_max_fs < 0 ) goto 213
      work%iter_max_fs = iter_max_fs
      if( iter_max_em < 0 ) goto 214
      work%iter_max_em = iter_max_em
      if( iter_max_mstep < 0 ) goto 215
      work%iter_max_mstep = iter_max_mstep
      if( work%method == "EM" ) then
         work%iter_max = work%iter_max_em
      else if( work%method == "NR" ) then
         work%iter_max = work%iter_max_nr
      else if( work%method == "FS" ) then
         work%iter_max = work%iter_max_fs
      else
         work%iter_max = 0 ! fix later
      end if
      if( crit_converge <= 0.D0 ) goto 216
      work%crit_converge = crit_converge
      if( crit_boundary <= 0.D0 ) goto 217
      work%crit_boundary = crit_boundary
      !
      if( size( model_matrix, 1, kind=our_int ) &
           /= work%n_cov_patt ) goto 220
      if( size( model_matrix, 2, kind=our_int ) &
           /= work%ncol_model_matrix ) goto 220
      allocate( work%model_matrix( work%n_cov_patt, work%ncol_model_matrix ), &
           stat=status )
      if( status /= 0 ) goto 100
      work%model_matrix(:,:) = model_matrix(:,:)
      if( size( fitweight_row_input_data, kind=our_int ) &
           /= work%nrow_input_data ) goto 230
      sum = 0.D0
      do i = 1, work%nrow_input_data
         if( fitweight_row_input_data(i) < 0.D0 ) goto 240
         sum = sum + fitweight_row_input_data(i)
      end do
      if( sum <= 0.D0 ) goto 245
      total_fitweight = sum
      allocate( work%fitweight_row_input_data( work%nrow_input_data ), &
           stat=status )
      if( status /= 0 ) goto 100
      work%fitweight_row_input_data(:) = fitweight_row_input_data(:)
      !
      if( size( fitweight_cov_patt, kind=our_int ) &
           /= work%n_cov_patt ) goto 250
      sum = 0.D0
      do i = 1, work%n_cov_patt
         if( fitweight_cov_patt(i) < 0.D0 ) goto 260
         sum = sum + fitweight_cov_patt(i)
      end do
      rtmp = abs( sum - total_fitweight ) / total_fitweight
      if( rtmp > tol ) goto 270
      allocate( work%fitweight_cov_patt( work%n_cov_patt ), &
           stat=status )
      if( status /= 0 ) goto 100
      work%fitweight_cov_patt(:) = fitweight_cov_patt(:)
      !
      if( size( fitweight_data_patt, kind=our_int ) &
           /= work%n_data_patt ) goto 280
      sum = 0.D0
      do i = 1, work%n_data_patt
         if( fitweight_data_patt(i) < 0.D0 ) goto 290
         sum = sum + fitweight_data_patt(i)
      end do
      rtmp = abs( sum - total_fitweight ) / total_fitweight
      if( rtmp > tol ) goto 300
      allocate( work%fitweight_data_patt( work%n_data_patt ), &
           stat=status )
      if( status /= 0 ) goto 100
      work%fitweight_data_patt(:) = fitweight_data_patt(:)
      !
      if( size( freq_for_data_patt, kind=our_int ) &
           /= work%n_data_patt ) goto 305
      !
      if( size( cov_patt, kind=our_int ) /= work%nrow_input_data ) goto 310
      do i = 1, work%nrow_input_data
         if( cov_patt(i) < 1 ) goto 320
         if( cov_patt(i) > work%n_cov_patt ) goto 320
      end do
      allocate( work%cov_patt( work%nrow_input_data ), &
           stat=status )
      if( status /= 0 ) goto 100
      work%cov_patt(:) = cov_patt(:)
      !
      if( work%wide_format ) then
         if( size( data_patt, kind=our_int ) /= 0 ) goto 330
         allocate( work%data_patt(0), stat=status )
         if( status /= 0 ) goto 100
      else
         if( size( data_patt, kind=our_int ) /= work%nrow_input_data ) goto 330
         do i = 1, work%nrow_input_data
            if( data_patt(i) < 1 ) goto 340
            if( data_patt(i) > work%n_data_patt ) goto 340
         end do
         allocate( work%data_patt( work%nrow_input_data ), stat=status )
         if( status /= 0 ) goto 100
         work%data_patt(:) = data_patt(:)
       end if
     !
       if( size( cov_patt_for_data_patt, kind=our_int ) &
            /= work%n_data_patt ) goto 350
      do i = 1, work%n_data_patt
         if( cov_patt_for_data_patt(i) < 1 ) goto 360 
         if( cov_patt_for_data_patt(i) > work%n_cov_patt ) goto 360 
      end do
      allocate( work%cov_patt_for_data_patt( work%n_data_patt ), &
           stat=status )
      if( status /= 0 ) goto 100
      work%cov_patt_for_data_patt(:) = cov_patt_for_data_patt(:)
      !
      if( size( response_for_data_patt, kind=our_int ) &
           /= work%n_data_patt ) goto 370
      do i = 1, work%n_data_patt
         if( response_for_data_patt(i) < 1 ) goto 380 
         if( response_for_data_patt(i) > work%n_levels ) goto 380 
      end do
      allocate( work%response_for_data_patt( work%n_data_patt ), &
           stat=status )
      if( status /= 0 ) goto 100
      work%response_for_data_patt(:) = response_for_data_patt(:)
      !
      allocate( work%data_patt_st( work%n_cov_patt ), &
           work%data_patt_fin( work%n_cov_patt ), stat=status )
      if( status /= 0 ) goto 100
      cur_cov_patt = 0
      posn = 0
      do dp = 1, work%n_data_patt
         cp = work%cov_patt_for_data_patt( dp )
         if( cp /= cur_cov_patt ) then
            posn = posn + 1
            cur_cov_patt = cp
            if( posn > work%n_cov_patt ) goto 502
            work%data_patt_st( posn ) = dp
         end if
      end do
      if( posn /= work%n_cov_patt ) goto 503
      do posn = 2, work%n_cov_patt
         work%data_patt_fin( posn - 1 ) = &
              work%data_patt_st( posn ) - 1
      end do
      work%data_patt_fin( work%n_cov_patt ) = work%n_data_patt
      !
      allocate( work%active_ystar( work%n_cov_patt, work%n_levels ), &
           stat=status )
      if( status /= 0 ) goto 100
      work%active_ystar(:,:) = .false.
      do cp = 1, work%n_cov_patt
         do dp = work%data_patt_st(cp), work%data_patt_fin(cp)
            ystar = work%response_for_data_patt(dp)
            if( work%fitweight_data_patt(dp) > 0.D0 ) &
                 work%active_ystar( cp, ystar ) = .true.
         end do
      end do
      !
      if( prior_int == 1 ) then
         work%prior = "none"
         work%prior_freq_tot = 0.D0
         work%prior_alloc_supplied = .false.
      else if( prior_int == 2 ) then
         work%prior = "DAP"
         if( prior_freq_tot < 0.D0 ) goto 430
         work%prior_freq_tot = prior_freq_tot
         if( prior_alloc_supplied_int == 0 ) then
            work%prior_alloc_supplied = .false.
         else if( prior_alloc_supplied_int == 1 ) then
            work%prior_alloc_supplied = .true.
         else
            goto 431
         end if
      else
         goto 420
      end if
      if( size( prior_alloc, kind=our_int ) /= work%n_levels ) goto 432
      allocate( work%prior_alloc( work%n_levels ), stat=status )
      if( status /= 0 ) goto 100
      if( work%prior_alloc_supplied ) then
         sum = 0.D0
         do j = 1, work%n_levels
            if( prior_alloc(j) < 0.D0 ) goto 433
            if( prior_alloc(j) > 1.D0 ) goto 433
            sum = sum + prior_alloc(j)
         end do
         if( abs( sum - 1.D0 ) > tol ) goto 434
         work%prior_alloc(:) = prior_alloc(:) / sum
      else
         work%prior_alloc(:) = 0.D0  ! handled later if prior="DAP"
      end if
      !
      allocate( work%empty_cov_patt( work%n_cov_patt ), stat=status )
      if( status /= 0 ) goto 100
      work%empty_cov_patt(:) = .true.
      work%n_cov_patt_empty = 0
      do cp = 1, work%n_cov_patt
         do ystar = 1, work%n_levels
            if( work%active_ystar(cp, ystar) ) then
               work%empty_cov_patt(cp) = .false.
               exit
            end if
         end do
         if( work%empty_cov_patt(cp) ) work%n_cov_patt_empty = &
              work%n_cov_patt_empty + 1
      end do
      !
      allocate( work%cov_patt_order( work%nrow_input_data ), &
           work%cov_patt_st( work%n_cov_patt ), &
           work%cov_patt_fin( work%n_cov_patt ), stat=status )
      if( status /= 0 ) goto 100
      if( qsort( work%cov_patt, work%cov_patt_order, work%nrow_input_data, &
           .false., err ) == RETURN_FAIL ) goto 800
      cur_cov_patt = 0
      posn = 0
      do i = 1, work%nrow_input_data
         cp = work%cov_patt( work%cov_patt_order(i) )
         if( cp /= cur_cov_patt ) then
            posn = posn + 1
            cur_cov_patt = cp
            if( posn > work%n_cov_patt ) goto 602
            work%cov_patt_st( posn ) = i
         end if
      end do
      if( posn /= work%n_cov_patt ) goto 603
      do posn = 2, work%n_cov_patt
         work%cov_patt_fin( posn - 1 ) = &
              work%cov_patt_st( posn ) - 1
      end do
      work%cov_patt_fin( work%n_cov_patt ) = work%nrow_input_data
      !
      !#### note: to cycle through all the data rows corresponding
      !#### to a given covariate pattern cp, do this:
      ! cp = something between 1 and work%n_cov_patt
      ! do i = work%cov_patt_st(cp), work%cov_patt_fin(cp)
      !    dr = work%cov_patt_order(i)  ( a row of the input data set )
      !    ( do something with dr )
      ! end do
      !
      if( size( pert_mat, 1, kind=our_int ) /= work%n_levels ) goto 390
      if( size( pert_mat, 2, kind=our_int ) /= work%n_levels ) goto 390
      allocate( work%pert_mat( work%n_levels, work%n_levels ), &
           stat=status )
      if( status /= 0 ) goto 100
      do j = 1, work%n_levels
         sum = 0.D0
         do i = 1, work%n_levels
            if( pert_mat(i,j) < 0.D0 ) goto 400
            if( pert_mat(i,j) > 1.D0 ) goto 400
            sum = sum + pert_mat(i,j)
         end do
         if( abs( sum - 1.D0 ) > tol ) goto 410
      end do
      work%pert_mat(:,:) = pert_mat(:,:)
      work%pert_mat_identity = .true.
      do i = 1, work%n_levels
         if( work%pert_mat(i,i) /= 1.D0 ) work%pert_mat_identity = .false.
         do j = 1, work%n_levels
            if( j == i ) cycle
            if( work%pert_mat(i,j) /= 0.D0 ) work%pert_mat_identity = .false.
         end do
      end do
      if( size( pert_mat_inv, 1, kind=our_int ) /= work%n_levels ) goto 415
      if( size( pert_mat_inv, 2, kind=our_int ) /= work%n_levels ) goto 415
      allocate( work%pert_mat_inv( work%n_levels, work%n_levels ), &
           stat=status )
      if( status /= 0 ) goto 100
      work%pert_mat_inv(:,:) = pert_mat_inv(:,:)
      do i = 1, work%n_levels
         do j = 1, work%n_levels
            sum = 0.D0
            do k = 1, work%n_levels
               sum = sum + work%pert_mat_inv(i,k) * work%pert_mat(k,j) 
            end do
            if( i == j ) then
               if( abs( sum - 1.D0 ) > tol ) goto 416
            else
               if( abs( sum ) > tol ) goto 416
            end if
         end do
      end do
      !
      work%n = work%nrow_input_data
      work%r = work%n_levels
      work%p = work%ncol_model_matrix
      work%d = work%nparam_this_model
      allocate( work%beta( work%p, work%r ), &
           work%beta_vec( work%d ), work%beta_vec_old( work%d ), &
           work%beta_vec_new( work%d ), &
           work%beta_mstep( work%d ), work%beta_mstep_old( work%d ), &
           work%beta_mstep_new( work%d ), &
           work%pi_mat( work%n_cov_patt, work%r), &
           work%pistar_mat( work%n_cov_patt, work%r), &
           work%score( work%d ), work%hess( work%d, work%d ), &
           work%pi_mat_old( work%n_cov_patt, work%r), &
           work%pistar_mat_old( work%n_cov_patt, work%r), &
           work%pi_mat_mstep( work%n_cov_patt, work%r), &
           work%pi_mat_mstep_old( work%n_cov_patt, work%r), &
           work%phi_mat( work%r, work%r ), &
           work%fhat_mat( work%n_cov_patt, work%r ), &
           work%fstar_mat( work%n_cov_patt, work%r ), &
           work%fitweight_mat( work%n_cov_patt, work%r ), &
           work%loglik_vec( work%iter_max ), &
           work%logP_vec( work%iter_max ), &
           work%iter_em_saturated( work%n_cov_patt ), &
           work%converged_em_saturated( work%n_cov_patt ), &
           work%boundary_em_saturated( work%n_cov_patt ), &
           work%scoreA( work%d ), work%hessA( work%d, work%d ), &
           work%prior_freq( work%n_cov_patt, work%r ), &
           work%prior_freq_tot_cov_patt( work%n_cov_patt ), &
           work%marg_proportions_ystar( work%r ), &
           work%marg_proportions_y( work%r ), &
           work%mini_em_prior_freq( work%r ), work%mini_em_result( work%r ), &
           work%pi_marg( work%r ), &
           work%vhat_coef( work%d, work%d ), &
           work%pi_vec( work%r ), work%pi_vec_old( work%r ), &
           work%pistar_vec( work%r ), work%fhat_vec( work%r ), &
           work%se_fit_mat(work%n_cov_patt, work%r), &
           work%vhat_fit_array(work%n_cov_patt, work%r, work%r), &
           work%fimp_vec( work%r ), &
           work%fimp_mat( work%n_cov_patt, work%r ), &
           stat=status )
      if( status /= 0 ) goto 100
      work%beta_vec(:) = 0.D0
      work%score(:) = 0.D0
      work%hess(:,:) = 0.D0
      work%vhat_coef(:,:) = 0.D0
      work%pi_mat(:,:) = 0.D0
      work%pistar_mat(:,:) = 0.D0
      work%marg_proportions_ystar(:) = 0.D0
      work%marg_proportions_y(:) = 0.D0
      work%fhat_mat(:,:) = 0.D0
      work%hessA(:,:) = 0.D0
      work%loglik_vec(:) = 0.D0
      work%logP_vec(:) = 0.D0
      allocate( work%dldpistar( work%r ), work%dldpi( work%r ), &
           work%d2ldpi2( work%r, work%r ), &
           work%Edldpistar( work%r ), work%Edldpi( work%r ), &
           work%dpidbeta( work%r, work%d ), &
           work%d2pidbeta2( work%r, work%d, work%d), &
           stat=status )
      if( status /= 0 ) goto 100
      allocate( work%wkrA( work%r ), work%wkrB( work%r ), &
           work%wkrdA( work%r, work%d ), work%wkrdB( work%r, work%d ), &
           work%wkrdC( work%r, work%d ), &
           work%wkddA( work%d, work%d ), work%wkddB( work%d, work%d ), &
           work%wkdA( work%d ), work%wkdB( work%d ), work%wkdC( work%d ), &
           work%wkppA( work%p, work%p ), work%wkppB( work%p, work%p ), &
           work%wkpA( work%p ), &
           work%wkpB( work%p ), work%wkprA( work%p, work%r ), &
           work%wkprB( work%p, work%r ), &
           stat=status )
      if( status /= 0 ) goto 100
      work%mvcode = mvcode
      work%nancode = nancode
      work%infcode = infcode
      work%neginfcode = neginfcode
      !
      work%fstar_mat(:,:) = 0.D0
      do cp = 1, work%n_cov_patt
         do dp = work%data_patt_st(cp), work%data_patt_fin(cp)
            ystar = work%response_for_data_patt(dp)
            work%fstar_mat( cp, ystar ) = freq_for_data_patt(dp)
            work%fitweight_mat( cp, ystar ) = work%fitweight_data_patt(dp)
         end do
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
100   call err_handle(err, 1, &
            comment = "Unable to allocate array" )
      goto 800
110   call err_handle(err, 1, &
            comment = "Negative value for nrow_input_data" )
      goto 800
120   call err_handle(err, 1, &
            comment = "Value for n_levels must be at least 2" )
      goto 800
130   call err_handle(err, 1, &
            comment = "Value for n_cov_patt is out of range" )
      goto 800
140   call err_handle(err, 1, &
            comment = "Value for n_data_patt is out of range" )
      goto 800
150   call err_handle(err, 1, &
            comment = "Value for ncol_model_matrix is out of range" )
      goto 800
160   call err_handle(err, 1, &
            comment = "Value for nparam_this_model is out of range" )
      goto 800
170   call err_handle(err, 1, &
            comment = "Value for nparam_sat_model is out of range" )
      goto 800
180   call err_handle(err, 1, &
            comment = "Model has more parameters than a saturated model" )
      goto 800
190   call err_handle(err, 1, &
            comment = "Value for baseline_int is out of range" )
      goto 800
200   call err_handle(err, 1, &
            comment = "Value for survey_mode_int is out of range" )
      goto 800
205   call err_handle(err, 1, &
            comment = "Value for saturated_int is out of range" )
      goto 800
210   call err_handle(err, 1, &
            comment = "Value for wide_format_int is out of range" )
      goto 800
211   call err_handle(err, 1, &
            comment = "Value for method_int is out of range" )
      goto 800
212   call err_handle(err, 1, &
            comment = "Value for iter_max_nr is out of range" )
      goto 800
213   call err_handle(err, 1, &
            comment = "Value for iter_max_fs is out of range" )
      goto 800
214   call err_handle(err, 1, &
            comment = "Value for iter_max_em is out of range" )
      goto 800
215   call err_handle(err, 1, &
            comment = "Value for iter_max_mstep is out of range" )
      goto 800
216   call err_handle(err, 1, &
            comment = "Value for crit_converge is out of range" )
      goto 800
217   call err_handle(err, 1, &
            comment = "Value for crit_boundary is out of range" )
      goto 800
220   call err_handle(err, 1, &
            comment = "Array model_matrix has incorrect size" )
      goto 800
230   call err_handle(err, 1, &
            comment = "Array fitweight_row_input_data has incorrect size" )
      goto 800
240   call err_handle(err, 1, &
            comment = "Negative value found in fitweight_row_input_data" )
      goto 800
245   call err_handle(err, 1, &
            comment = "Sum of fitweight_row_input_data is not positive" )
      goto 800
250   call err_handle(err, 1, &
            comment = "Array fitweight_cov_patt has incorrect size" )
      goto 800
260   call err_handle(err, 1, &
            comment = "Negative value found in fitweight_cov_patt" )
      goto 800
270   call err_handle(err, 1, &
            comment = "Sum of fitweight_cov_patt is not correct" )
      goto 800
280   call err_handle(err, 1, &
            comment = "Array fitweight_data_patt has incorrect size" )
      goto 800
290   call err_handle(err, 1, &
            comment = "Negative value found in fitweight_data_patt" )
      goto 800
300   call err_handle(err, 1, &
            comment = "Sum of fitweight_data_patt is not correct" )
      goto 800
305   call err_handle(err, 1, &
            comment = "Array freq_for_data_patt has incorrect size" )
      goto 800
310   call err_handle(err, 1, &
            comment = "Array cov_patt has incorrect size" )
      goto 800
320   call err_handle(err, 1, &
            comment = "Element of cov_patt is out of range" )
      goto 800
330   call err_handle(err, 1, &
            comment = "Array data_patt has incorrect size" )
      goto 800
340   call err_handle(err, 1, &
            comment = "Element of data_patt is out of range" )
      goto 800
350   call err_handle(err, 1, &
            comment = "Array cov_patt_for_data_patt has incorrect size" )
      goto 800
360   call err_handle(err, 1, &
            comment = "Element of cov_patt_for_data_patt is out of range" )
      goto 800
370   call err_handle(err, 1, &
            comment = "Array response_for_data_patt has incorrect size" )
      goto 800
380   call err_handle(err, 1, &
            comment = "Element of response_for_data_patt is out of range" )
      goto 800
390   call err_handle(err, 1, &
            comment = "Array pert_mat has incorrect size" )
      goto 800
400   call err_handle(err, 1, &
            comment = "Element of pert_mat is out of range" )
      goto 800
410   call err_handle(err, 1, &
            comment = "Column of pert_mat does not sum to one" )
      goto 800
415   call err_handle(err, 1, &
            comment = "Array pert_mat_inv has incorrect size" )
      goto 800
416   call err_handle(err, 1, &
            comment = "Array pert_mat_inv is not the inverse of pert_mat" )
      goto 800
420   call err_handle(err, 1, &
            comment = "Value of prior_int not recognized" )
      goto 800
430   call err_handle(err, 1, &
            comment = "Value of prior_freq_tot out of range" )
      goto 800
431   call err_handle(err, 1, &
            comment = "Value of prior_alloc_supplied_int out of range" )
      goto 800
432   call err_handle(err, 1, &
            comment = "Array prior_alloc has incorrect size" )
      goto 800
433   call err_handle(err, 1, &
            comment = "Element of prior_alloc is out of range" )
      goto 800
434   call err_handle(err, 1, &
            comment = "Elements of prior_alloc do not sum to 1.0" )
      goto 800
502   call err_handle(err, 1, &
           comment = "Array bounds exceeded for data_patt_st" ) 
      goto 800
503   call err_handle(err, 1, &
           comment = "Bad data in cov_patt_for_data_patt" ) 
      goto 800
602   call err_handle(err, 1, &
           comment = "Array bounds exceeded for cov_patt_st" ) 
      goto 800
603   call err_handle(err, 1, &
           comment = "Bad data in cov_patt_order" ) 
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
   end function put_data_into_workspace_rrlogit
   !###################################################################
   integer(kind=our_int) function create_mcmc_objects_rrlogit( &
        iter_approx_bayes, impute_approx_bayes_int, &
        iter_mcmc, burn_mcmc, thin_mcmc, &
        impute_every, type_mcmc_int, stuck_limit, &
        df_da, step_size_da, scale_fac_da, df_rwm, scale_fac_rwm, &
        series_length, n_impute, vhat_coef_rwm, &
        micro_data_int, freq_row_input_data_int, &
        f_mat_int, response_row_input_data, freq_for_data_patt_int, &
        coef_mcmc, coef_vec_mcmc, vhat_coef_mcmc, fitted_pi_mcmc, &
        fitted_pistar_mcmc, &
        n_iter_actual, n_sample_actual, n_imp_actual, &
        mh_accept_rate, start_logP, logP_series, coef_vec_series, &
        imp_mat_series, imp_vec_series, pi_marg_mcmc, pi_marg_series, &
        approx_bayes_log_imp_ratios, work, err ) result(answer)
      implicit none
      ! inputs
      integer(kind=our_int), intent(in) :: iter_approx_bayes
      integer(kind=our_int), intent(in) :: impute_approx_bayes_int
      integer(kind=our_int), intent(in) :: iter_mcmc
      integer(kind=our_int), intent(in) :: burn_mcmc
      integer(kind=our_int), intent(in) :: thin_mcmc
      integer(kind=our_int), intent(in) :: impute_every
      integer(kind=our_int), intent(in) :: type_mcmc_int
      integer(kind=our_int), intent(in) :: stuck_limit
      real(kind=our_dble), intent(in) :: df_da
      real(kind=our_dble), intent(in) :: step_size_da
      real(kind=our_dble), intent(in) :: scale_fac_da
      real(kind=our_dble), intent(in) :: df_rwm
      real(kind=our_dble), intent(in) :: scale_fac_rwm
      integer(kind=our_int), intent(in) :: series_length
      integer(kind=our_int), intent(in) :: n_impute
      real(kind=our_dble), intent(in) :: vhat_coef_rwm(:,:)
      integer(kind=our_int), intent(in) :: micro_data_int
      integer(kind=our_int), intent(in) :: freq_row_input_data_int(:)
      integer(kind=our_int), intent(in) :: f_mat_int(:,:)
      integer(kind=our_int), intent(in) :: response_row_input_data(:)
      integer(kind=our_int), intent(in) :: freq_for_data_patt_int(:)
      real(kind=our_dble), intent(out) :: coef_mcmc(:,:)
      real(kind=our_dble), intent(out) :: coef_vec_mcmc(:)
      real(kind=our_dble), intent(out) :: vhat_coef_mcmc(:,:)
      real(kind=our_dble), intent(out) :: fitted_pi_mcmc(:,:)
      real(kind=our_dble), intent(out) :: fitted_pistar_mcmc(:,:)
      integer(kind=our_int), intent(out) :: n_iter_actual
      integer(kind=our_int), intent(out) :: n_sample_actual
      integer(kind=our_int), intent(out) :: n_imp_actual
      real(kind=our_dble), intent(out) :: mh_accept_rate
      real(kind=our_dble), intent(out) :: start_logP
      real(kind=our_dble), intent(out) :: logP_series(:)
      real(kind=our_dble), intent(out) :: coef_vec_series(:,:)
      integer(kind=our_int), intent(out) :: imp_mat_series(:,:,:)
      integer(kind=our_int), intent(out) :: imp_vec_series(:,:)
      real(kind=our_dble), intent(out) :: pi_marg_mcmc(:)
      real(kind=our_dble), intent(out) :: pi_marg_series(:,:)
      real(kind=our_dble), intent(out) :: approx_bayes_log_imp_ratios(:)
      ! workspaces
      type(workspace_type_rrlogit), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! locals
      integer(kind=our_int) :: itmp
      integer(kind=our_int) :: status
      character(len=*), parameter :: &
           subname = "create_mcmc_objects_rrlogit"
      ! begin
      answer = RETURN_FAIL
      !####
      if( .not. ( work%method == "MCMC" .or. &
          work%method == "approxBayes" ) ) goto 10    ! nothing to do
      !
      if( iter_approx_bayes < 0 ) goto 156
      work%iter_approx_bayes = iter_approx_bayes
      if( impute_approx_bayes_int == 0 ) then
         work%impute_approx_bayes = .false.
      else if( impute_approx_bayes_int == 1 ) then
         work%impute_approx_bayes = .true.
      else
         goto 157
      end if
      if( iter_mcmc < 0 ) goto 200
      work%iter_mcmc_nominal = iter_mcmc
      if( burn_mcmc < 0 ) goto 210
      work%burn_mcmc = burn_mcmc
      if( thin_mcmc < 0 ) goto 220
      work%thin_mcmc = thin_mcmc
      if( impute_every < 0 ) goto 230
      work%impute_every = impute_every
      work%iter_mcmc = floor( real(iter_mcmc, our_dble) / &
         real(thin_mcmc, our_dble) ) * thin_mcmc
      if( mod( iter_mcmc, thin_mcmc ) > 0 ) work%iter_mcmc = &
           work%iter_mcmc + work%thin_mcmc
      if( type_mcmc_int == 1 ) then
         work%type_mcmc = "DA"
      else if( type_mcmc_int == 2 ) then
         work%type_mcmc = "RWM"
      else
         goto 250
      end if
      if( work%saturated .and. ( work%method == "approxBayes" ) ) goto 251
      if( work%saturated .and. ( work%type_mcmc == "RWM" ) ) then
         work%type_mcmc = "DA"
      end if
      if( stuck_limit <= 0 ) goto 255
      work%stuck_limit = stuck_limit
      if( df_da <= 0.D0 ) goto 260
      work%df_da = df_da
      work%step_size_da = step_size_da
      if( scale_fac_da <= 0.D0 ) goto 270
      work%scale_fac_da = scale_fac_da
      if( df_rwm <= 0.D0 ) goto 260
      work%df_rwm = df_rwm
      if( scale_fac_rwm <= 0.D0 ) goto 270
      work%scale_fac_rwm = scale_fac_rwm
      !
      if( series_length < 0 ) goto 280
      itmp = work%iter_mcmc / work%thin_mcmc
      if( ( work%method == "MCMC" ) .and. ( series_length /= itmp ) ) goto 280
      if( ( work%method == "approxBayes" ) .and. &
           ( series_length /= work%iter_approx_Bayes ) ) goto 280
      if( ( work%method /= "MCMC" ) .and. ( work%method /= "approxBayes" ) &
           .and. ( series_length /= 0 ) ) goto 280
      work%series_length = series_length
      !
      if( n_impute < 0 ) goto 290
      if( work%method == "MCMC" ) then
         if( work%impute_every == 0 ) then
            if( n_impute /= 0 ) goto 290
         else
            itmp = floor( real( work%iter_mcmc, our_dble ) / &
                 real( work%impute_every, our_dble ) )
            if( n_impute /= itmp ) goto 290
         end if
      else if( ( work%method == "approxBayes" ) .and. &
           work%impute_approx_bayes ) then
         if( n_impute /= work%iter_approx_bayes ) goto 290
      else
         if( n_impute /= 0 ) goto 290
      end if
      work%n_impute = n_impute
      !
      if( work%saturated ) then
         if( size( vhat_coef_rwm, 1, kind=our_int ) /= 0 ) goto 300
         if( size( vhat_coef_rwm, 2, kind=our_int ) /= 0 ) goto 300
         if( size( coef_mcmc, 1, kind=our_int ) /= 0 ) goto 310
         if( size( coef_mcmc, 2, kind=our_int ) /= 0 ) goto 310
         if( size( coef_vec_mcmc, kind=our_int ) /= 0 ) goto 320
         if( size( vhat_coef_mcmc, 1, kind=our_int ) /= 0 ) goto 330
         if( size( vhat_coef_mcmc, 2, kind=our_int ) /= 0 ) goto 330
         if( size( coef_vec_series, 1, kind=our_int ) &
              /= work%series_length ) goto 340
         if( size( coef_vec_series, 2, kind=our_int ) /= 0 ) goto 340
      else
         if( size( vhat_coef_rwm, 1, kind=our_int ) /= work%d ) goto 300
         if( size( vhat_coef_rwm, 2, kind=our_int ) /= work%d ) goto 300
         allocate( work%vhat_beta_rwm(work%d, work%d), stat=status )
         if( status /= 0 ) goto 100
         work%vhat_beta_rwm(:,:) = vhat_coef_rwm(:,:)
         if( size( coef_mcmc, 1, kind=our_int ) /= work%p ) goto 310
         if( size( coef_mcmc, 2, kind=our_int ) /= work%r ) goto 310
         coef_mcmc(:,:) = work%mvcode
         if( size( coef_vec_mcmc, kind=our_int ) /= work%d ) goto 320
         coef_vec_mcmc(:) = work%mvcode
         if( size( vhat_coef_mcmc, 1, kind=our_int ) /= work%d ) goto 330
         if( size( vhat_coef_mcmc, 2, kind=our_int ) /= work%d ) goto 330
         vhat_coef_mcmc(:,:) = work%mvcode
         if( size( coef_vec_series, 1, kind=our_int ) &
              /= work%series_length ) goto 340
         if( size( coef_vec_series, 2, kind=our_int ) &
              /= work%d ) goto 340
         coef_vec_series(:,:) = work%mvcode
      end if
      if( size( fitted_pi_mcmc, 1, kind=our_int ) /= work%n_cov_patt ) goto 350
      if( size( fitted_pi_mcmc, 2, kind=our_int ) /= work%r ) goto 350
      fitted_pi_mcmc(:,:) = work%mvcode
      if( size( fitted_pistar_mcmc, 1, kind=our_int ) &
           /= work%n_cov_patt ) goto 360
      if( size( fitted_pistar_mcmc, 2, kind=our_int ) &
           /= work%r ) goto 360
      fitted_pistar_mcmc(:,:) = work%mvcode
      n_iter_actual = 0
      n_sample_actual = 0
      n_imp_actual = 0
      mh_accept_rate = 0.D0
      start_logP = 0.D0
      if( size( logP_series, kind=our_int ) /= work%series_length ) goto 370
      logP_series(:) = work%mvcode
      if( size( imp_mat_series, 1, kind=our_int ) &
           /= work%nrow_input_data ) goto 380
      if( size( imp_mat_series, 2, kind=our_int ) /= work%r ) goto 380
      if( size( imp_mat_series, 3, kind=our_int ) /= work%n_impute ) goto 380
      imp_mat_series(:,:,:) = -1 ! integer missing value code
      if( size( imp_vec_series, 1, kind=our_int ) &
           /= work%nrow_input_data ) goto 390
      if( size( imp_vec_series, 2, kind=our_int ) &
           /= work%n_impute ) goto 390
      imp_vec_series(:,:) = -1 ! integer missing value code
      if( size( pi_marg_mcmc, kind=our_int ) /= work%r ) goto 391
      pi_marg_mcmc(:) = work%mvcode
      if( size( pi_marg_series, 1, kind=our_int ) &
           /= work%series_length ) goto 392
      if( size( pi_marg_series, 2, kind=our_int ) /= work%r ) goto 392
      pi_marg_series(:,:) = work%mvcode
      if( size( approx_bayes_log_imp_ratios, kind=our_int ) &
           /= work%series_length ) goto 393
      approx_bayes_log_imp_ratios(:) = work%mvcode
      !
      allocate( work%mh_ratios_beta(work%iter_mcmc + work%burn_mcmc), &
           work%mh_accept_beta(work%iter_mcmc + work%burn_mcmc), &
           work%beta_can(work%d), work%beta_center(work%d), &
           work%beta_scale(work%d,work%d), & 
           work%beta_scale_inv(work%d,work%d), & 
           work%beta_scale_inv_sqrt(work%d,work%d), & 
           work%beta_scale_sqrt(work%d,work%d), & 
           work%beta_mean(work%d), work%beta_cov_mat(work%d, work%d), &
           work%pi_mat_mean(work%n_cov_patt, work%r), &
           work%pistar_mat_mean(work%n_cov_patt, work%r), &
           work%pi_marg_mean( work%r ), &
           stat=status )
      if( status /= 0 ) goto 100
      !
      if( micro_data_int == 0 ) then
         work%micro_data = .false.
      else if( micro_data_int == 1 ) then
         work%micro_data = .true.
      else
         goto 400
      end if
      !
      if( size(freq_row_input_data_int, kind=our_int) &
           /= work%nrow_input_data ) goto 410
      if( size(f_mat_int, 1, kind=our_int)  /= work%nrow_input_data ) goto 420
      if( size(f_mat_int, 2, kind=our_int)  /= work%r ) goto 420
      allocate( work%freq_row_input_data_int(work%nrow_input_data), &
           work%f_mat_int(work%nrow_input_data,work%r), &
           stat = status )
      if( status /= 0 ) goto 100
      work%freq_row_input_data_int(:) = freq_row_input_data_int(:)
      work%f_mat_int(:,:) = f_mat_int(:,:)
      if( work%wide_format ) then
         if( size(response_row_input_data, kind=our_int) /= 0 ) goto 430
          allocate( work%response_row_input_data(0), stat = status )
         if( status /= 0 ) goto 100
      else
         if( size(response_row_input_data, kind=our_int) &
              /= work%nrow_input_data ) goto 430
         allocate( work%response_row_input_data(work%nrow_input_data), &
              stat = status )
         if( status /= 0 ) goto 100
         work%response_row_input_data(:) = response_row_input_data(:)
      end if
      if( size(freq_for_data_patt_int, kind=our_int) &
           /= work%n_data_patt ) goto 440
      allocate( work%freq_for_data_patt_int(work%n_data_patt), &
           stat = status )
      if( status /= 0 ) goto 100
      work%freq_for_data_patt_int(:) = freq_for_data_patt_int(:)
10    continue
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
100   call err_handle(err, 1, &
            comment = "Unable to allocate array" )
      goto 800
156   call err_handle(err, 1, &
            comment = "Negative value for iter_approx_bayes" )
      goto 800
157   call err_handle(err, 1, &
            comment = "Value for impute_approx_bayes_int not recognized" )
      goto 800
200   call err_handle(err, 1, &
            comment = "Negative value for iter_mcmc" )
      goto 800
210   call err_handle(err, 1, &
            comment = "Negative value for burn_mcmc" )
      goto 800
220   call err_handle(err, 1, &
            comment = "Negative value for thin_mcmc" )
      goto 800
230   call err_handle(err, 1, &
            comment = "Negative value for impute_every" )
      goto 800
250   call err_handle(err, 1, &
            comment = "Value for type_mcmc_int not recognized" )
      goto 800
251   call err_handle(err, 1, &
            comment = "Approx Bayes not implemented for saturated model" )
      goto 800
255   call err_handle(err, 1, &
            comment = "Value for stuck_limit not positive" )
      goto 800
260   call err_handle(err, 1, &
            comment = "Non-positive value for mcmc proposal df" )
      goto 800
270   call err_handle(err, 1, &
            comment = "Non-positive value for mcmc proposal scale_fac" )
      goto 800
280   call err_handle(err, 1, &
            comment = "Inconsistent value for series_length" )
      goto 800
290   call err_handle(err, 1, &
            comment = "Inconsistent value for n_impute" )
      goto 800
300   call err_handle(err, 1, &
            comment = "Array vhat_coef_rwm has incorrect size" )
      goto 800
310   call err_handle(err, 1, &
            comment = "Array coef_mcmc has incorrect size" )
      goto 800
320   call err_handle(err, 1, &
            comment = "Array coef_vec_mcmc has incorrect size" )
      goto 800
330   call err_handle(err, 1, &
            comment = "Array vhat_coef_mcmc has incorrect size" )
      goto 800
340   call err_handle(err, 1, &
            comment = "Array coef_vec_series has incorrect size" )
      goto 800
350   call err_handle(err, 1, &
            comment = "Array fitted_pi_mcmc has incorrect size" )
      goto 800
360   call err_handle(err, 1, &
            comment = "Array fitted_pistar_mcmc has incorrect size" )
      goto 800
370   call err_handle(err, 1, &
            comment = "Array logP_series has incorrect size" )
      goto 800
380   call err_handle(err, 1, &
            comment = "Array imp_mat_series has incorrect size" )
      goto 800
390   call err_handle(err, 1, &
            comment = "Array imp_vec_series has incorrect size" )
      goto 800
391   call err_handle(err, 1, &
            comment = "Array pi_marg_mcmc has incorrect size" )
      goto 800
392   call err_handle(err, 1, &
            comment = "Array pi_marg_series has incorrect size" )
      goto 800
393   call err_handle(err, 1, &
            comment = "Array approx_bayes_log_imp_ratios has incorrect size" )
      goto 800
400   call err_handle(err, 1, &
            comment = "Value of micro_data_int out of range" )
      goto 800
410   call err_handle(err, 1, &
            comment = "Array freq_row_input_data_int has incorrect size" )
      goto 800
420   call err_handle(err, 1, &
            comment = "Array f_mat_int has incorrect size" )
      goto 800
430   call err_handle(err, 1, &
            comment = "Array response_row_input_data has incorrect size" )
      goto 800
440   call err_handle(err, 1, &
            comment = "Array freq_for_data_patt_int has incorrect size" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
   end function create_mcmc_objects_rrlogit
   !###################################################################
   integer(kind=our_int) function compute_prior_freq_rrlogit( &
        work, err ) result(answer)
      implicit none
      ! workspaces
      type(workspace_type_rrlogit), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! locals
      integer(kind=our_int) :: cp, dp, ystar, y
      real(kind=our_dble) :: sum, rtmp
      character(len=*), parameter :: &
           subname = "compute_prior_freq_rrlogit"
      ! begin
      answer = RETURN_FAIL
      ! estimate marginal proportions for ystar in observed data
      work%wkrA(:) = 0.D0
      do cp = 1, work%n_cov_patt
         if( work%empty_cov_patt(cp) ) cycle
         do dp = work%data_patt_st(cp), work%data_patt_fin(cp)
            ystar = work%response_for_data_patt(dp)
            work%wkrA(ystar) = work%wkrA(ystar) + &
                 work%fitweight_data_patt(dp)
         end do
      end do
      sum = 0.D0
      do ystar = 1, work%r
         sum = sum + work%wkrA(ystar)
      end do
      do ystar = 1, work%r
         work%marg_proportions_ystar(ystar) = work%wkrA(ystar) / sum
      end do
      ! estimate marginal proportions for y
      work%neg_prob = .false.
      do y = 1, work%r
         sum = 0.D0
         do ystar = 1, work%r
            sum = sum + work%pert_mat_inv(y,ystar) * &
                 work%marg_proportions_ystar(ystar)
         end do
         if( sum < 0.D0 ) work%neg_prob = .true.
         work%marg_proportions_y(y) = sum
      end do
      if( ( work%prior == "DAP" ) .and. ( .not. work%prior_alloc_supplied ) &
           .and. work%neg_prob ) then
         call err_handle(err, 1, &
           comment = "Estimated marginal prob. for true response is negative" ) 
      end if
      ! fill in matrix of prior frequencies
      work%dap_failed = .false.
      work%prior_freq(:,:) = 0.D0
      if( work%prior == "DAP" ) then
         if( .not. work%prior_alloc_supplied ) then
            if( work%neg_prob ) then
               work%dap_failed = .true.
               call err_handle(err, 1, &
                    comment = "Default prior allocation for DAP failed" )
            else
               work%prior_alloc(:) = work%marg_proportions_y(:)
            end if
         end if
         if( .not. work%dap_failed ) then
            rtmp = work%prior_freq_tot / &
                 real( work%n_cov_patt - work%n_cov_patt_empty, our_dble )
            do cp = 1, work%n_cov_patt
               if( work%empty_cov_patt(cp) ) then
                  work%prior_freq_tot_cov_patt(cp) = 0.D0
                  work%prior_freq(cp,:) = 0.D0
               else
                  work%prior_freq_tot_cov_patt(cp) = rtmp
                  work%prior_freq(cp,:) = work%prior_alloc(:) * rtmp
               end if
            end do
         end if
      end if
      !####
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function compute_prior_freq_rrlogit
   !###################################################################
    integer(kind=our_int) function compute_start_val_rrlogit( &
         coefficients, fitted_pi, start_val_source_int, start_val_jitter, &
         work, err ) result(answer)
      ! starting values procedure
      implicit none
      ! inputs
      real(kind=our_dble), intent(in) :: coefficients(:,:)
      real(kind=our_dble), intent(in) :: fitted_pi(:,:)
      integer(kind=our_int), intent(in) :: start_val_source_int
      real(kind=our_dble), intent(in) :: start_val_jitter
      ! workspaces
      type(workspace_type_rrlogit), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! locals
      integer(kind=our_int) :: cp, y, j, k
      real(kind=our_dble) :: sum, rtmp
      character(len=*), parameter :: &
           subname = "compute_start_val_rrlogit"
      ! begin
      answer = RETURN_FAIL
      !####
      sum = 0.D0 ! to suppress spurious warning message
      if( size( coefficients, 1, kind=our_int ) /= work%p ) goto 100
      if( size( coefficients, 2, kind=our_int ) /= work%r ) goto 100
      if( size( fitted_pi, 1, kind=our_int ) /= work%n_cov_patt ) goto 130
      if( size( fitted_pi, 2, kind=our_int ) /= work%r ) goto 130
      if( start_val_source_int == 1 ) then
         work%start_val_source = "default"
      else if( start_val_source_int == 2 ) then
         work%start_val_source = "supplied by user"
      else if( start_val_source_int == 3 ) then
         work%start_val_source = "from rrLogit object"
      else
         goto 140
      end if
      if( start_val_jitter < 0.D0 ) goto 150
      work%start_val_jitter = start_val_jitter
      !
!      if( run_mini_em_rrlogit( work, err ) == RETURN_FAIL ) goto 800
      !
      if( work%saturated ) then
         ! set starting values for fitted probabilities
         if( work%start_val_source == "default" ) then
            if( work%neg_prob ) then
               ! set to uniform values
               rtmp = 1.D0 / real(work%r, our_dble)
               work%pi_mat(:,:) = rtmp
            else
               ! set to estimated marginal proportions
               do cp = 1, work%n_cov_patt
                  work%pi_mat(cp,:) = work%marg_proportions_y(:)
               end do
            end if
         else
            ! normalize the supplied probs
            do cp = 1, work%n_cov_patt
               sum = 0.D0
               do y = 1, work%r
                  sum = sum + fitted_pi(cp,y)
               end do
               work%pi_mat(cp,:) = fitted_pi(cp,:) / sum
            end do
         end if
         if( work%start_val_jitter > 0.D0 ) then
            do cp = 1, work%n_cov_patt
               do y = 1, work%r
                  sum = 0.D0
                  if( rnorm_R( rtmp, err ) == RETURN_FAIL ) goto 800
                  work%pi_mat(cp,y) = work%pi_mat(cp,y) * &
                       exp( rtmp * work%start_val_jitter )
                  sum = sum + work%pi_mat(cp,y)
               end do
               work%pi_mat(cp,:) = work%pi_mat(cp,:) / sum
            end do
         end if
      else
         if( work%start_val_source == "default" ) then
            ! fill pi_mat with default values
            if( work%neg_prob ) then
               rtmp = 1.D0 / real(work%r, our_dble)
               work%pi_mat(:,:) = rtmp
            else
               do cp = 1, work%n_cov_patt
                  work%pi_mat(cp,:) = work%marg_proportions_y(:)
               end do
            end if
            ! accumulate xtx (lower tri) and xty
            work%wkppA(:,:) = 0.D0
            work%wkprA(:,:) = 0.D0
            do cp = 1, work%n_cov_patt
               do j = 1, work%p
                  do k = 1, j
                     work%wkppA(j,k) = work%wkppA(j,k) + &
                          work%model_matrix(cp,j) * work%model_matrix(cp,k)
                  end do
                  ! convert to logits
                  do k = 1, work%r
                     if( k == work%baseline ) cycle
                     if( work%pi_mat(cp, work%baseline ) == 0.D0 ) goto 200
                     rtmp = work%pi_mat(cp,k) / work%pi_mat(cp, work%baseline)
                     if( rtmp <= 0.D0 ) goto 250
                     work%wkrA(k) = log(rtmp)
                  end do
                  work%wkrA(work%baseline) = 0.D0
                  do k = 1, work%r
                     work%wkprA(j,k) = work%wkprA(j,k) + &
                          work%model_matrix(cp,j) * work%wkrA(k)
                  end do
               end do
            end do
            if( cholesky_in_place(work%wkppA, err) == RETURN_FAIL ) goto 800
            if( invert_lower(work%wkppA, err) == RETURN_FAIL ) goto 800
            if( premult_lower_by_transpose(work%wkppA, &
                 work%wkppB, err) == RETURN_FAIL ) goto 800
            work%wkprB = matmul( work%wkppB, work%wkprA )
            if( stack_coef( work%wkprB, work%beta_vec, work ) &
                 == RETURN_FAIL ) goto 800
         else
            if( stack_coef( coefficients, work%beta_vec, work ) &
                 == RETURN_FAIL ) goto 800
         end if
         if( work%start_val_jitter > 0.D0 ) then
            do j = 1, work%d
               if( rnorm_R( rtmp, err ) == RETURN_FAIL ) goto 800
               work%beta_vec(j) = work%beta_vec(j) + &
                    rtmp * work%start_val_jitter
            end do
         end if
      end if
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
100   call err_handle(err, 1, &
            comment = "Array coefficients has incorrect size" )
      goto 800
130   call err_handle(err, 1, &
            comment = "Array fitted_pi has incorrect size" )
      goto 800
140   call err_handle(err, 1, &
            comment = "Value of start_val_source_int not recognized" )
      goto 800
150   call err_handle(err, 1, &
            comment = "Value of start_val_jitter is negative" )
      goto 800
200   call err_handle(err, 1, &
            comment = "Attempted division by zero" )
      goto 800
250   call err_handle(err, 1, &
            comment = "Attempted logarithm of non-positive number" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function compute_start_val_rrlogit
   !###################################################################
   integer(kind=our_int) function run_mini_em_rrlogit( &
        work, err ) result(answer)
      ! runs a mini EM algorithm to estimate the marginal proportions
      ! for the noise-free response, taking into account the DAP
      ! prior if prior allocation was supplied by user, and smoothing
      ! toward uniform values by adding 1/2 observation to each
      ! category of y
      implicit none
      ! workspaces
      type(workspace_type_rrlogit), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! locals
      integer(kind=our_int) :: cp, dp, ystar, y, iter
      integer(kind=our_int), parameter :: iter_max = 1000
      real(kind=our_dble), parameter :: crit = 1D-06
      logical :: converged
      real(kind=our_dble) :: sum, rtmp, loglik, logprior, logP, max_abs_diff
      character(len=*), parameter :: &
           subname = "run_mini_em"
      ! begin
      answer = RETURN_FAIL
      ! initial flattening
      work%mini_em_prior_freq_tot = 1.D0 * real( work%r, kind=our_dble )
      work%mini_em_prior_freq(:) = 1.D0
      ! add DAP if supplied
      if( work%prior == "DAP" ) then
         work%mini_em_prior_freq_tot = work%mini_em_prior_freq_tot + &
              work%prior_freq_tot
         if( work%dap_failed ) then
            rtmp =  work%prior_freq_tot / real( work%r, kind=our_dble )
            do y = 1, work%r
               work%mini_em_prior_freq(y) = work%mini_em_prior_freq(y) + rtmp
            end do
         else
            do y = 1, work%r
               rtmp = work%prior_freq_tot * work%prior_alloc(y)
               work%mini_em_prior_freq(y) = work%mini_em_prior_freq(y) + rtmp
            end do
         end if
      end if
      ! put observed marginal total frequencies for ystar into wkrA
      work%wkrA(:) = 0.D0
      do cp = 1, work%n_cov_patt
         if( work%empty_cov_patt(cp) ) cycle
         do dp = work%data_patt_st(cp), work%data_patt_fin(cp)
            ystar = work%response_for_data_patt(dp)
            work%wkrA(ystar) = work%wkrA(ystar) + &
                 work%fitweight_data_patt(dp)
         end do
      end do
      ! put starting values for marginal pi into mini_em_result
      sum = 0.D0
      do y = 1, work%r
         sum = sum + work%mini_em_prior_freq(y)
      end do
      work%mini_em_result(:) = work%mini_em_prior_freq(:) / sum
      ! main iteration
      iter = 0
      converged = .false.
      do
         if( iter >= iter_max ) exit
         if( converged ) exit
         iter = iter + 1
         if( iter > 1 ) work%mini_em_result(:) = work%wkrB(:)
         ! compute logprior
         logprior = 0.D0
         do y = 1, work%r
            if( work%mini_em_result(y) <= 0.D0 ) goto 100
            logprior = logprior + work%mini_em_prior_freq(y) * &
                 log( work%mini_em_result(y) )
         end do
         ! compute loglik and logP
         loglik = 0.D0
         do ystar = 1, work%r
            sum = 0.D0
            do y = 1, work%r
               sum = sum + work%pert_mat(ystar,y) * work%mini_em_result(y)
            end do
            if( sum <= 0.D0 ) goto 100
            loglik = loglik + work%wkrA(ystar) * log(sum)
         end do
         logP = logprior + loglik
         ! compute phi_mat assuming pi = mini_em_result
         do ystar = 1, work%r
            sum = 0.D0
            do y = 1, work%r
               rtmp = work%pert_mat( ystar, y ) * work%mini_em_result(y)
               work%phi_mat( ystar, y ) =  rtmp
               sum = sum + rtmp
            end do
            if( sum == 0.D0 ) goto 200
            do y = 1, work%r
               work%phi_mat( ystar, y ) = work%phi_mat( ystar, y ) / sum
            end do
         end do
         ! pre-multiply wkrA by t(phi_mat), put result into wkrB
         do y = 1, work%r
            sum = 0.D0
            do ystar = 1, work%r
               sum = sum + work%phi_mat(ystar,y) * work%wkrA(ystar)
            end do
            work%wkrB(y) = sum
         end do
         ! add in prior freqs, then normalize
         sum = 0.D0
         do y = 1, work%r
            work%wkrB(y) = work%wkrB(y) + work%mini_em_prior_freq(y)
            sum = sum + work%wkrB(y)
         end do
         if( sum == 0.D0 ) goto 200
         do y = 1, work%r
            work%wkrB(y) = work%wkrB(y) / sum
         end do
         ! check convergence
         max_abs_diff = 0.D0
         do y = 1, work%r
            max_abs_diff = max( max_abs_diff, &
                 abs( work%wkrB(y) - work%mini_em_result(y) ) )
         end do
         if( max_abs_diff <= crit ) converged = .true.
      end do
      if( iter > 0 ) work%mini_em_result(:) = work%wkrB(:)
      !####
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
100   call err_handle(err, 1, &
            comment = "Attempted logarithm of non-positive number" )
      goto 800
200   call err_handle(err, 1, &
            comment = "Attempted division by zero" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function run_mini_em_rrlogit
   !###################################################################
   integer(kind=our_int) function stack_coef( beta, beta_vec, &
        work ) result(answer)
      ! stacks columns of beta into beta_vec, omitting the column 
      ! corresponding to the baseline level
      implicit none
      ! inputs
      real(kind=our_dble), intent(in) :: beta(:,:)
      ! outputs
      real(kind=our_dble), intent(out) :: beta_vec(:)
      ! workspaces
      type(workspace_type_rrlogit), intent(inout) :: work
!      type(error_type), intent(inout) :: err
      ! locals
      integer(kind=our_int) :: j, k, posn
      character(len=*), parameter :: &
           subname = "stack_coef"
      ! begin
      answer = RETURN_FAIL
      posn = 0
      do k = 1, work%r
         if( k == work%baseline ) cycle
         do j = 1, work%p
            posn = posn + 1
            beta_vec(posn) = beta(j,k)
         end do
      end do
      !####
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
!800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
!      goto 999
      ! cleanup
999   continue
   end function stack_coef
   !###################################################################
   integer(kind=our_int) function unstack_coef( beta_vec, beta, &
        work ) result(answer)
      ! unstacks beta_vec into columns of beta
      implicit none
      ! inputs
      real(kind=our_dble), intent(in) :: beta_vec(:)
      ! outputs
      real(kind=our_dble), intent(out) :: beta(:,:)
      ! workspaces
      type(workspace_type_rrlogit), intent(inout) :: work
!      type(error_type), intent(inout) :: err
      ! locals
      integer(kind=our_int) :: j, k, posn
      character(len=*), parameter :: &
           subname = "unstack_coef"
      ! begin
      answer = RETURN_FAIL
      posn = 0
      do k = 1, work%r
         if( k == work%baseline ) then
            do j = 1, work%p
               beta(j,k) = 0.D0
            end do
         else
            do j = 1, work%p
               posn = posn + 1
               beta(j,k) = beta_vec(posn)
            end do
         end if
      end do
      !####
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
!800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
!      goto 999
      ! cleanup
999   continue
   end function unstack_coef
   !###################################################################
   integer(kind=our_int) function compute_logP_score_hess_rrlogit( &
        work, err, skip_compute_pi_mat ) result(answer)
      ! computes logP and its derivatives at current value of beta_vec
      implicit none
      ! workspaces
      type(workspace_type_rrlogit), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! optionals
      logical, intent(in), optional :: skip_compute_pi_mat
      ! locals
      logical :: skip_compute_pi_mat_local, done
      integer(kind=our_int) :: cp, dp, ystar, y, j, k, l, m, posn, posnA
      real(kind=our_dble) :: sum, rtmp, rtmpA, rtmpB
      character(len=*), parameter :: &
           subname = "compute_logP_score_hess_rrlogit"
      ! begin
      answer = RETURN_FAIL
      if( present( skip_compute_pi_mat ) ) then
         skip_compute_pi_mat_local = skip_compute_pi_mat
      else
         skip_compute_pi_mat_local = .false.
      end if
      if( .not. skip_compute_pi_mat_local ) then
         if( compute_pi_mat_rrlogit( work, err ) &
              == RETURN_FAIL ) goto 800
      end if
      !
      work%loglik = 0.D0
      work%score(:) = 0.D0
      work%hess(:,:) = 0.D0
      work%logprior = 0.D0
      do cp = 1, work%n_cov_patt
         if( work%empty_cov_patt(cp) ) cycle
         if( work%prior == "DAP" ) then
            ! accumulate logprior
            do y = 1, work%r
               if( work%pi_mat(cp, y) <= 0.D0 ) goto 200
               work%logprior = work%logprior + work%prior_freq(cp,y) * &
                    log( work%pi_mat(cp,y) )
            end do
         end if
         ! accumulate loglik and compute dldpistar
         do dp = work%data_patt_st(cp), work%data_patt_fin(cp)
            ystar = work%response_for_data_patt(dp)
            if( .not. work%active_ystar(cp, ystar) ) cycle
            if( work%pistar_mat(cp, ystar) <= 0.D0 ) goto 200
            work%loglik = work%loglik + work%fitweight_data_patt(dp) * &
                 log( work%pistar_mat(cp, ystar) )
            ! note: zeros in dldpistar are not zeroed out; we will skip
            ! over these elements in all later computations
            work%dldpistar( ystar ) = work%fitweight_data_patt(dp) / &
                 work%pistar_mat(cp, ystar) 
         end do
         ! compute dldpi
         if( work%pert_mat_identity ) then
            ! copy dldpistar into dldpi, inserting zeros where necessary
            do y = 1, work%r
               if( .not. work%active_ystar(cp, y) ) then
                  work%dldpi(y) = 0.D0
               else
                  work%dldpi(y) = work%dldpistar(y)
               end if
            end do
         else
            ! pre-multiply dldpistar by pert_mat, put result into dldpi
            do y = 1, work%r
               sum = 0.D0
               do ystar = 1, work%r
                  if( .not. work%active_ystar(cp, ystar) ) cycle
                  sum = sum + work%pert_mat(ystar, y) * work%dldpistar(ystar)
               end do
               work%dldpi(y) = sum
            end do
         end if
         ! compute wkrdA = pert_mat %*% dpidbeta, skipping the rows that
         ! correspond to non-ystar-active positions
         if( work%pert_mat_identity .and. ( work%prior == "none" )  ) then
            ! compute ystar-active rows of dpidbeta, putting the results
            ! into wkrdA; the remaining rows will not be needed
            do y = 1, work%r
               if( .not. work%active_ystar(cp, y) ) cycle
               posn = 0               ! indexes element of beta_vec
               do k = 1, work%r       ! indexes column of beta
                  if( k == work%baseline ) cycle
                  do j = 1, work%p    ! indexes row of beta
                     posn = posn + 1
                     rtmp = 0.D0
                     if( k == y ) rtmp = 1.D0
                     work%wkrdA(y, posn) = work%pi_mat(cp, y) * &
                          ( rtmp - work%pi_mat(cp, k) ) * &
                          work%model_matrix(cp, j)
                  end do
               end do
            end do
         else
            ! compute the entire dpidbeta matrix, then compute the
            ! ystar-active rows of wkrdA
            do y = 1, work%r
               posn = 0               ! indexes element of beta_vec
               do k = 1, work%r       ! indexes column of beta
                  if( k == work%baseline ) cycle
                  do j = 1, work%p    ! indexes row of beta
                     posn = posn + 1
                     rtmp = 0.D0
                     if( k == y ) rtmp = 1.D0
                     work%dpidbeta(y, posn) = work%pi_mat(cp, y) * &
                          ( rtmp - work%pi_mat(cp, k) ) * &
                          work%model_matrix(cp, j)
                  end do
               end do
            end do
            do ystar = 1, work%r
               if( .not. work%active_ystar(cp, ystar) ) cycle
               do posn = 1, work%d
                  sum = 0.D0
                  do y = 1, work%r
                     sum = sum + work%pert_mat(ystar, y) * &
                          work%dpidbeta(y, posn)
                  end do
                  work%wkrdA(ystar, posn) = sum
               end do
            end do
         end if
         ! pre-multiply wkrdA by sqrt( - d2ldpistar2 ), put result into
         ! wkrdB. Skip over the non-ystar-active rows, because these
         ! are known to be zero, even though they are not actually
         ! set to zero in wkrdB.
         do dp = work%data_patt_st(cp), work%data_patt_fin(cp)
            ystar = work%response_for_data_patt(dp)
            if( .not. work%active_ystar(cp, ystar) ) cycle
            if( work%pistar_mat(cp, ystar) == 0.D0 ) goto 210
            rtmp = sqrt( work%fitweight_data_patt(dp) ) / &
                 work%pistar_mat(cp, ystar)
            do posn = 1, work%d
               work%wkrdB(ystar, posn) = rtmp * work%wkrdA(ystar, posn)
            end do
         end do
         ! if prior = "DAP", pre-multiply dpidbeta by 
         ! sqrt( - d2logpriordpi2 ), put the result into wkrdC
         if( work%prior == "DAP" ) then
            do y = 1, work%r
               if( work%pi_mat(cp,y) == 0.D0 ) goto 210
               rtmp = sqrt( work%prior_freq(cp,y) ) / work%pi_mat(cp,y)
               do posn = 1, work%d     
                  work%wkrdC(y, posn) = rtmp * work%dpidbeta(y,posn)
               end do
            end do
         end if
         ! compute layers of d2pidbeta2, skipping over those
         ! corresponding to elements of dldpi that happen to be zero
         ! (that will happen a lot if pert_mat is the identity)
         ! but don't skip any if prior == "DAP"
         do y = 1, work%r
            if( ( work%dldpi(y) == 0.D0 ) .and. &
                 ( work%prior == "none"  ) ) cycle
            posn = 0
            do k = 1, work%r
               if( k == work%baseline ) cycle
               rtmp = 0.D0
               if( k == y ) rtmp = 1.D0
               do j = 1, work%p
                  posn = posn + 1
                  done = .false.
                  posnA = 0
                  do m = 1, work%r
                     if( m == work%baseline ) cycle
                     rtmpA = 0.D0
                     if( k == m ) rtmpA = 1.D0
                     rtmpB = 0.D0
                     if( m == y ) rtmpB = 1.D0
                     do l = 1, work%p
                        posnA = posnA + 1
                        if( posnA > posn ) then  ! only compute the
                           done = .true.         ! lower triangle
                           exit
                        end if
                        work%d2pidbeta2(y, posn, posnA) = &
                              - work%pi_mat(cp, y) * ( &
                              work%pi_mat(cp, k ) * &
                              ( rtmpA - work%pi_mat(cp, m ) ) &
                              - ( rtmp - work%pi_mat(cp, k ) ) &
                              * ( rtmpB - work%pi_mat(cp, m ) ) &
                              ) * work%model_matrix(cp, j) * &
                              work%model_matrix(cp, l)
                     end do
                     if( done ) exit
                  end do
               end do
            end do
         end do
         ! put derivatives of log-prior wrt pi into wkrA
         if( work%prior == "DAP" ) then
            do y = 1, work%r
               if( work%pi_mat(cp,y) == 0.D0 ) goto 210
               work%wkrA(y) = work%prior_freq(cp,y) / work%pi_mat(cp,y)
            end do
         end if
         ! accumulate score = t(wkrdA) %*% dldpistar
         ! and lower-triangle of first term of hessian
         do posn = 1, work%d
            sum = 0.D0
            do ystar = 1, work%r
               if( .not. work%active_ystar(cp, ystar) ) cycle
               sum = sum + work%wkrdA( ystar, posn ) * &
                    work%dldpistar( ystar )
            end do
            work%score(posn) = work%score(posn) + sum
            ! add in first derivative of log prior
            if( work%prior == "DAP" ) then
               sum = 0.D0
               do y = 1, work%r
                  sum = sum + work%dpidbeta(y,posn) * work%wkrA(y)
               end do
               work%score(posn) = work%score(posn) + sum
            end if
            ! first term of hessian
            do posnA = 1, posn
               do y = 1, work%r
                  if( work%dldpi(y) == 0.D0 ) cycle
                  work%hess( posn, posnA ) = work%hess( posn, posnA ) + &
                       work%dldpi(y) * work%d2pidbeta2(y, posn, posnA)
               end do
            end do
            ! corresponding term for log-prior
            if( work%prior == "DAP" ) then
               do posnA = 1, posn
                  do y = 1, work%r
                     work%hess( posn, posnA ) = work%hess( posn, posnA ) + &
                          work%wkrA(y) * work%d2pidbeta2(y, posn, posnA)
                  end do
               end do
            end if
         end do
         ! accumulate second term of hessian, lower-triangle only
         do posn = 1, work%d
            do posnA = 1, posn
               sum = 0.D0
               do ystar = 1, work%r
                  if( .not. work%active_ystar(cp, ystar) ) cycle
                  sum = sum + work%wkrdB( ystar, posn ) * &
                       work%wkrdB( ystar, posnA )
               end do
               work%hess( posn, posnA ) = work%hess( posn, posnA ) - sum
            end do
         end do
         ! add in corresponding term for log-prior
         if( work%prior == "DAP" ) then
            do posn = 1, work%d
               do posnA = 1, posn
                  sum = 0.D0
                  do y = 1, work%r
                     sum = sum + work%wkrdC(y,posn) * work%wkrdC(y,posnA) 
                  end do
                  work%hess( posn, posnA ) = work%hess( posn, posnA ) - sum
               end do
            end do
         end if
      end do
      ! fill in upper triangle of hess
      do posn = 1, work%d
         do posnA = posn + 1, work%d
            work%hess( posn, posnA ) = work%hess( posnA, posn )
         end do
      end do
      !####
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
200   call err_handle(err, 1, &
            comment = "Attempted logarithm of non-positive number" )
      goto 800
210   call err_handle(err, 1, &
            comment = "Attempted division by zero" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
   end function compute_logP_score_hess_rrlogit
   !###################################################################
   integer(kind=our_int) function compute_logP_score_Ehess_rrlogit( &
        work, err, skip_compute_pi_mat ) result(answer)
      ! computes logP, first derivative and expected Hessian at 
      ! current value of beta_vec
      implicit none
      ! workspaces
      type(workspace_type_rrlogit), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! optionals
      logical, intent(in), optional :: skip_compute_pi_mat
      ! locals
      logical :: skip_compute_pi_mat_local, done
      integer(kind=our_int) :: cp, dp, ystar, y, j, k, l, m, posn, posnA
      real(kind=our_dble) :: sum, sumA, rtmp, rtmpA, rtmpB
      character(len=*), parameter :: &
           subname = "compute_logP_score_Ehess_rrlogit"
      ! begin
      answer = RETURN_FAIL
      if( present( skip_compute_pi_mat ) ) then
         skip_compute_pi_mat_local = skip_compute_pi_mat
      else
         skip_compute_pi_mat_local = .false.
      end if
      if( .not. skip_compute_pi_mat_local ) then
         if( compute_pi_mat_rrlogit( work, err ) &
              == RETURN_FAIL ) goto 800
      end if
      !
      work%loglik = 0.D0
      work%score(:) = 0.D0
      work%hess(:,:) = 0.D0
      work%logprior = 0.D0
      do cp = 1, work%n_cov_patt
         if( work%empty_cov_patt(cp) ) cycle
         if( work%prior == "DAP" ) then
            ! accumulate logprior
            do y = 1, work%r
               if( work%pi_mat(cp, y) <= 0.D0 ) goto 200
               work%logprior = work%logprior + work%prior_freq(cp,y) * &
                    log( work%pi_mat(cp,y) )
            end do
         end if
         ! accumulate loglik and compute dldpistar, Edldpistar
         do dp = work%data_patt_st(cp), work%data_patt_fin(cp)
            ystar = work%response_for_data_patt(dp)
            if( .not. work%active_ystar(cp, ystar) ) cycle
            if( work%pistar_mat(cp, ystar) <= 0.D0 ) goto 200
            work%loglik = work%loglik + work%fitweight_data_patt(dp) * &
                 log( work%pistar_mat(cp, ystar) )
            ! note: zeros in dldpistar are not zeroed out; we will skip
            ! over these elements in all later computations
            work%dldpistar( ystar ) = work%fitweight_data_patt(dp) / &
                 work%pistar_mat(cp, ystar) 
         end do
         work%Edldpistar(:) = work%fitweight_cov_patt(cp)
         ! compute dldpi, Edldpi
         if( work%pert_mat_identity ) then
            work%Edldpi(:) = work%Edldpistar(:)
            ! copy dldpistar into dldpi, inserting zeros where necessary
            do y = 1, work%r
               if( .not. work%active_ystar(cp, y) ) then
                  work%dldpi(y) = 0.D0
               else
                  work%dldpi(y) = work%dldpistar(y)
               end if
            end do
         else
            ! pre-multiply dldpistar by pert_mat, put result into dldpi
            do y = 1, work%r
               sum = 0.D0
               sumA = 0.D0
               do ystar = 1, work%r
                  sumA = sumA + work%pert_mat(ystar, y) * &
                       work%Edldpistar(ystar)
                  if( .not. work%active_ystar(cp, ystar) ) cycle
                  sum = sum + work%pert_mat(ystar, y) * work%dldpistar(ystar)
               end do
               work%dldpi(y) = sum
               work%Edldpi(y) = sumA
            end do
         end if
         ! compute wkrdA = pert_mat %*% dpidbeta, all rows 
         if( work%pert_mat_identity ) then
            ! compute all rows of dpidbeta, putting the results into wkrdA
            do y = 1, work%r
               posn = 0               ! indexes element of beta_vec
               do k = 1, work%r       ! indexes column of beta
                  if( k == work%baseline ) cycle
                  do j = 1, work%p    ! indexes row of beta
                     posn = posn + 1
                     rtmp = 0.D0
                     if( k == y ) rtmp = 1.D0
                     work%wkrdA(y, posn) = work%pi_mat(cp, y) * &
                          ( rtmp - work%pi_mat(cp, k) ) * &
                          work%model_matrix(cp, j)
                     work%dpidbeta(y, posn) = work%wkrdA(y, posn)
                  end do
               end do
            end do
         else
            ! compute the entire dpidbeta matrix, then compute all
            ! rows of wkrdA
            do y = 1, work%r
               posn = 0               ! indexes element of beta_vec
               do k = 1, work%r       ! indexes column of beta
                  if( k == work%baseline ) cycle
                  do j = 1, work%p    ! indexes row of beta
                     posn = posn + 1
                     rtmp = 0.D0
                     if( k == y ) rtmp = 1.D0
                     work%dpidbeta(y, posn) = work%pi_mat(cp, y) * &
                          ( rtmp - work%pi_mat(cp, k) ) * &
                          work%model_matrix(cp, j)
                  end do
               end do
            end do
            do ystar = 1, work%r
               do posn = 1, work%d
                  sum = 0.D0
                  do y = 1, work%r
                     sum = sum + work%pert_mat(ystar, y) * &
                          work%dpidbeta(y, posn)
                  end do
                  work%wkrdA(ystar, posn) = sum
               end do
            end do
         end if
         ! pre-multiply wkrdA by E(d2ldpistar2), put result into wkrdB.
         do ystar = 1, work%r
            rtmp =  - work%fitweight_cov_patt(cp) / &
                 work%pistar_mat(cp, ystar)
            do posn = 1, work%d
               work%wkrdB(ystar, posn) = rtmp * work%wkrdA(ystar, posn)
            end do
         end do

         ! if prior = "DAP", pre-multiply dpidbeta by 
         ! sqrt( - d2logpriordpi2 ), put the result into wkrdC
         if( work%prior == "DAP" ) then
            do y = 1, work%r
               if( work%pi_mat(cp,y) == 0.D0 ) goto 210
               rtmp = sqrt( work%prior_freq(cp,y) ) / work%pi_mat(cp,y)
               do posn = 1, work%d     
                  work%wkrdC(y, posn) = rtmp * work%dpidbeta(y,posn)
               end do
            end do
         end if
         ! compute all layers of d2pidbeta2
         do y = 1, work%r
            posn = 0
            do k = 1, work%r
               if( k == work%baseline ) cycle
               rtmp = 0.D0
               if( k == y ) rtmp = 1.D0
               do j = 1, work%p
                  posn = posn + 1
                  done = .false.
                  posnA = 0
                  do m = 1, work%r
                     if( m == work%baseline ) cycle
                     rtmpA = 0.D0
                     if( k == m ) rtmpA = 1.D0
                     rtmpB = 0.D0
                     if( m == y ) rtmpB = 1.D0
                     do l = 1, work%p
                        posnA = posnA + 1
                        if( posnA > posn ) then  ! only compute the
                           done = .true.         ! lower triangle
                           exit
                        end if
                        work%d2pidbeta2(y, posn, posnA) = &
                              - work%pi_mat(cp, y) * ( &
                              work%pi_mat(cp, k ) * &
                              ( rtmpA - work%pi_mat(cp, m ) ) &
                              - ( rtmp - work%pi_mat(cp, k ) ) &
                              * ( rtmpB - work%pi_mat(cp, m ) ) &
                              ) * work%model_matrix(cp, j) * &
                              work%model_matrix(cp, l)
                     end do
                     if( done ) exit
                  end do
               end do
            end do
         end do
         ! put derivatives of log-prior wrt pi into wkrA
         if( work%prior == "DAP" ) then
            do y = 1, work%r
               if( work%pi_mat(cp,y) == 0.D0 ) goto 210
               work%wkrA(y) = work%prior_freq(cp,y) / work%pi_mat(cp,y)
            end do
         end if
         ! accumulate score = t(wkrdA) %*% dldpistar
         ! and lower-triangle of first term of Ehessian
         do posn = 1, work%d
            sum = 0.D0
            do ystar = 1, work%r
               if( .not. work%active_ystar(cp, ystar) ) cycle
               sum = sum + work%wkrdA( ystar, posn ) * &
                    work%dldpistar( ystar )
            end do
            work%score(posn) = work%score(posn) + sum
            ! add in first derivative of log prior
            if( work%prior == "DAP" ) then
               sum = 0.D0
               do y = 1, work%r
                  sum = sum + work%dpidbeta(y,posn) * work%wkrA(y)
               end do
               work%score(posn) = work%score(posn) + sum
            end if
            ! first term of hessian
            do posnA = 1, posn
               do y = 1, work%r
                  work%hess( posn, posnA ) = work%hess( posn, posnA ) + &
                       work%Edldpi(y) * work%d2pidbeta2(y, posn, posnA)
               end do
            end do
            ! corresponding term for log-prior
            if( work%prior == "DAP" ) then
               do posnA = 1, posn
                  do y = 1, work%r
                     work%hess( posn, posnA ) = work%hess( posn, posnA ) + &
                          work%wkrA(y) * work%d2pidbeta2(y, posn, posnA)
                  end do
               end do
            end if
         end do
         ! accumulate second term of Ehessian, lower-triangle only
         do posn = 1, work%d
            do posnA = 1, posn
               sum = 0.D0
               do ystar = 1, work%r
                  sum = sum + work%wkrdA( ystar, posn ) * &
                       work%wkrdB( ystar, posnA )
               end do
               work%hess( posn, posnA ) = work%hess( posn, posnA ) + sum
            end do
         end do
         ! add in corresponding term for log-prior
         if( work%prior == "DAP" ) then
            do posn = 1, work%d
               do posnA = 1, posn
                  sum = 0.D0
                  do y = 1, work%r
                     sum = sum + work%wkrdC(y,posn) * work%wkrdC(y,posnA) 
                  end do
                  work%hess( posn, posnA ) = work%hess( posn, posnA ) - sum
               end do
            end do
         end if
      end do
      ! fill in upper triangle of Ehess
      do posn = 1, work%d
         do posnA = posn + 1, work%d
            work%hess( posn, posnA ) = work%hess( posnA, posn )
         end do
      end do
      !####
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
200   call err_handle(err, 1, &
            comment = "Attempted logarithm of non-positive number" )
      goto 800
210   call err_handle(err, 1, &
            comment = "Attempted division by zero" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
   end function compute_logP_score_Ehess_rrlogit
   !###################################################################
   integer(kind=our_int) function compute_pi_mat_rrlogit( &
        work, err, skip_pistar_mat, logit_scale, logP ) result(answer)
      ! computes pi_mat and pistar_mat given current value of beta_vec
      ! Also computes logprior and loglik, putting them into logprior_new
      ! and loglik_new (if logP is .true.)
      implicit none
      ! workspaces
      type(workspace_type_rrlogit), intent(inout) :: work
      type(error_type), intent(inout) :: err
      logical, intent(in), optional :: skip_pistar_mat
      logical, intent(in), optional :: logit_scale
      logical, intent(in), optional :: logP
      ! locals
      integer(kind=our_int) :: cov_patt, j, posn, k, dp, ystar, cp
      real(kind=our_dble) :: sum, eta_max
      logical :: skip_pistar_mat_local, logit_scale_local, logP_local
      character(len=*), parameter :: &
           subname = "compute_pi_mat_rrlogit"
      ! begin
      answer = RETURN_FAIL
      if( present(skip_pistar_mat) ) then
         skip_pistar_mat_local = skip_pistar_mat
      else
         skip_pistar_mat_local = .false.
      end if
      if( present(logit_scale) ) then
         logit_scale_local = logit_scale
      else
         logit_scale_local = .false.
      end if
      if( present(logP) ) then
         logP_local = logP
      else
         logP_local = .false.
      end if
      ! compute pi_mat from beta_vec
      work%any_zero_pi = .false.
      if( logP_local ) work%logprior_new = 0.D0
      do cov_patt = 1, work%n_cov_patt
         ! put linear predictors (logits) into wkrA
         posn = 0
         do j = 1, work%r
            if( j == work%baseline ) then
               work%wkrA(j) = 0.D0
               cycle
            end if
            sum = 0.D0
            do k = 1, work%p
               posn = posn + 1
               sum = sum + work%model_matrix(cov_patt,k) * &
                    work%beta_vec(posn)
            end do
            work%wkrA(j) = sum
         end do
         if( logit_scale_local ) then
            work%pi_mat(cov_patt,:) = work%wkrA(:)
            cycle
         end if
         ! exponentiate and normalize, being careful to avoid overflow
         ! and underflow; put result into wkrB
         eta_max = log_tiny
         do j = 1, work%r
            eta_max = max(eta_max, work%wkrA(j) )
         end do
         work%wkrA(:) = work%wkrA(:) - eta_max
         sum = 0.D0
         do j = 1, work%r
            if( work%wkrA(j) < log_tiny ) then
               work%wkrB(j) = 0.D0
            else if( work%wkrA(j) > log_huge ) then
               goto 100
            else
               work%wkrB(j) = exp( work%wkrA(j) )
            end if
            sum = sum + work%wkrB(j)
         end do
         if( sum == 0.D0 ) goto 200
         do j = 1, work%r
            work%pi_mat(cov_patt,j) = work%wkrB(j) / sum
            if( work%pi_mat(cov_patt,j) == 0.D0 ) work%any_zero_pi = .true.
         end do
         if( logP_local .and. ( work%prior == "DAP" ) ) then
            do j = 1, work%r
               if( work%logprior_new == work%neginfcode ) exit
               if( work%prior_freq(cov_patt, j) == 0.D0 ) cycle
               if( work%pi_mat(cov_patt, j) <= 0.D0 ) then
                  work%logprior_new = work%neginfcode
               else
                  work%logprior_new = work%logprior_new + &
                    work%prior_freq(cov_patt, j) * &
                    log( work%pi_mat(cov_patt, j) )
               end if
            end do
         end if
      end do
      ! compute pistar_mat from pi_mat
      if( logP_local ) work%loglik_new = 0.D0
      if( .not. skip_pistar_mat_local ) then
         if( work%pert_mat_identity ) then
            work%pistar_mat(:,:) = work%pi_mat(:,:)
         else
            do cov_patt = 1, work%n_cov_patt
               do j = 1, work%r
                  sum = 0.D0
                  do k = 1, work%r
                     sum = sum + work%pert_mat(j,k) * work%pi_mat(cov_patt,k)
                  end do
                  work%pistar_mat(cov_patt,j) = sum
               end do
            end do
         end if
         if( logP_local ) then
            do cp = 1, work%n_cov_patt
               if( work%loglik_new == work%neginfcode ) exit
               do dp = work%data_patt_st(cp), work%data_patt_fin(cp)
                  ystar = work%response_for_data_patt(dp)
                  if( .not. work%active_ystar(cp, ystar) ) cycle
                  if( work%pistar_mat(cp, ystar) <= 0.D0 ) then
                     work%loglik_new = work%neginfcode
                     exit
                  else
                     work%loglik_new = work%loglik_new + &
                          work%fitweight_data_patt(dp) * &
                          log( work%pistar_mat(cp, ystar) )
                  end if
               end do
            end do
         end if
      end if
      !####
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
100   call err_handle(err, 1, &
            comment = "Overflow; fitted value became too large" )
      goto 800
200   call err_handle(err, 1, &
            comment = "Attempted division by zero" )
      goto 800
300   call err_handle(err, 1, &
            comment = "Attempted logarithm of non-positive number (prior)" )
      goto 800
301   call err_handle(err, 1, &
            comment = "Attempted logarithm of non-positive number" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
   end function compute_pi_mat_rrlogit
   !###################################################################
   integer(kind=our_int) function compute_pi_marg( work, err ) &
        result(answer)
      ! computes pi_marg from pi_mat
      implicit none
      ! workspaces
      type(workspace_type_rrlogit), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! locals
      logical :: skip
      integer(kind=our_int) :: cp, y
      real(kind=our_dble) :: sum
      character(len=*), parameter :: &
           subname = "compute_pi_marg"
      ! begin
      answer = RETURN_FAIL
      work%pi_marg(:) = 0.D0
      sum = 0.D0
      do cp = 1, work%n_cov_patt
         skip = .false.
         do y = 1, work%r
            if( work%pi_mat(cp,y) == work%mvcode ) then
               skip = .true.
               exit
            end if
         end do
         if( skip ) cycle
         do y = 1, work%r
            work%pi_marg(y) = work%pi_marg(y) + &
                 work%fitweight_cov_patt(cp) * work%pi_mat(cp,y)
         end do
      end do
      sum = 0.D0
      do y = 1, work%r
         sum = sum + work%pi_marg(y)
      end do
      if( sum == 0.D0 ) goto 200
      work%pi_marg(:) = work%pi_marg(:) / sum
      !####
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
200   call err_handle(err, 1, &
            comment = "Attempted division by zero" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
   end function compute_pi_marg
   !###################################################################
   integer(kind=our_int) function run_rrlogit_nr( work, err ) &
        result(answer)
      ! Newton-Raphson procedure
      implicit none
      ! workspaces
      type(workspace_type_rrlogit), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! locals
      real(kind=our_dble) :: logP_new
      logical :: skip, do_step_halving
      character(len=*), parameter :: &
           subname = "run_rrlogit_nr"
      ! begin
      answer = RETURN_FAIL
      ! prepare for main iteration
      work%iter = 0
      work%converged = .false.
      work%aborted = .false.
      work%loglik_vec(:) = 0.D0
      work%logP_vec(:) = 0.D0
      work%vhat_failed = .true.
      work%score(:) = 0.D0
      work%vhat_coef(:,:) = 0.D0
      do
         if( work%converged .or. ( work%iter >= work%iter_max_nr ) ) exit
         work%iter = work%iter + 1
         work%aborted = .true. ! set to .false. at end of iteration
         skip = .false.
         if( work%iter > 1 ) skip = .true.
         if( compute_logP_score_hess_rrlogit( work, err, &
              skip_compute_pi_mat = skip ) == RETURN_FAIL ) exit
         work%loglik_vec( work%iter ) = work%loglik
         work%logP_vec( work%iter ) = work%loglik + work%logprior
         ! put negative hessian into wkddA and inverse into wkddB
         work%wkddA(:,:) = - work%hess(:,:)
         if( cholesky_in_place(work%wkddA, err, suppress_msg=.false. ) &
              == RETURN_FAIL ) exit
         if( invert_lower(work%wkddA, err, suppress_msg = .false. ) &
              == RETURN_FAIL ) exit
         if( premult_lower_by_transpose( work%wkddA, &
              work%wkddB, err) == RETURN_FAIL ) exit
         work%wkdC = matmul( work%wkddB, work%score ) 
         ! save current params
         work%beta_vec_old(:) = work%beta_vec(:)
         work%pi_mat_old(:,:) = work%pi_mat(:,:)
         ! attempt full step
         do_step_halving = .false.
         work%beta_vec = work%beta_vec + work%wkdC
         if( compute_pi_mat_rrlogit( work, err, logP=.true. )  &
              == RETURN_FAIL ) goto 800
         logP_new = work%logprior_new + work%loglik_new
         if( logP_new < work%logP_vec( work%iter ) ) do_step_halving = .true.
         if( work%any_zero_pi ) do_step_halving = .true.
         if( do_step_halving ) then
            work%step_halving = 0
            do
               work%step_halving = work%step_halving + 1
               work%beta_vec = work%beta_vec_old + work%wkdC * &
                    0.5D0**work%step_halving
               if( compute_pi_mat_rrlogit( work, err, logP=.true. )  &
                    == RETURN_FAIL ) goto 800
               logP_new = work%logprior_new + work%loglik_new
               if( logP_new >= work%logP_vec( work%iter ) ) exit
            end do
         end if
         if( check_convergence_rrlogit( work )  == RETURN_FAIL ) goto 800
         work%aborted = .false.
      end do
      !      
      if( work%aborted ) then
         call err_handle(err, 1, &
              comment = "Newton-Raphson aborted" )
         call err_handle(err, 5, iiter = int(work%iter, kind=our_int) )
         work%beta_vec(:) = work%beta_vec_old(:)
         work%pi_mat(:,:) = work%pi_mat_old(:,:)
         work%pistar_mat(:,:) = work%pistar_mat_old(:,:)
      end if
      ! update loglik, logP, score, vhat_coef
      if( work%iter == 0 ) then
         skip = .false.
      else
         skip = .true.
      end if
      if( compute_logP_score_hess_rrlogit( work, err, &
           skip_compute_pi_mat = skip ) == RETURN_FAIL ) goto 600
      work%wkddA(:,:) = - work%hess(:,:)
      if( cholesky_in_place(work%wkddA, err, suppress_msg=.false. ) &
           == RETURN_FAIL ) goto 600
      if( invert_lower(work%wkddA, err, suppress_msg = .false. ) &
           == RETURN_FAIL ) goto 600
      if( premult_lower_by_transpose( work%wkddA, &
           work%vhat_coef, err) == RETURN_FAIL ) goto 600
      work%vhat_failed = .false.
600   continue
      if( work%vhat_failed ) then
         call err_handle(err, 1, &
              comment = "Hessian-based SEs not available" )
      end if
      !####
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
   end function run_rrlogit_nr
   !###################################################################
   integer(kind=our_int) function run_rrlogit_fs( work, err ) &
        result(answer)
      ! Fisher scoring procedure
      implicit none
      ! workspaces
      type(workspace_type_rrlogit), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! locals
      real(kind=our_dble) :: logP_new
      logical :: skip, do_step_halving
      character(len=*), parameter :: &
           subname = "run_rrlogit_fs"
      ! begin
      answer = RETURN_FAIL
      ! prepare for main iteration
      work%iter = 0
      work%converged = .false.
      work%aborted = .false.
      work%loglik_vec(:) = 0.D0
      work%logP_vec(:) = 0.D0
      work%vhat_failed = .true.
      work%score(:) = 0.D0
      work%vhat_coef(:,:) = 0.D0
      do
         if( work%converged .or. ( work%iter >= work%iter_max_fs ) ) exit
         work%iter = work%iter + 1
         work%aborted = .true. ! set to .false. at end of iteration
         skip = .false.
         if( work%iter > 1 ) skip = .true.
         if( compute_logP_score_Ehess_rrlogit( work, err, &
              skip_compute_pi_mat = skip ) == RETURN_FAIL ) exit
         work%loglik_vec( work%iter ) = work%loglik
         work%logP_vec( work%iter ) = work%loglik + work%logprior
         ! put negative hessian into wkddA and inverse into wkddB
         work%wkddA(:,:) = - work%hess(:,:)
         if( cholesky_in_place(work%wkddA, err, suppress_msg=.false. ) &
              == RETURN_FAIL ) exit
         if( invert_lower(work%wkddA, err, suppress_msg = .false. ) &
              == RETURN_FAIL ) exit
         if( premult_lower_by_transpose( work%wkddA, &
              work%wkddB, err) == RETURN_FAIL ) exit
         work%wkdC = matmul( work%wkddB, work%score ) 
         ! save current params
         work%beta_vec_old(:) = work%beta_vec(:)
         work%pi_mat_old(:,:) = work%pi_mat(:,:)
         ! attempt full step
         do_step_halving = .false.
         work%beta_vec = work%beta_vec + work%wkdC
         if( compute_pi_mat_rrlogit( work, err, logP=.true. ) &
              == RETURN_FAIL ) goto 800
         logP_new = work%logprior_new + work%loglik_new
         if( logP_new < work%logP_vec( work%iter ) ) do_step_halving = .true.
         if( work%any_zero_pi ) do_step_halving = .true.
         if( do_step_halving ) then
            work%step_halving = 0
            do
               work%step_halving = work%step_halving + 1
               work%beta_vec = work%beta_vec_old + work%wkdC * &
                    0.5D0**work%step_halving
               if( compute_pi_mat_rrlogit( work, err, logP=.true. )  &
                    == RETURN_FAIL ) goto 800
               logP_new = work%logprior_new + work%loglik_new
               if( logP_new >= work%logP_vec( work%iter ) ) exit
            end do
         end if
         if( compute_pi_mat_rrlogit( work, err )  == RETURN_FAIL ) goto 800
         if( check_convergence_rrlogit( work )  == RETURN_FAIL ) goto 800
         work%aborted = .false.
      end do
      !
      if( work%aborted ) then
         call err_handle(err, 1, &
              comment = "Fisher scoring aborted" )
         call err_handle(err, 5, iiter = int(work%iter, kind=our_int) )
         work%beta_vec(:) = work%beta_vec_old(:)
         work%pi_mat(:,:) = work%pi_mat_old(:,:)
         work%pistar_mat(:,:) = work%pistar_mat_old(:,:)
      end if
      ! update loglik, logP, score, vhat_coef
      if( work%iter == 0 ) then
         skip = .false.
      else
         skip = .true.
      end if
      if( compute_logP_score_hess_rrlogit( work, err, &
           skip_compute_pi_mat = skip ) == RETURN_FAIL ) goto 600
      work%wkddA(:,:) = - work%hess(:,:)
      if( cholesky_in_place(work%wkddA, err, suppress_msg=.false. ) &
           == RETURN_FAIL ) goto 600
      if( invert_lower(work%wkddA, err, suppress_msg = .false. ) &
           == RETURN_FAIL ) goto 600
      if( premult_lower_by_transpose( work%wkddA, &
           work%vhat_coef, err) == RETURN_FAIL ) goto 600
      work%vhat_failed = .false.
600   continue
      if( work%vhat_failed ) then
         call err_handle(err, 1, &
              comment = "Hessian-based SEs not available" )
      end if
      !####
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
   end function run_rrlogit_fs
   !###################################################################
   integer(kind=our_int) function run_rrlogit_em( work, err ) &
        result(answer)
      ! EM algorithm
      implicit none
      ! workspaces
      type(workspace_type_rrlogit), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! locals
      logical :: skip
      character(len=*), parameter :: &
           subname = "run_rrlogit_em"
      ! begin
      answer = RETURN_FAIL
      ! prepare for main iteration
      work%iter = 0
      work%converged = .false.
      work%aborted = .false.
      work%loglik_vec(:) = 0.D0
      work%logP_vec(:) = 0.D0
      work%vhat_failed = .true.
      work%score(:) = 0.D0
      work%vhat_coef(:,:) = 0.D0
      do
         if( work%converged .or. ( work%iter >= work%iter_max_em ) ) exit
         work%iter = work%iter + 1
         work%aborted = .true. ! set to .false. at end of iteration
         skip = .false.
         if( work%iter > 1 ) skip = .true.
         if( estep_rrlogit( work, err, skip ) == RETURN_FAIL ) exit
         if( work%aborted_estep ) exit
         work%loglik_vec( work%iter ) = work%loglik
         if( ( work%loglik == work%mvcode ) .or. &
              ( work%logprior == work%mvcode ) ) then
            work%logP_vec( work%iter ) = work%loglik + work%logprior
         else
            work%logP_vec( work%iter ) = work%mvcode
         end if
         if( mstep_rrlogit( work, err ) == RETURN_FAIL ) exit
         if( work%aborted_mstep ) exit
         work%beta_vec_old(:) = work%beta_vec(:)
         work%beta_vec(:) = work%beta_mstep(:)
         ! update pi_mat and pistar_mat
         work%pi_mat_old(:,:) = work%pi_mat(:,:)
         work%pistar_mat_old(:,:) = work%pistar_mat(:,:)
         if( compute_pi_mat_rrlogit( work, err )  == RETURN_FAIL ) exit
         if( check_convergence_rrlogit( work )  == RETURN_FAIL ) exit
         work%aborted = .false.
      end do
      !
      if( work%aborted ) then
         call err_handle(err, 1, &
              comment = "EM algorithm aborted" )
         call err_handle(err, 5, iiter = int(work%iter, kind=our_int) )
         work%beta_vec(:) = work%beta_vec_old(:)
         work%pi_mat(:,:) = work%pi_mat_old(:,:)
         work%pistar_mat(:,:) = work%pistar_mat_old(:,:)
      end if
      ! update loglik, logP, score, vhat_coef
      if( work%iter == 0 ) then
         skip = .false.
      else
         skip = .true.
      end if
      if( compute_logP_score_hess_rrlogit( work, err, &
           skip_compute_pi_mat = skip ) == RETURN_FAIL ) goto 600
      work%wkddA(:,:) = - work%hess(:,:)
      if( cholesky_in_place(work%wkddA, err, suppress_msg=.false. ) &
           == RETURN_FAIL ) goto 600
      if( invert_lower(work%wkddA, err, suppress_msg = .false. ) &
           == RETURN_FAIL ) goto 600
      if( premult_lower_by_transpose( work%wkddA, &
           work%vhat_coef, err) == RETURN_FAIL ) goto 600
      work%vhat_failed = .false.
600   continue
      if( work%vhat_failed ) then
         call err_handle(err, 1, &
              comment = "Hessian-based SEs not available" )
      end if
      !####
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
!800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
!      goto 999
      ! cleanup
999   continue
   end function run_rrlogit_em
   !###################################################################
   integer(kind=our_int) function estep_rrlogit( &
        work, err, skip_compute_pi_mat ) result(answer)
      ! computes expected counts for the true response variable in
      ! each covariate pattern given the current value of beta_vec,
      ! putting the results into fhat_mat 
      ! also computes logprior and loglik
      ! If model is saturated, make sure that pi_mat is current and
      ! call this function with skip_compute_pi_mat = .true.
      implicit none
      ! workspaces
      type(workspace_type_rrlogit), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! optionals
      logical, intent(in), optional :: skip_compute_pi_mat
      ! locals
      logical :: skip_compute_pi_mat_local
      integer(kind=our_int) :: dp, cp, y, ystar
      real(kind=our_dble) :: sum
      character(len=*), parameter :: &
           subname = "estep_rrlogit"
      ! begin
      answer = RETURN_FAIL
      work%aborted_estep = .true.
      if( present( skip_compute_pi_mat ) ) then
         skip_compute_pi_mat_local = skip_compute_pi_mat
      else
         skip_compute_pi_mat_local = .false.
      end if
      if( .not. skip_compute_pi_mat_local ) then
         if( compute_pi_mat_rrlogit( work, err ) &
              == RETURN_FAIL ) goto 800
      end if
      !
      work%loglik = 0.D0
      work%logprior = 0.D0
      do cp = 1, work%n_cov_patt
         if( work%saturated .and. work%empty_cov_patt(cp) ) then
            work%fhat_mat(cp,:) = work%mvcode
            cycle
         end if
         ! accumulate logprior
         if( work%prior == "DAP" ) then
            do y = 1, work%r
               if( work%prior_freq(cp,y) == 0.D0 ) cycle
               if( work%logprior == work%mvcode ) cycle
               if( work%pi_mat(cp, y) <= 0.D0 ) then
                  work%logprior = work%mvcode
               else
                  work%logprior = work%logprior + work%prior_freq(cp,y) * &
                       log( work%pi_mat(cp,y) )
               end if
            end do
         end if
         ! compute phi_mat for current covariate pattern, skipping
         ! over the rows that have no observations of ystar; also
         ! accumulate loglik
         do dp = work%data_patt_st(cp), work%data_patt_fin(cp)
            ystar = work%response_for_data_patt(dp)
            if( .not. work%active_ystar(cp, ystar) ) cycle
            if( work%loglik /= work%mvcode ) then 
               if( work%pistar_mat(cp, ystar) <= 0.D0 ) then
                  work%loglik = work%mvcode
               else
                  work%loglik = work%loglik + work%fitweight_data_patt(dp) * &
                       log( work%pistar_mat(cp, ystar) )
               end if
            end if
            if( .not. work%pert_mat_identity ) then
               sum = 0.D0
               do y = 1, work%r
                  work%phi_mat(ystar,y) = work%pert_mat(ystar,y) * &
                       work%pi_mat(cp,y)
                  sum = sum + work%phi_mat(ystar,y)
               end do
               if( sum == 0.D0 ) goto 300
               do y = 1, work%r
                  work%phi_mat(ystar,y) = work%phi_mat(ystar,y) / sum
               end do
            end if
         end do
         ! fill in the row of fhat_mat
         if( work%pert_mat_identity ) then
            work%fhat_mat(cp,:) = 0.D0
            do dp = work%data_patt_st(cp), work%data_patt_fin(cp)
               y = work%response_for_data_patt(dp)
               if( .not. work%active_ystar(cp, y) ) cycle
               work%fhat_mat(cp,y) = work%fitweight_data_patt(dp)
            end do
         else
            do y = 1, work%r
               sum = 0.D0
               do dp = work%data_patt_st(cp), work%data_patt_fin(cp)
                  ystar = work%response_for_data_patt(dp)
                  if( .not. work%active_ystar(cp, ystar) ) cycle
                  sum = sum + work%fitweight_data_patt(dp) * &
                       work%phi_mat(ystar,y)
               end do
               work%fhat_mat(cp,y) = sum
            end do
         end if
      end do
      !####
      work%aborted_estep = .false.
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
200   call err_handle(err, 1, &
            comment = "Attempted logarithm of non-positive number" )
      goto 800
300   call err_handle(err, 1, &
            comment = "Attempted division by zero" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
   end function estep_rrlogit
   !###################################################################
   integer(kind=our_int) function mstep_rrlogit( &
        work, err ) result(answer)
      ! maximized expected augmented-data loglikelihood based on the
      ! expected true-response frequencies in fhat_mat
      implicit none
      ! workspaces
      type(workspace_type_rrlogit), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! locals
      real(kind=our_dble) :: max_abs_diff, logP, logP_new
      logical :: do_step_halving
      integer(kind=our_int) :: cp, y
      character(len=*), parameter :: &
           subname = "mstep_rrlogit"
      ! begin
      answer = RETURN_FAIL
      work%beta_mstep(:) = work%beta_vec(:)
      work%iter_mstep = 0
      work%converged_mstep = .false.
      work%aborted_mstep = .false.
      do
         if( work%converged_mstep .or. &
              ( work%iter_mstep > work%iter_max_mstep) ) exit
         work%iter_mstep = work%iter_mstep + 1
         work%aborted_mstep = .true. ! set to .false. at end of iteration
         if( work%iter_mstep == 1 ) then
            if( compute_logP_score_hess_mstep_rrlogit( &
                 work, err, .false. ) == RETURN_FAIL ) goto 600
         else
            if( compute_logP_score_hess_mstep_rrlogit( &
                 work, err, .true. ) == RETURN_FAIL ) goto 600
         end if
         logP = work%loglikA + work%logpriorA
         ! put negative hessA into wkddA and inverse into wkddB
         work%wkddA(:,:) = - work%hessA(:,:)
         if( cholesky_in_place(work%wkddA, err ) &
              == RETURN_FAIL ) exit
         if( invert_lower(work%wkddA, err ) == RETURN_FAIL ) exit
         if( premult_lower_by_transpose( work%wkddA, &
              work%wkddB, err) == RETURN_FAIL ) exit
         work%wkdC = matmul( work%wkddB, work%scoreA )
         ! save current params
         work%pi_mat_mstep_old(:,:) = work%pi_mat_mstep(:,:)
         work%beta_mstep_old(:) = work%beta_mstep(:)
         ! attempt full step
         do_step_halving = .false.
         work%beta_mstep = work%beta_mstep + work%wkdC
         if( compute_pi_mat_mstep_rrlogit( work, err, logP=.true. ) &
              == RETURN_FAIL ) goto 800
         logP_new = work%loglikA_new + work%logpriorA_new
         if( logP_new < logP ) do_step_halving = .true.
         if( work%any_zero_pi_mstep ) do_step_halving = .true.
         if( do_step_halving ) then
            work%step_halving = 0
            do
               work%step_halving = work%step_halving + 1
               work%beta_mstep = work%beta_mstep_old + work%wkdC * &
                    0.5D0**work%step_halving
               if( compute_pi_mat_mstep_rrlogit( work, err, logP=.true. ) &
                    == RETURN_FAIL ) goto 800
               logP_new = work%loglikA_new + work%logpriorA_new
               if( logP_new >= logP ) exit
            end do
         end if
         max_abs_diff = 0.D0
         do cp = 1, work%n_cov_patt
            do y = 1, work%r
               max_abs_diff = max( max_abs_diff, &
                    abs( work%pi_mat_mstep(cp,y) - &
                    work%pi_mat_mstep_old(cp,y) ) )
            end do
         end do
         if( max_abs_diff < work%crit_converge ) work%converged_mstep = .true.
         work%aborted_mstep = .false.
      end do
600   continue
      if( work%aborted_mstep ) then
         call err_handle(err, 1, &
              comment = "M-step aborted" )
         call err_handle(err, 5, iiter = int(work%iter_mstep, kind=our_int) )
      end if
      if( .not. work%converged_mstep ) then
         call err_handle(err, 1, &
              comment = "M-step failed to converge by" )
         call err_handle(err, 5, iiter = int(work%iter_mstep, kind=our_int) )
      end if
      !####
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function mstep_rrlogit
   !###################################################################
   integer(kind=our_int) function compute_logP_score_hess_mstep_rrlogit( &
        work, err, skip_compute_pi_mat_mstep ) result(answer)
      ! computes loglikA, first derivative and Hessian at 
      ! current value of beta_mstep
      ! stores derivatives in scoreA and hessA
      implicit none
      ! workspaces
      type(workspace_type_rrlogit), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! optionals
      logical, intent(in), optional :: skip_compute_pi_mat_mstep
      ! locals
      logical :: skip_compute_pi_mat_mstep_local, done
      integer(kind=our_int) :: cp, y, j, k, l, m, posn, posnA
      real(kind=our_dble) :: q1, q2, q3, rtmp
      character(len=*), parameter :: &
           subname = "compute_logP_score_hess_mstep_rrlogit"
      ! begin
      answer = RETURN_FAIL
      if( present( skip_compute_pi_mat_mstep ) ) then
         skip_compute_pi_mat_mstep_local = skip_compute_pi_mat_mstep
      else
         skip_compute_pi_mat_mstep_local = .false.
      end if
      if( .not. skip_compute_pi_mat_mstep_local ) then
         if( compute_pi_mat_mstep_rrlogit( work, err ) &
              == RETURN_FAIL ) goto 800
      end if
      !
      work%loglikA = 0.D0
      work%logpriorA = 0.D0
      work%scoreA(:) = 0.D0
      work%hessA(:,:) = 0.D0
      do cp = 1, work%n_cov_patt
         if( work%empty_cov_patt(cp) ) cycle
         ! accumulate loglikA and logpriorA
         do y = 1, work%r
            if( ( work%fhat_mat(cp,y) > 0.D0 ) .or. &
                 (work%prior == "DAP" )  ) then
               if( work%pi_mat_mstep(cp,y) <= 0.D0 ) goto 200
               rtmp = log( work%pi_mat_mstep(cp,y) )
               if( work%prior == "DAP" ) work%logpriorA = &
                    work%logpriorA + work%prior_freq(cp,y) * rtmp
               if( work%fhat_mat(cp,y) > 0.D0 ) work%loglikA = &
                    work%loglikA + work%fhat_mat(cp,y) * rtmp
            end if
         end do
         ! accumulate scoreA and lower triangle of hess
         posn = 0
         do k = 1, work%r
            if( k == work%baseline ) cycle
            if( work%prior == "DAP" ) then
               q1 = ( work%prior_freq(cp,k) + work%fhat_mat(cp,k) ) &
                    - ( work%prior_freq_tot_cov_patt(cp) + &
                    work%fitweight_cov_patt(cp) ) * &
                    work%pi_mat_mstep(cp,k)
               q2 = ( work%prior_freq_tot_cov_patt(cp) + &
                    work%fitweight_cov_patt(cp) ) &
                    * work%pi_mat_mstep(cp,k) * &
                    ( 1.D0 - work%pi_mat_mstep(cp,k) )
            else
               q1 = work%fhat_mat(cp,k) - work%fitweight_cov_patt(cp) * &
                    work%pi_mat_mstep(cp,k)
               q2 = work%fitweight_cov_patt(cp) * work%pi_mat_mstep(cp,k) * &
                    ( 1.D0 - work%pi_mat_mstep(cp,k) )
            end if
            do j = 1, work%p
               posn = posn + 1
               work%scoreA(posn) = work%scoreA(posn) &
                    + q1 * work%model_matrix(cp,j)
               posnA = 0
               done = .false.
               do m = 1, work%r
                  if( m == work%baseline ) cycle
                  if( work%prior == "DAP" ) then
                     q3 = ( work%prior_freq_tot_cov_patt(cp) + &
                          work%fitweight_cov_patt(cp) ) * &
                          work%pi_mat_mstep(cp,k) * work%pi_mat_mstep(cp,m)
                  else
                     q3 = work%fitweight_cov_patt(cp) * &
                          work%pi_mat_mstep(cp,k) * work%pi_mat_mstep(cp,m)
                  end if
                  do l = 1, work%p
                     posnA = posnA + 1
                     if( posnA > posn ) then
                        done = .true.
                        exit
                     end if
                     if( k == m ) then
                        work%hessA(posn,posnA) = work%hessA(posn,posnA) &
                             - q2 * work%model_matrix(cp,j) * &
                             work%model_matrix(cp,l)
                     else
                        work%hessA(posn,posnA) = work%hessA(posn,posnA) &
                             + q3 * work%model_matrix(cp,j) * &
                             work%model_matrix(cp,l)
                     end if
                  end do
                  if( done ) exit
               end do
            end do
         end do
      end do
      ! fill in upper triangle of hessA
      do posn = 1, work%d
         do posnA = posn + 1, work%d
            work%hessA( posn, posnA ) = work%hessA( posnA, posn )
         end do
      end do
      !####
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
200   call err_handle(err, 1, &
            comment = "Attempted logarithm of non-positive number" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function compute_logP_score_hess_mstep_rrlogit
   !###################################################################
   integer(kind=our_int) function compute_pi_mat_mstep_rrlogit( &
        work, err, logP ) result(answer)
      ! computes pi_mat_mstep from beta_mstep
      ! If logP is .true., computes logprior and augmented-data
      ! loglikelihood, putting them into logpriorA_new and loglikA_new
      implicit none
      ! workspaces
      type(workspace_type_rrlogit), intent(inout) :: work
      type(error_type), intent(inout) :: err
      logical, intent(in), optional :: logP
      ! locals
      logical :: logP_local
      integer(kind=our_int) :: cov_patt, j, posn, k
      real(kind=our_dble) :: sum, eta_max
      character(len=*), parameter :: &
           subname = "compute_pi_mat_mstep_rrlogit"
      ! begin
      answer = RETURN_FAIL
      if( present(logP) ) then
         logP_local = logP
      else
         logP_local = .false.
      end if
      work%any_zero_pi_mstep = .false.
      if( logP_local ) then
         work%logpriorA_new = 0.D0
         work%loglikA_new = 0.D0
      end if
      do cov_patt = 1, work%n_cov_patt
         ! put linear predictors (logits) into wkrA
         posn = 0
         do j = 1, work%r
            if( j == work%baseline ) then
               work%wkrA(j) = 0.D0
               cycle
            end if
            sum = 0.D0
            do k = 1, work%p
               posn = posn + 1
               sum = sum + work%model_matrix(cov_patt,k) * &
                    work%beta_mstep(posn)
            end do
            work%wkrA(j) = sum
         end do
         ! exponentiate and normalize, being careful to avoid overflow
         ! and underflow; put result into wkrB
         eta_max = log_tiny
         do j = 1, work%r
            eta_max = max(eta_max, work%wkrA(j) )
         end do
         work%wkrA(:) = work%wkrA(:) - eta_max
         sum = 0.D0
         do j = 1, work%r
            if( work%wkrA(j) < log_tiny ) then
               work%wkrB(j) = 0.D0
            else if( work%wkrA(j) > log_huge ) then
               goto 100
            else
               work%wkrB(j) = exp( work%wkrA(j) )
            end if
            sum = sum + work%wkrB(j)
         end do
         if( sum == 0.D0 ) goto 200
         do j = 1, work%r
            work%pi_mat_mstep(cov_patt,j) = work%wkrB(j) / sum
            if( work%pi_mat_mstep(cov_patt,j) == 0.D0 ) &
                 work%any_zero_pi_mstep = .true.
         end do
         !
         if( logP_local ) then
            ! accumulate logpriorA_new and loglikA_new
            if( work%prior == "DAP" ) then
               do j = 1, work%r
                  if( work%logpriorA_new == work%neginfcode ) exit
                  if( work%prior_freq(cov_patt, j) == 0.D0 ) cycle
                  if( work%pi_mat_mstep(cov_patt, j) <= 0.D0 ) then
                     work%logpriorA_new = work%neginfcode
                  else
                     work%logpriorA_new = work%logpriorA_new + &
                          work%prior_freq(cov_patt, j) * &
                          log( work%pi_mat_mstep(cov_patt, j) )
                  end if
               end do
               do j = 1, work%r
                  if( work%loglikA_new == work%neginfcode ) exit
                  if( work%fhat_mat(cov_patt, j) == 0.D0 ) cycle
                  if( work%pi_mat_mstep(cov_patt, j) <= 0.D0 ) then
                     work%loglikA_new = work%neginfcode
                  else
                     work%loglikA_new = work%loglikA_new + &
                          work%fhat_mat(cov_patt, j) * &
                          log( work%pi_mat_mstep(cov_patt, j) )
                  end if
               end do
            end if
         end if
      end do
      !####
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
100   call err_handle(err, 1, &
            comment = "Overflow; fitted value became too large" )
      goto 800
200   call err_handle(err, 1, &
            comment = "Attempted division by zero" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
   end function compute_pi_mat_mstep_rrlogit
   !###################################################################
   integer(kind=our_int) function check_convergence_rrlogit( &
        work ) result(answer)
      ! compares pi_mat and pi_mat_old
      implicit none
      ! workspaces
      type(workspace_type_rrlogit), intent(inout) :: work
!      type(error_type), intent(inout) :: err
      ! locals
      integer(kind=our_int) :: cp, y
      real(kind=our_dble) :: max_abs_diff, smallest_pi
      character(len=*), parameter :: &
           subname = "check_convergence_rrlogit"
      ! begin
      answer = RETURN_FAIL
      ! compute pi_mat from beta_vec
      max_abs_diff = 0.D0
      smallest_pi = 1.D0
      do cp = 1, work%n_cov_patt
         do y = 1, work%r
            max_abs_diff = max( max_abs_diff, &
                 abs( work%pi_mat(cp,y) - work%pi_mat_old(cp,y) ) )
            smallest_pi = min( smallest_pi, work%pi_mat(cp,y) )
         end do
      end do
      work%converged = .false.
      if( max_abs_diff < work%crit_converge ) work%converged = .true.
      work%boundary = .false.
      if( smallest_pi < work%crit_boundary ) work%boundary = .true.
      !####
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
!800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
!      goto 999
      ! cleanup
999   continue
   end function check_convergence_rrlogit
   !###################################################################
   integer(kind=our_int) function run_rrlogit_em_saturated( work, err ) &
        result(answer)
     ! fits the saturated model, which estimates pi separately within
     ! each covariate pattern
     ! skips the empty covariate patterns, filling the corresponding
     ! rows of pi_mat with missing value codes
     implicit none
      ! workspaces
      type(workspace_type_rrlogit), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! locals
      real(kind=our_dble) :: loglik, logprior, sum, max_abs_diff, &
           smallest_pi
      integer(kind=our_int) :: cp, dp, ystar, y
      character(len=*), parameter :: &
           subname = "run_rrlogit_em_saturated"
      ! begin
      answer = RETURN_FAIL
      work%loglik = 0.D0    ! for all covariate patterns
      work%logprior = 0.D0  ! for all covariate patterns
      work%converged = .true. ! .false. if not converged for any cov patt
      work%boundary = .false. ! set to .true. if boundary for any cov patt
      work%iter = 0 ! set to max across cov patts
      do cp = 1, work%n_cov_patt
         if( work%empty_cov_patt(cp) ) then
            work%pi_mat(cp,:) = work%mvcode
            work%pistar_mat(cp,:) = work%mvcode
            cycle
         end if
         work%pi_vec(:) = work%pi_mat(cp,:) ! starting values
         work%iter_em_saturated(cp) = 0
         work%converged_em_saturated(cp) = .false.
         work%boundary_em_saturated(cp) = .false.
         do
            if( work%converged_em_saturated(cp) .or. &
                 ( work%iter_em_saturated(cp) >= work%iter_max_em ) ) exit
            work%iter_em_saturated(cp) = work%iter_em_saturated(cp) + 1
            work%aborted = .true. ! set to .false. at end of iteration
            ! compute pistar_vec, logprior, loglik for this cov patt at pi_vec
            ! along with the active rows of phi_mat
            if( work%pert_mat_identity ) then
               work%pistar_vec(:) = work%pi_vec(:)
            else
               do ystar = 1, work%r
                  sum = 0.D0
                  do y = 1, work%r
                     sum = sum + work%pert_mat(ystar,y) * work%pi_vec(y)
                  end do
                  work%pistar_vec(ystar) = sum
               end do
            end if
            logprior = 0.D0
            if( work%prior == "DAP" ) then
               do y = 1, work%r
                  if( work%prior_freq(cp,y) == 0.D0 ) cycle
                  if( work%pi_vec(y) <= 0.D0 ) then
                     call err_handle(err, 1, &
                       comment = "Attempted logarithm of non-positive number" )
                     goto 10
               end if
                  logprior = logprior + work%prior_freq(cp,y) * &
                       log( work%pi_vec(y) )
               end do
            end if
            loglik = 0.D0
            do dp = work%data_patt_st(cp), work%data_patt_fin(cp)
               ystar = work%response_for_data_patt(dp)
               if( .not. work%active_ystar(cp, ystar) ) cycle
               if( work%pistar_vec(ystar) <= 0.D0 ) then
                  call err_handle(err, 1, &
                       comment = "Attempted logarithm of non-positive number" )
                  goto 10
               end if
               loglik = loglik + work%fitweight_data_patt(dp) * &
                    log( work%pistar_vec(ystar) )
               if( .not. work%pert_mat_identity ) then
                  sum = 0.D0
                  do y = 1, work%r
                     work%phi_mat(ystar,y) = work%pert_mat(ystar,y) * &
                          work%pi_vec(y)
                     sum = sum + work%phi_mat(ystar,y)
                  end do
                  do y = 1, work%r
                     work%phi_mat(ystar,y) = work%phi_mat(ystar,y) / sum
                  end do
               end if
            end do
            ! do estep, putting expected true frequencies into fhat_vec
            if( work%pert_mat_identity ) then
               work%fhat_vec(:) = 0.D0
               do dp = work%data_patt_st(cp), work%data_patt_fin(cp)
                  y = work%response_for_data_patt(dp)
                  if( .not. work%active_ystar(cp, y) ) cycle
                  work%fhat_vec(y) = work%fitweight_data_patt(dp)
               end do
            else
               do y = 1, work%r
                  sum = 0.D0
                  do dp = work%data_patt_st(cp), work%data_patt_fin(cp)
                     ystar = work%response_for_data_patt(dp)
                     if( .not. work%active_ystar(cp, ystar) ) cycle
                     sum = sum + work%fitweight_data_patt(dp) * &
                          work%phi_mat(ystar,y)
                  end do
                  work%fhat_vec(y) = sum
               end do
            end if
            ! do m-step, updating pi_vec
            work%pi_vec_old(:) = work%pi_vec(:)
            sum = 0.D0
            do y = 1, work%r
               work%pi_vec(y) = work%fhat_vec(y)
               if( work%prior == "DAP" ) work%pi_vec(y) = work%pi_vec(y) &
                    + work%prior_freq(cp,y)
               sum = sum + work%pi_vec(y)
            end do
            work%pi_vec(:) = work%pi_vec(:) / sum
            ! check convergence
            max_abs_diff = 0.D0
            smallest_pi = 1.D0
            do y = 1, work%r
               max_abs_diff = max( max_abs_diff, &
                    abs( work%pi_vec(y) - work%pi_vec_old(y) ) )
               smallest_pi = min( smallest_pi, work%pi_vec(y) )
            end do
            if( max_abs_diff < work%crit_converge ) &
                 work%converged_em_saturated(cp) = .true.
            if( smallest_pi < work%crit_boundary ) &
                 work%boundary_em_saturated(cp) = .true.
            work%aborted = .false.
         end do
         if( .not. work%converged_em_saturated(cp) ) work%converged = .false.
         if( work%boundary_em_saturated(cp) ) work%boundary = .true.
         ! update pistar, logpri, loglik
         if( work%pert_mat_identity ) then
            work%pistar_vec(:) = work%pi_vec(:)
         else
            do ystar = 1, work%r
               sum = 0.D0
               do y = 1, work%r
                  sum = sum + work%pert_mat(ystar,y) * work%pi_vec(y)
               end do
               work%pistar_vec(ystar) = sum
            end do
         end if
         logprior = 0.D0
         if( work%prior == "DAP" ) then
            do y = 1, work%r
               if( work%prior_freq(cp,y) == 0.D0 ) cycle
               if( work%pi_vec(y) <= 0.D0 ) then
                  call err_handle(err, 1, &
                       comment = "Attempted logarithm of non-positive number" )
                  goto 10
               end if
               logprior = logprior + work%prior_freq(cp,y) * &
                    log( work%pi_vec(y) )
            end do
         end if
         loglik = 0.D0
         do dp = work%data_patt_st(cp), work%data_patt_fin(cp)
            ystar = work%response_for_data_patt(dp)
            if( .not. work%active_ystar(cp, ystar) ) cycle
            if( work%pistar_vec(ystar) <= 0.D0 ) then
                call err_handle(err, 1, &
                     comment = "Attempted logarithm of non-positive number" )
               goto 10
            end if
            loglik = loglik + work%fitweight_data_patt(dp) * &
                 log( work%pistar_vec(ystar) )
         end do
         work%logprior = work%logprior + logprior
         work%loglik = work%loglik + loglik
         work%pi_mat(cp,:) = work%pi_vec(:)
         work%pistar_mat(cp,:) = work%pistar_vec(:)
         work%iter = max( work%iter, work%iter_em_saturated(cp) )
      end do
10    continue
      if( work%aborted ) then
         call err_handle(err, 1, &
              comment = "EM algorithm aborted" )
         call err_handle(err, 5, iiter = &
              int(work%iter_em_saturated(cp), kind=our_int) )
         call err_handle(err, 2, whichsub = subname, whichmod = modname )
      end if
      work%logP = work%loglik + work%logprior
      !####
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
!800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      !goto 999
      ! cleanup
999   continue
   end function run_rrlogit_em_saturated
   !###################################################################
   integer(kind=our_int) function run_rrlogit_da_saturated( &
        start_logP, logP_series, imp_mat_series, imp_vec_series, &
        pi_marg_series, &
        n_iter_actual, n_sample_actual, n_imp_actual, &
        work, err ) result(answer)
      ! runs mcmc on the saturated model, which estimates pi separately
      !   within each covariate pattern
      ! running means are stored in work%pi_mat_mean and
      !   work%pistar_mat_mean 
      ! skips the empty covariate patterns, filling the corresponding
      ! rows of work%pi_mat_mean and work%pistar_mat_mean with missing
      !   value codes
      implicit none
      ! outputs
      real(kind=our_dble), intent(out) :: start_logP
      real(kind=our_dble), intent(out) :: logP_series(:)
      integer(kind=our_int), intent(out) :: imp_mat_series(:,:,:)
      integer(kind=our_int), intent(out) :: imp_vec_series(:,:)
      real(kind=our_dble), intent(out) :: pi_marg_series(:,:)
      integer(kind=our_int), intent(out) :: n_iter_actual
      integer(kind=our_int), intent(out) :: n_sample_actual
      integer(kind=our_int), intent(out) :: n_imp_actual
      ! workspaces
      type(workspace_type_rrlogit), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! locals
      integer(kind=our_int) :: iter, thin_step, imp_step, itmp
      integer(kind=our_int) :: status, ijunk
      integer(kind=our_int), allocatable :: harvest(:)
      real(kind=our_dble) :: loglik, logprior, sum, rtmpA, rtmpB, nfloat, rtmp
      integer(kind=our_int) :: cp, dp, ystar, y
      character(len=*), parameter :: &
           subname = "run_rrlogit_da_saturated"
      ! begin
      answer = RETURN_FAIL
      if( .not. work%saturated ) goto 20
      if( work%method /= "MCMC" ) goto 30
      allocate( harvest( work%r ), stat=status )
      if( status /= 0 ) goto 100
      ! initialize iteration counters
      work%iter = 0
      work%iter_past_burn_in = 0
      work%store_count = 0
      work%imp_count = 0
      thin_step = 0
      imp_step = 0
      n_iter_actual = 0
      n_sample_actual = 0
      n_imp_actual = 0
      ! initialize running means and returned series
      start_logP = 0.D0
      do cp = 1, work%n_cov_patt
         if( work%empty_cov_patt(cp) ) then
            work%pi_mat(cp,:) = work%mvcode
            work%pistar_mat(cp,:) = work%mvcode
            work%pi_mat_mean(cp,:) = work%mvcode
            work%pistar_mat_mean(cp,:) = work%mvcode
         else
            work%pi_mat_mean(cp,:) = 0.D0
            work%pistar_mat_mean(cp,:) = 0.D0
         end if
      end do
      work%pi_marg_mean(:) = 0.D0
      logP_series(:) = work%mvcode
      imp_mat_series(:,:,:) = -1
      imp_vec_series(:,:) = -1
      work%fimp_mat(:,:) = 0
      work%aborted = .false.
      ! main iteration
      do iter = 1, work%iter_mcmc + work%burn_mcmc
         work%aborted = .true.  ! will be set to .false. at end of cycle
         work%iter = iter
         work%iter_past_burn_in = work%iter - work%burn_mcmc ! could be neg
         work%store_this_iter = .false.
         work%imp_this_iter = .false.
         if( work%iter_past_burn_in > 0 ) then
            thin_step = thin_step + 1
            if( thin_step == work%thin_mcmc ) then
               work%store_this_iter = .true.
               work%store_count = work%store_count + 1
               thin_step = 0
            end if
            imp_step = imp_step + 1
            if( ( imp_step == work%impute_every ) .and. &
                 ( work%impute_every /= 0 ) ) then
               work%imp_this_iter = .true.
               work%imp_count = work%imp_count + 1
               imp_step = 0
            end if
         end if
         !#######################
         ! run I-step for each covariate pattern
         work%loglik = 0.D0   ! for all cov patts
         work%logprior = 0.D0 ! for all cov patts
         do cp = 1, work%n_cov_patt
            if( work%empty_cov_patt(cp) ) cycle
            work%pi_vec(:) = work%pi_mat(cp,:) ! current params
            ! compute pistar_vec, logprior, loglik for this cov patt at pi_vec
            ! along with the active rows of phi_mat
            if( work%pert_mat_identity ) then
               work%pistar_vec(:) = work%pi_vec(:)
            else
               do ystar = 1, work%r
                  sum = 0.D0
                  do y = 1, work%r
                     sum = sum + work%pert_mat(ystar,y) * work%pi_vec(y)
                  end do
                  work%pistar_vec(ystar) = sum
               end do
            end if
            work%pistar_mat(cp,:) = work%pistar_vec(:)
            logprior = 0.D0
            if( work%prior == "DAP" ) then
               do y = 1, work%r
                  if( work%prior_freq(cp,y) == 0.D0 ) cycle
                  if( work%pi_vec(y) <= 0.D0 ) then
                     call err_handle(err, 1, &
                       comment = "Attempted logarithm of non-positive number" )
                     goto 10
                  end if
                  logprior = logprior + work%prior_freq(cp,y) * &
                       log( work%pi_vec(y) )
               end do
            end if
            work%logprior = work%logprior + logprior
            loglik = 0.D0
            do dp = work%data_patt_st(cp), work%data_patt_fin(cp)
               ystar = work%response_for_data_patt(dp)
               if( .not. work%active_ystar(cp, ystar) ) cycle
               if( work%pistar_vec(ystar) <= 0.D0 ) then
                  call err_handle(err, 1, &
                       comment = "Attempted logarithm of non-positive number" )
                  goto 10
               end if
               loglik = loglik + &
                    real(work%freq_for_data_patt_int(dp), our_dble) * &
                    log( work%pistar_vec(ystar) )
               if( .not. work%pert_mat_identity ) then
                  sum = 0.D0
                  do y = 1, work%r
                     work%phi_mat(ystar,y) = work%pert_mat(ystar,y) * &
                          work%pi_vec(y)
                     sum = sum + work%phi_mat(ystar,y)
                  end do
                  do y = 1, work%r
                     work%phi_mat(ystar,y) = work%phi_mat(ystar,y) / sum
                  end do
               end if
            end do
            work%loglik = work%loglik + loglik
            ! impute true frequencies, putting result into fimp_vec
            if( work%pert_mat_identity ) then
               work%fimp_vec(:) = 0
               do dp = work%data_patt_st(cp), work%data_patt_fin(cp)
                  y = work%response_for_data_patt(dp)
                  if( .not. work%active_ystar(cp, y) ) cycle
                  work%fimp_vec(y) = work%freq_for_data_patt_int(dp)
               end do
            else
               work%fimp_vec(:) = 0
               do dp = work%data_patt_st(cp), work%data_patt_fin(cp)
                  ystar = work%response_for_data_patt(dp)
                  if( .not. work%active_ystar(cp, ystar) ) cycle
                  itmp = work%freq_for_data_patt_int(dp)
                  work%wkrA(:) = work%phi_mat(ystar,:)
                  if( rmult_n( itmp, work%wkrA, harvest, err ) &
                       == RETURN_FAIL ) goto 10
                  work%fimp_vec(:) = work%fimp_vec(:) + harvest(:)
               end do
            end if
            work%fimp_mat(cp,:) = work%fimp_vec(:)
         end do  ! end of I-step cycle
         work%logP = work%logprior + work%loglik
         if( work%iter == 1 ) start_logP = work%logP
         !#######################
         ! run P-step for each covariate pattern
         do cp = 1, work%n_cov_patt
            if( work%empty_cov_patt(cp) ) cycle
            sum = 0.D0
            do y = 1, work%r
               rtmpA = work%prior_freq(cp,y) + &
                    real( work%fimp_mat(cp,y), our_dble )
               if( rgamma_R( rtmpA + 1.D0, 1.D0, rtmpB, err ) &
                    == RETURN_FAIL ) goto 10
               work%pi_vec(y) = rtmpB
               sum = sum + rtmpB
            end do
            if( sum == 0.D0 ) then
               call err_handle(err, 1, &
                    comment = "Attempted division by zero" )
               goto 10
            end if
            work%pi_vec(:) = work%pi_vec(:) / sum
            work%pi_mat(cp,:) = work%pi_vec(:)
         end do  ! end of P-step cycle
         if( compute_pi_marg( work, err ) == RETURN_FAIL ) goto 10
         !#######################
         if( work%store_this_iter ) then
            ! store results in series
            if( work%store_count < 0 ) goto 300
            if( work%store_count > size(logp_series, kind=our_int) ) goto 300
            logP_series( work%store_count ) = work%logP
            pi_marg_series( work%store_count, : ) = work%pi_marg(:)
         end if
         !#######################
         if( work%iter_past_burn_in > 0 ) then
            ! update running means in work%pi_mat_mean and work%pi_marg_mean
            nfloat = real( work%iter_past_burn_in, our_dble )
            if( nfloat == 0.D0 ) goto 250
            do cp = 1, work%n_cov_patt
               if( work%empty_cov_patt(cp) ) cycle
               do y = 1, work%r
                  rtmp = work%pi_mat(cp,y) - work%pi_mat_mean(cp,y)
                  work%pi_mat_mean(cp,y) = work%pi_mat_mean(cp,y) &
                       + rtmp / nfloat
               end do
            end do
            do y = 1, work%r
               rtmp = work%pi_marg(y) - work%pi_marg_mean(y)
               work%pi_marg_mean(y) = work%pi_marg_mean(y) &
                    + rtmp / nfloat
            end do
         end if
         !#######################
         if( work%imp_this_iter ) then
            if( work%imp_count < 0 ) goto 300
            if( work%imp_count > size(imp_vec_series, 2, kind=our_int) ) &
                 goto 300
            if( work%imp_count > size(imp_mat_series, 3, kind=our_int) ) &
                 goto 300
            if( run_rrlogit_impute_random( work%wide_format, &
                 work%micro_data, work%freq_row_input_data_int, &
                 work%response_row_input_data, work%f_mat_int, &
                 imp_mat_series(:,:,work%imp_count), &
                 imp_vec_series(:,work%imp_count), &
                 work, err ) == RETURN_FAIL ) goto 800
         end if
         !#######################
         work%aborted = .false.
      end do
10    continue
      if( work%aborted ) then
         !#### issue warning message and continue
         call err_handle(err, 1, &
              comment = "MCMC procedure aborted" )
         call err_handle(err, 5, iiter = int(work%iter, kind=our_int) )
         call err_handle(err, 2, whichsub = subname, whichmod = modname )
      end if
      !#####################
      ! compute pistar_mat_mean and pistar_mat
      if( work%pert_mat_identity ) then
         work%pistar_mat_mean(:,:) = work%pi_mat_mean(:,:)
      else
         do cp = 1, work%n_cov_patt
            if( work%empty_cov_patt(cp) ) cycle
            do ystar = 1, work%r
               sum = 0.D0
               do y = 1, work%r
                  sum = sum + work%pert_mat(ystar,y) * work%pi_mat_mean(cp,y)
               end do
               work%pistar_mat_mean(cp,ystar) = sum
               sum = 0.D0
               do y = 1, work%r
                  sum = sum + work%pert_mat(ystar,y) * work%pi_mat(cp,y)
               end do
               work%pistar_mat(cp,ystar) = sum
            end do
        end do
      end if
      !####
      if( compute_loglik_logprior_rrlogit( work, err, &
           logprior = logprior, &
           loglik = loglik ) == RETURN_FAIL ) goto 10
      work%loglik = loglik
      work%logprior = logprior
      work%logP = loglik + logprior
      !####
      n_iter_actual = work%iter
      n_sample_actual = work%store_count
      n_imp_actual = work%imp_count
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
20    call err_handle(err, 1, &
            comment = "The model is not saturated" )
      goto 800
30    call err_handle(err, 1, &
            comment = "Not prepared for MCMC" )
      goto 800
100   call err_handle(err, 1, &
            comment = "Unable to allocate array" )
      goto 800
250   call err_handle(err, 1, &
            comment = "Attempted division by zero" )
      goto 800
300   call err_handle(err, 1, &
            comment = "Array bounds exceeded" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
      if( allocated( harvest ) ) deallocate( harvest, stat=ijunk )
   end function run_rrlogit_da_saturated
   !###################################################################
   integer(kind=our_int) function run_rrlogit_da_nonsaturated( &
        start_logP, logP_series, coef_vec_series, &
        imp_mat_series, imp_vec_series, pi_marg_series, &
        n_iter_actual, n_sample_actual, n_imp_actual, &
        work, err ) result(answer)
      ! runs mcmc on the nonsaturated model
      ! running means are stored in work%beta_mean, work%pi_mat_mean,
      !   work%pistar_mat_mean; running cov matrix in work%beta_cov_mat
      implicit none
      ! outputs
      real(kind=our_dble), intent(out) :: start_logP
      real(kind=our_dble), intent(out) :: logP_series(:)
      real(kind=our_dble), intent(out) :: coef_vec_series(:,:)
      integer(kind=our_int), intent(out) :: imp_mat_series(:,:,:)
      integer(kind=our_int), intent(out) :: imp_vec_series(:,:)
      real(kind=our_dble), intent(out) :: pi_marg_series(:,:)
      integer(kind=our_int), intent(out) :: n_iter_actual
      integer(kind=our_int), intent(out) :: n_sample_actual
      integer(kind=our_int), intent(out) :: n_imp_actual
      ! workspaces
      type(workspace_type_rrlogit), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! locals
      integer(kind=our_int), allocatable :: harvest(:)
      integer(kind=our_int) :: iter, thin_step, imp_step, itmp
      integer(kind=our_int) :: status, ijunk
      real(kind=our_dble) :: loglik, logprior, sum, nfloat, &
           rtmp, log_mh_ratio, mh_ratio, runif, loglik_tmp, logprior_tmp
      logical :: accept
      integer(kind=our_int) :: cp, dp, ystar, y, j, k
      character(len=*), parameter :: &
           subname = "run_rrlogit_da_nonsaturated"
      ! begin
      answer = RETURN_FAIL
      if( work%saturated ) goto 20
      if( work%method /= "MCMC" ) goto 30
      allocate( harvest( work%r ), stat=status )
      if( status /= 0 ) goto 100
      ! initialize iteration counters
      work%iter = 0
      work%iter_past_burn_in = 0
      work%store_count = 0
      work%imp_count = 0
      thin_step = 0
      imp_step = 0
      n_iter_actual = 0
      n_sample_actual = 0
      n_imp_actual = 0
      ! initialize running means and returned series
      start_logP = 0.D0
      work%beta_mean(:) = 0.D0
      work%beta_cov_mat(:,:) = 0.D0
      work%pi_mat_mean(:,:) = 0.D0
      work%pistar_mat_mean(:,:) = 0.D0
      work%pi_marg_mean(:) = 0.D0
      logP_series(:) = work%mvcode
      coef_vec_series(:,:) = work%mvcode
      imp_mat_series(:,:,:) = -1
      imp_vec_series(:,:) = -1
      work%fimp_mat(:,:) = 0
      ! initialize MH diagnostics
      work%beta_accept_count = 0
      work%beta_current_reject_run = 0
      work%beta_accept_rate = 0.D0
      !
      work%aborted = .false.
      accept = .false.
      ! main iteration
      do iter = 1, work%iter_mcmc + work%burn_mcmc
         work%aborted = .true.  ! will be set to .false. at end of cycle
         work%iter = iter
         work%iter_past_burn_in = work%iter - work%burn_mcmc ! could be neg
         work%store_this_iter = .false.
         work%imp_this_iter = .false.
         if( work%iter_past_burn_in > 0 ) then
            thin_step = thin_step + 1
            if( thin_step == work%thin_mcmc ) then
               work%store_this_iter = .true.
               work%store_count = work%store_count + 1
               thin_step = 0
            end if
            imp_step = imp_step + 1
            if( ( imp_step == work%impute_every ) .and. &
                 ( work%impute_every /= 0 ) ) then
               work%imp_this_iter = .true.
               work%imp_count = work%imp_count + 1
               imp_step = 0
            end if
         end if
         !#######################
         ! run I-step for each covariate pattern
         if( ( iter == 1 ) .or. ( ( iter > 1 ) .and. accept ) ) then 
            if( compute_pi_mat_rrlogit( work, err ) == RETURN_FAIL ) goto 10
            if( compute_pi_marg( work, err ) == RETURN_FAIL ) goto 10
         end if
         work%loglik = 0.D0   ! for all cov patts
         work%logprior = 0.D0 ! for all cov patts
         do cp = 1, work%n_cov_patt
            if( work%empty_cov_patt(cp) ) cycle
            work%pi_vec(:) = work%pi_mat(cp,:) ! current params
            ! compute pistar_vec, logprior, loglik for this cov patt at pi_vec
            ! along with the active rows of phi_mat
            if( work%pert_mat_identity ) then
               work%pistar_vec(:) = work%pi_vec(:)
            else
               do ystar = 1, work%r
                  sum = 0.D0
                  do y = 1, work%r
                     sum = sum + work%pert_mat(ystar,y) * work%pi_vec(y)
                  end do
                  work%pistar_vec(ystar) = sum
               end do
            end if
            work%pistar_mat(cp,:) = work%pistar_vec(:)
            logprior = 0.D0
            if( work%prior == "DAP" ) then
               do y = 1, work%r
                  if( work%prior_freq(cp,y) == 0.D0 ) cycle
                  if( work%pi_vec(y) <= 0.D0 ) then
                     call err_handle(err, 1, &
                       comment = "Attempted logarithm of non-positive number" )
                     goto 10
                  end if
                  logprior = logprior + work%prior_freq(cp,y) * &
                       log( work%pi_vec(y) )
               end do
            end if
            work%logprior = work%logprior + logprior
            loglik = 0.D0
            do dp = work%data_patt_st(cp), work%data_patt_fin(cp)
               ystar = work%response_for_data_patt(dp)
               if( .not. work%active_ystar(cp, ystar) ) cycle
               if( work%pistar_vec(ystar) <= 0.D0 ) then
                  call err_handle(err, 1, &
                       comment = "Attempted logarithm of non-positive number" )
                  goto 10
               end if
               loglik = loglik + &
                    real(work%freq_for_data_patt_int(dp), our_dble) * &
                    log( work%pistar_vec(ystar) )
               if( .not. work%pert_mat_identity ) then
                  sum = 0.D0
                  do y = 1, work%r
                     work%phi_mat(ystar,y) = work%pert_mat(ystar,y) * &
                          work%pi_vec(y)
                     sum = sum + work%phi_mat(ystar,y)
                  end do
                  do y = 1, work%r
                     work%phi_mat(ystar,y) = work%phi_mat(ystar,y) / sum
                  end do
               end if
            end do
            work%loglik = work%loglik + loglik
            ! impute true frequencies, putting result into fimp_vec
            if( work%pert_mat_identity ) then
               work%fimp_vec(:) = 0
               do dp = work%data_patt_st(cp), work%data_patt_fin(cp)
                  y = work%response_for_data_patt(dp)
                  if( .not. work%active_ystar(cp, y) ) cycle
                  work%fimp_vec(y) = work%freq_for_data_patt_int(dp)
               end do
            else
               work%fimp_vec(:) = 0
               do dp = work%data_patt_st(cp), work%data_patt_fin(cp)
                  ystar = work%response_for_data_patt(dp)
                  if( .not. work%active_ystar(cp, ystar) ) cycle
                  itmp = work%freq_for_data_patt_int(dp)
                  work%wkrA(:) = work%phi_mat(ystar,:)
                  if( rmult_n( itmp, work%wkrA, harvest, err ) &
                       == RETURN_FAIL ) goto 10
                  work%fimp_vec(:) = work%fimp_vec(:) + harvest(:)
               end do
            end if
            work%fimp_mat(cp,:) = work%fimp_vec(:)
         end do  ! end of I-step cycle
         work%logP = work%logprior + work%loglik
         if( work%iter == 1 ) start_logP = work%logP
         !#######################
         ! run P-step Metropolis-Hastings
         do cp = 1, work%n_cov_patt
            if( work%empty_cov_patt(cp) ) cycle
            do y = 1, work%r
               work%fhat_mat(cp,y) = real( work%fimp_mat(cp,y), our_dble )
            end do
         end do
         work%beta_mstep(:) = work%beta_vec(:)
         work%pi_mat_mstep(:,:) = work%pi_mat(:,:)
         if( compute_logP_score_hess_mstep_rrlogit( &
              work, err, skip_compute_pi_mat_mstep = .true. ) &
              == RETURN_FAIL ) goto 10
         log_mh_ratio = - ( work%loglikA + work%logpriorA )
         if( compute_center_and_scale_beta_da( work%beta_vec, &
              work%df_da, work%step_size_da, work%scale_fac_da, &
              work, err ) == RETURN_FAIL ) goto 10
         if( draw_candidate_beta( work%df_da, &
              work, err ) == RETURN_FAIL ) goto 10
         if( compute_log_proposal_beta( work%beta_can, work%df_da, &
              rtmp, work, err ) == RETURN_FAIL ) goto 10
         log_mh_ratio = log_mh_ratio - rtmp
         work%beta_mstep(:) = work%beta_can(:)
         if( compute_logP_score_hess_mstep_rrlogit( &
              work, err, skip_compute_pi_mat_mstep = .false. ) &
              == RETURN_FAIL ) goto 10
         log_mh_ratio = log_mh_ratio + ( work%loglikA + work%logpriorA )
         if( compute_center_and_scale_beta_da( work%beta_can, &
              work%df_da, work%step_size_da, work%scale_fac_da, &
              work, err ) == RETURN_FAIL ) goto 10
         if( compute_log_proposal_beta( work%beta_vec, work%df_da, &
              rtmp, work, err ) == RETURN_FAIL ) goto 800
         log_mh_ratio = log_mh_ratio + rtmp
         ! to prevent over/underflow
         if( log_mh_ratio > log_huge ) then
            mh_ratio = huge(0.D0)
         else if( log_mh_ratio < log_tiny ) then
            mh_ratio = 0.D0
         else
            mh_ratio = exp( log_mh_ratio )
         end if
         ! compare to uniform variate
         if( runif_R( runif, err ) == RETURN_FAIL ) goto 10
         if( runif <= mh_ratio ) then
            accept = .true.
            work%beta_vec(:) = work%beta_can(:)
            work%pi_mat(:,:) = work%pi_mat_mstep(:,:)
            if( compute_pi_marg( work, err ) == RETURN_FAIL ) goto 10
         else
            accept = .false.
         end if
         work%mh_ratios_beta( work%iter ) = mh_ratio
         work%mh_accept_beta( work%iter ) = accept
         !#######################
         ! update counters and running estimates
         if( accept ) then
            work%beta_accept_count = work%beta_accept_count + 1
            work%beta_current_reject_run = 0
            rtmp = 1.D0
         else
            work%beta_current_reject_run = &
                 work%beta_current_reject_run + 1
            rtmp = 0.D0
         end if
         work%beta_accept_rate = work%beta_accept_rate + &
              ( rtmp - work%beta_accept_rate ) / real( work%iter, our_dble )
         !#######################
         if( work%store_this_iter ) then
            ! store results in series
            if( work%store_count < 0 ) goto 300
            if( work%store_count > size(logp_series, kind=our_int) ) goto 300
            logP_series( work%store_count ) = work%logP
            if( work%store_count > size(coef_vec_series, 1, kind=our_int) ) &
                 goto 300
            coef_vec_series( work%store_count, : ) = work%beta_vec(:)
            if( work%store_count > size(pi_marg_series, 1, kind=our_int) ) &
                 goto 300
            pi_marg_series( work%store_count, : ) = work%pi_marg(:)
         end if
         !#######################
         if( work%iter_past_burn_in > 0 ) then
            ! update running means in work%pi_mat_mean, work%beta_mean,
            ! and work%beta_cov_mat
            nfloat = real( work%iter_past_burn_in, our_dble )
            if( nfloat == 0.D0 ) goto 250
            do cp = 1, work%n_cov_patt
               do y = 1, work%r
                  rtmp = work%pi_mat(cp,y) - work%pi_mat_mean(cp,y)
                  work%pi_mat_mean(cp,y) = work%pi_mat_mean(cp,y) &
                       + rtmp / nfloat
               end do
            end do
            do y = 1, work%r
               rtmp = work%pi_marg(y) - work%pi_marg_mean(y)
               work%pi_marg_mean(y) = work%pi_marg_mean(y) &
                    + rtmp / nfloat
            end do
            do j = 1, work%d
               work%wkdA(j) = work%beta_vec(j) - work%beta_mean(j)
            end do
            do j = 1, work%d
               do k = j, work%d
                  work%wkddA(j,k) = work%wkdA(j) * work%wkdA(k)
               end do
            end do
            do j = 1, work%d
               work%beta_mean(j) = work%beta_mean(j) + work%wkdA(j) / nfloat
               do k = j, work%d
                  work%beta_cov_mat(j,k) = &
                       ( nfloat - 1.D0 ) * work%beta_cov_mat(j,k) + &
                       ( ( nfloat - 1.D0 ) / nfloat ) * work%wkddA(j,k)
                  work%beta_cov_mat(j,k) = work%beta_cov_mat(j,k) / nfloat
                  work%beta_cov_mat(k,j) = work%beta_cov_mat(j,k)
               end do
            end do
         end if
         !#######################
         if( work%imp_this_iter ) then
            if( work%imp_count < 0 ) goto 300
            if( work%imp_count > size(imp_vec_series, 2, kind=our_int) ) &
                 goto 300
            if( work%imp_count > size(imp_mat_series, 3, kind=our_int) ) &
                 goto 300
            if( run_rrlogit_impute_random( work%wide_format, &
                 work%micro_data, work%freq_row_input_data_int, &
                 work%response_row_input_data, work%f_mat_int, &
                 imp_mat_series(:,:,work%imp_count), &
                 imp_vec_series(:,work%imp_count), &
                 work, err, skip_compute_pi_mat = .true. ) &
                 == RETURN_FAIL ) goto 800
         end if
         !#######################
         if( work%beta_current_reject_run > work%stuck_limit ) then
            call err_handle(err, 1, &
                 comment = "Metropolis-Hastings got stuck" )
            goto 10
         end if
         work%aborted = .false.
      end do
10    continue
      if( work%aborted ) then
         !#### issue warning message and continue
         call err_handle(err, 1, &
              comment = "MCMC procedure aborted" )
         call err_handle(err, 5, iiter = int(work%iter, kind=our_int) )
         call err_handle(err, 2, whichsub = subname, whichmod = modname )
      end if
      !#####################
      ! compute pistar_mat_mean and pistar_mat
      if( work%pert_mat_identity ) then
         work%pistar_mat_mean(:,:) = work%pi_mat_mean(:,:)
      else
         do cp = 1, work%n_cov_patt
            do ystar = 1, work%r
               sum = 0.D0
               do y = 1, work%r
                  sum = sum + work%pert_mat(ystar,y) * work%pi_mat_mean(cp,y)
               end do
               work%pistar_mat_mean(cp,ystar) = sum
               sum = 0.D0
               do y = 1, work%r
                  sum = sum + work%pert_mat(ystar,y) * work%pi_mat(cp,y)
               end do
               work%pistar_mat(cp,ystar) = sum
            end do
        end do
      end if
      !####
      if( compute_pi_mat_rrlogit( work, err ) == RETURN_FAIL ) goto 10
      if( compute_loglik_logprior_rrlogit( work, err, &
           logprior = logprior_tmp, &
           loglik = loglik_tmp ) == RETURN_FAIL ) goto 10
      work%loglik = loglik_tmp
      work%logprior = logprior_tmp
      work%logP = loglik_tmp + logprior_tmp
      !####
      n_iter_actual = work%iter
      n_sample_actual = work%store_count
      n_imp_actual = work%imp_count
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
20    call err_handle(err, 1, &
            comment = "The model is saturated" )
      goto 800
30    call err_handle(err, 1, &
            comment = "Not prepared for MCMC" )
      goto 800
100   call err_handle(err, 1, &
            comment = "Unable to allocate array" )
      goto 800
250   call err_handle(err, 1, &
            comment = "Attempted division by zero" )
      goto 800
300   call err_handle(err, 1, &
            comment = "Array bounds exceeded" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
      if( allocated( harvest ) ) deallocate( harvest, stat=ijunk )
   end function run_rrlogit_da_nonsaturated
   !###################################################################
    integer(kind=our_int) function compute_center_and_scale_beta_da( &
        beta, df, step_size, scale_fac, work, err ) result(answer)
      ! computes beta_center and beta_scale for MH proposal given the first
      ! two derivatives in scoreA and hessA
      ! also computes beta_scale_inv, beta_scale_inv_sqrt, and
      ! beta_scale_sqrt
      implicit none
      real(kind=our_dble), intent(in) :: beta(:)
      real(kind=our_dble), intent(in) :: df
      real(kind=our_dble), intent(in) :: step_size
      real(kind=our_dble), intent(in) :: scale_fac
      ! declare workspaces
      type(workspace_type_rrlogit), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      integer(kind=our_int) :: j, k
      real(kind=our_dble) :: rtmp, sum
      character(len=*), parameter :: &
           subname = "compute_center_and_scale_beta_da"
      ! begin
      answer = RETURN_FAIL
      if( work%saturated ) goto 20 
      if( size(beta, kind=our_int) /= work%d ) goto 30
      rtmp = - ( df / &
           ( df + real( work%d, our_dble) ) ) &
                 * ( 1.D0 / ( scale_fac**2 ) )
      work%beta_scale_inv(:,:) = rtmp * work%hessA(:,:)
      work%beta_scale_inv_sqrt(:,:) = work%beta_scale_inv(:,:)
      if( cholesky_in_place( work%beta_scale_inv_sqrt, err ) &
           == RETURN_FAIL ) goto 200
      work%beta_scale_sqrt(:,:) = work%beta_scale_inv_sqrt(:,:)
      if( invert_lower(work%beta_scale_sqrt, err ) == RETURN_FAIL ) goto 200
      if( premult_lower_by_transpose( work%beta_scale_sqrt, &
           work%beta_scale, err) == RETURN_FAIL ) goto 800
      ! put inverse of negative hessian into wkddA
      work%wkddA(:,:) = - rtmp * work%beta_scale(:,:)
      ! find center of proposal
      do j = 1, work%d
         sum = 0.D0
         do k = 1, work%d
            sum = sum + work%wkddA(j,k) * work%scoreA(k)
         end do
         work%beta_center(j) = beta(j) + step_size * sum
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
20    call err_handle(err, 1, &
            comment = "The model is saturated" )
      goto 800
30    call err_handle(err, 1, &
            comment = "Argument beta has incorrect size" )
      goto 800
200   call err_handle(err, 1, &
           comment = "Hessian not neg-def" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function compute_center_and_scale_beta_da
   !###################################################################
    integer(kind=our_int) function draw_candidate_beta( df, &
         work, err ) result(answer)
      ! Draws candidate beta from multivariate t proposal, storing
      ! the result in mmw%beta_can 
      ! Depends on beta_center and beta_scale, so should be run after
      ! compute_center_and_scale_beta_da
      implicit none
      real(kind=our_dble), intent(in) :: df
      ! declare workspaces
      type(workspace_type_rrlogit), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      integer(kind=our_int) :: j, k
      real(kind=our_dble) :: rtmp, sum, rnorm, rchisq
      character(len=*), parameter :: &
           subname = "draw_candidate_beta"
      ! begin
      answer = RETURN_FAIL
      if( work%saturated ) goto 20
      if( df <= 0.D0 ) goto 30
      !####
      ! draw independent t variates, put into wkpA
      if( rchisq_R( df, rchisq, err ) == RETURN_FAIL ) goto 800
      rtmp = sqrt( df / rchisq )
      do j = 1, work%d
         if( rnorm_R( rnorm, err ) == RETURN_FAIL ) goto 800
         work%wkdA(j) = rnorm * rtmp
      end do
      ! premultiply by t(scale_sqrt), put into wkdB
      do j = 1, work%d
         sum = 0.D0
         do k = j, work%d
            sum = sum + work%beta_scale_sqrt(k,j) * work%wkdA(k)
         end do
         work%wkdB(j) = sum
      end do
      ! add center to obtain candidate
      do j = 1, work%d
         work%beta_can(j) = work%wkdB(j) + work%beta_center(j)
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
20    call err_handle(err, 1, &
            comment = "The model is saturated" )
      goto 800
30    call err_handle(err, 1, &
            comment = "Degrees of freedom are not positive" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function draw_candidate_beta
   !###################################################################
    integer(kind=our_int) function compute_log_proposal_beta( beta, df, &
         ans, work, err ) result(answer)
      ! computes log-proposal density at beta, assuming that the
      ! proposal center and scale are in beta_center and beta_scale
      implicit none
      real(kind=our_dble), intent(in) :: beta(:)
      real(kind=our_dble), intent(in) :: df
      real(kind=our_dble), intent(out) :: ans
      ! declare workspaces
      type(workspace_type_rrlogit), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      integer(kind=our_int) :: j, k
      real(kind=our_dble) :: sum
      character(len=*), parameter :: &
           subname = "compute_log_proposal_beta"
      ! begin
      answer = RETURN_FAIL
      if( work%saturated ) goto 20
      if( df <= 0.D0 ) goto 30
      if( size(beta, kind=our_int) /= work%d ) goto 40
      !####
      ! wkdB = t( beta - center) %*% (lower tri scale_inv_sqrt)
      work%wkdA(:) = beta(:) - work%beta_center(:)
      do j = 1, work%d
         sum = 0.D0
         do k = j, work%d
            sum = sum + work%wkdA(k) * work%beta_scale_inv_sqrt(k,j)
         end do
         work%wkdB(j) = sum
      end do
      ! sum = t( beta - center) %*% scale_inv %*% ( beta - center )
      sum = 0.D0
      do j = 1, work%d
         sum = sum + work%wkdB(j)**2
      end do
      ans = - (( df + real(work%d, our_dble) ) &
           / 2.D0 ) * log( 1.D0 + sum / df )
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
20    call err_handle(err, 1, &
            comment = "The model is saturated" )
      goto 800
30    call err_handle(err, 1, &
            comment = "Degrees of freedom are not positive" )
      goto 800
40    call err_handle(err, 1, &
            comment = "Argument beta has incorrect size" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function compute_log_proposal_beta
   !###################################################################
   integer(kind=our_int) function run_rrlogit_rwm_nonsaturated( &
        start_logP, logP_series, coef_vec_series, &
        imp_mat_series, imp_vec_series, pi_marg_series, &
        n_iter_actual, n_sample_actual, n_imp_actual, &
        work, err ) result(answer)
      ! runs mcmc on the nonsaturated model
      ! running means are stored in work%beta_mean, work%pi_mat_mean,
      !   work%pistar_mat_mean; running cov matrix in work%beta_cov_mat
      implicit none
      ! outputs
      real(kind=our_dble), intent(out) :: start_logP
      real(kind=our_dble), intent(out) :: logP_series(:)
      real(kind=our_dble), intent(out) :: coef_vec_series(:,:)
      integer(kind=our_int), intent(out) :: imp_mat_series(:,:,:)
      integer(kind=our_int), intent(out) :: imp_vec_series(:,:)
      real(kind=our_dble), intent(out) :: pi_marg_series(:,:)
      integer(kind=our_int), intent(out) :: n_iter_actual
      integer(kind=our_int), intent(out) :: n_sample_actual
      integer(kind=our_int), intent(out) :: n_imp_actual
      ! workspaces
      type(workspace_type_rrlogit), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! locals
      integer(kind=our_int) :: iter, thin_step, imp_step
      real(kind=our_dble) :: loglik_tmp, logprior_tmp, &
           sum, nfloat, rtmp, log_mh_ratio, mh_ratio, runif
      logical :: accept
      integer(kind=our_int) :: cp, ystar, y, j, k
      character(len=*), parameter :: &
           subname = "run_rrlogit_rwm_nonsaturated"
      ! begin
      answer = RETURN_FAIL
      if( work%saturated ) goto 20
      if( work%method /= "MCMC" ) goto 30
      ! initialize iteration counters
      work%iter = 0
      work%iter_past_burn_in = 0
      work%store_count = 0
      work%imp_count = 0
      thin_step = 0
      imp_step = 0
      n_iter_actual = 0
      n_sample_actual = 0
      n_imp_actual = 0
      ! initialize running means and returned series
      start_logP = 0.D0
      work%beta_mean(:) = 0.D0
      work%beta_cov_mat(:,:) = 0.D0
      work%pi_mat_mean(:,:) = 0.D0
      work%pistar_mat_mean(:,:) = 0.D0
      work%pi_marg_mean(:) = 0.D0
      logP_series(:) = work%mvcode
      coef_vec_series(:,:) = work%mvcode
      imp_mat_series(:,:,:) = -1
      imp_vec_series(:,:) = -1
      work%fimp_mat(:,:) = 0
      !
      if( compute_scale_rwm( work%df_rwm, work%scale_fac_rwm, work, err ) &
           == RETURN_FAIL ) goto 800
      ! initialize MH diagnostics
      work%beta_accept_count = 0
      work%beta_current_reject_run = 0
      work%beta_accept_rate = 0.D0
      !
      work%aborted = .false.
      accept = .false.
      ! main iteration
      do iter = 1, work%iter_mcmc + work%burn_mcmc
         work%aborted = .true.  ! will be set to .false. at end of cycle
         work%iter = iter
         work%iter_past_burn_in = work%iter - work%burn_mcmc ! could be neg
         work%store_this_iter = .false.
         work%imp_this_iter = .false.
         if( work%iter_past_burn_in > 0 ) then
            thin_step = thin_step + 1
            if( thin_step == work%thin_mcmc ) then
               work%store_this_iter = .true.
               work%store_count = work%store_count + 1
               thin_step = 0
            end if
            imp_step = imp_step + 1
            if( ( imp_step == work%impute_every ) .and. &
                 ( work%impute_every /= 0 ) ) then
               work%imp_this_iter = .true.
               work%imp_count = work%imp_count + 1
               imp_step = 0
            end if
         end if
         !#######################
         if( work%iter == 1 ) then
            if( compute_pi_mat_rrlogit( work, err ) == RETURN_FAIL ) goto 10
            if( compute_pi_marg( work, err ) == RETURN_FAIL ) goto 10
            if( compute_loglik_logprior_rrlogit( work, err, &
                 logprior = logprior_tmp, &
                 loglik = loglik_tmp ) == RETURN_FAIL ) goto 10
            work%loglik = loglik_tmp
            work%logprior = logprior_tmp
            work%logP = loglik_tmp + logprior_tmp
            start_logP = work%logP
         else
            if( accept ) then
               work%loglik = loglik_tmp
               work%logprior = logprior_tmp
               work%logP = loglik_tmp + logprior_tmp
            end if
         end if
         log_mh_ratio = - work%logP
         if( draw_candidate_beta_rwm( work%df_da, &
              work, err ) == RETURN_FAIL ) goto 10
         work%beta_vec_old(:) = work%beta_vec(:)
         work%pi_mat_old(:,:) = work%pi_mat(:,:)
         work%pistar_mat_old(:,:) = work%pistar_mat(:,:)
         work%beta_vec(:) = work%beta_can(:)
         if( compute_pi_mat_rrlogit( work, err ) == RETURN_FAIL ) goto 10
         if( compute_loglik_logprior_rrlogit( work, err, &
              logprior = logprior_tmp, &
              loglik = loglik_tmp ) == RETURN_FAIL ) goto 10
         log_mh_ratio = log_mh_ratio + logprior_tmp + loglik_tmp
         ! to prevent over/underflow
         if( log_mh_ratio > log_huge ) then
            mh_ratio = huge(0.D0)
         else if( log_mh_ratio < log_tiny ) then
            mh_ratio = 0.D0
         else
            mh_ratio = exp( log_mh_ratio )
         end if
         ! compare to uniform variate
         if( runif_R( runif, err ) == RETURN_FAIL ) goto 10
         if( runif <= mh_ratio ) then
            accept = .true.
            if( compute_pi_marg( work, err ) == RETURN_FAIL ) goto 10
         else
            accept = .false.
            work%beta_vec(:) = work%beta_vec_old(:)
            work%pi_mat(:,:) = work%pi_mat_old(:,:) 
            work%pistar_mat(:,:) = work%pistar_mat_old(:,:) 
         end if
         work%mh_ratios_beta( work%iter ) = mh_ratio
         work%mh_accept_beta( work%iter ) = accept
         !#######################
         ! update counters and running estimates
         if( accept ) then
            work%beta_accept_count = work%beta_accept_count + 1
            work%beta_current_reject_run = 0
            rtmp = 1.D0
         else
            work%beta_current_reject_run = &
                 work%beta_current_reject_run + 1
            rtmp = 0.D0
         end if
         work%beta_accept_rate = work%beta_accept_rate + &
              ( rtmp - work%beta_accept_rate ) / real( work%iter, our_dble )
         !#######################
         if( work%store_this_iter ) then
            ! store results in series
            if( work%store_count < 0 ) goto 300
            if( work%store_count > size(logp_series, kind=our_int) ) goto 300
            logP_series( work%store_count ) = work%logP
            if( work%store_count > size(coef_vec_series, 1, kind=our_int) ) &
                 goto 300
            coef_vec_series( work%store_count, : ) = work%beta_vec(:)
            if( work%store_count > size(pi_marg_series, 1, kind=our_int) ) &
                 goto 300
            pi_marg_series( work%store_count, : ) = work%pi_marg(:)
         end if
         !#######################
         if( work%iter_past_burn_in > 0 ) then
            ! update running means in work%pi_mat_mean, work%beta_mean,
            ! and work%beta_cov_mat
            nfloat = real( work%iter_past_burn_in, our_dble )
            if( nfloat == 0.D0 ) goto 250
            do cp = 1, work%n_cov_patt
               do y = 1, work%r
                  rtmp = work%pi_mat(cp,y) - work%pi_mat_mean(cp,y)
                  work%pi_mat_mean(cp,y) = work%pi_mat_mean(cp,y) &
                       + rtmp / nfloat
               end do
            end do
            do y = 1, work%r
               rtmp = work%pi_marg(y) - work%pi_marg_mean(y)
               work%pi_marg_mean(y) = work%pi_marg_mean(y) &
                    + rtmp / nfloat
            end do
            do j = 1, work%d
               work%wkdA(j) = work%beta_vec(j) - work%beta_mean(j)
            end do
            do j = 1, work%d
               do k = j, work%d
                  work%wkddA(j,k) = work%wkdA(j) * work%wkdA(k)
               end do
            end do
            do j = 1, work%d
               work%beta_mean(j) = work%beta_mean(j) + work%wkdA(j) / nfloat
               do k = j, work%d
                  work%beta_cov_mat(j,k) = &
                       ( nfloat - 1.D0 ) * work%beta_cov_mat(j,k) + &
                       ( ( nfloat - 1.D0 ) / nfloat ) * work%wkddA(j,k)
                  work%beta_cov_mat(j,k) = work%beta_cov_mat(j,k) / nfloat
                  work%beta_cov_mat(k,j) = work%beta_cov_mat(j,k)
               end do
            end do
         end if
         !#######################
         if( work%imp_this_iter ) then
            if( work%imp_count < 0 ) goto 300
            if( work%imp_count > size(imp_vec_series, 2, kind=our_int) ) &
                 goto 300
            if( work%imp_count > size(imp_mat_series, 3, kind=our_int) ) &
                 goto 300
            if( run_rrlogit_impute_random( work%wide_format, &
                 work%micro_data, work%freq_row_input_data_int, &
                 work%response_row_input_data, work%f_mat_int, &
                 imp_mat_series(:,:,work%imp_count), &
                 imp_vec_series(:,work%imp_count), &
                 work, err, skip_compute_pi_mat = .true. ) &
                 == RETURN_FAIL ) goto 800
         end if
         !#######################
         if( work%beta_current_reject_run > work%stuck_limit ) then
            call err_handle(err, 1, &
                 comment = "Random-walk Metropolis got stuck" )
            goto 10
         end if
         work%aborted = .false.
      end do
10    continue
      if( work%aborted ) then
         !#### issue warning message and continue
         call err_handle(err, 1, &
              comment = "MCMC procedure aborted" )
         call err_handle(err, 5, iiter = int(work%iter, kind=our_int) )
         call err_handle(err, 2, whichsub = subname, whichmod = modname )
      end if
      !#####################
      if( accept ) then
         work%loglik = loglik_tmp
         work%logprior = logprior_tmp
         work%logP = loglik_tmp + logprior_tmp
      end if
      if( work%iter == 0 ) then
         if( compute_pi_mat_rrlogit( work, err ) == RETURN_FAIL ) goto 10
         if( compute_loglik_logprior_rrlogit( work, err, &
              logprior = logprior_tmp, &
              loglik = loglik_tmp ) == RETURN_FAIL ) goto 10
         work%loglik = loglik_tmp
         work%logprior = logprior_tmp
         work%logP = loglik_tmp + logprior_tmp
      else
         ! compute pistar_mat_mean
         if( work%pert_mat_identity ) then
            work%pistar_mat_mean(:,:) = work%pi_mat_mean(:,:)
            work%pistar_mat(:,:) = work%pi_mat(:,:)
         else
            do cp = 1, work%n_cov_patt
               do ystar = 1, work%r
                  sum = 0.D0
                  do y = 1, work%r
                     sum = sum + work%pert_mat(ystar,y) * work%pi_mat_mean(cp,y)
                  end do
                  work%pistar_mat_mean(cp,ystar) = sum
                  sum = 0.D0
                  do y = 1, work%r
                     sum = sum + work%pert_mat(ystar,y) * work%pi_mat(cp,y)
                  end do
                  work%pistar_mat(cp,ystar) = sum
               end do
            end do
         end if
      end if
      !####
      n_iter_actual = work%iter
      n_sample_actual = work%store_count
      n_imp_actual = work%imp_count
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
20    call err_handle(err, 1, &
            comment = "The model is saturated" )
      goto 800
30    call err_handle(err, 1, &
            comment = "Not prepared for MCMC" )
      goto 800
250   call err_handle(err, 1, &
            comment = "Attempted division by zero" )
      goto 800
300   call err_handle(err, 1, &
            comment = "Array bounds exceeded" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
   end function run_rrlogit_rwm_nonsaturated
   !###################################################################
    integer(kind=our_int) function compute_scale_rwm( &
        df, scale_fac, work, err ) result(answer)
      ! beta_scale for RWM proposal from work%vhat_beta_rwm
      ! also computes beta_scale_sqrt
      implicit none
      real(kind=our_dble), intent(in) :: df
      real(kind=our_dble), intent(in) :: scale_fac
      ! declare workspaces
      type(workspace_type_rrlogit), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      real(kind=our_dble) :: rtmp
      character(len=*), parameter :: &
           subname = "compute_scale_rwm"
      ! begin
      answer = RETURN_FAIL
      if( work%saturated ) goto 20 
      rtmp = ( df / &
           ( df + real( work%d, our_dble) ) ) &
                 * ( 1.D0 / ( scale_fac**2 ) )
      rtmp = 1.D0 / rtmp
      work%beta_scale(:,:) = rtmp * work%vhat_beta_rwm
      work%beta_scale_sqrt(:,:) = work%beta_scale(:,:)
      if( cholesky_in_place( work%beta_scale_sqrt, err ) &
           == RETURN_FAIL ) goto 200
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
20    call err_handle(err, 1, &
            comment = "The model is saturated" )
      goto 800
200   call err_handle(err, 1, &
            comment = "Matrix vhat_beta_rm not positive definite" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function compute_scale_rwm
    !##################################################################
    integer(kind=our_int) function compute_loglik_logprior_rrlogit( &
         work, err, logprior, loglik ) result(answer)
      ! Given the current pi_mat and pistar_mat,
      ! accumulates the loglik and logprior
      implicit none
      ! declare workspaces
      type(workspace_type_rrlogit), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare outputs
      real(kind=our_dble), intent(out) :: loglik, logprior
      ! declare locals
      integer(kind=our_int) :: cp, ystar, y, dp
      character(len=*), parameter :: &
           subname = "compute_loglik_logprior_rrlogit"
      ! begin
      answer = RETURN_FAIL
      logprior = 0.D0
      loglik = 0.D0
      do cp = 1, work%n_cov_patt
         if( work%empty_cov_patt(cp) ) cycle
         if( work%prior == "DAP" ) then
            do y = 1, work%r
               if( work%prior_freq(cp,y) == 0.D0 ) cycle
               if( work%pi_mat(cp,y) <= 0.D0 ) goto 200
               logprior = logprior + work%prior_freq(cp,y) * &
                    log( work%pi_mat(cp,y) )
            end do
         end if
         do dp = work%data_patt_st(cp), work%data_patt_fin(cp)
            ystar = work%response_for_data_patt(dp)
            if( .not. work%active_ystar(cp, ystar) ) cycle
            if( work%pistar_mat(cp,ystar) <= 0.D0 ) goto 200
            loglik = loglik + work%fitweight_data_patt(dp) * &
                    log( work%pistar_mat(cp,ystar) )
         end do
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
200   call err_handle(err, 1, &
            comment = "Attempted logarithm of non-positive number" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
    end function compute_loglik_logprior_rrlogit
   !##################################################################
    integer(kind=our_int) function draw_candidate_beta_rwm( df, &
         work, err ) result(answer)
      ! Draws candidate beta from multivariate t proposal, storing
      ! the result in work%beta_can 
      ! Depends on beta_scale_sqrt, so should be run after 
      ! compute_scale_rwm
      implicit none
      real(kind=our_dble), intent(in) :: df
      ! declare workspaces
      type(workspace_type_rrlogit), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      integer(kind=our_int) :: j, k
      real(kind=our_dble) :: rtmp, sum, rnorm, rchisq
      character(len=*), parameter :: &
           subname = "draw_candidate_beta_rwm"
      ! begin
      answer = RETURN_FAIL
      if( work%saturated ) goto 20
      if( df <= 0.D0 ) goto 30
      !####
      ! draw independent t variates, put into wkpA
      if( rchisq_R( df, rchisq, err ) == RETURN_FAIL ) goto 800
      rtmp = sqrt( df / rchisq )
      do j = 1, work%d
         if( rnorm_R( rnorm, err ) == RETURN_FAIL ) goto 800
         work%wkdA(j) = rnorm * rtmp
      end do
      ! premultiply by beta_scale_sqrt, put into wkdB
      do j = 1, work%d
         sum = 0.D0
         do k = 1, j
            sum = sum + work%beta_scale_sqrt(j,k) * work%wkdA(k)
         end do
         work%wkdB(j) = sum
      end do
      ! add center to obtain candidate
      do j = 1, work%d
         work%beta_can(j) = work%wkdB(j) + work%beta_vec(j)
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
20    call err_handle(err, 1, &
            comment = "The model is saturated" )
      goto 800
30    call err_handle(err, 1, &
            comment = "Degrees of freedom are not positive" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function draw_candidate_beta_rwm
    !##################################################################
   integer(kind=our_int) function run_rrlogit_approx_bayes( &
        start_logP, logP_series, coef_vec_series, &
        imp_mat_series, imp_vec_series, pi_marg_series, &
        approx_bayes_log_imp_ratios, &
        n_iter_actual, n_sample_actual, n_imp_actual, &
        work, err ) result(answer)
      ! runs approximate bayes on the nonsaturated model
      ! running means are stored in work%beta_mean, work%pi_mat_mean,
      !   work%pistar_mat_mean; running cov matrix in work%beta_cov_mat
      implicit none
      ! outputs
      real(kind=our_dble), intent(out) :: start_logP
      real(kind=our_dble), intent(out) :: logP_series(:)
      real(kind=our_dble), intent(out) :: coef_vec_series(:,:)
      integer(kind=our_int), intent(out) :: imp_mat_series(:,:,:)
      integer(kind=our_int), intent(out) :: imp_vec_series(:,:)
      real(kind=our_dble), intent(out) :: pi_marg_series(:,:)
      real(kind=our_dble), intent(out) :: approx_bayes_log_imp_ratios(:)
      integer(kind=our_int), intent(out) :: n_iter_actual
      integer(kind=our_int), intent(out) :: n_sample_actual
      integer(kind=our_int), intent(out) :: n_imp_actual
      ! workspaces
      type(workspace_type_rrlogit), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! locals
      integer(kind=our_int) :: iter
      real(kind=our_dble) :: loglik_tmp, logprior_tmp, sum, nfloat, rtmp
      integer(kind=our_int) :: cp, ystar, y, j, k
      character(len=*), parameter :: &
           subname = "run_rrlogit_approx_bayes"
      ! begin
      answer = RETURN_FAIL
      if( work%saturated ) goto 20
      if( work%method /= "approxBayes" ) goto 30
      ! initialize iteration counters
      work%iter = 0
      work%store_count = 0
      work%imp_count = 0
      n_iter_actual = 0
      n_sample_actual = 0
      n_imp_actual = 0
      ! initialize running means and returned series
      start_logP = 0.D0
      work%beta_mean(:) = 0.D0
      work%beta_cov_mat(:,:) = 0.D0
      work%pi_mat_mean(:,:) = 0.D0
      work%pistar_mat_mean(:,:) = 0.D0
      work%pi_marg_mean(:) = 0.D0
      logP_series(:) = work%mvcode
      coef_vec_series(:,:) = work%mvcode
      imp_mat_series(:,:,:) = -1
      imp_vec_series(:,:) = -1
      approx_bayes_log_imp_ratios(:) = work%mvcode
      ! get center and scale; assume we are starting at mode
      work%beta_center(:) = work%beta_vec(:)
      work%beta_scale(:,:) = work%vhat_beta_rwm(:,:)
      work%beta_scale_sqrt(:,:) = work%beta_scale(:,:)
      if( cholesky_in_place( work%beta_scale_sqrt, err ) &
           == RETURN_FAIL ) goto 100
      work%beta_scale_inv_sqrt(:,:) = work%beta_scale_sqrt(:,:)
      if( invert_lower( work%beta_scale_inv_sqrt, err ) &
           == RETURN_FAIL ) goto 100
      ! initialize MH diagnostics
      work%beta_accept_count = 0
      work%beta_accept_rate = 0.D0
      work%store_this_iter = .true.
      if( work%impute_approx_bayes ) then
         work%imp_this_iter = .true.
      else
         work%imp_this_iter = .false.
      end if
      ! main iteration
      do iter = 1, work%iter_approx_bayes
         work%aborted = .true.  ! will be set to .false. at end of cycle
         work%iter = iter
         work%iter_past_burn_in = iter
         work%store_count = work%store_count + 1
         if( work%imp_this_iter ) work%imp_count = work%imp_count + 1
         !#######################
         if( work%iter == 1 ) then
            if( compute_pi_mat_rrlogit( work, err ) == RETURN_FAIL ) goto 10
            if( compute_loglik_logprior_rrlogit( work, err, &
                 logprior = logprior_tmp, &
                 loglik = loglik_tmp ) == RETURN_FAIL ) goto 10
         end if
         work%loglik = loglik_tmp
         work%logprior = logprior_tmp
         work%logP = loglik_tmp + logprior_tmp
         if( work%iter == 1 ) start_logP = work%logP
         if( draw_approx_bayes_beta( work, err ) == RETURN_FAIL ) goto 10
         if( compute_pi_mat_rrlogit( work, err ) == RETURN_FAIL ) goto 10
         if( compute_pi_marg( work, err ) == RETURN_FAIL ) goto 10
         if( compute_loglik_logprior_rrlogit( work, err, &
              logprior = logprior_tmp, &
              loglik = loglik_tmp ) == RETURN_FAIL ) goto 10
         !#######################
         ! update counters and running estimates
         work%beta_accept_rate = 1.D0
         work%beta_accept_count = work%beta_accept_count + 1
         !#######################
         if( work%store_this_iter ) then
            ! store results in series
            if( work%store_count < 0 ) goto 300
            if( work%store_count > size(logp_series, kind=our_int) ) goto 300
            logP_series( work%store_count ) = work%logP
            if( work%store_count > size(coef_vec_series, 1, kind=our_int) ) &
                 goto 300
            coef_vec_series( work%store_count, : ) = work%beta_vec(:)
            if( work%store_count > size(pi_marg_series, 1, kind=our_int) ) &
                 goto 300
            pi_marg_series( work%store_count, : ) = work%pi_marg(:)
            if( work%store_count > size(approx_bayes_log_imp_ratios, &
                 kind=our_int ) ) goto 300
            rtmp = loglik_tmp + logprior_tmp - work%beta_vec_log_dens
            approx_bayes_log_imp_ratios( work%store_count ) = rtmp
         end if
         !#######################
         if( work%iter_past_burn_in > 0 ) then
            ! update running means in work%pi_mat_mean, work%beta_mean,
            ! and work%beta_cov_mat
            nfloat = real( work%iter_past_burn_in, our_dble )
            if( nfloat == 0.D0 ) goto 250
            do cp = 1, work%n_cov_patt
               do y = 1, work%r
                  rtmp = work%pi_mat(cp,y) - work%pi_mat_mean(cp,y)
                  work%pi_mat_mean(cp,y) = work%pi_mat_mean(cp,y) &
                       + rtmp / nfloat
               end do
            end do
            do y = 1, work%r
               rtmp = work%pi_marg(y) - work%pi_marg_mean(y)
               work%pi_marg_mean(y) = work%pi_marg_mean(y) &
                    + rtmp / nfloat
            end do
            do j = 1, work%d
               work%wkdA(j) = work%beta_vec(j) - work%beta_mean(j)
            end do
            do j = 1, work%d
               do k = j, work%d
                  work%wkddA(j,k) = work%wkdA(j) * work%wkdA(k)
               end do
            end do
            do j = 1, work%d
               work%beta_mean(j) = work%beta_mean(j) + work%wkdA(j) / nfloat
               do k = j, work%d
                  work%beta_cov_mat(j,k) = &
                       ( nfloat - 1.D0 ) * work%beta_cov_mat(j,k) + &
                       ( ( nfloat - 1.D0 ) / nfloat ) * work%wkddA(j,k)
                  work%beta_cov_mat(j,k) = work%beta_cov_mat(j,k) / nfloat
                  work%beta_cov_mat(k,j) = work%beta_cov_mat(j,k)
               end do
            end do
         end if
         !#######################
         if( work%imp_this_iter ) then
            if( work%imp_count < 0 ) goto 300
            if( work%imp_count > size(imp_vec_series, 2, kind=our_int) ) &
                 goto 300
            if( work%imp_count > size(imp_mat_series, 3, kind=our_int) ) &
                 goto 300
            if( run_rrlogit_impute_random( work%wide_format, &
                 work%micro_data, work%freq_row_input_data_int, &
                 work%response_row_input_data, work%f_mat_int, &
                 imp_mat_series(:,:,work%imp_count), &
                 imp_vec_series(:,work%imp_count), &
                 work, err, skip_compute_pi_mat = .true. ) &
                 == RETURN_FAIL ) goto 800
         end if
         !#######################
         work%aborted = .false.
      end do
10    continue
      if( work%aborted ) then
         !#### issue warning message and continue
         call err_handle(err, 1, &
              comment = "Approx Bayes procedure aborted" )
         call err_handle(err, 5, iiter = int(work%iter, kind=our_int) )
         call err_handle(err, 2, whichsub = subname, whichmod = modname )
      end if
      !#####################
      if( work%iter == 0 ) then
         if( compute_pi_mat_rrlogit( work, err ) == RETURN_FAIL ) goto 10
         if( compute_loglik_logprior_rrlogit( work, err, &
              logprior = logprior_tmp, &
              loglik = loglik_tmp ) == RETURN_FAIL ) goto 10
         work%loglik = loglik_tmp
         work%logprior = logprior_tmp
         work%logP = loglik_tmp + logprior_tmp
      else
         work%loglik = loglik_tmp
         work%logprior = logprior_tmp
         work%logP = loglik_tmp + logprior_tmp
         ! compute pistar_mat_mean
         if( work%pert_mat_identity ) then
            work%pistar_mat_mean(:,:) = work%pi_mat_mean(:,:)
         else
            do cp = 1, work%n_cov_patt
               do ystar = 1, work%r
                  sum = 0.D0
                  do y = 1, work%r
                     sum = sum + work%pert_mat(ystar,y) * work%pi_mat_mean(cp,y)
                  end do
                  work%pistar_mat_mean(cp,ystar) = sum
                  sum = 0.D0
                  do y = 1, work%r
                     sum = sum + work%pert_mat(ystar,y) * work%pi_mat(cp,y)
                  end do
                  work%pistar_mat(cp,ystar) = sum
               end do
            end do
         end if
      end if
      !####
      n_iter_actual = work%iter
      n_sample_actual = work%store_count
      n_imp_actual = work%imp_count
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
20    call err_handle(err, 1, &
            comment = "The model is saturated" )
      goto 800
30    call err_handle(err, 1, &
            comment = "Not prepared for approx Bayes" )
      goto 800
100   call err_handle(err, 1, &
            comment = "Matrix vhat_beta_rwm not positive definite" )
      goto 800
250   call err_handle(err, 1, &
            comment = "Attempted division by zero" )
      goto 800
300   call err_handle(err, 1, &
            comment = "Array bounds exceeded" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
   end function run_rrlogit_approx_bayes
   !###################################################################
    integer(kind=our_int) function draw_approx_bayes_beta( &
         work, err ) result(answer)
      ! Draws beta from multivariate normal 
      ! puts result in work%beta_vec
      ! Depends on beta_center, beta_scale_sqrt, beta_scale_inv_sqrt
      ! Also evaluates the multivariate normal log-density at
      ! the drawn value, putting the result in work%beta_vec_log_dens 
      implicit none
      ! declare workspaces
      type(workspace_type_rrlogit), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! declare locals
      integer(kind=our_int) :: j, k
      real(kind=our_dble) :: sum, rnorm
      character(len=*), parameter :: &
           subname = "draw_approx_bayes_beta"
      ! begin
      answer = RETURN_FAIL
      if( work%saturated ) goto 20
      !####
      do j = 1, work%d
         if( rnorm_R( rnorm, err ) == RETURN_FAIL ) goto 800
         work%wkdA(j) = rnorm
      end do
      ! premultiply by beta_scale_sqrt, put into wkdB
      do j = 1, work%d
         sum = 0.D0
         do k = 1, j
            sum = sum + work%beta_scale_sqrt(j,k) * work%wkdA(k)
         end do
         work%wkdB(j) = sum
      end do
      ! add center to obtain beta
      do j = 1, work%d
         work%beta_vec(j) = work%wkdB(j) + work%beta_center(j)
      end do
      ! compute log-density, omitting additive constant
!     work%beta_vec_log_dens = 0.D0
!      do j = 1, work%d
!         sum = 0.D0
!         do k = 1, j
!            sum = sum + work%beta_scale_inv_sqrt(j,k) * work%wkdB(k)
!         end do
!         work%wkdA(j) = sum
!      end do
!      sum = 0.D0
!      do j = 1, work%d
!         sum = sum + work%wkdA(j)**2
!      end do
!      work%beta_vec_log_dens = -0.5D0 * sum

      do j = 1, work%d
         sum = 0.D0
         do k = 1, work%d
            sum = sum + work%beta_scale_inv(j,k) * work%wkdB(k)
         end do
         work%wkdA(j) = sum
      end do
      sum = 0.D0
      do j = 1, work%d
         sum = sum + work%wkdA(j) * work%wkdB(j)
      end do
      work%beta_vec_log_dens = -0.5D0 * sum

      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
20    call err_handle(err, 1, &
            comment = "The model is saturated" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function draw_approx_bayes_beta
   !###################################################################
   integer(kind=our_int) function run_rrlogit_predict( &
        nrow_input_data, n_levels, n_cov_patt, n_data_patt, &
        ncol_model_matrix, nparam_this_model, nparam_sat_model, &
        wide_format_int, &
        model_matrix, fitweight_row_input_data, fitweight_cov_patt, &
        fitweight_data_patt, freq_for_data_patt, cov_patt, data_patt, &
        cov_patt_for_data_patt, response_for_data_patt, &
        survey_mode_int, baseline_int, pert_mat, pert_mat_inv, &
        prior_int, prior_freq_tot, prior_alloc_supplied_int, prior_alloc, &
        saturated_int, method_int, &
        iter_max_nr, iter_max_fs, iter_max_em, iter_max_mstep, &
        crit_converge, crit_boundary, &
        mvcode, nancode, infcode, neginfcode, &
        !
        type_int, noisy_int, se_fit_int, freq_row_input_data, &
        coefficients, vhat_coef, &
        fitted_pi, fitted_pistar, &
        fitted, se_mat, vhat_fitted_array, &
        !
        work, err ) result(answer)
      ! fits randomized response logistic model using Newton-Raphson
      ! or Fisher scoring
      implicit none
      ! inputs
      integer(kind=our_int), intent(in) :: nrow_input_data
      integer(kind=our_int), intent(in) :: n_levels
      integer(kind=our_int), intent(in) :: n_cov_patt
      integer(kind=our_int), intent(in) :: n_data_patt
      integer(kind=our_int), intent(in) :: ncol_model_matrix
      integer(kind=our_int), intent(in) :: nparam_this_model
      integer(kind=our_int), intent(in) :: nparam_sat_model
      integer(kind=our_int), intent(in) :: wide_format_int
      real(kind=our_dble), intent(in) :: model_matrix(:,:)
      real(kind=our_dble), intent(in) :: fitweight_row_input_data(:)
      real(kind=our_dble), intent(in) :: fitweight_cov_patt(:)
      real(kind=our_dble), intent(in) :: fitweight_data_patt(:)
      real(kind=our_dble), intent(in) :: freq_for_data_patt(:)
      integer(kind=our_int), intent(in) :: cov_patt(:)
      integer(kind=our_int), intent(in) :: data_patt(:)
      integer(kind=our_int), intent(in) :: cov_patt_for_data_patt(:)
      integer(kind=our_int), intent(in) :: response_for_data_patt(:)
      integer(kind=our_int), intent(in) :: survey_mode_int
      integer(kind=our_int), intent(in) :: baseline_int
      real(kind=our_dble), intent(in) :: pert_mat(:,:)
      real(kind=our_dble), intent(in) :: pert_mat_inv(:,:)
      integer(kind=our_int), intent(in) :: prior_int
      real(kind=our_dble), intent(in) :: prior_freq_tot
      integer(kind=our_int), intent(in) :: prior_alloc_supplied_int
      real(kind=our_dble), intent(inout) :: prior_alloc(:)
      integer(kind=our_int), intent(in) :: saturated_int
      integer(kind=our_int), intent(in) :: method_int
      integer(kind=our_int), intent(in) :: iter_max_nr
      integer(kind=our_int), intent(in) :: iter_max_fs
      integer(kind=our_int), intent(in) :: iter_max_em
      integer(kind=our_int), intent(in) :: iter_max_mstep
      real(kind=our_dble), intent(in) :: crit_converge
      real(kind=our_dble), intent(in) :: crit_boundary
      real(kind=our_dble), intent(in) :: mvcode
      real(kind=our_dble), intent(in) :: nancode
      real(kind=our_dble), intent(in) :: infcode
      real(kind=our_dble), intent(in) :: neginfcode
      integer(kind=our_int), intent(in) :: type_int
      integer(kind=our_int), intent(in) :: noisy_int
      integer(kind=our_int), intent(in) :: se_fit_int
      real(kind=our_dble), intent(in) :: freq_row_input_data(:)
      real(kind=our_dble), intent(in) :: coefficients(:,:) 
      real(kind=our_dble), intent(in) :: vhat_coef(:,:)
      ! inouts
      real(kind=our_dble), intent(inout) :: fitted_pi(:,:)
      real(kind=our_dble), intent(inout) :: fitted_pistar(:,:)
      ! outputs
      real(kind=our_dble), intent(out) :: fitted(:,:)
      real(kind=our_dble), intent(out) :: se_mat(:,:)
      real(kind=our_dble), intent(out) :: vhat_fitted_array(:,:,:)
      ! workspaces
      type(workspace_type_rrlogit), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! locals
      integer(kind=our_int) :: ijunk
      character(len=*), parameter :: &
           subname = "run_rrlogit_predict"
      ! begin
      answer = RETURN_FAIL
      !####
      if( put_data_into_workspace_rrlogit( &
           nrow_input_data, n_levels, n_cov_patt, n_data_patt, &
           ncol_model_matrix, nparam_this_model, nparam_sat_model, &
           wide_format_int, &
           model_matrix, fitweight_row_input_data, fitweight_cov_patt, &
           fitweight_data_patt, freq_for_data_patt, cov_patt, data_patt, &
           cov_patt_for_data_patt, response_for_data_patt, &
           survey_mode_int, baseline_int, pert_mat, pert_mat_inv, &
           prior_int, prior_freq_tot, prior_alloc_supplied_int, prior_alloc, & 
           saturated_int, method_int, &
           iter_max_nr, iter_max_fs, iter_max_em, iter_max_mstep, &
           crit_converge, crit_boundary, &
           mvcode, nancode, infcode, neginfcode, &
           work, err ) == RETURN_FAIL ) goto 800
      if( work%saturated ) then
         if( size( fitted_pi, 1, kind=our_int ) /= work%n_cov_patt ) goto 100
         if( size( fitted_pi, 2, kind=our_int ) /= work%r ) goto 100
         work%pi_mat(:,:) = fitted_pi(:,:)
         if( size( fitted_pistar, 1, kind=our_int ) &
              /= work%n_cov_patt ) goto 105
         if( size( fitted_pistar, 2, kind=our_int ) /= work%r ) goto 105
         work%pistar_mat(:,:) = fitted_pistar(:,:)
         if( run_rrlogit_predict_saturated( &
              type_int, noisy_int, freq_row_input_data, fitted, &
              work, err ) == RETURN_FAIL ) goto 800
      else
         if( size( coefficients, 1, kind=our_int ) /= work%p ) goto 110
         if( size( coefficients, 2, kind=our_int ) /= work%r ) goto 110
         if( stack_coef( coefficients, work%beta_vec, work ) &
              == RETURN_FAIL ) goto 800
         if( size( vhat_coef, 1, kind=our_int ) /= work%d ) goto 115
         if( size( vhat_coef, 2, kind=our_int ) /= work%d ) goto 115
         work%vhat_coef(:,:) = vhat_coef(:,:)
         if( run_rrlogit_predict_nonsaturated( &
              type_int, noisy_int, se_fit_int, freq_row_input_data, fitted, &
              se_mat, vhat_fitted_array, work, err ) == RETURN_FAIL ) goto 800
      end if
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
100   call err_handle(err, 1, &
            comment = "Array pi_mat has incorrect size" )
      goto 800
105   call err_handle(err, 1, &
            comment = "Array pistar_mat has incorrect size" )
      goto 800
110   call err_handle(err, 1, &
            comment = "Array coefficients has incorrect size" )
      goto 800
115   call err_handle(err, 1, &
            comment = "Array vhat_coef has incorrect size" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
      ijunk = nullify_workspace_type_rrlogit( work, err )
   end function run_rrlogit_predict
   !###################################################################
   integer(kind=our_int) function run_rrlogit_predict_saturated( &
        type_int, noisy_int, freq_row_input_data, fitted, &
        work, err ) result(answer)
      implicit none
      ! inputs
      integer(kind=our_int), intent(in) :: type_int
      integer(kind=our_int), intent(in) :: noisy_int
      real(kind=our_dble), intent(in) :: freq_row_input_data(:)
      ! outputs
      real(kind=our_dble), intent(out) :: fitted(:,:)
      ! workspaces
      type(workspace_type_rrlogit), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! locals
      character(len=method_str_len) :: typeStr = ""
      logical :: noisy
      integer(kind=our_int) :: dr, cp, y, ystar, baseline
      real(kind=our_dble) :: rtmp, num, den
      character(len=*), parameter :: &
           subname = "run_logit_predict_saturated"
      ! begin
      answer = RETURN_FAIL
      ! check all args
      if( type_int == 1 ) then
         typeStr = "prob"
      else if( type_int == 2 ) then
         typeStr = "link"
      else if( type_int == 3 ) then
         typeStr = "mean"
      else
         goto 100
      end if
      if( noisy_int == 0 ) then
         noisy = .false.
      else if( noisy_int == 1 ) then
         noisy = .true.
      else
         goto 105
      end if
      if( size(freq_row_input_data, kind=our_int) /= work%nrow_input_data ) goto 110
      do dr = 1, work%nrow_input_data
         if( freq_row_input_data(dr) < 0.D0 ) goto 115
      end do
      if( size(fitted, 1, kind=our_int) /= work%nrow_input_data ) goto 120
      if( size(fitted, 2, kind=our_int) /= work%r ) goto 120
      !####
      if( noisy ) then
         if( typeStr == "prob" ) then
            ! copy rows of pistar_mat into fitted, including mvcode values
            do dr = 1, work%nrow_input_data
               cp = work%cov_patt(dr)
               fitted(dr,:) = work%pistar_mat(cp,:)
            end do
         else if( typeStr == "link" ) then
            goto 200
         else if( typeStr == "mean" ) then
            ! copy rows of pistar_mat into fitted, including mvcode values,
            ! and multiply by the frequency
            do dr = 1, work%nrow_input_data
               cp = work%cov_patt(dr)
               fitted(dr,:) = work%pistar_mat(cp,:)
               rtmp = freq_row_input_data(dr)
               do ystar = 1, work%r
                  if( fitted(dr,ystar) == work%mvcode ) cycle
                  fitted(dr,ystar) = rtmp * fitted(dr,ystar)
               end do
            end do
         else
            goto 100
         end if
      else
         if( typeStr == "prob" ) then
            ! copy rows of pi_mat into fitted, including mvcode values
            do dr = 1, work%nrow_input_data
               cp = work%cov_patt(dr)
               fitted(dr,:) = work%pi_mat(cp,:)
            end do
         else if( typeStr == "link" ) then
            ! convert probs in pi_mat to logit scale
            baseline = work%baseline
            do dr = 1, work%nrow_input_data
               cp = work%cov_patt(dr)
               den = work%pi_mat(cp,baseline)
               if( den == work%mvcode ) then
                  fitted(dr,:) = work%mvcode
                  cycle
               end if
               do y = 1, work%r
                  if( y == baseline ) cycle
                  num = work%pi_mat(cp,y)
                  if( num == work%mvcode ) then
                     fitted(dr,y) = work%mvcode
                     cycle
                  end if
                  if( den == 0.D0 ) then
                     if( num == 0.D0 ) then
                        fitted(dr,y) = work%nancode
                     else
                        fitted(dr,y) = work%infcode
                     end if
                  else
                     if( num == 0.D0 ) then
                        fitted(dr,y) = work%neginfcode
                     else
                        rtmp = num / den
                        if( rtmp <= 0.D0 ) then
                           fitted(dr,y) = work%neginfcode
                        else
                           fitted(dr,y) = log(rtmp)
                        end if
                     end if
                  end if
               end do
            end do
         else if( typeStr == "mean" ) then
            ! copy rows of pi_mat into fitted, including mvcode values,
            ! and multiply by the frequency
            do dr = 1, work%nrow_input_data
               cp = work%cov_patt(dr)
               fitted(dr,:) = work%pi_mat(cp,:)
               rtmp = freq_row_input_data(dr)
               do y = 1, work%r
                  if( fitted(dr,y) == work%mvcode ) cycle
                  fitted(dr,y) = rtmp * fitted(dr,y)
               end do
            end do
         else
            goto 100
         end if
      end if
      !####
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
100   call err_handle(err, 1, &
            comment = "Value of type_int not recognized" )
      goto 800
105   call err_handle(err, 1, &
            comment = "Value of noisy_int not recognized" )
      goto 800
110   call err_handle(err, 1, &
            comment = "Array freq_row_input_data has incorrect size" )
      goto 800
115   call err_handle(err, 1, &
            comment = "Negative frequency encountered" )
      goto 800
120   call err_handle(err, 1, &
            comment = "Array fitted has incorrect size" )
      goto 800
200   call err_handle(err, 1, &
            comment = "type=link not allowed when noisy is .true." )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function run_rrlogit_predict_saturated
   !###################################################################
   integer(kind=our_int) function run_rrlogit_predict_nonsaturated( &
        type_int, noisy_int, se_fit_int, freq_row_input_data, fitted, &
        se_mat, vhat_fitted_array, work, err ) result(answer)
      implicit none
      ! inputs
      integer(kind=our_int), intent(in) :: type_int
      integer(kind=our_int), intent(in) :: noisy_int
      integer(kind=our_int), intent(in) :: se_fit_int
      real(kind=our_dble), intent(in) :: freq_row_input_data(:)
      ! outputs
      real(kind=our_dble), intent(out) :: fitted(:,:)
      real(kind=our_dble), intent(out) :: se_mat(:,:)
      real(kind=our_dble), intent(out) :: vhat_fitted_array(:,:,:)
      ! workspaces
      type(workspace_type_rrlogit), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! locals
      character(len=method_str_len) :: typeStr = ""
      logical :: noisy, se_fit, skip_pistar, logit_scale
      integer(kind=our_int) :: dr, cp, y, ystar
      real(kind=our_dble) :: rtmp
      character(len=*), parameter :: &
           subname = "run_logit_predict_nonsaturated"
      ! begin
      answer = RETURN_FAIL
      ! check all args
      if( type_int == 1 ) then
         typeStr = "prob"
      else if( type_int == 2 ) then
         typeStr = "link"
      else if( type_int == 3 ) then
         typeStr = "mean"
      else
         goto 100
      end if
      if( noisy_int == 0 ) then
         noisy = .false.
      else if( noisy_int == 1 ) then
         noisy = .true.
      else
         goto 105
      end if
      if( size(freq_row_input_data, kind=our_int) &
           /= work%nrow_input_data ) goto 110
      do dr = 1, work%nrow_input_data
         if( freq_row_input_data(dr) < 0.D0 ) goto 115
      end do
      if( size(fitted, 1, kind=our_int) /= work%nrow_input_data ) goto 120
      if( size(fitted, 2, kind=our_int) /= work%r ) goto 120
      if( se_fit_int == 0 ) then
         se_fit = .false.
      else if( se_fit_int == 1 ) then
         se_fit = .true.
      else
         goto 125
      end if
      if( se_fit ) then
         if( size(se_mat, 1, kind=our_int) /= work%nrow_input_data ) goto 130
         if( size(se_mat, 2, kind=our_int) /= work%r ) goto 130
         if( size(vhat_fitted_array, 1, kind=our_int) &
              /= work%nrow_input_data ) goto 135
         if( size(vhat_fitted_array, 2, kind=our_int) /= work%r ) goto 135
         if( size(vhat_fitted_array, 3, kind=our_int) /= work%r ) goto 135
      else
         if( size(se_mat, 1, kind=our_int) /= 0 ) goto 130
         if( size(se_mat, 2, kind=our_int) /= 0 ) goto 130
         if( size(vhat_fitted_array, 1, kind=our_int) /= 0 ) goto 135
         if( size(vhat_fitted_array, 2, kind=our_int) /= 0 ) goto 135
         if( size(vhat_fitted_array, 3, kind=our_int) /= 0 ) goto 135
      end if
      !####
      if( noisy ) then
         skip_pistar = .false.
      else
         skip_pistar = .true.
      end if
      if( typeStr == "link" ) then
         logit_scale = .true.
      else
         logit_scale = .false.
      end if
      if( compute_pi_mat_rrlogit( work, err, skip_pistar, logit_scale ) &
           == RETURN_FAIL ) goto 800
      !####
      if( noisy ) then
         if( typeStr == "prob" ) then
            ! copy rows of pistar_mat into fitted, including mvcode values
            do dr = 1, work%nrow_input_data
               cp = work%cov_patt(dr)
               fitted(dr,:) = work%pistar_mat(cp,:)
            end do
         else if( typeStr == "link" ) then
            goto 200
         else if( typeStr == "mean" ) then
            ! copy rows of pistar_mat into fitted, including mvcode values,
            ! and multiply by the frequency
            do dr = 1, work%nrow_input_data
               cp = work%cov_patt(dr)
               fitted(dr,:) = work%pistar_mat(cp,:)
               rtmp = freq_row_input_data(dr)
               do ystar = 1, work%r
                  if( fitted(dr,ystar) == work%mvcode ) cycle
                  fitted(dr,ystar) = rtmp * fitted(dr,ystar)
               end do
            end do
         else
            goto 100
         end if
      else
         if( typeStr == "prob" ) then
            ! copy rows of pi_mat into fitted, including mvcode values
            do dr = 1, work%nrow_input_data
               cp = work%cov_patt(dr)
               fitted(dr,:) = work%pi_mat(cp,:)
            end do
         else if( typeStr == "link" ) then
            ! copy linear predictors from pi_mat into fitted
            do dr = 1, work%nrow_input_data
               cp = work%cov_patt(dr)
               fitted(dr,:) = work%pi_mat(cp,:)
            end do
         else if( typeStr == "mean" ) then
            ! copy rows of pi_mat into fitted, including mvcode values,
            ! and multiply by the frequency
            do dr = 1, work%nrow_input_data
               cp = work%cov_patt(dr)
               fitted(dr,:) = work%pi_mat(cp,:)
               rtmp = freq_row_input_data(dr)
               do y = 1, work%r
                  if( fitted(dr,y) == work%mvcode ) cycle
                  fitted(dr,y) = rtmp * fitted(dr,y)
               end do
            end do
         else
            goto 100
         end if
      end if
      !####
      if( se_fit ) then
         if( compute_rrlogit_predict_se(typeStr, noisy, freq_row_input_data, &
              se_mat, vhat_fitted_array, work, err ) &
              == RETURN_FAIL ) goto 800
      end if
      !####
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
100   call err_handle(err, 1, &
            comment = "Value of type_int not recognized" )
      goto 800
105   call err_handle(err, 1, &
            comment = "Value of noisy_int not recognized" )
      goto 800
110   call err_handle(err, 1, &
            comment = "Array freq_row_input_data has incorrect size" )
      goto 800
115   call err_handle(err, 1, &
            comment = "Negative frequency encountered" )
      goto 800
120   call err_handle(err, 1, &
            comment = "Array fitted has incorrect size" )
      goto 800
125   call err_handle(err, 1, &
            comment = "Value of se_fit_int not recognized" )
      goto 800
130   call err_handle(err, 1, &
            comment = "Array se_mat has incorrect size" )
      goto 800
135   call err_handle(err, 1, &
            comment = "Array vhat_fitted_array has incorrect size" )
      goto 800
200   call err_handle(err, 1, &
            comment = "type=link not allowed when noisy is .true." )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function run_rrlogit_predict_nonsaturated
   !###################################################################
   integer(kind=our_int) function compute_rrlogit_predict_se( &
        typeStr, noisy, freq_row_input_data, &
        se_mat, vhat_fitted_array, work, err ) result(answer)
      implicit none
      ! inputs
      character(len=method_str_len), intent(in) :: typeStr
      logical, intent(in) :: noisy
      real(kind=our_dble), intent(in) :: freq_row_input_data(:)
      ! outputs
      real(kind=our_dble), intent(out) :: se_mat(:,:)
      real(kind=our_dble), intent(out) :: vhat_fitted_array(:,:,:)
      ! workspaces
      type(workspace_type_rrlogit), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! locals
      integer(kind=our_int) :: dr, cp, y, baseline, posn, j, k, l
      real(kind=our_dble) :: rtmp, sum
      character(len=*), parameter :: &
           subname = "compute_rrlogit_predict_se"
      ! begin
      answer = RETURN_FAIL
      ! cycle through covariate patterns to compute vhat_fit_array, se_fit_mat
      baseline = work%baseline
      do cp = 1, work%n_cov_patt
         if( noisy ) then
            if( ( typeStr == "prob" ) .or. ( typeStr == "mean" ) ) then
               ! compute derivs of pi wrt beta, store result in wkrdA
               do y = 1, work%r
                  posn = 0
                  do k = 1, work%r
                     if( k == baseline ) cycle
                     do j = 1, work%p
                        posn = posn + 1
                        if( k == y ) then
                           rtmp = 1.D0
                        else
                           rtmp = 0.D0
                        end if
                        work%wkrdA(y, posn) = work%pi_mat(cp,y) * &
                             ( rtmp - work%pi_mat(cp,k) ) * &
                             work%model_matrix(cp,j)
                     end do
                  end do
               end do
               ! premultiply by pert_mat to get derivs of pistar wrt beta,
               ! store result in dpidbeta
               do j = 1, work%r
                  do k = 1, work%d
                     sum = 0.D0
                     do l = 1, work%r
                        sum = sum + work%pert_mat(j,l) * work%wkrdA(l,k)
                     end do
                     work%dpidbeta(j,k) = sum
                  end do
               end do
            else if( typeStr == "link" ) then
               goto 200
            else
               goto 100
            end if
         else
            if( ( typeStr == "prob" ) .or. ( typeStr == "mean" ) ) then
               ! compute derivs of pi wrt beta, store result in dpidbeta
               do y = 1, work%r
                  posn = 0
                  do k = 1, work%r
                     if( k == baseline ) cycle
                     do j = 1, work%p
                        posn = posn + 1
                        if( k == y ) then
                           rtmp = 1.D0
                        else
                           rtmp = 0.D0
                        end if
                        work%dpidbeta(y, posn) = work%pi_mat(cp,y) * &
                             ( rtmp - work%pi_mat(cp,k) ) * &
                             work%model_matrix(cp,j)
                     end do
                  end do
               end do
            else if( typeStr == "link" ) then
               ! compute derivs of eta wrt beta, store result in dpidbeta
               do y = 1, work%r
                  if( y == baseline ) then
                     work%dpidbeta(y,:) = 0.D0
                     cycle
                  end if
                  posn = 0
                  do k = 1, work%r
                     if( k == baseline ) cycle
                     do j = 1, work%p
                        posn = posn + 1
                        if( k == y ) then
                           work%dpidbeta(y,posn) = work%model_matrix(cp,j)
                        else
                           work%dpidbeta(y,posn) = 0.D0
                        end if
                     end do
                  end do
               end do
            else
               goto 100
            end if
         end if
         ! premultiply vhat_beta by dpidbeta
         do j = 1, work%r
            do k = 1, work%d
               sum = 0.D0
               do l = 1, work%d
                  sum = sum + work%dpidbeta(j,l) * work%vhat_coef(l,k)
               end do
               work%wkrdA(j,k) = sum
            end do
         end do
         ! post-multiply by t(dpidbeta), taking advantage of symmetry
         do j = 1, work%r
            do k = 1, j
               sum = 0.D0
               do l = 1, work%d
                  sum = sum + work%wkrdA(j,l) * work%dpidbeta(k,l)
               end do
               work%vhat_fit_array(cp,j,k) = sum
            end do
         end do
         do j = 1, ( work%r - 1 )
            do k = ( j + 1 ), work%r
               work%vhat_fit_array(cp,j,k) = work%vhat_fit_array(cp,k,j)
            end do
         end do
         ! square roots of diagonals
         do j = 1, work%r
            if( work%vhat_fit_array(cp,j,j) < 0.D0 ) goto 300
            work%se_fit_mat(cp,j) = sqrt( work%vhat_fit_array(cp,j,j) )
         end do
      end do
      ! copy results into se_mat and vhat_fitted_array
      if( typeStr == "mean" ) then
         do dr = 1, work%nrow_input_data
            cp = work%cov_patt(dr)
            rtmp = freq_row_input_data(dr)
            se_mat(dr,:) = rtmp * work%se_fit_mat(cp,:)
            rtmp = freq_row_input_data(dr)**2
            vhat_fitted_array(dr,:,:) = rtmp * work%vhat_fit_array(cp,:,:)
         end do
      else
         do dr = 1, work%nrow_input_data
            cp = work%cov_patt(dr)
            se_mat(dr,:) = work%se_fit_mat(cp,:)
            vhat_fitted_array(dr,:,:) = work%vhat_fit_array(cp,:,:)
         end do
      end if
      !####
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
100   call err_handle(err, 1, &
            comment = "Value of typeStr not recognized" )
      goto 800
200   call err_handle(err, 1, &
            comment = "type=link not allowed when noisy is .true." )
      goto 800
300    call err_handle(err, 1, &
            comment = "Attempted square root of negative number" )
       goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function compute_rrlogit_predict_se
   !###################################################################
   integer(kind=our_int) function run_rrlogit_impute( &
        nrow_input_data, n_levels, n_cov_patt, n_data_patt, &
        ncol_model_matrix, nparam_this_model, nparam_sat_model, &
        wide_format_int, &
        model_matrix, fitweight_row_input_data, fitweight_cov_patt, &
        fitweight_data_patt, freq_for_data_patt, cov_patt, data_patt, &
        cov_patt_for_data_patt, response_for_data_patt, &
        survey_mode_int, baseline_int, pert_mat, pert_mat_inv, &
        prior_int, prior_freq_tot, prior_alloc_supplied_int, prior_alloc, &
        saturated_int, method_int, &
        iter_max_nr, iter_max_fs, iter_max_em, iter_max_mstep, &
        crit_converge, crit_boundary, &
        mvcode, nancode, infcode, neginfcode, &
        type_int, micro_data_int, &
        freq_row_input_data, freq_row_input_data_int, f_mat, f_mat_int, &
        response_row_input_data, &
        coefficients, &
        fitted_pi, fitted_pistar, &
        cond_means, imp_mat, imp_vec, &
        work, err ) result(answer)
      ! imputes true responses under randomized response logistic model
      implicit none
      ! inputs
      integer(kind=our_int), intent(in) :: nrow_input_data
      integer(kind=our_int), intent(in) :: n_levels
      integer(kind=our_int), intent(in) :: n_cov_patt
      integer(kind=our_int), intent(in) :: n_data_patt
      integer(kind=our_int), intent(in) :: ncol_model_matrix
      integer(kind=our_int), intent(in) :: nparam_this_model
      integer(kind=our_int), intent(in) :: nparam_sat_model
      integer(kind=our_int), intent(in) :: wide_format_int
      real(kind=our_dble), intent(in) :: model_matrix(:,:)
      real(kind=our_dble), intent(in) :: fitweight_row_input_data(:)
      real(kind=our_dble), intent(in) :: fitweight_cov_patt(:)
      real(kind=our_dble), intent(in) :: fitweight_data_patt(:)
      real(kind=our_dble), intent(in) :: freq_for_data_patt(:)
      integer(kind=our_int), intent(in) :: cov_patt(:)
      integer(kind=our_int), intent(in) :: data_patt(:)
      integer(kind=our_int), intent(in) :: cov_patt_for_data_patt(:)
      integer(kind=our_int), intent(in) :: response_for_data_patt(:)
      integer(kind=our_int), intent(in) :: survey_mode_int
      integer(kind=our_int), intent(in) :: baseline_int
      real(kind=our_dble), intent(in) :: pert_mat(:,:)
      real(kind=our_dble), intent(in) :: pert_mat_inv(:,:)
      integer(kind=our_int), intent(in) :: prior_int
      real(kind=our_dble), intent(in) :: prior_freq_tot
      integer(kind=our_int), intent(in) :: prior_alloc_supplied_int
      real(kind=our_dble), intent(inout) :: prior_alloc(:)
      integer(kind=our_int), intent(in) :: saturated_int
      integer(kind=our_int), intent(in) :: method_int
      integer(kind=our_int), intent(in) :: iter_max_nr
      integer(kind=our_int), intent(in) :: iter_max_fs
      integer(kind=our_int), intent(in) :: iter_max_em
      integer(kind=our_int), intent(in) :: iter_max_mstep
      real(kind=our_dble), intent(in) :: crit_converge
      real(kind=our_dble), intent(in) :: crit_boundary
      real(kind=our_dble), intent(in) :: mvcode
      real(kind=our_dble), intent(in) :: nancode
      real(kind=our_dble), intent(in) :: infcode
      real(kind=our_dble), intent(in) :: neginfcode
      integer(kind=our_int), intent(in) :: type_int
      integer(kind=our_int), intent(in) :: micro_data_int
      real(kind=our_dble), intent(in) :: freq_row_input_data(:)
      integer(kind=our_int), intent(in) :: freq_row_input_data_int(:)
      real(kind=our_dble), intent(in) :: f_mat(:,:)
      integer(kind=our_int), intent(in) :: f_mat_int(:,:)
      integer(kind=our_int), intent(in) :: response_row_input_data(:)
      real(kind=our_dble), intent(in) :: coefficients(:,:) 
      ! inouts
      real(kind=our_dble), intent(inout) :: fitted_pi(:,:)
      real(kind=our_dble), intent(inout) :: fitted_pistar(:,:)
      ! outputs
      real(kind=our_dble), intent(out) :: cond_means(:,:)
      integer(kind=our_int), intent(out) :: imp_mat(:,:)
      integer(kind=our_int), intent(out) :: imp_vec(:)
      ! workspaces
      type(workspace_type_rrlogit), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! locals
      integer(kind=our_int) :: dr, y
      integer(kind=our_int) :: ijunk
      character(len=method_str_len) :: typeStr = ""
      logical :: wide_format, micro_data
      character(len=*), parameter :: &
           subname = "run_rrlogit_impute"
      ! begin
      answer = RETURN_FAIL
      !####
      if( put_data_into_workspace_rrlogit( &
           nrow_input_data, n_levels, n_cov_patt, n_data_patt, &
           ncol_model_matrix, nparam_this_model, nparam_sat_model, &
           wide_format_int, &
           model_matrix, fitweight_row_input_data, fitweight_cov_patt, &
           fitweight_data_patt, freq_for_data_patt, cov_patt, data_patt, &
           cov_patt_for_data_patt, response_for_data_patt, &
           survey_mode_int, baseline_int, pert_mat, pert_mat_inv, &
           prior_int, prior_freq_tot, prior_alloc_supplied_int, prior_alloc, & 
           saturated_int, method_int, &
           iter_max_nr, iter_max_fs, iter_max_em, iter_max_mstep, &
           crit_converge, crit_boundary, &
           mvcode, nancode, infcode, neginfcode, &
           work, err ) == RETURN_FAIL ) goto 800
      !####
      if( type_int == 1 ) then
         typeStr = "random"
      else if( type_int == 2 ) then
         typeStr = "condMean"
      else
         goto 100
      end if
      if( wide_format_int == 0 ) then
         wide_format = .false.
      else if( wide_format_int == 1 ) then
         wide_format = .true.
      else
         goto 110
      end if
      if( micro_data_int == 0 ) then
         micro_data = .false.
      else if( micro_data_int == 1 ) then
         micro_data = .true.
      else
         goto 120
      end if

      if( typeStr == "condMean" ) then
         if( size(freq_row_input_data, kind=our_int) &
              /= work%nrow_input_data ) goto 140
         do dr = 1, work%nrow_input_data
            if( freq_row_input_data(dr) < 0.D0 ) goto 145
         end do
         if( size(freq_row_input_data_int, kind=our_int) /= 0 ) goto 141
         if( size(f_mat, 1, kind=our_int) /= work%nrow_input_data ) goto 142
         if( size(f_mat, 2, kind=our_int) /= work%r ) goto 142
         do dr = 1, work%nrow_input_data
            do y = 1, work%r
               if( f_mat(dr,y) < 0.D0 ) goto 145
            end do
         end do
         if( size(f_mat_int, 1, kind=our_int) /= 0 ) goto 143
         if( size(f_mat_int, 2, kind=our_int) /= 0 ) goto 143
      else
         if( size(freq_row_input_data, kind=our_int) /= 0 ) goto 140
         if( size(freq_row_input_data_int, kind=our_int) &
              /= work%nrow_input_data ) goto 141
         do dr = 1, work%nrow_input_data
            if( freq_row_input_data_int(dr) < 0 ) goto 145
         end do
         if( size(f_mat, 1, kind=our_int) /= 0 ) goto 142
         if( size(f_mat, 2, kind=our_int) /= 0 ) goto 142
         if( size(f_mat_int, 1, kind=our_int) &
              /= work%nrow_input_data ) goto 143
         if( size(f_mat_int, 2, kind=our_int) /= work%r ) goto 143
         do dr = 1, work%nrow_input_data
            do y = 1, work%r
               if( f_mat_int(dr,y) < 0 ) goto 145
            end do
         end do
      end if
      if( wide_format ) then
         if( size(response_row_input_data, kind=our_int) /= 0 ) goto 150
      else
         if( size(response_row_input_data, kind=our_int) &
              /= work%nrow_input_data ) goto 150
         ! add check for elements
      end if
      if( work%saturated ) then
         if( size( fitted_pi, 1, kind=our_int ) /= work%n_cov_patt ) goto 300
         if( size( fitted_pi, 2, kind=our_int ) /= work%r ) goto 300
         work%pi_mat(:,:) = fitted_pi(:,:)
         if( size( fitted_pistar, 1, kind=our_int ) &
              /= work%n_cov_patt ) goto 305
         if( size( fitted_pistar, 2, kind=our_int ) /= work%r ) goto 305
         work%pistar_mat(:,:) = fitted_pistar(:,:)
      else
         if( size( coefficients, 1, kind=our_int ) /= work%p ) goto 310
         if( size( coefficients, 2, kind=our_int ) /= work%r ) goto 310
         if( stack_coef( coefficients, work%beta_vec, work ) &
              == RETURN_FAIL ) goto 800
      end if
      if( typeStr == "condMean" ) then
         if( size( cond_means, 1, kind=our_int ) &
              /= work%nrow_input_data ) goto 320
         if( size( cond_means, 2, kind=our_int ) /= work%r ) goto 320
         if( size( imp_mat, 1, kind=our_int ) /= 0 ) goto 330
         if( size( imp_mat, 2, kind=our_int ) /= 0 ) goto 330
         if( size( imp_vec, kind=our_int ) /= 0 ) goto 340
      else
         if( size( cond_means, 1, kind=our_int ) /= 0 ) goto 320
         if( size( cond_means, 2, kind=our_int ) /= 0 ) goto 320
         if( size( imp_mat, 1, kind=our_int ) &
              /= work%nrow_input_data ) goto 330
         if( size( imp_mat, 2, kind=our_int ) /= work%r ) goto 330
         if( size( imp_vec, kind=our_int ) /= work%nrow_input_data ) goto 340
      end if
      !####
      if( typeStr == "condMean" ) then
         if( run_rrlogit_impute_cond_mean( wide_format, micro_data, &
              freq_row_input_data, response_row_input_data, f_mat, &
              cond_means, work, err ) == RETURN_FAIL ) goto 800
      else if( typeStr == "random" ) then
         if( run_rrlogit_impute_random( wide_format, micro_data, &
              freq_row_input_data_int, response_row_input_data, f_mat_int, &
              imp_mat, imp_vec, work, err ) == RETURN_FAIL ) goto 800
      end if
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
100   call err_handle(err, 1, &
            comment = "Value of type_int not recognized" )
      goto 800
110   call err_handle(err, 1, &
            comment = "Value of wide_format_int not recognized" )
      goto 800
120   call err_handle(err, 1, &
            comment = "Value of micro_data_int not recognized" )
      goto 800
140   call err_handle(err, 1, &
            comment = "Array freq_row_input_data has incorrect size" )
      goto 800
141   call err_handle(err, 1, &
            comment = "Array freq_row_input_data_int has incorrect size" )
      goto 800
142   call err_handle(err, 1, &
            comment = "Array f_mat has incorrect size" )
      goto 800
143   call err_handle(err, 1, &
            comment = "Array f_mat_int has incorrect size" )
      goto 800
145   call err_handle(err, 1, &
            comment = "Negative frequency encountered" )
      goto 800
150   call err_handle(err, 1, &
            comment = "Array response_row_input_data has incorrect size" )
      goto 800
300   call err_handle(err, 1, &
            comment = "Array pi_mat has incorrect size" )
      goto 800
305   call err_handle(err, 1, &
            comment = "Array pistar_mat has incorrect size" )
      goto 800
310   call err_handle(err, 1, &
            comment = "Array coefficients has incorrect size" )
      goto 800
320   call err_handle(err, 1, &
            comment = "Array cond_means has incorrect size" )
      goto 800
330   call err_handle(err, 1, &
            comment = "Array imp_mat has incorrect size" )
      goto 800
340   call err_handle(err, 1, &
            comment = "Array imp_vec has incorrect size" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
      ijunk = nullify_workspace_type_rrlogit( work, err )
    end function run_rrlogit_impute
   !###################################################################
   integer(kind=our_int) function run_rrlogit_impute_cond_mean( &
        wide_format, micro_data, &
        freq_row_input_data, response_row_input_data, f_mat, &
        cond_means,  &
        work, err ) result(answer)
      implicit none
      ! inputs
      logical, intent(in) :: wide_format
      logical, intent(in) :: micro_data
      real(kind=our_dble), intent(in) :: freq_row_input_data(:)
      integer(kind=our_int), intent(in) :: response_row_input_data(:)
      real(kind=our_dble), intent(in) :: f_mat(:,:)
      ! outputs
      real(kind=our_dble), intent(out) :: cond_means(:,:)
      ! workspaces
      type(workspace_type_rrlogit), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! locals
      integer(kind=our_int) :: i, cp, dr, y, ystar
      integer(kind=our_int) :: status, ijunk
      real(kind=our_dble) :: sum
      logical, allocatable :: done(:)
      logical :: any_miss
      character(len=*), parameter :: &
           subname = "run_rrlogit_impute_cond_mean"
      ! begin
      answer = RETURN_FAIL
      !####
      allocate( done( work%r ), stat=status )
      if( status /= 0 ) goto 100
      if( .not. work%saturated ) then
         if( compute_pi_mat_rrlogit( work, err ) == RETURN_FAIL ) goto 800
      end if
      !####
      if( wide_format ) then
         do cp = 1, work%n_cov_patt
            any_miss = .false.
            do y = 1, work%r
               if( work%pi_mat(cp,y) == work%mvcode ) then
                  any_miss = .true.
                  exit
               end if
            end do
            if( any_miss ) then
               do i = work%cov_patt_st(cp), work%cov_patt_fin(cp)
                  dr = work%cov_patt_order(i)
                  cond_means(dr,:) = work%mvcode
               end do
               cycle
            end if
            done(:) = .false.
            do i = work%cov_patt_st(cp), work%cov_patt_fin(cp)
               dr = work%cov_patt_order(i)
               cond_means(dr,:) = 0.D0
               do ystar = 1, work%r
                  if( f_mat(dr,ystar) == 0.D0 ) cycle
                  if( .not. done(ystar) ) then
                     ! fill in row of phi_mat
                     sum = 0.D0
                     do y = 1, work%r
                        work%phi_mat(ystar,y) = work%pert_mat(ystar,y) * &
                             work%pi_mat(cp,y)
                        sum = sum + work%phi_mat(ystar,y)
                     end do
                     if( sum == 0.D0 ) goto 200
                     do y = 1, work%r
                        work%phi_mat(ystar,y) = work%phi_mat(ystar,y) / sum
                     end do
                     done(ystar) = .true.
                  end if
                  do y = 1, work%r
                     cond_means(dr,y) = cond_means(dr,y) + &
                          f_mat(dr,ystar) * work%phi_mat(ystar,y)
                  end do
               end do
            end do
         end do
      else
         if( micro_data ) then
            ! all frequencies assumed to be one
            do cp = 1, work%n_cov_patt
               any_miss = .false.
               do y = 1, work%r
                  if( work%pi_mat(cp,y) == work%mvcode ) then
                     any_miss = .true.
                     exit
                  end if
               end do
               if( any_miss ) then
                  do i = work%cov_patt_st(cp), work%cov_patt_fin(cp)
                     dr = work%cov_patt_order(i)
                     cond_means(dr,:) = work%mvcode
                  end do
                  cycle
               end if
               done(:) = .false.
               do i = work%cov_patt_st(cp), work%cov_patt_fin(cp)
                  dr = work%cov_patt_order(i)
                  ystar = response_row_input_data(dr)
                  if( .not. done(ystar) ) then
                     ! fill in row of phi_mat
                     sum = 0.D0
                     do y = 1, work%r
                        work%phi_mat(ystar,y) = work%pert_mat(ystar,y) * &
                             work%pi_mat(cp,y)
                        sum = sum + work%phi_mat(ystar,y)
                     end do
                     if( sum == 0.D0 ) goto 200
                     do y = 1, work%r
                        work%phi_mat(ystar,y) = work%phi_mat(ystar,y) / sum
                     end do
                     done(ystar) = .true.
                  end if
                  do y = 1, work%r
                     cond_means(dr,y) = work%phi_mat(ystar,y)
                  end do
               end do
            end do
         else
            ! narrow format with frequencies
            do cp = 1, work%n_cov_patt
               any_miss = .false.
               do y = 1, work%r
                  if( work%pi_mat(cp,y) == work%mvcode ) then
                     any_miss = .true.
                     exit
                  end if
               end do
               if( any_miss ) then
                  do i = work%cov_patt_st(cp), work%cov_patt_fin(cp)
                     dr = work%cov_patt_order(i)
                     cond_means(dr,:) = work%mvcode
                  end do
                  cycle
               end if
               done(:) = .false.
               do i = work%cov_patt_st(cp), work%cov_patt_fin(cp)
                  dr = work%cov_patt_order(i)
                  if( freq_row_input_data(dr) == 0.D0 ) then
                     cond_means(dr,:) = 0.D0
                     cycle
                  end if
                  ystar = response_row_input_data(dr)
                  if( .not. done(ystar) ) then
                     ! fill in row of phi_mat
                     sum = 0.D0
                     do y = 1, work%r
                        work%phi_mat(ystar,y) = work%pert_mat(ystar,y) * &
                             work%pi_mat(cp,y)
                        sum = sum + work%phi_mat(ystar,y)
                     end do
                     if( sum == 0.D0 ) goto 200
                     do y = 1, work%r
                        work%phi_mat(ystar,y) = work%phi_mat(ystar,y) / sum
                     end do
                     done(ystar) = .true.
                  end if
                  do y = 1, work%r
                     cond_means(dr,y) = freq_row_input_data(dr) * &
                          work%phi_mat(ystar,y)
                  end do
               end do
            end do
         end if
      end if
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
100   call err_handle(err, 1, &
            comment = "Unable to allocate array" )
      goto 800
200   call err_handle(err, 1, &
            comment = "Attempted division by zero" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
      if( allocated( done ) ) deallocate( done, stat=ijunk )
    end function run_rrlogit_impute_cond_mean
    !###################################################################
    integer(kind=our_int) function run_rrlogit_impute_random( &
        wide_format, micro_data, &
        freq_row_input_data_int, response_row_input_data, f_mat_int, &
        imp_mat, imp_vec,  &
        work, err, skip_compute_pi_mat ) result(answer)
      implicit none
      ! inputs
      logical, intent(in) :: wide_format
      logical, intent(in) :: micro_data
      integer(kind=our_int), intent(in) :: freq_row_input_data_int(:)
      integer(kind=our_int), intent(in) :: response_row_input_data(:)
      integer(kind=our_int), intent(in) :: f_mat_int(:,:)
      ! outputs
      integer(kind=our_int), intent(out) :: imp_mat(:,:)
      integer(kind=our_int), intent(out) :: imp_vec(:)
      ! workspaces
      type(workspace_type_rrlogit), intent(inout) :: work
      type(error_type), intent(inout) :: err
      ! optionals
      logical, intent(in), optional :: skip_compute_pi_mat
      ! locals
      integer(kind=our_int) :: i, cp, dr, y, ystar
      integer(kind=our_int) :: status, ijunk
      real(kind=our_dble) :: sum
      logical, allocatable :: done(:)
      logical :: any_miss, skip_compute_pi_mat_local
      integer(kind=our_int), allocatable :: harvest(:)
      character(len=*), parameter :: &
           subname = "run_rrlogit_impute_random"
      ! begin
      answer = RETURN_FAIL
      !####
      if( present(skip_compute_pi_mat) ) then
         skip_compute_pi_mat_local = skip_compute_pi_mat
      else
         skip_compute_pi_mat_local = .false.
      end if
      allocate( done( work%r ), stat=status )
      if( status /= 0 ) goto 100
      allocate( harvest( work%r ), stat=status )
      if( status /= 0 ) goto 100
      if( ( .not. work%saturated ) .and. &
           ( .not. skip_compute_pi_mat_local ) ) then
         if( compute_pi_mat_rrlogit( work, err ) == RETURN_FAIL ) goto 800
      end if
      !####
      if( wide_format ) then
         imp_vec(:) = 0
         do cp = 1, work%n_cov_patt
            any_miss = .false.
            do y = 1, work%r
               if( work%pi_mat(cp,y) == work%mvcode ) then
                  any_miss = .true.
                  exit
               end if
            end do
            if( any_miss ) then
               do i = work%cov_patt_st(cp), work%cov_patt_fin(cp)
                  dr = work%cov_patt_order(i)
                  imp_mat(dr,:) = -1 ! integer missing value code 
               end do
               cycle
            end if
            done(:) = .false.
            do i = work%cov_patt_st(cp), work%cov_patt_fin(cp)
               dr = work%cov_patt_order(i)
               imp_mat(dr,:) = 0
               do ystar = 1, work%r
                  if( f_mat_int(dr,ystar) == 0 ) cycle
                  if( .not. done(ystar) ) then
                     ! fill in row of phi_mat
                     sum = 0.D0
                     do y = 1, work%r
                        work%phi_mat(ystar,y) = work%pert_mat(ystar,y) * &
                             work%pi_mat(cp,y)
                        sum = sum + work%phi_mat(ystar,y)
                     end do
                     if( sum == 0.D0 ) goto 200
                     do y = 1, work%r
                        work%phi_mat(ystar,y) = work%phi_mat(ystar,y) / sum
                     end do
                     done(ystar) = .true.
                  end if
                  work%wkrA(:) = work%phi_mat(ystar,:)
                  if( rmult_n( f_mat_int(dr,ystar), work%wkrA, harvest, err) &
                       == RETURN_FAIL ) goto 800
                  imp_mat(dr,:) = imp_mat(dr,:) + harvest
               end do
            end do
         end do
      else
         if( micro_data ) then
            ! all frequencies assumed to be one
            imp_mat(:,:) = 0
            do cp = 1, work%n_cov_patt
               any_miss = .false.
               do y = 1, work%r
                  if( work%pi_mat(cp,y) == work%mvcode ) then
                     any_miss = .true.
                     exit
                  end if
               end do
               if( any_miss ) then
                  do i = work%cov_patt_st(cp), work%cov_patt_fin(cp)
                     dr = work%cov_patt_order(i)
                     imp_vec(dr) = -1 ! integer missing value code 
                  end do
                  cycle
               end if
               done(:) = .false.
               do i = work%cov_patt_st(cp), work%cov_patt_fin(cp)
                  dr = work%cov_patt_order(i)
                  ystar = response_row_input_data(dr)
                  if( .not. done(ystar) ) then
                     ! fill in row of phi_mat
                     sum = 0.D0
                     do y = 1, work%r
                        work%phi_mat(ystar,y) = work%pert_mat(ystar,y) * &
                             work%pi_mat(cp,y)
                        sum = sum + work%phi_mat(ystar,y)
                     end do
                     if( sum == 0.D0 ) goto 200
                     do y = 1, work%r
                        work%phi_mat(ystar,y) = work%phi_mat(ystar,y) / sum
                     end do
                     done(ystar) = .true.
                  end if
                  work%wkrA(:) = work%phi_mat(ystar,:)
                  if( rmult_one( work%wkrA, y, err ) == RETURN_FAIL ) goto 800
                  imp_vec(dr) = y
               end do
            end do
         else
            ! narrow format with frequencies
            imp_vec(:) = 0
            do cp = 1, work%n_cov_patt
               any_miss = .false.
               do y = 1, work%r
                  if( work%pi_mat(cp,y) == work%mvcode ) then
                     any_miss = .true.
                     exit
                  end if
               end do
               if( any_miss ) then
                  do i = work%cov_patt_st(cp), work%cov_patt_fin(cp)
                     dr = work%cov_patt_order(i)
                     imp_mat(dr,:) = -1 ! integer missing value code 
                  end do
                  cycle
               end if
               done(:) = .false.
               do i = work%cov_patt_st(cp), work%cov_patt_fin(cp)
                  dr = work%cov_patt_order(i)
                  if( freq_row_input_data_int(dr) == 0 ) then
                     imp_mat(dr,:) = 0
                     cycle
                  end if
                  ystar = response_row_input_data(dr)
                  if( .not. done(ystar) ) then
                     ! fill in row of phi_mat
                     sum = 0.D0
                     do y = 1, work%r
                        work%phi_mat(ystar,y) = work%pert_mat(ystar,y) * &
                             work%pi_mat(cp,y)
                        sum = sum + work%phi_mat(ystar,y)
                     end do
                     if( sum == 0.D0 ) goto 200
                     do y = 1, work%r
                        work%phi_mat(ystar,y) = work%phi_mat(ystar,y) / sum
                     end do
                     done(ystar) = .true.
                  end if
                  work%wkrA(:) = work%phi_mat(ystar,:)
                  if( rmult_n( freq_row_input_data_int(dr), &
                       work%wkrA, harvest, err) == RETURN_FAIL ) goto 800
                  imp_mat(dr,:) = harvest
               end do
            end do
         end if
      end if
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
100   call err_handle(err, 1, &
            comment = "Unable to allocate array" )
      goto 800
200   call err_handle(err, 1, &
            comment = "Attempted division by zero" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
      if( allocated( done ) ) deallocate( done, stat=ijunk )
      if( allocated( harvest ) ) deallocate( harvest, stat=ijunk )
    end function run_rrlogit_impute_random
    !###################################################################
    integer(kind=our_int) function rmult_n(n, p, harvest, err) &
         result(answer)
      ! draws a random vector from Mult(n,p) by repeatedly drawing
      ! from binomials
      ! p is automatically scaled to sum to one
      implicit none
      ! inputs
      integer(kind=our_int), intent(in) :: n
      real(kind=our_dble), intent(in) :: p(:)
      ! output
      integer(kind=our_int), intent(out) :: harvest(:)
      ! workspaces
      type(error_type), intent(inout) :: err
      ! locals
      integer(kind=our_int) :: i, n_remaining, rb
      real(kind=our_dble) :: sum, ptmp, ntmp, rtmp
      character(len=*), parameter :: &
           subname = "rmult_n"
      ! begin
      answer = RETURN_FAIL
      !####
      if( n < 0 ) goto 50
      if( size(p, kind=our_int) /= size(harvest, kind=our_int) ) goto 100
      sum = 0.D0
      do i = 1, size(p, kind=our_int)
         if( p(i) < 0.D0 ) goto 150
         sum = sum + p(i)
      end do
      if( ( size(p, kind=our_int) > 0 ) .and. ( sum == 0.D0 ) ) goto 200
      !####
      n_remaining = n
      harvest(:) = 0
      do i = 1, size(p, kind=our_int)
         if( i == size(p, kind=our_int) ) then
            harvest(i) = n_remaining
            exit
         end if
         if( sum == 0.D0 ) goto 250
         ptmp = p(i) / sum
         if( ptmp < 0.D0 ) ptmp = 0.D0  ! correction for rounding error
         if( ptmp > 1.D0 ) ptmp = 1.D0  ! correction for rounding error
         ntmp = real( n_remaining, our_dble )
         if( rbinom_R( ntmp, ptmp, rtmp, err ) == RETURN_FAIL ) goto 800
         rb = int(rtmp, our_int)
         harvest(i) = rb
         n_remaining = n_remaining - rb
         sum = sum - p(i)
         if( n_remaining < 0 ) goto 300
         if( n_remaining == 0 ) exit
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
50    call err_handle(err, 1, &
           comment = "Negative value for n" )
      goto 800
100   call err_handle(err, 1, &
            comment = "Sizes of p and harvest mismatched" )
      goto 800
150   call err_handle(err, 1, &
            comment = "Negative probability encountered" )
      goto 800
200   call err_handle(err, 1, &
            comment = "All probabilities are zero" )
      goto 800
250   call err_handle(err, 1, &
            comment = "Attempted division by zero" )
      goto 800
300   call err_handle(err, 1, &
            comment = "Drawn binomials exceeded limit" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function rmult_n
    !###################################################################
    integer(kind=our_int) function rmult_one(p, harvest, err) &
         result(answer)
      ! draw from Mult(1,p) expressed as a single integer
      ! p is automatically scaled to sum to one
      implicit none
      ! inputs
      real(kind=our_dble), intent(in) :: p(:)
      ! output
      integer(kind=our_int), intent(out) :: harvest
      ! workspaces
      type(error_type), intent(inout) :: err
      ! locals
      integer(kind=our_int) :: i
      real(kind=our_dble) :: sum, rtmp, num
      character(len=*), parameter :: &
           subname = "rmult_one"
      ! begin
      answer = RETURN_FAIL
      !####
      sum = 0.D0
      do i = 1, size(p, kind=our_int)
         if( p(i) < 0.D0 ) goto 150
         sum = sum + p(i)
      end do
      if( ( size(p, kind=our_int) > 0 ) .and. ( sum == 0.D0 ) ) goto 200
      !####
      if( runif_R( rtmp, err ) == RETURN_FAIL ) goto 800
      harvest = 0
      num = 0.D0
      do i = 1, size(p, kind=our_int)
         if( i == size(p, kind=our_int) ) then
            harvest = i
            exit
         end if
         num = num + p(i)
         if( rtmp <= ( num / sum ) ) then
            harvest = i
            exit
         end if
      end do
      ! normal exit
      answer = RETURN_SUCCESS
      goto 999
      ! error traps
150   call err_handle(err, 1, &
            comment = "Negative probability encountered" )
      goto 800
200   call err_handle(err, 1, &
            comment = "All probabilities are zero" )
      goto 800
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      goto 999
      ! cleanup
999   continue
    end function rmult_one
    !###################################################################
end module rrlogit_engine
!#####################################################################
