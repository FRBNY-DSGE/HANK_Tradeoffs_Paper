
ELB_frequencies_final = zeros(num_simul_noexplosive,1);

for jj = 1:num_simul_noexplosive

	ELB_frequencies_final(jj,1) = sum(X_vec(jj).ELB_Taylor_hits);

end

ELB_large_idx = nonzeros((ELB_frequencies_final>ELB_threshold).*[1:1:num_simul_noexplosive]');

X_Taylor_series_trim     = zeros(grid.num_endo,length(ELB_large_idx)*refresh_size);
ELB_Taylor_hits_trim     = zeros(1,length(ELB_large_idx)*refresh_size);
ELB_Taylor_expected_trim = zeros(1,length(ELB_large_idx)*refresh_size);

% X_AIT_series_trim     = zeros(grid.num_endo,length(ELB_large_idx)*refresh_size);
% ELB_AIT_hits_trim     = zeros(1,length(ELB_large_idx)*refresh_size);
% ELB_AIT_expected_trim = zeros(1,length(ELB_large_idx)*refresh_size);

shock_series_trim        = zeros(grid.numstates_shocks,length(ELB_large_idx)*refresh_size);
shock_series_raw_trim    = zeros(grid.numstates_shocks,length(ELB_large_idx)*refresh_size);


for jj = 1:length(ELB_large_idx)

	idx_aux_ELB = ELB_large_idx(jj);

	X_Taylor_series_trim(:,(jj-1)*refresh_size+1:jj*refresh_size)     = X_vec(idx_aux_ELB).X_Taylor_series;
    ELB_Taylor_hits_trim(:,(jj-1)*refresh_size+1:jj*refresh_size)     = X_vec(idx_aux_ELB).ELB_Taylor_hits;
    ELB_Taylor_expected_trim(:,(jj-1)*refresh_size+1:jj*refresh_size) = X_vec(idx_aux_ELB).ELB_Taylor_expected;
    % X_AIT_series_trim(:,(jj-1)*refresh_size+1:jj*refresh_size)     = X_vec(idx_aux_ELB).X_AIT_series;
    % ELB_AIT_hits_trim(:,(jj-1)*refresh_size+1:jj*refresh_size)     = X_vec(idx_aux_ELB).ELB_AIT_hits;
    % ELB_AIT_expected_trim(:,(jj-1)*refresh_size+1:jj*refresh_size) = X_vec(idx_aux_ELB).ELB_AIT_expected;
    shock_series_trim(:,(jj-1)*refresh_size+1:jj*refresh_size)        = simul_shock_series(idx_aux_ELB).shock_series;
    shock_series_raw_trim(:,(jj-1)*refresh_size+1:jj*refresh_size)    = simul_shock_series_raw(idx_aux_ELB).shock_series;
    X_vec_trim(jj).X_Taylor_series      = X_vec(idx_aux_ELB).X_Taylor_series;
    X_vec_trim(jj).ELB_Taylor_hits      = X_vec(idx_aux_ELB).ELB_Taylor_hits;
    X_vec_trim(jj).ELB_Taylor_expected  = X_vec(idx_aux_ELB).ELB_Taylor_expected;
    % X_vec_trim(jj).X_AIT_series      = X_vec(idx_aux_ELB).X_AIT_series;
    % X_vec_trim(jj).ELB_AIT_hits      = X_vec(idx_aux_ELB).ELB_AIT_hits;
    % X_vec_trim(jj).ELB_AIT_expected  = X_vec(idx_aux_ELB).ELB_AIT_expected;
    simul_shock_series_trim(jj).shock_series     = simul_shock_series(idx_aux_ELB).shock_series;
    simul_shock_series_raw_trim(jj).shock_series = simul_shock_series_raw(idx_aux_ELB).shock_series;


end


results_trim_cd4_extra_main.X_Taylor_series_trim = X_Taylor_series_trim;
results_trim_cd4_extra_main.ELB_Taylor_hits_trim = ELB_Taylor_hits_trim;
results_trim_cd4_extra_main.ELB_Taylor_expected_trim = ELB_Taylor_expected_trim;
% results_trim_cd4_extra_main.X_AIT_series_trim = X_AIT_series_trim;
% results_trim_cd4_extra_main.ELB_AIT_hits_trim = ELB_AIT_hits_trim;
% results_trim_cd4_extra_main.ELB_AIT_expected_trim = ELB_AIT_expected_trim;
results_trim_cd4_extra_main.shock_series_trim = shock_series_trim;
results_trim_cd4_extra_main.shock_series_raw_trim = shock_series_raw_trim;
results_trim_cd4_extra_main.X_vec_trim = X_vec_trim;
results_trim_cd4_extra_main.X_vec_trim = X_vec_trim;
results_trim_cd4_extra_main.X_vec_trim = X_vec_trim;
results_trim_cd4_extra_main.X_vec_trim = X_vec_trim;
results_trim_cd4_extra_main.X_vec_trim = X_vec_trim;
results_trim_cd4_extra_main.X_vec_trim = X_vec_trim;
results_trim_cd4_extra_main.simul_shock_series_trim = simul_shock_series_trim;
results_trim_cd4_extra_main.simul_shock_series_raw_trim = simul_shock_series_raw_trim;