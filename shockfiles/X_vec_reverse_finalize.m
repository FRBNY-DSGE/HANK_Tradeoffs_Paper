
% (1) Finds shocks in Trim 2 that are not explosive in Trim 3

idx_common_noexplosive = nonzeros((1-results_trim_reverse.Explosive_check).*[1:1:length(results_trim_reverse.Explosive_check)]');

% idx_common_noexplosive = 

for jj = 1:length(idx_common_noexplosive)

	kk = idx_common_noexplosive(jj);

	X_vec_tmp(jj).X_Taylor_series               = results_trim.X_vec_trim(kk).X_Taylor_series;
	X_vec_tmp(jj).ELB_Taylor_hits               = results_trim.X_vec_trim(kk).ELB_Taylor_hits;
	X_vec_tmp(jj).ELB_Taylor_expected           = results_trim.X_vec_trim(kk).ELB_Taylor_expected;
	simul_shock_series_tmp(jj).shock_series     = results_trim.simul_shock_series_trim(kk).shock_series;
	simul_shock_series_raw_tmp(jj).shock_series = results_trim.simul_shock_series_raw_trim(kk).shock_series;

	X_vec_reverse_tmp(jj).X_Taylor_series               = results_trim_reverse.X_vec_trim(kk).X_Taylor_series;
	X_vec_reverse_tmp(jj).ELB_Taylor_hits               = results_trim_reverse.X_vec_trim(kk).ELB_Taylor_hits;
	X_vec_reverse_tmp(jj).ELB_Taylor_expected           = results_trim_reverse.X_vec_trim(kk).ELB_Taylor_expected;
	simul_shock_series_reverse_tmp(jj).shock_series     = results_trim_reverse.simul_shock_series_trim(kk).shock_series;
	simul_shock_series_raw_reverse_tmp(jj).shock_series = results_trim_reverse.simul_shock_series_raw_trim(kk).shock_series;

end

% (2) Find ZLB threshold in Trim 2 that ensures at least 30k period sequence of
% shocks

ELB_frequencies_final = zeros(num_simul_noexplosive,1);

for jj = 1:num_simul_noexplosive

	% ELB_frequencies_final(jj,1) = sum(X_vec_tmp(jj).ELB_Taylor_hits);
	ELB_frequencies_final(jj,1) = sum(X_vec_reverse_tmp(jj).ELB_Taylor_hits);

end

num_simul_trimmed_aux = 0;

while num_simul_trimmed_aux < num_simul_thresholds 

	ELB_threshold = ELB_threshold - 1;

	ELB_large_idx_aux = nonzeros((ELB_frequencies_final>ELB_threshold).*[1:1:num_simul_noexplosive]');	

	num_simul_trimmed_aux = length(ELB_large_idx_aux)*refresh_size;

end

% (3) Obtain blocks that meet or exceed threshold

ELB_large_idx = nonzeros((ELB_frequencies_final>ELB_threshold).*[1:1:num_simul_noexplosive]');
ELB_large_idx = ELB_large_idx(1:100);

for jj = 1:length(ELB_large_idx)

	idx_aux_ELB = ELB_large_idx(jj);

    X_vec(jj).X_Taylor_series     = X_vec_tmp(idx_aux_ELB).X_Taylor_series;
    X_vec(jj).ELB_Taylor_hits     = X_vec_tmp(idx_aux_ELB).ELB_Taylor_hits;
    X_vec(jj).ELB_Taylor_expected = X_vec_tmp(idx_aux_ELB).ELB_Taylor_expected;
    simul_shock_series(jj).shock_series     = simul_shock_series_tmp(idx_aux_ELB).shock_series;
    simul_shock_series_raw(jj).shock_series = simul_shock_series_raw_tmp(idx_aux_ELB).shock_series;   

    X_vec_reverse(jj).X_Taylor_series     = X_vec_reverse_tmp(idx_aux_ELB).X_Taylor_series;
    X_vec_reverse(jj).ELB_Taylor_hits     = X_vec_reverse_tmp(idx_aux_ELB).ELB_Taylor_hits;
    X_vec_reverse(jj).ELB_Taylor_expected = X_vec_reverse_tmp(idx_aux_ELB).ELB_Taylor_expected;
    simul_shock_series_reverse(jj).shock_series     = simul_shock_series_reverse_tmp(idx_aux_ELB).shock_series;
    simul_shock_series_raw_reverse(jj).shock_series = simul_shock_series_raw_reverse_tmp(idx_aux_ELB).shock_series;    

end
%% 

% (4) Turn 300-period blocks into continuous 60k period sequence and append
% sign-reversed shocks


X_Taylor_series_trim     = zeros(grid.num_endo,length(ELB_large_idx)*refresh_size);
ELB_Taylor_hits_trim     = zeros(1,length(ELB_large_idx)*refresh_size);
ELB_Taylor_expected_trim = zeros(1,length(ELB_large_idx)*refresh_size);
shock_series_trim        = zeros(grid.numstates_shocks,length(ELB_large_idx)*refresh_size);
shock_series_raw_trim    = zeros(grid.numstates_shocks,length(ELB_large_idx)*refresh_size);

X_Taylor_series_trim_reverse     = zeros(grid.num_endo,length(ELB_large_idx)*refresh_size);
ELB_Taylor_hits_trim_reverse     = zeros(1,length(ELB_large_idx)*refresh_size);
ELB_Taylor_expected_trim_reverse = zeros(1,length(ELB_large_idx)*refresh_size);
shock_series_trim_reverse        = zeros(grid.numstates_shocks,length(ELB_large_idx)*refresh_size);
shock_series_raw_trim_reverse    = zeros(grid.numstates_shocks,length(ELB_large_idx)*refresh_size);

for jj = 1:length(ELB_large_idx)

	X_Taylor_series_trim(:,(jj-1)*refresh_size+1:jj*refresh_size)     = X_vec(jj).X_Taylor_series;
	ELB_Taylor_hits_trim(:,(jj-1)*refresh_size+1:jj*refresh_size)     = X_vec(jj).ELB_Taylor_hits;
	ELB_Taylor_expected_trim(:,(jj-1)*refresh_size+1:jj*refresh_size) = X_vec(jj).ELB_Taylor_expected;
	shock_series_trim(:,(jj-1)*refresh_size+1:jj*refresh_size)        = simul_shock_series(jj).shock_series;
	shock_series_raw_trim(:,(jj-1)*refresh_size+1:jj*refresh_size)    = simul_shock_series_raw(jj).shock_series;

	X_Taylor_series_trim_reverse(:,(jj-1)*refresh_size+1:jj*refresh_size)     = X_vec_reverse(jj).X_Taylor_series;
	ELB_Taylor_hits_trim_reverse(:,(jj-1)*refresh_size+1:jj*refresh_size)     = X_vec_reverse(jj).ELB_Taylor_hits;
	ELB_Taylor_expected_trim_reverse(:,(jj-1)*refresh_size+1:jj*refresh_size) = X_vec_reverse(jj).ELB_Taylor_expected;
	shock_series_trim_reverse(:,(jj-1)*refresh_size+1:jj*refresh_size)        = simul_shock_series_reverse(jj).shock_series;
	shock_series_raw_trim_reverse(:,(jj-1)*refresh_size+1:jj*refresh_size)    = simul_shock_series_raw_reverse(jj).shock_series;

end

% X_Taylor_series_final     = [X_Taylor_series_trim,X_Taylor_series_trim_reverse];
% ELB_Taylor_hits_final     = [ELB_Taylor_hits_trim,ELB_Taylor_hits_trim_reverse];
% ELB_Taylor_expected_final = [ELB_Taylor_expected_trim,ELB_Taylor_expected_trim_reverse];
shock_series_final        = [shock_series_trim,shock_series_trim_reverse];
shock_series_raw_final    = [shock_series_raw_trim,shock_series_raw_trim_reverse];