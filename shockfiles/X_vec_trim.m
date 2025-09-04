ELB_frequencies_final = zeros(num_simul_noexplosive,1);

for jj = 1:num_simul_noexplosive

	ELB_frequencies_final(jj,1) = sum(X_vec(jj).ELB_Taylor_hits);

end

ELB_large_idx = nonzeros((ELB_frequencies_final>ELB_threshold).*[1:1:num_simul_noexplosive]');

X_Taylor_series_trim     = zeros(grid.num_endo,length(ELB_large_idx)*refresh_size);
ELB_Taylor_hits_trim     = zeros(1,length(ELB_large_idx)*refresh_size);
ELB_Taylor_expected_trim = zeros(1,length(ELB_large_idx)*refresh_size);

shock_series_trim        = zeros(grid.numstates_shocks,length(ELB_large_idx)*refresh_size);
shock_series_raw_trim    = zeros(grid.numstates_shocks,length(ELB_large_idx)*refresh_size);


for jj = 1:length(ELB_large_idx)

	idx_aux_ELB = ELB_large_idx(jj);

	X_Taylor_series_trim(:,(jj-1)*refresh_size+1:jj*refresh_size)     = X_vec(idx_aux_ELB).X_Taylor_series;
    ELB_Taylor_hits_trim(:,(jj-1)*refresh_size+1:jj*refresh_size)     = X_vec(idx_aux_ELB).ELB_Taylor_hits;
    ELB_Taylor_expected_trim(:,(jj-1)*refresh_size+1:jj*refresh_size) = X_vec(idx_aux_ELB).ELB_Taylor_expected;
    shock_series_trim(:,(jj-1)*refresh_size+1:jj*refresh_size)        = simul_shock_series(idx_aux_ELB).shock_series;
    shock_series_raw_trim(:,(jj-1)*refresh_size+1:jj*refresh_size)    = simul_shock_series_raw(idx_aux_ELB).shock_series;
    simul_shock_series_trim(jj).shock_series                          = simul_shock_series(idx_aux_ELB).shock_series;
    simul_shock_series_raw_trim(jj).shock_series                      = simul_shock_series_raw(idx_aux_ELB).shock_series;

end