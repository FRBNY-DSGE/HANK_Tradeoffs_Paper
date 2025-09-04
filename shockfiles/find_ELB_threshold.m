
ELB_threshold = ELB_threshold;

ELB_frequencies_final = zeros(num_simul_noexplosive,1);

for jj = 1:num_simul_noexplosive

	ELB_frequencies_final(jj,1) = sum(X_vec(jj).ELB_Taylor_hits);

end

num_simul_trimmed_aux = 0;

while num_simul_trimmed_aux < num_simul_thresholds 

	ELB_threshold = ELB_threshold - 1;

	ELB_large_idx_aux = nonzeros((ELB_frequencies_final>ELB_threshold).*[1:1:num_simul_noexplosive]');	

	num_simul_trimmed_aux = length(ELB_large_idx_aux)*refresh_size;

end