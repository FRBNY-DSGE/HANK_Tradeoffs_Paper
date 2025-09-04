

mean_GR_shocks = mean(gr_shocks,2);
sum_GR_shocks  = sum(gr_shocks,2);

min_r_aux = 0;

while min_r_aux < 0.3

	r_aux_GR = unifrnd(-3,3,9,num_cooldown_periods)+1;

	mean_r_aux_GR = mean(r_aux_GR,2);

	min_r_aux = min(mean_r_aux_GR);

end

% r_aux_GR = r_aux_GR./repmat(mean_r_aux_GR*num_cooldown_periods/6,[1 num_cooldown_periods]);

% gr_cooldown_shocks = -r_aux_GR.*(repmat(mean_GR_shocks(2:end),[1 num_cooldown_periods]));

% r_aux_GR = r_aux_GR./repmat(mean_r_aux_GR*num_cooldown_periods/6,[1 num_cooldown_periods]);

r_aux_GR = r_aux_GR./repmat(mean_r_aux_GR,[1 num_cooldown_periods]);

gr_cooldown_shocks = -r_aux_GR.*(repmat(mean_GR_shocks(2:end)*6/num_cooldown_periods,[1 num_cooldown_periods]));

% gr_cooldown_shocks = -r_aux_GR.*(repmat(sum_GR_shocks(2:end),[1 num_cooldown_periods]));

% gr_cooldown_shocks/num_cooldown_periods/6

new_gr_shocks = [gr_shocks,[zeros(1,num_cooldown_periods);gr_cooldown_shocks]];