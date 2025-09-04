

mean_covid_shocks = mean(covid_shocks,2);
sum_covid_shocks  = sum(covid_shocks,2);

min_r_aux = 0;

while min_r_aux < 0.3

	r_aux_covid = unifrnd(-3,3,9,num_cooldown_periods)+1;

	mean_r_aux_covid = mean(r_aux_covid,2);

	min_r_aux = min(mean_r_aux_covid);

end

% r_aux_covid = r_aux_covid./repmat(mean_r_aux_covid*num_cooldown_periods/6,[1 num_cooldown_periods]);

% covid_cooldown_shocks = -r_aux_covid.*(repmat(mean_covid_shocks(2:end),[1 num_cooldown_periods]));

% r_aux_covid = r_aux_covid./repmat(mean_r_aux_covid*num_cooldown_periods/6,[1 num_cooldown_periods]);

r_aux_covid = r_aux_covid./repmat(mean_r_aux_covid,[1 num_cooldown_periods]);

covid_cooldown_shocks = -r_aux_covid.*(repmat(mean_covid_shocks(2:end)*6/num_cooldown_periods,[1 num_cooldown_periods]));

% covid_cooldown_shocks = -r_aux_covid.*(repmat(sum_covid_shocks(2:end),[1 num_cooldown_periods]));

% covid_cooldown_shocks/num_cooldown_periods/6

new_covid_shocks = [covid_shocks,[zeros(1,num_cooldown_periods);covid_cooldown_shocks]];