function [regime regimestart] = map_regime(violvecbool)

nperiods = length(violvecbool)-1;

% analyse violvec and isolate contiguous periods in the other regime

regime(1)      = violvecbool(1);
regimeindx     = 1;
regimestart(1) = 1;

for i = 2:nperiods

	if violvecbool(i)~=regime(regimeindx)
		regimeindx = regimeindx + 1;
		regime(regimeindx) = violvecbool(i);
		regimestart(regimeindx) = i;
	end
		
end

% if (regime(1)==1 && length(regimestart)==1) % simulation starts at the alternative regime and does not converge in the currently guessed number of periods
% 	warning('Increase nperiods');
% end

% if (regime(end)==1)                         % simulation starts at the reference regime and moves to the alternative regime at the end of the guessed number of periods
% 	warning('Increase nperiods');
% end
