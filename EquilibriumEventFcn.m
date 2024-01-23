function [position,isterminal,direction] = EquilibriumEventFcn(t,x,DE_params_equil,Init,t_scale)

%%%%%%%%%%%%%%%%%%%%
% Specifically, this function is used to find equilibrium PD1, PD-L1 and 
% PD1:PD-L1 complex levels in the untreated case.
%%%%%%%%%%%%%%%%%%%%

isterminal = 1;
direction = [];
threshold = 0.005;

%--------------------------------------------------------------------------
dxdt = integratingfunction_equilibrium(t,x,DE_params_equil,Init);
%--------------------------------------------------------------------------
Init = [Init(1); Init(2); Init(1)];
derivatives_vec=Init.*dxdt/t_scale;
norm_vec=norm(derivatives_vec);

% Logical(cond) evaluates to 1 if cond is true and to 0 if cond is false.
% So we want to make the condition the opposite of what we actually want,
% e.g. check if norm_vec is above the threshold. Then when it is below,
% position will drop to zero and the integration will stop.
if norm_vec > threshold
   position = 1; 
else
   position = 0;
end

end

