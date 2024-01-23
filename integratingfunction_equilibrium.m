function dxdt = integratingfunction_equilibrium(t,x,params,Init)

% Specifically, this function is used to simulate equilibration of PD-1 and 
% PD-L1 levels on the PD1-PDL1 pathway.
%--------------------------------------------------------------------------
% Variable definitions:
%--------------------------------------------------------------------------
%Proteins:
%PD1=x(1); PD-L1=x(2); PD-1:PD-L1=x(3); 

%--------------------------------------------------------------------------
% Parameters
%--------------------------------------------------------------------------
% Read in parameters
%--------------------------------------------------------------------------
% complex association/disassociation

beta_plus_tilde = params(1); beta_minus_tilde = params(2); 
%--------------------------------------------------------------------------
PD1_0 = Init(1); PDL1_0 = Init(2);

%--------------------------------------------------------------------------
dxdt = zeros(max(size(x)),1);
%--------------------------------------------------------------------------
% Equations defining output:  
%--------------------------------------------------------------------------

% Normalized to PD1_0               
%(1) PD-1: 
dxdt(1) = -beta_plus_tilde*PDL1_0*x(1)*x(2) + beta_minus_tilde*x(3);
            
% Normalized to PDL1_0
%(2) PD-L1: 
dxdt(2) = -beta_plus_tilde*PD1_0*x(1)*x(2) + beta_minus_tilde*PD1_0*x(3)/PDL1_0;
 
% Normalized to PD1_0
%(3) PD-1:PD-L1: 
dxdt(3) = beta_plus_tilde*PDL1_0*x(1)*x(2) - beta_minus_tilde*x(3);
            
end
