% Name: Main.m

function Main
t_scale = 12*60;  %Time scale
total_cells_flow = 1e8;
nivo_0 = 0.17e11; %pg/ml
%--------------------------------------------------------------------------
% Baseline Parameter Values
%--------------------------------------------------------------------------

pars = [0.45, 1.54, 3339.16, 11,...
        1000, 0.0175, 0.0175, 0.024, 0.024,...
        0.0205, 0.035, 0.019, 0.0175,...
        0.018, 0.01, 1e-10, 1e-10,...
        0.095, 0.018, 0.01047, 0.000125,...
        0.01134, 0.000005470, 0.000005470, 1.9e-6,...
        0.00065, 0.00008, 0.00001, 0.0006046,...
        1.5e-8, 0.5, 0.5, 8.5,...
        0.3465, 0.02405, 0.006, 0.252,...
        159.05, 0.4049, 0.6307, 4.0322,...
        0.8436, 130.04, 0.006307, 0.03447,...
        4, 0.7, 0.08895, 0.7467, 141.94,...
        1000, 7127.55, 7198.54, 0.1266,...
        0.0008437, 700, 0.0007, 48.515,...
        2.1451, 2, 19.433, 0.97322,...
        0.00486, 0.00973, 0.00121, 1, 1.4, 1.84, 0.05,...
        0.5e-4, 0.32, 5.9, 2, 2, 2,...
        1.8, 50, 0.5, 1.9, 0.5, 0.66, 2.425e-8,...
        6500, 1600, 650, 23, 160, 4.5e8,...
        0.75, 0.75, 0.75, 0.75,...
        10000, 0.00523191, 0.00523191, 0.005, 0.005, 0.005, 0.005, 0.005, 0.05, 0.085, 0.07, 10000];
%--------------------------------------------------------------------------
% Initial cytokine levels
%--------------------------------------------------------------------------
IFNg_0 = pars(1);  
IL12_0 = pars(2); 
IL6_0 = pars(3); 
IL4_0 = pars(4);    
TGF_beta_0 = pars(5);
Ec_0 = pars(93);
m214_0 = pars(66); 
m21_0 = pars(67);
m203_0 = pars(68);
HSP70_0 = pars(69);
%--------------------------------------------------------------------------
% Optimization parameters
%--------------------------------------------------------------------------
% net proliferation

n_4_tilde = pars(6); 
n_8_tilde = pars(7); 
n_1_tilde = pars(8); 
n_reg_tilde = pars(9); 
n_C_tilde = pars(10); 
n_cancer_tilde = pars(11);
n_reg_214_tilde = pars(70);
n_cancer_21_tilde = pars(71);
%---------------------------
% growth of cells

g_2_tilde = pars(12); 
g_2_4_tilde = pars(13); 
g_C_12_tilde = pars(14); 
g_NK_12_tilde = pars(15);
g_DC_cancer_tilde = pars(16); 
g_DC_th1_tilde = pars(17);
%---------------------------
% differentiation

d_1_IFN_tilde = pars(18); 
d_1_12_tilde = pars(19); 
d_2_tilde = pars(20); 
d_Treg_TGF_tilde = pars(21); 
d_C_tilde = pars(22);
%---------------------------
% cell killing

k_C_tilde = pars(23); 
k_NK_tilde = pars(24);
%---------------------------
% production

p_Treg_TGF_tilde = pars(25)/TGF_beta_0;
p_1_IFN_tilde = pars(26)/IFNg_0;
p_2_4_6_tilde = pars(27)/IL4_0;
p_DC_6_tilde = pars(28)/IL6_0;
p_DC_12_tilde = pars(29)/IL12_0;
p_NK_IFN_tilde = pars(30)/IFNg_0;
%--------------------------
% decay

delta_IFN_tilde = pars(31); 
delta_IL4_tilde = pars(32); 
delta_IL6_tilde = pars(33); 
delta_IL12_tilde = pars(34); 
delta_A_tilde = pars(35); 
delta_2_tilde = pars(36);
delta_TGF_tilde = pars(37);
d_deg_tilde = pars(72);
delta_21_tilde = pars(89);
delta_214_tilde = pars(90);
delta_203_tilde = pars(91);
delta_HSP70_tilde = pars(92);
%---------------------------
% induction/upregulation

q_1 = pars(38); 
q_IFN_tilde = pars(39)/IFNg_0; 
q_IFN_PDL1_tilde = pars(40)/IFNg_0; 
q_gIL4_tilde = pars(41)/IL4_0; 
q_dIL4_tilde = pars(42)/IL4_0; 
q_IL6_tilde = pars(43)/IL6_0; 
q_dIL12_tilde = pars(44)/IL12_0; 
q_gIL12_tilde = pars(45)/IL12_0;
q_TGF_tilde = pars(46)/TGF_beta_0;
q_214_tilde = pars(73)/m214_0; 
q_HSP70_tilde = pars(74)/HSP70_0; 
q_21_tilde = pars(75)/m21_0;
%---------------------------
% inhibition

r_Treg_tilde = pars(47);
r_IFN_tilde = pars(48)/IFNg_0; 
r_IL4_tilde = pars(49)/IL4_0; 
r_IL6_tilde = pars(50)/IL6_0; 
r_TGF_tilde = pars(51)/TGF_beta_0;
r_203_tilde = pars(76)/m203_0;
%---------------------------
% protein expression

rho = pars(52)/t_scale; 
lambda = pars(53)/t_scale; 
lambda_Can_IFN = pars(54)/t_scale;
lambda_exo_PDL1 = pars(77)/t_scale;
lambda_exo_HSP70 = pars(78)/t_scale; 
lambda_exo_m21 = pars(79)/t_scale;
lambda_exo_m214 = pars(80)/t_scale; 
lambda_exo_m203 = pars(81)/t_scale; 
%---------------------------
% Cells saturation

k_4_tilde = pars(83);
k_Treg_tilde = pars(84);
k_8_tilde = pars(85);
k_DC_tilde = pars(86);
k_NKS_tilde = pars(87);
k_can_tilde = pars(88);
%---------------------------
% complex association/disassociation

beta_plus_tilde = pars(55); 
beta_minus_tilde = pars(56); 
alpha_minus_tilde = pars(57);
alpha_plus_tilde = pars(103);
%---------------------------
% inhibition by PD-1:PD-L1

s_1 = pars(58);
s_2 = pars(59);
s_3 = pars(60);
s_C = pars(61);
%---------------------------
% cell population fractions

C_frac = pars(62); 
Th1_frac = pars(63); 
Th2_frac = pars(64);
Treg_frac = pars(65); 
%---------------------------
% Cells natural decay

dec_n4_tilde = pars(94);
dec_th1_tilde = pars(95);
dec_th2_tilde = pars(96);
dec_treg_tilde = pars(97);
dec_n8_tilde = pars(98);
dec_tc_tilde = pars(99);
dec_nk_tilde = pars(100);
dec_dc_tilde = pars(101);
dec_c_tilde = pars(102);
delta_exo_tilde = pars(104);

%--------------------------------------------------------------------------
% Initialize remaning variables
%--------------------------------------------------------------------------
% T-cell fractions at t=0

TN8_frac = 0.00194;
Tc_frac = 0.000486;
CD4_frac = 0.00292;
NK_frac = 0.00121; 
DC_frac = 1e-5;

%1.Initialize cell populations
%--------------------------------------------------------------------------
C_0 = C_frac*total_cells_flow;
T_cell_0 = (1.0-C_frac)*total_cells_flow;
TN4_0 = CD4_frac*T_cell_0;
Th1_0 = Th1_frac*T_cell_0;
Th2_0 = Th2_frac*T_cell_0;
Treg_0 = Treg_frac*T_cell_0;
TN8_0 = TN8_frac*T_cell_0;
Tc_0 = Tc_frac*T_cell_0;
NK_0 = NK_frac*T_cell_0;
DC_0 = DC_frac*T_cell_0;
%--------------------------------------------------------------------------
% 2. Calculate remaining protein parameters:
%--------------------------------------------------------------------------
p_2_4_tilde = (delta_IL4_tilde-p_2_4_6_tilde*Th2_0/(q_IL6_tilde+1.0))/(Th2_0*(r_Treg_tilde/(r_Treg_tilde+Treg_0)));
p_C_IFN_tilde = (delta_IFN_tilde - p_1_IFN_tilde*Th1_0*(r_IL4_tilde/(r_IL4_tilde+1.0))*(r_IL6_tilde/(r_IL6_tilde+1.0))*(r_Treg_tilde/(r_Treg_tilde+1))-p_NK_IFN_tilde*NK_0*((r_TGF_tilde/(r_TGF_tilde+1.0))+(1/(q_HSP70_tilde+1))))/(Tc_0*((r_TGF_tilde/(r_TGF_tilde+1.0)))); 
p_2_6_tilde = 0.055; 
p_1_12_tilde = delta_IL12_tilde/Th1_0; 
p_NK_12_tilde = delta_IL12_tilde/NK_0;

lambda_E_c = Ec_0*(d_deg_tilde*(((Treg_0/(k_Treg_tilde+Treg_0))+(DC_0/(k_DC_tilde+DC_0))+(C_0/(k_can_tilde+C_0))))+delta_exo_tilde)/(C_0);


PDL1_0 = (lambda*(Th1_0 + Th2_0 + Tc_0 + C_0) + lambda_exo_PDL1*Ec_0);

PD1_0 = rho*(Th1_0 + Th2_0 + Tc_0 + DC_0 + NK_0);
%==========================================================================
% 1. Simulations to equilibrate PD-1:PD-L1, PD-1, and PD-L1
%==========================================================================
IC_0 = [1.0, 1.0, 0];   
IC_0 = IC_0.';

DE_params_equil = [beta_plus_tilde, beta_minus_tilde];

Init = [PD1_0, PDL1_0];

% An event function is included here to ensure the system reaches an 
% equilibrium state.  
%--------------------------------------------------------------------------
tmin = 0.0;
tmax_equil = 1000000/t_scale;
num_equil_points = 50000;
tspan_equil=linspace(tmin,tmax_equil,num_equil_points);

% Pass the parameters to the event function first:
%--------------------------------------------------------------------------
ParamEventHdl = @(t_equil, x_equil)EquilibriumEventFcn(t_equil,x_equil,DE_params_equil,Init,t_scale); 
opts_equil = odeset('RelTol',1e-5,'AbsTol',1e-12,'Events',ParamEventHdl); % Add the events function here. 

% Run the equilibration simulation
%--------------------------------------------------------------------------
[t_equil,x_equil,te,ye,ie] = ode15s(@(t_equil,x_equil)...
    integratingfunction_equilibrium(t_equil,x_equil,DE_params_equil,Init),tspan_equil,IC_0,opts_equil);

% Check to make sure that the levels did actually equilibrate.
%--------------------------------------------------------------------------
if max(t_equil)==tmax_equil
   disp('ERROR: PD-1:PD-L1 levels did not equilibrate.');
   return
end 
 
%--------------------------------------------------------------------------
% Calculate new SS initial protein levels based on cell populations, then
% calculate new scaled parameters:
%--------------------------------------------------------------------------
p_DC_12 = p_DC_12_tilde*IL12_0;
p_1_12 = p_1_12_tilde*IL12_0;
p_NK_12 = p_NK_12_tilde*IL12_0;

IL12_0_FC = (p_DC_12*DC_0*r_203_tilde/(r_203_tilde+m203_0) + p_1_12*Th1_0 + p_NK_12*NK_0 )/delta_IL12_tilde;
%---------------------------------
p_2_6 = p_2_6_tilde*IL6_0;
p_DC_6 = p_DC_6_tilde*IL6_0;

IL6_0_FC = (p_2_6*Th2_0 + p_DC_6*C_0)/delta_IL6_tilde;
%---------------------------------
p_2_4 = p_2_4_tilde*IL4_0;
p_2_4_6 = p_2_4_6_tilde*IL4_0;
q_IL6 = q_IL6_tilde*IL6_0;
r_Treg = r_Treg_tilde*Treg_0;

IL4_0_FC = (p_2_4*Th2_0*(r_Treg_tilde/(r_Treg_tilde + Treg_0)) + p_2_4_6*Th2_0*IL6_0_FC/(q_IL6+IL6_0_FC))/delta_IL4_tilde;
%---------------------------------
p_Treg_TGF = p_Treg_TGF_tilde*TGF_beta_0;

TGF_beta_0_FC = (p_Treg_TGF*Treg_0)/delta_TGF_tilde;
%---------------------------------
p_1_IFN = p_1_IFN_tilde*IFNg_0;
r_IL4 = r_IL4_tilde*IL4_0;
r_IL6 = r_IL6_tilde*IL6_0;
r_TGF = r_TGF_tilde*TGF_beta_0;
p_C_IFN = p_C_IFN_tilde*IFNg_0;
p_NK_IFN = p_NK_IFN_tilde*IFNg_0;


IFNg_0_FC = ((p_1_IFN*Th1_0)*(r_IL4/(r_IL4 + IL4_0_FC))*(r_IL6/(r_IL6 + IL6_0_FC))*...
        (r_Treg/(r_Treg + Treg_0))+...
        (p_C_IFN*Tc_0)*(r_TGF/(TGF_beta_0+r_TGF))+...
        (p_NK_IFN*NK_0)*((r_TGF/(TGF_beta_0+r_TGF))+(HSP70_0/(HSP70_0+q_HSP70_tilde))))/delta_IFN_tilde;
%---------------------------------
PD1_0_FC = rho*(Th1_0 + Th2_0 + Tc_0);
rho_tilde_FC = rho/PD1_0_FC;
s_1_tilde_FC = s_1/PD1_0_FC;
s_2_tilde_FC = s_2/PD1_0_FC;
s_3_tilde_FC = s_3/PD1_0_FC;
s_C_tilde_FC = s_C/PD1_0_FC;
%---------------------------------
%HSP70

HSP70_0_FC = lambda_exo_HSP70*Ec_0;
%---------------------------------
%miRNAs

m214_0_FC = lambda_exo_m214*Ec_0*(Treg_0/Treg_0+k_Treg_tilde)/delta_214_tilde;
m203_0_FC = lambda_exo_m203*Ec_0*(DC_0/(DC_0+k_DC_tilde))/delta_203_tilde;
m21_0_FC = lambda_exo_m21*Ec_0*(C_0/(C_0+k_can_tilde))/delta_21_tilde;
lambda_exo_HSP70_tilde_FC = lambda_exo_HSP70/HSP70_0_FC;
lambda_exo_m21_tilde_FC = lambda_exo_m21/m21_0_FC;
lambda_exo_m214_tilde_FC = lambda_exo_m214/m214_0_FC;
lambda_exo_m203_tilde_FC = lambda_exo_m203/m203_0_FC;
%---------------------------------
q_IFN_PDL1 = q_IFN_PDL1_tilde*IFNg_0;

PDL1_0_FC = (lambda*(Th1_0 + Th2_0 + Tc_0 + C_0) + ...
            lambda_Can_IFN*C_0*IFNg_0_FC/(q_IFN_PDL1 + IFNg_0_FC) + lambda_exo_PDL1*Ec_0);
%---------------------------------
lambda_tilde_FC = lambda/PDL1_0_FC;
lambda_Can_IFN_tilde_FC = lambda_Can_IFN/PDL1_0_FC;
lambda_exo_PDL1_tilde_FC = lambda_exo_PDL1/PDL1_0_FC;
%------------------------------------------------------
%Now rescale the parameters appropriately

p_1_IFN_tilde_FC = p_1_IFN_tilde*IFNg_0/IFNg_0_FC;
p_1_12_tilde_FC = p_1_12_tilde*IL12_0/IL12_0_FC;
p_2_4_tilde_FC = p_2_4_tilde*IL4_0/IL4_0_FC;
p_2_6_tilde_FC = p_2_6_tilde*IL6_0/IL6_0_FC;
p_2_4_6_tilde_FC = p_2_4_6_tilde*IL4_0/IL4_0_FC; 
p_C_IFN_tilde_FC = p_C_IFN_tilde*IFNg_0/IFNg_0_FC;
p_Treg_TGF_tilde_FC = p_Treg_TGF_tilde*TGF_beta_0/TGF_beta_0_FC;  
p_NK_IFN_tilde_FC = p_NK_IFN_tilde*IFNg_0/IFNg_0_FC;
p_DC_6_tilde_FC = p_DC_6_tilde*IL6_0/IL6_0_FC;
p_DC_12_tilde_FC = p_DC_12_tilde*IL12_0/IL12_0_FC;
p_NK_12_tilde_FC = p_NK_12_tilde*IL12_0/IL12_0_FC;
q_IFN_tilde_FC = q_IFN_tilde*IFNg_0/IFNg_0_FC; 
q_IFN_PDL1_tilde_FC = q_IFN_PDL1_tilde*IFNg_0/IFNg_0_FC; 
q_gIL4_tilde_FC = q_gIL4_tilde*IL4_0/IL4_0_FC; 
q_dIL4_tilde_FC = q_dIL4_tilde*IL4_0/IL4_0_FC; 
q_IL6_tilde_FC = q_IL6_tilde*IL6_0/IL6_0_FC; 
q_dIL12_tilde_FC = q_dIL12_tilde*IL12_0/IL12_0_FC; 
q_gIL12_tilde_FC = q_gIL12_tilde*IL12_0/IL12_0_FC;
q_TGF_tilde_FC = q_TGF_tilde*TGF_beta_0/TGF_beta_0_FC;
q_214_tilde_FC = q_214_tilde*m214_0/m214_0_FC; 
q_HSP70_tilde_FC = q_HSP70_tilde*HSP70_0/HSP70_0_FC; 
q_21_tilde_FC = q_21_tilde*m21_0/m21_0_FC;
r_IFN_tilde_FC = r_IFN_tilde*IFNg_0/IFNg_0_FC; 
r_IL4_tilde_FC = r_IL4_tilde*IL4_0/IL4_0_FC; 
r_IL6_tilde_FC = r_IL6_tilde*IL6_0/IL6_0_FC;
r_TGF_tilde_FC = r_TGF_tilde*TGF_beta_0/TGF_beta_0_FC;
r_203_tilde_FC = r_203_tilde*m203_0/m203_0_FC;
%==========================================================================
% 1. Simulations to equilibrate PD-1:PD-L1, PD-1, and PD-L1
%==========================================================================
IC_0 = [1.0, 1.0, 0];   
IC_0 = IC_0.';

DE_params_equil = [beta_plus_tilde, beta_minus_tilde];

Init = [PD1_0_FC, PDL1_0_FC];

% An event function is included here to ensure the system reaches an 
% equilibrium state.  
%--------------------------------------------------------------------------
tmin = 0.0;
tmax_equil = 1000000/t_scale;
num_equil_points = 50000;
tspan_equil=linspace(tmin,tmax_equil,num_equil_points);

% Pass the parameters to the event function first:
%--------------------------------------------------------------------------
ParamEventHdl = @(t_equil, x_equil)EquilibriumEventFcn(t_equil,x_equil,DE_params_equil,Init,t_scale); 
opts_equil = odeset('RelTol',1e-5,'AbsTol',1e-12,'Events',ParamEventHdl); % Add the events function here. 

% Run the equilibration simulation: (Call ode15s solver while passing the 
% parameters and ODE options to the integrating function)
%--------------------------------------------------------------------------
[t_equil,x_equil,te,ye,ie] = ode15s(@(t_equil,x_equil)...
    integratingfunction_equilibrium(t_equil,x_equil,DE_params_equil,Init),tspan_equil,IC_0,opts_equil);


%==========================================================================
% Parameter array to pass to integrating function:


%--------------------------------------------------------------------------
DE_params = [n_4_tilde, n_8_tilde, n_1_tilde, n_reg_tilde,... 
            n_cancer_tilde,  n_C_tilde,...
            g_2_tilde, g_2_4_tilde, g_C_12_tilde,... 
            g_NK_12_tilde, g_DC_cancer_tilde, g_DC_th1_tilde,... 
            d_1_12_tilde, d_1_IFN_tilde, d_2_tilde, d_Treg_TGF_tilde, d_C_tilde,...
            k_C_tilde, k_NK_tilde,...
            s_1_tilde_FC, s_2_tilde_FC, s_3_tilde_FC, s_C_tilde_FC,...
            q_dIL12_tilde_FC, q_gIL12_tilde_FC, q_IFN_tilde_FC, q_dIL4_tilde_FC,...
            q_gIL4_tilde_FC, q_TGF_tilde_FC,...
            q_1, q_IL6_tilde_FC, q_IFN_PDL1_tilde_FC,...
            r_IFN_tilde_FC, r_IL4_tilde_FC, r_IL6_tilde_FC,...
            r_TGF_tilde_FC, r_Treg_tilde,... 
            p_Treg_TGF_tilde_FC, p_1_IFN_tilde_FC, p_C_IFN_tilde_FC, p_NK_IFN_tilde_FC,...
            p_2_4_tilde_FC, p_2_4_6_tilde_FC, p_2_6_tilde_FC, p_DC_6_tilde_FC,...
            p_1_12_tilde_FC, p_DC_12_tilde_FC, p_NK_12_tilde_FC,...
            delta_TGF_tilde, delta_IFN_tilde, delta_IL6_tilde, delta_IL4_tilde, ...
            delta_IL12_tilde, delta_A_tilde, delta_2_tilde, 0,...
            rho_tilde_FC, lambda_tilde_FC, lambda_Can_IFN_tilde_FC,...
            beta_plus_tilde, beta_minus_tilde,...  
            alpha_plus_tilde, alpha_minus_tilde, n_reg_214_tilde, n_cancer_21_tilde, q_214_tilde_FC, q_HSP70_tilde_FC, q_21_tilde_FC,...
            r_203_tilde_FC, d_deg_tilde, lambda_exo_PDL1_tilde_FC, lambda_exo_HSP70_tilde_FC,...
            lambda_exo_m21_tilde_FC, lambda_exo_m214_tilde_FC, lambda_exo_m203_tilde_FC, lambda_E_c,...
            k_4_tilde, k_Treg_tilde, k_8_tilde, k_DC_tilde, k_NKS_tilde, k_can_tilde,...
            delta_21_tilde, delta_214_tilde, delta_203_tilde, delta_HSP70_tilde, dec_n4_tilde, dec_th1_tilde, dec_th2_tilde, dec_treg_tilde, dec_n8_tilde, dec_tc_tilde, dec_nk_tilde, dec_dc_tilde, dec_c_tilde, delta_exo_tilde];
             
Init = [PD1_0_FC, PDL1_0_FC, nivo_0];

tmin = 0.0;
tmax = 24*60/t_scale;
num_points = 4000;
tspan=linspace(tmin,tmax,num_points);

%--------------------------------------------------------------------------         
% i) First 24 hours:         
%--------------------------------------------------------------------------
% Initial conditions array to pass to integrating function:
%--------------------------------------------------------------------------
IC = [TN4_0, Th1_0, Th2_0, Treg_0, TN8_0, Tc_0, NK_0, DC_0, C_0,...
      1.0, 1.0, 1.0, 1.0, ...
      x_equil(end,1), x_equil(end,2), x_equil(end,3),...
      1, 1, 0,...
      Ec_0, 1, 1, 1, 1]; 
  
IC = IC.';

% Run the treatment simulations
%--------------------------------------------------------------------------
opts = odeset('RelTol',1e-5,'AbsTol',1e-12); 

[t_treatment_1d,x_treatment_1d] = ode15s(@(t_treatment_1d,x_treatment_1d)...
       integratingfunction_flow_cytometry_treatment(t_treatment_1d,x_treatment_1d,DE_params,Init),tspan,IC,opts);

%--------------------------------------------------------------------------         
% ii) Second dosage:         
%--------------------------------------------------------------------------
% Initial conditions array to pass to integrating function:
%--------------------------------------------------------------------------
IC = x_treatment_1d(end,:);
IC(18) = IC(18) + 1.0;   %Drug remains from the previous day
IC = IC.';

% Run the treatment simulations
%--------------------------------------------------------------------------
opts = odeset('RelTol',1e-7,'AbsTol',1e-12); 

[t_treatment_2d,x_treatment_2d] = ode15s(@(t_treatment_2d,x_treatment_2d)...
       integratingfunction_flow_cytometry_treatment(t_treatment_2d,x_treatment_2d,DE_params,Init),tspan,IC,opts);
   
t_treatment_2d = t_treatment_2d + 24*60/t_scale;
%--------------------------------------------------------------------------         
% iii) Third dosage:         
%--------------------------------------------------------------------------
% Initial conditions array to pass to integrating function:
%--------------------------------------------------------------------------
IC = x_treatment_2d(end,:);
IC(18) = IC(18) + 1.0;    %Drug remains from the previous day
IC = IC.';

% Run the treatment simulations
%--------------------------------------------------------------------------
opts = odeset('RelTol',1e-7,'AbsTol',1e-12); 

[t_treatment_3d,x_treatment_3d] = ode15s(@(t_treatment_3d,x_treatment_3d)...
       integratingfunction_flow_cytometry_treatment(t_treatment_3d,x_treatment_3d,DE_params,Init),tspan,IC,opts);
   
t_treatment_3d = t_treatment_2d + 24*60/t_scale;
%--------------------------------------------------------------------------         
% iii) fourth dosage:         
%--------------------------------------------------------------------------
% Initial conditions array to pass to integrating function:
%--------------------------------------------------------------------------
IC = x_treatment_3d(end,:);
IC(18) = IC(18) + 1.0;    %Drug remains from the previous day
IC = IC.';

% Run the treatment simulations
%--------------------------------------------------------------------------
opts = odeset('RelTol',1e-7,'AbsTol',1e-12); 

[t_treatment_4d,x_treatment_4d] = ode15s(@(t_treatment_4d,x_treatment_4d)...
       integratingfunction_flow_cytometry_treatment(t_treatment_4d,x_treatment_4d,DE_params,Init),tspan,IC,opts);
   
t_treatment_4d = t_treatment_3d + 24*60/t_scale;
%--------------------------------------------------------------------------         
% iii) fifth dosage:         
%--------------------------------------------------------------------------
% Initial conditions array to pass to integrating function:
%--------------------------------------------------------------------------
IC = x_treatment_4d(end,:);
IC(18) = IC(18) + 1.0;    %Drug remains from the previous day
IC = IC.';

% Run the treatment simulations
%--------------------------------------------------------------------------
opts = odeset('RelTol',1e-7,'AbsTol',1e-12); 

[t_treatment_5d,x_treatment_5d] = ode15s(@(t_treatment_5d,x_treatment_5d)...
       integratingfunction_flow_cytometry_treatment(t_treatment_5d,x_treatment_5d,DE_params,Init),tspan,IC,opts);
   
t_treatment_5d = t_treatment_4d + 24*60/t_scale;
%--------------------------------------------------------------------------         
% iii) sixth dosage:         
%--------------------------------------------------------------------------
% Initial conditions array to pass to integrating function:
%--------------------------------------------------------------------------
IC = x_treatment_5d(end,:);
IC(18) = IC(18) + 1.0;   %Drug remains from the previous day
IC = IC.';

% Run the treatment simulations
%--------------------------------------------------------------------------
opts = odeset('RelTol',1e-7,'AbsTol',1e-12); 

[t_treatment_6d,x_treatment_6d] = ode15s(@(t_treatment_6d,x_treatment_6d)...
       integratingfunction_flow_cytometry_treatment(t_treatment_6d,x_treatment_6d,DE_params,Init),tspan,IC,opts);
   
t_treatment_6d = t_treatment_5d + 24*60/t_scale;
%--------------------------------------------------------------------------         
% iii) seventh dosage:         
%--------------------------------------------------------------------------
% Initial conditions array to pass to integrating function:
%--------------------------------------------------------------------------
IC = x_treatment_6d(end,:);
IC(18) = IC(18) + 1.0;    %Drug remains from the previous day
IC = IC.';

% Run the treatment simulations
%--------------------------------------------------------------------------
opts = odeset('RelTol',1e-7,'AbsTol',1e-12); 

[t_treatment_7d,x_treatment_7d] = ode15s(@(t_treatment_7d,x_treatment_7d)...
       integratingfunction_flow_cytometry_treatment(t_treatment_7d,x_treatment_7d,DE_params,Init),tspan,IC,opts);
   
t_treatment_7d = t_treatment_6d + 24*60/t_scale;
%--------------------------------------------------------------------------         
% iii) eigths dosage:         
%--------------------------------------------------------------------------
% Initial conditions array to pass to integrating function:
%--------------------------------------------------------------------------
IC = x_treatment_7d(end,:);
IC(18) = IC(18) + 1.0;    %Drug remains from the previous day
IC = IC.';

% Run the treatment simulations
%--------------------------------------------------------------------------
opts = odeset('RelTol',1e-7,'AbsTol',1e-12); 

[t_treatment_8d,x_treatment_8d] = ode15s(@(t_treatment_8d,x_treatment_8d)...
       integratingfunction_flow_cytometry_treatment(t_treatment_8d,x_treatment_8d,DE_params,Init),tspan,IC,opts);
   
t_treatment_8d = t_treatment_7d + 24*60/t_scale;
%--------------------------------------------------------------------------         
% iii) ninth dosage:         
%--------------------------------------------------------------------------
% Initial conditions array to pass to integrating function:
%--------------------------------------------------------------------------
IC = x_treatment_8d(end,:);
IC(18) = IC(18) + 1.0;    %Drug remains from the previous day
IC = IC.';

% Run the treatment simulations
%--------------------------------------------------------------------------
opts = odeset('RelTol',1e-7,'AbsTol',1e-12); 

[t_treatment_9d,x_treatment_9d] = ode15s(@(t_treatment_9d,x_treatment_9d)...
       integratingfunction_flow_cytometry_treatment(t_treatment_9d,x_treatment_9d,DE_params,Init),tspan,IC,opts);
   
t_treatment_9d = t_treatment_8d + 24*60/t_scale;
%--------------------------------------------------------------------------         
% iii) tenth dosage:         
%--------------------------------------------------------------------------
% Initial conditions array to pass to integrating function:
%--------------------------------------------------------------------------
IC = x_treatment_9d(end,:);
IC(18) = IC(18) + 1.0;    %Drug remains from the previous day
IC = IC.';

% Run the treatment simulations
%--------------------------------------------------------------------------
opts = odeset('RelTol',1e-7,'AbsTol',1e-12); 

[t_treatment_10d,x_treatment_10d] = ode15s(@(t_treatment_10d,x_treatment_10d)...
       integratingfunction_flow_cytometry_treatment(t_treatment_10d,x_treatment_10d,DE_params,Init),tspan,IC,opts);
   
t_treatment_10d = t_treatment_9d + 24*60/t_scale;
%--------------------------------------------------------------------------         
% iii) eleventh dosage:         
%--------------------------------------------------------------------------
% Initial conditions array to pass to integrating function:
%--------------------------------------------------------------------------
IC = x_treatment_10d(end,:);
IC(18) = IC(18) + 1.0;    %Drug remains from the previous day
IC = IC.';

% Run the treatment simulations
%--------------------------------------------------------------------------
opts = odeset('RelTol',1e-7,'AbsTol',1e-12); 

[t_treatment_11d,x_treatment_11d] = ode15s(@(t_treatment_11d,x_treatment_11d)...
       integratingfunction_flow_cytometry_treatment(t_treatment_11d,x_treatment_11d,DE_params,Init),tspan,IC,opts);
   
t_treatment_11d = t_treatment_10d + 24*60/t_scale;
%--------------------------------------------------------------------------         
% iii) twelwth dosage:         
%--------------------------------------------------------------------------
% Initial conditions array to pass to integrating function:
%--------------------------------------------------------------------------
IC = x_treatment_11d(end,:);
IC(18) = IC(18) + 1.0;    %Drug remains from the previous day
IC = IC.';

% Run the treatment simulations
%--------------------------------------------------------------------------
opts = odeset('RelTol',1e-7,'AbsTol',1e-12); 

[t_treatment_12d,x_treatment_12d] = ode15s(@(t_treatment_12d,x_treatment_12d)...
       integratingfunction_flow_cytometry_treatment(t_treatment_12d,x_treatment_12d,DE_params,Init),tspan,IC,opts);
   
t_treatment_12d = t_treatment_11d + 24*60/t_scale;
%--------------------------------------------------------------------------         
% iii) thirtheenth dosage:         
%--------------------------------------------------------------------------
% Initial conditions array to pass to integrating function:
%--------------------------------------------------------------------------
IC = x_treatment_12d(end,:);
IC(18) = IC(18) + 1.0;    %Drug remains from the previous day
IC = IC.';

% Run the treatment simulations
%--------------------------------------------------------------------------
opts = odeset('RelTol',1e-7,'AbsTol',1e-12); 

[t_treatment_13d,x_treatment_13d] = ode15s(@(t_treatment_13d,x_treatment_13d)...
       integratingfunction_flow_cytometry_treatment(t_treatment_13d,x_treatment_13d,DE_params,Init),tspan,IC,opts);
   
t_treatment_13d = t_treatment_12d + 24*60/t_scale;
%--------------------------------------------------------------------------         
% iii) fourteenth dosage:         
%--------------------------------------------------------------------------
% Initial conditions array to pass to integrating function:
%--------------------------------------------------------------------------
IC = x_treatment_13d(end,:);
IC(18) = IC(18) + 1.0;    %Drug remains from the previous day
IC = IC.';

% Run the treatment simulations
%--------------------------------------------------------------------------
opts = odeset('RelTol',1e-7,'AbsTol',1e-12); 

[t_treatment_14d,x_treatment_14d] = ode15s(@(t_treatment_14d,x_treatment_14d)...
       integratingfunction_flow_cytometry_treatment(t_treatment_14d,x_treatment_14d,DE_params,Init),tspan,IC,opts);
   
t_treatment_14d = t_treatment_13d + 24*60/t_scale;

%-----------DOSAGE STOPS FOR NEXT 2 WEEKS-----------

tmin = 14*24*60/t_scale;
tmax = 28*24*60/t_scale;
num_points = 56000;
tspan=linspace(tmin,tmax,num_points);


IC = x_treatment_14d(end,:);
%IC(18) = IC(18) + 1.0;    %No Drug Addition
IC = IC.';

% Run the treatment simulations
%--------------------------------------------------------------------------
opts = odeset('RelTol',1e-3,'AbsTol',1e-12); 

[t_treatment_15d,x_treatment_15d] = ode15s(@(t_treatment_15d,x_treatment_15d)...
       integratingfunction_flow_cytometry_treatment2(t_treatment_15d,x_treatment_15d,DE_params,Init),tspan,IC,opts);


% Append all the treatment arrays
%--------------------------------------------------------------------------

t_treatment_FC1 = [t_treatment_1d; t_treatment_2d; t_treatment_3d; t_treatment_4d; t_treatment_5d; t_treatment_6d; t_treatment_7d; t_treatment_8d; t_treatment_9d; t_treatment_10d; t_treatment_11d; t_treatment_12d; t_treatment_13d; t_treatment_14d; t_treatment_15d];
x_treatment_FC1 = [x_treatment_1d; x_treatment_2d; x_treatment_3d; x_treatment_4d; x_treatment_5d; x_treatment_6d; x_treatment_7d; x_treatment_8d; x_treatment_9d; x_treatment_10d; x_treatment_11d; x_treatment_12d; x_treatment_13d; x_treatment_14d; x_treatment_15d];

t_treatment_FC2 = linspace(0,tmax/2,112000);

semilogy(t_treatment_FC2, x_treatment_FC1(:,9), 'color','black')

xlabel('Time (days)')
ylabel('Cancer cells number')

L = legend('show');
title(L,'Cancer cells population')


grid off
hold on


end