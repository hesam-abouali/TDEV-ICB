function dxdt = integratingfunction_flow_cytometry_treatment2(t,x,params,Init)

% The integrating function for simulations without treatment and drugs disabled
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Variable definitions:
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Cells:

%TN4=x(1); 
%Th1=x(2); 
%Th2=x(3); 
%Treg=x(4); 
%TN8=x(5); 
%Tc=x(6);
%NK=x(7); 
%DC=x(8); 
%C=x(9);
%--------------------------------------------------------------------------
%Proteins:

%IFN_gamma=x(10); 
%IL-4=x(11); 
%IL-6=x(12); 
%IL-12=x(13); 
%PD-1=x(14); 
%PD-L1=x(15); 
%PD-1:PD-L1=x(16); 
%TGF_beta=x(17);
%--------------------------------------------------------------------------
%Drug-related:

%A=x(18); 
%A:PD-1=x(19);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Exosomal Content

%Ec=x(20); 
%m21=x(21); 
%m203=x(22); 
%m214=x(23); 
%HSP70=x(24);
% Parameters
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Read in all parameters
%--------------------------------------------------------------------------
% net proliferation

n_4_tilde = params(1); 
n_8_tilde = params(2); 
n_1_tilde = params(3);
n_reg_tilde = params(4);  
n_cancer_tilde = params(5);
n_C_tilde = params(6);
n_reg_214_tilde = params(64);
n_cancer_21_tilde = params(65);
%--------------------------------------------------------------------------
% growth of cells

g_2_tilde = params(7); 
g_2_4_tilde = params(8); 
g_C_12_tilde = params(9); 
g_NK_12_tilde = params(10); 
g_DC_cancer_tilde = params(11); 
g_DC_th1_tilde = params(12); 
%--------------------------------------------------------------------------
% differentiation

d_1_12_tilde = params(13); 
d_1_IFN_tilde = params(14); 
d_2_tilde = params(15); 
d_Treg_TGF_tilde = params(16); 
d_C_tilde = params(17);
%--------------------------------------------------------------------------
% cell killing

k_C_tilde = params(18); 
k_NK_tilde = params(19);
%--------------------------------------------------------------------------
% inhibition by PD-1:PD-L1

s_1_tilde = params(20); 
s_2_tilde = params(21); 
s_3_tilde = params(22); 
s_C_tilde = params(23);
%--------------------------------------------------------------------------
% induction/upregulation

q_dIL12_tilde = params(24); 
q_gIL12_tilde = params(25);
q_IFN_tilde = params(26); 
q_dIL4_tilde = params(27);
q_gIL4_tilde = params(28); 
q_TGF_tilde = params(29);
q_1 = params(30);
q_IL6_tilde = params(31); 
q_IFN_PDL1_tilde = params(32);
q_214_tilde = params(66); 
q_HSP70_tilde = params(67);
q_21_tilde = params(68); 
%--------------------------------------------------------------------------
% inhibition

r_IFN_tilde = params(33); 
r_IL4_tilde = params(34); 
r_IL6_tilde = params(35); 
r_TGF_tilde = params(36);
r_Treg_tilde = params(37);
r_203_tilde = params(69); 
%--------------------------------------------------------------------------
% production

p_Treg_TGF_tilde = params(38); 
p_1_IFN_tilde = params(39); 
p_C_IFN_tilde = params(40); 
p_NK_IFN_tilde = params(41); 
p_2_4_tilde = params(42); 
p_2_4_6_tilde = params(43);
p_2_6_tilde = params(44); 
p_DC_6_tilde = params(45); 
p_1_12_tilde = params(46); 
p_DC_12_tilde = params(47);
p_NK_12_tilde = params(48); 
%--------------------------------------------------------------------------
% decay

delta_TGF_tilde = params(49); 
delta_IFN_tilde = params(50); 
delta_IL6_tilde = params(51); 
delta_IL4_tilde = params(52); 
delta_IL12_tilde = params(53); 
delta_A_tilde = params(54); 
delta_2_tilde = params(55);
d_deg_tilde = params(70);
delta_21_tilde = params(83);
delta_214_tilde = params(84);
delta_203_tilde = params(85);
delta_HSP70_tilde = params(86);
%--------------------------------------------------------------------------
% protein expression

rho_tilde = params(57); 
lambda_tilde = params(58); 
lambda_Can_IFN_tilde = params(59);
lambda_exo_PDL1_tilde = params(71);
lambda_exo_HSP70_tilde = params(72); 
lambda_exo_m21_tilde = params(73);
lambda_exo_m214_tilde = params(74); 
lambda_exo_m203_tilde = params(75); 
lambda_E_c = params(76);
%--------------------------------------------------------------------------
% complex association/disassociation

beta_plus_tilde = params(60); 
beta_minus_tilde = params(61); 
alpha_plus_tilde = params(62); 
alpha_minus_tilde = params(63);
%--------------------------------------------------------------------------
% Cells saturation

k_4_tilde = params(77);
k_Treg_tilde = params(78);
k_8_tilde = params(79);
k_DC_tilde = params(80);
k_NKS_tilde = params(81);
k_can_tilde = params(82);
%--------------------------------------------------------------------------
% Cells natural decay

dec_n4_tilde = params(87);
dec_th1_tilde = params(88);
dec_treg_tilde = params(90);
dec_n8_tilde = params(91);
dec_tc_tilde = params(92);
dec_nk_tilde = params(93);
dec_dc_tilde = params(94);
dec_c_tilde = params(95);
delta_exo = params(96);
%--------------------------------------------------------------------------
PD1_0 = Init(1); 
PDL1_0 = Init(2); 
A_0 = Init(3);
%--------------------------------------------------------------------------
dxdt = zeros(max(size(x)),1);
%--------------------------------------------------------------------------
% Equations defining output:  
%--------------------------------------------------------------------------   
%(1) T_N4: 
dxdt(1) = n_4_tilde*x(1)...
        - (d_1_12_tilde*x(1)*(x(13)/(q_dIL12_tilde + x(13)))*(r_TGF_tilde/r_TGF_tilde+x(17))...
        +  d_1_IFN_tilde*x(1)*(x(10)/(q_IFN_tilde + x(10))))*(s_1_tilde/(s_1_tilde+x(16)))...
        - (d_2_tilde*x(1)*(x(11)/(q_dIL4_tilde + x(11))))*(s_2_tilde/(s_2_tilde+x(16)))...
        - (d_Treg_TGF_tilde*x(1)*(x(17)/(q_TGF_tilde + x(17))) + ...
        d_Treg_TGF_tilde*x(1)*(x(17)/(q_TGF_tilde + x(17)))*(x(15)/x(15)+s_3_tilde))*(r_IL6_tilde/(x(12)+r_IL6_tilde)) - dec_n4_tilde*x(1);
            
            
%(2) Th1: 
dxdt(2) = n_1_tilde*x(2) + (d_1_12_tilde*x(1)*(x(13)/(q_dIL12_tilde + x(13)))*(r_TGF_tilde/(r_TGF_tilde+x(17))) + ...
          d_1_IFN_tilde*x(1)*(x(10)/(q_IFN_tilde + x(10))))*(s_1_tilde/(s_1_tilde+x(16))) - dec_th1_tilde*x(2);

      
%(3) Th2: 
dxdt(3) = (g_2_tilde*x(3) + g_2_4_tilde*x(3) * (x(11)/(q_gIL4_tilde + x(11))))* ...
          (r_IFN_tilde/(r_IFN_tilde + x(10))) + ...
          (d_2_tilde*x(1)*(x(11)/(q_dIL4_tilde + x(11))))*(s_2_tilde/(s_2_tilde+x(16))) - ...
           delta_2_tilde*x(3);
              
%(4) Treg:
dxdt(4) = n_reg_tilde*x(4) + (n_reg_214_tilde*x(4)*(x(23)/(x(23)+q_214_tilde))) + ...
        (d_Treg_TGF_tilde*x(1)*(x(17)/(q_TGF_tilde + x(17))) + ...
        d_Treg_TGF_tilde*x(1)*(x(17)/(q_TGF_tilde + x(17)))*(x(15)/(x(15)+s_3_tilde)))*(r_IL6_tilde/(x(12)+r_IL6_tilde))-dec_treg_tilde*x(4);
   
%(5) T_N8: 
dxdt(5) = n_8_tilde*x(5) - d_C_tilde*x(5)*(x(2)/(q_1+x(2)))*(s_C_tilde/(s_C_tilde + x(16))) - dec_n8_tilde*x(5);


%(6) Tc: 
dxdt(6) = n_C_tilde*x(6) + g_C_12_tilde*x(6)*(x(13)/(q_gIL12_tilde + x(13)))...
        + d_C_tilde*x(5)*(x(2)/(q_1+x(2)))*(s_C_tilde/(s_C_tilde + x(16))) - dec_tc_tilde*x(6);

%(7) NK:
dxdt(7) = g_NK_12_tilde*x(7)*(x(13)/(q_gIL12_tilde+x(13))) - dec_nk_tilde*x(7);

%(8) DC:
dxdt(8) = g_DC_cancer_tilde*x(9)*(r_203_tilde/(x(22)+r_203_tilde)) + ...
        g_DC_th1_tilde*x(2) - dec_dc_tilde*x(8);

%(9) C:
dxdt(9) = n_cancer_tilde*x(9) + n_cancer_21_tilde*x(9)*(x(21)/(x(21)+q_21_tilde)) - ...
        k_C_tilde*x(9)*x(6) - k_NK_tilde*x(9)*x(7) - dec_c_tilde*x(9);
    

% Normalized to IFNg_0    
%(10) IFN_gamma: 

dxdt(10) = p_1_IFN_tilde*x(2)*(r_IL4_tilde/(r_IL4_tilde + x(11)))*(r_IL6_tilde/(r_IL6_tilde + x(12)))*...
        (r_Treg_tilde/(r_Treg_tilde + x(4))) +  p_C_IFN_tilde*x(6)*((r_TGF_tilde/(x(17)+r_TGF_tilde))) + ... 
        p_NK_IFN_tilde*x(7)*((r_TGF_tilde/(x(17)+r_TGF_tilde))+(x(24)/(x(24)+q_HSP70_tilde))) - delta_IFN_tilde*x(10);    
            
% Normalized to IL4_0
%(11) IL-4: 
dxdt(11) = p_2_4_tilde*x(3)*(r_Treg_tilde/(r_Treg_tilde + x(4))) + p_2_4_6_tilde*x(3)*(x(12)/(q_IL6_tilde + x(12))) - ...
          delta_IL4_tilde*x(11);
            
% Normalized to IL6_0
%(12) IL-6: 
dxdt(12) = p_2_6_tilde*x(3) + p_DC_6_tilde*x(9) - delta_IL6_tilde*x(12);


% Normalized to IL12_0            
%(13) IL-12: 
dxdt(13) = p_DC_12_tilde*x(8)*(r_203_tilde/(r_203_tilde+x(22))) + p_NK_12_tilde*x(7) + p_1_12_tilde*x(2) - ... 
delta_IL12_tilde*x(13);
    
% Normalized to PD1_0        
%(14) PD-1: 
dxdt(14) =  rho_tilde*(dxdt(2) + dxdt(3) + dxdt(6) + dxdt(7) + dxdt(8))...  
            - beta_plus_tilde*PDL1_0*x(14)*x(15) + beta_minus_tilde*x(16)...
            - alpha_plus_tilde*A_0*x(14)*x(18) + alpha_minus_tilde*A_0*x(19)/PD1_0;
        
% Normalized to PDL1_0   
%(15) PD-L1: 

dxdt(15) =  lambda_tilde*(dxdt(2) + dxdt(3) + dxdt(6) + dxdt(9))...
            + lambda_Can_IFN_tilde*dxdt(9)*(x(10)/(x(10)+q_IFN_PDL1_tilde)) + lambda_exo_PDL1_tilde*dxdt(20) - beta_plus_tilde*PD1_0*x(14)*x(15) + beta_minus_tilde*PD1_0*x(16)/PDL1_0;

        
        
% Normalized to PD1_0 
%(16) PD-1:PD-L1: 
dxdt(16) =  beta_plus_tilde*PDL1_0*x(14)*x(15) - beta_minus_tilde*x(16);

% Normalized to TGF_beta 
%(17) TGF
dxdt(17) = p_Treg_TGF_tilde*x(4) - delta_TGF_tilde*x(17);


% Normalized to A_0
%(18) A (free drug): 
dxdt(18) = 0;            

% Normalized to A_0
%(19) A:PD-1: 
dxdt(19) = 0;

% TDEV
%(20)

dxdt(20)= lambda_E_c*x(9)-d_deg_tilde*x(20)*((x(4)/(k_Treg_tilde+x(4)))+(x(8)/(k_DC_tilde+x(8)))+(x(9)/(k_can_tilde+x(9))))-delta_exo*x(20);

    
% Normalized to miRNA-21
%(21)
dxdt(21) = lambda_exo_m21_tilde*x(20)*(x(9)/(x(9)+k_can_tilde)) - delta_21_tilde*x(21);

% Normalized to miRNA-203
%(22)
dxdt(22) = lambda_exo_m203_tilde*x(20)*(x(8)/(x(8)+k_DC_tilde)) - delta_203_tilde*x(22);

% Normalized to miRNA-214
%(23)
dxdt(23) = lambda_exo_m214_tilde*x(20)*(x(4)/(x(4)+k_Treg_tilde)) - delta_214_tilde*x(23);

% Normalized to HSP70
%(24)
dxdt(24) = lambda_exo_HSP70_tilde*x(20) - delta_HSP70_tilde*x(24);


end