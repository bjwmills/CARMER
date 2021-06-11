function dy = CARMER_equations(t,y)
%%%%%% CAR-MER
%%%%%% Coupled C-Hg 4-box ocean-atmosphere model
%%%%%% As used in Dal Corso et al. 2020 Nat. Comm.
%%%%%% Coded by Benjamin JW Mills // b.mills@leeds.ac.uk
%%%%%% Model equations file: do not run this code directly

%%%%%%% setup dy array
dy = zeros(19,1);  

%%%%%%% set up global parameters
global stepnumber
global pars
global workingstate

%%%%%%%%%%%%% get variables from Y 
CO2_a = y(1) ;
DIC_s = y(2) ;
DIC_h = y(3) ;
DIC_d = y(4) ;
ALK_s = y(5) ;
ALK_h = y(6) ;
ALK_d = y(7) ;

d13c_atm = y(8) / y(1) ;
d13c_DIC_s = y(9) / y(2) ;
d13c_DIC_h = y(10) / y(3) ;
d13c_DIC_d = y(11) / y(4) ;

Hg_a = y(12) ;
Hg_s = y(13) ;
Hg_h = y(14) ;
Hg_d = y(15) ;

Atmospheric_CO2_ppm = ( CO2_a / pars.CO2_a_0 ) * 280 ;

%%%% pCO2 in PAL
pCO2_a = ( CO2_a / pars.CO2_a_0 ) ;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%   SCENARIO 1: VOLC ONLY %%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%% SIBERIAN TRAPS INTRUSIVE MAGMATISM
% 
% %%% LOW INPUT SCENARIO
% CO2_ramp = 4e12 ;
% Hg_ramp = 25 ;
% 
% %%% HIGH INPUT SCENARIO
% CO2_ramp = 8e12 ;
% Hg_ramp = 43 ;
% 
% %%% TH slowdown
% TH_ramp = 1 ; 
% 
% %%% interpolaiton functions
% CO2_input_force = interp1([-253e6 -251.99e6 -251.98e6 -251.56e6 -251.55e6 -251e6 ],[ 0 0 CO2_ramp CO2_ramp 0 0 ],t) ;
% Hg_volc_force = interp1([-253e6 -251.99e6 -251.98e6 -251.56e6 -251.55e6 -251e6 ],[ 1 1 Hg_ramp Hg_ramp 1 1 ],t) ;
% thermo_force = interp1([-253e6 -251.99e6 -251.98e6 -251.56e6 -251.55e6 -251e6 ],[ 10 10 TH_ramp TH_ramp 10 10 ],t) ;
% 
% %%% LAND BIOTA OPTIONS
% 
% %%% carbon burial from land biota
% C_burial_force = 1.2 ;
% OXIDW_FORCE = 1 ;
% Hg_runoff_force = 1 ;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% SCENARIO 2: VOLC+BIOSPHERE %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% SIBERIAN TRAPS INTRUSIVE MAGMATISM

%%%% LOW INPUT SCENARIO
% CO2_ramp = 4e12 ;
% Hg_ramp = 25 ;

%%%% HIGH INPUT SCENARIO
CO2_ramp = 8e12 ;
Hg_ramp = 43 ;

%%%% TH slowdown
TH_ramp = 1 ; 

%%%% interpolaiton functions
CO2_input_force = interp1([-253e6 -251.99e6 -251.98e6 -251.56e6 -251.55e6 -251e6 ],[ 0 0 CO2_ramp CO2_ramp 0 0 ],t) ;
Hg_volc_force = interp1([-253e6 -251.99e6 -251.98e6 -251.56e6 -251.55e6 -251e6 ],[ 1 1 Hg_ramp Hg_ramp 1 1 ],t) ;
thermo_force = interp1([-253e6 -251.99e6 -251.98e6 -251.56e6 -251.55e6 -251e6 ],[ 10 10 TH_ramp TH_ramp 10 10 ],t) ;

%%%% LAND BIOTA COLLAPSE (1kyr)
C_burial_force = interp1([-253e6 -251.951e6 -251.950e6 -251.949e6 -251.948e6 -251e6 ],[ 1.2 1.2 0 0 1.2 1.2],t) ;
OXIDW_FORCE = interp1([-253e6 -251.951e6 -251.950e6 -251.949e6 -251.948e6 -251e6 ],[ 1 1 30 30 1 1],t) ;
Hg_runoff_force = interp1([-253e6 -251.951e6 -251.950e6 -251.949e6 -251.948e6 -251e6  ],[ 1 1 100 100 1 1],t) ;

CO2_input_force_spike = 0 ;




% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%  SI PLOT: VOLC + VOLC SPIKE %%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % %%% SIBERIAN TRAPS INTRUSIVE MAGMATISM
% % 
% %%% LOW INPUT SCENARIO
% % CO2_ramp = 4e12 ;
% % CO2_ramp_pulse = CO2_ramp * 5 ;
% % Hg_ramp = 25 ;
% % Hg_ramp_pulse = Hg_ramp * 5 ;
% 
% % %%%% HIGH INPUT SCENARIO
% CO2_ramp = 8e12 ;
% CO2_ramp_pulse = CO2_ramp * 5 ;
% Hg_ramp = 43 ;
% Hg_ramp_pulse = Hg_ramp * 5 ;
% 
% %%%% TH slowdown
% TH_ramp = 1 ; 
% 
% %%%% interpolaiton functions
% CO2_input_force = interp1([-253e6 -251.99e6 -251.98e6 -251.56e6 -251.55e6 -251e6 ],[ 0 0 CO2_ramp CO2_ramp 0 0 ],t) ; 
% 
% CO2_input_force_spike = interp1([-253e6 -251.951e6 -251.950e6 -251.949e6 -251.948e6 -251e6],[ 0 0 CO2_ramp_pulse CO2_ramp_pulse 0 0 ],t);
% 
% Hg_volc_force = interp1([-253e6 -251.99e6 -251.98e6 -251.56e6 -251.55e6 -251e6 ],[ 1 1 Hg_ramp Hg_ramp 1 1 ],t) ...
%        + interp1([-253e6 -251.951e6 -251.950e6 -251.949e6 -251.948e6 -251e6],[ 0 0 Hg_ramp_pulse Hg_ramp_pulse 0 0 ],t) ;
% 
% thermo_force = interp1([-253e6 -251.99e6 -251.98e6 -251.56e6 -251.55e6 -251e6 ],[ 10 10 TH_ramp TH_ramp 10 10 ],t) ;
% 
% %%%% LAND BIOTA OPTIONS
% 
% %%%% carbon burial from land biota
% C_burial_force = 1.2 ;
% OXIDW_FORCE = 1 ;
% Hg_runoff_force = 1 ;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Flux calculations   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% Thermohaline speed (Sv)
circ_TH_Sv = thermo_force ;

%%%% Thermohaline speed (m3/yr)
circ_TH_m3_yr = ( circ_TH_Sv * 1e6 ) * 3.15e7 ; 

%%%% DIC concentration in mol/m3
DIC_conc_s = DIC_s / pars.vol_s ;
DIC_conc_h = DIC_h / pars.vol_h ;
DIC_conc_d = DIC_d / pars.vol_d ;
ALK_conc_s = ALK_s / pars.vol_s ;
ALK_conc_h = ALK_h / pars.vol_h ;
ALK_conc_d = ALK_d / pars.vol_d ;
Hg_conc_s = Hg_s / pars.vol_s ;
Hg_conc_h = Hg_h / pars.vol_h ;
Hg_conc_d = Hg_d / pars.vol_d ;

%%%% Transport fluxes in mol/yr
Tran_DIC_s_h = DIC_conc_s * circ_TH_m3_yr ;
Tran_DIC_h_d = DIC_conc_h * circ_TH_m3_yr ;
Tran_DIC_d_s = DIC_conc_d * circ_TH_m3_yr ;
Tran_ALK_s_h = ALK_conc_s * circ_TH_m3_yr ;
Tran_ALK_h_d = ALK_conc_h * circ_TH_m3_yr ;
Tran_ALK_d_s = ALK_conc_d * circ_TH_m3_yr ;
f_Tran_Hg_s_h = Hg_conc_s * circ_TH_m3_yr ;
f_Tran_Hg_h_d = Hg_conc_h * circ_TH_m3_yr ;
f_Tran_Hg_d_s = Hg_conc_d * circ_TH_m3_yr ;

%%%% Global average surface temperature
Climate_Sensitivity = 3 ;
t_geol = -250 ;
GAST = 288 + Climate_Sensitivity * ( log( Atmospheric_CO2_ppm / 280 ) / log(2) )  - 7.4*(t_geol/-570)  ; 
T_s = 298 + (GAST - 288)*0.66 ;
T_h = max( 275.5 + (GAST - 288) , 271 ) ;
T_d = max( 275.5 + (GAST - 288) , 271 ) ;
T_cont = GAST ;


%%%% Carbonate chemistry parameters
k_2  = 7.4e-10 ;
k_carb_s = 0.000575 + 0.000006 * ( T_s - 278 ) ;
KCO2_s = 0.035 + 0.0019 * ( T_s - 278 ) ;
k_carb_h = 0.000575 + 0.000006 * ( T_h - 278 ) ;
KCO2_h = 0.035 + 0.0019 * ( T_h - 278 ) ;
k_carb_d = 0.000575 + 0.000006 * ( T_d - 278 ) ;


%%%% Carbonate speciation lowlat
HCO3_s = ( DIC_conc_s - ( DIC_conc_s^2 - ALK_conc_s * ( 2 * DIC_conc_s - ALK_conc_s ) * ( 1 - 4 * k_carb_s ) )^0.5  ) / ( 1 - 4 * k_carb_s ) ;
CO3_s = ( ALK_conc_s - HCO3_s ) / 2 ;
H_s =  k_2 * HCO3_s / CO3_s ;
pH_s = -1 * log10(H_s) ;
%%%% Carbonate speciation hilat
HCO3_h = ( DIC_conc_h - ( DIC_conc_h^2 - ALK_conc_h * ( 2 * DIC_conc_h - ALK_conc_h ) * ( 1 - 4 * k_carb_h ) )^0.5  ) / ( 1 - 4 * k_carb_h ) ;
CO3_h = ( ALK_conc_h - HCO3_h ) / 2 ;
H_h =  k_2 * HCO3_h / CO3_h ;
pH_h = -1 * log10(H_h) ;
%%%% Carbonate speciation deep
HCO3_d = ( DIC_conc_d - ( DIC_conc_d^2 - ALK_conc_d * ( 2 * DIC_conc_d - ALK_conc_d ) * ( 1 - 4 * k_carb_d ) )^0.5  ) / ( 1 - 4 * k_carb_d ) ;
CO3_d = ( ALK_conc_d - HCO3_d ) / 2 ;
H_d = k_2 * HCO3_d / CO3_d ;
pH_d = -1 * log10(H_d) ;

%%%% Air-sea exchange (mol/yr)  
%%%% pCO2 in PAL
pCO2_s =    KCO2_s * ( ( HCO3_s^2 ) / CO3_s ) ;
AirSea_s = 5e16 * 0.85 * 0.1 * ( pCO2_a - pCO2_s ) ;
pCO2_h =   KCO2_h * ( ( HCO3_h^2 ) / CO3_h ) ;
AirSea_h = 5e16 * 0.15 * 0.1 * ( pCO2_a - pCO2_h ) ;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Continental fluxes  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% weathering relationships
silconst = 0.33 ;
carbconst = 0.9 ;
UPLIFT = 1 ;
CARB_AREA = 1 ;
BAS_AREA = 1 ;
GRAN_AREA = 1 ;
ORG_AREA = 1 ;
PG = 1 ;
f_biota = 1 ;
O = pars.O0 ;

%%%%%% basalt and granite temp dependency - direct and runoff
f_T_bas =  exp(0.0608*(T_cont-288)) * ( (1 + 0.038*(T_cont - 288))^0.65 ) ; %%% 42KJ/mol
f_T_gran =  exp(0.0724*(T_cont-288)) * ( (1 + 0.038*(T_cont - 288))^0.65 ) ; %%% 50 KJ/mol
g_T = 1 + 0.087*(T_cont - 288) ;

%%%% basalt and granite weathering
basw = pars.k_basw * BAS_AREA * PG * f_biota * f_T_bas ;
granw = pars.k_granw * UPLIFT^silconst * GRAN_AREA * PG * f_biota * f_T_gran ;

%%% silicate weathering
silw = basw + granw ;

%%%% carbonate weathering
carbw = pars.k_carbw * CARB_AREA * UPLIFT^carbconst * PG * f_biota * g_T ;

%%%% oxidative weathering 
oxidw = pars.k_oxidw*UPLIFT^silconst*ORG_AREA*((O/pars.O0)^0.5) * OXIDW_FORCE ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Degassing fluxes  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEGASS = 1.5 ;
Bforcing = 1 ;
ccdeg = pars.k_ccdeg*DEGASS*Bforcing ;
ocdeg = pars.k_ocdeg*DEGASS ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Burial fluxes  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%% CaCO3 saturation
Ca_conc_s = 1.397e19 / 1.35e18 ;

ksp_s = 0.8 ; %%% in mM^2

sat = ( Ca_conc_s * CO3_s ) / ksp_s ;
satpresent = 3 ;

%%%% shallow carbonate burial
sat_minus_1 = max( sat - 1 , 0 ) ;
mccb = pars.k_mccb * (1 / satpresent) * ( (sat_minus_1)^1.7 ) ;

%%%% land organic C burial
locb = pars.k_locb * C_burial_force ;

%%%% marine organic C burial
mocb = pars.k_mocb  ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Mercury fluxes   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% surface ocean Hg conc at present (mM)
Hg_conc_s_0 = pars.Hg_s_0 / pars.vol_s ; 
Hg_conc_h_0 = pars.Hg_h_0 / pars.vol_h ; 

%%%% volcanic input
f_Hg_volc = pars.k_Hg_volc * Hg_volc_force ;

%%%% wildfire
f_Hg_wildfire = pars.k_Hg_wildfire ;

%%%% runoff to rivers
f_Hg_runoff = pars.k_Hg_runoff * Hg_runoff_force ;

%%%% marine burial rate
f_Hgb = pars.k_Hgb * ( Hg_conc_s / Hg_conc_s_0) ;

%%%% vegetation uptake and evasion
f_Hg_veg_dep = pars.k_Hg_vegdep * ( Hg_a / pars.Hg_a_0 ) ;
f_Hg_veg_eva = pars.k_Hg_vegevasion ;

%%%% ocean deposition and evasion
f_Hg_ocean_dep_s = pars.k_Hg_oceandep_s * ( Hg_a / pars.Hg_a_0 ) ;
f_Hg_ocean_dep_h = pars.k_Hg_oceandep_h * ( Hg_a / pars.Hg_a_0 ) ;
f_Hg_ocean_eva_s = pars.k_Hg_oceanevasion_s * ( Hg_conc_s / Hg_conc_s_0 ) ;
f_Hg_ocean_eva_h = pars.k_Hg_oceanevasion_h * ( Hg_conc_h / Hg_conc_h_0 ) ;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Reservoir calculations   %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% CO2_a
dy(1) = - AirSea_s - AirSea_h + ccdeg + ocdeg + oxidw - locb - carbw - 2*silw + CO2_input_force + CO2_input_force_spike ;

%%%% DIC_s
dy(2) = AirSea_s + Tran_DIC_d_s - Tran_DIC_s_h + 2*carbw + 2*silw - mccb - mocb ;

%%%% DIC_h
dy(3) = AirSea_h + Tran_DIC_s_h - Tran_DIC_h_d ;

%%%% DIC_d
dy(4) =  Tran_DIC_h_d - Tran_DIC_d_s  ;

%%%% ALK_s
dy(5) =  Tran_ALK_d_s - Tran_ALK_s_h + 2*carbw + 2*silw - 2*mccb ;

%%%% ALK_h
dy(6) =  Tran_ALK_s_h - Tran_ALK_h_d ;

%%%% ALK_d
dy(7) =  Tran_ALK_h_d - Tran_ALK_d_s  ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Carbon isotopes   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d13c_C = 3 ;
d13c_G = -24 ;
d13c_CO2_input_force = -25 ;
d13c_CO2_input_force_spike = -5 ;
capdelB = 27 ;

%%%% CO2_a * d13c_a
dy(8) = - AirSea_s*d13c_atm - AirSea_h*d13c_atm + ccdeg*d13c_C + ocdeg*d13c_G + oxidw*d13c_G - locb*(d13c_atm - capdelB) - carbw*d13c_atm - 2*silw*d13c_atm + CO2_input_force*d13c_CO2_input_force + CO2_input_force_spike*d13c_CO2_input_force_spike;

%%%% DIC_s * d13c_s
dy(9) = AirSea_s*d13c_atm + Tran_DIC_d_s*d13c_DIC_d - Tran_DIC_s_h*d13c_DIC_s + carbw*d13c_C + carbw*d13c_atm + 2*silw*d13c_atm - mccb*d13c_DIC_s - mocb*(d13c_DIC_s - capdelB) ;

%%%% DIC_h * d13c_h
dy(10) = AirSea_h*d13c_atm + Tran_DIC_s_h*d13c_DIC_s - Tran_DIC_h_d*d13c_DIC_h ;

%%%% DIC_d * d13c_d
dy(11) =  Tran_DIC_h_d*d13c_DIC_h - Tran_DIC_d_s*d13c_DIC_d  ;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Hg reservoir calculations %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Hg_a
dy(12) = f_Hg_volc + f_Hg_wildfire + f_Hg_ocean_eva_s + f_Hg_ocean_eva_h + f_Hg_veg_eva - f_Hg_veg_dep - f_Hg_ocean_dep_s - f_Hg_ocean_dep_h  ;

%%%% Hg_s
dy(13) = f_Hg_ocean_dep_s - f_Hg_ocean_eva_s + f_Tran_Hg_d_s - f_Tran_Hg_s_h + f_Hg_runoff - f_Hgb ;

%%%% Hg_h
dy(14) = f_Hg_ocean_dep_h - f_Hg_ocean_eva_h + f_Tran_Hg_s_h - f_Tran_Hg_h_d ;

%%%% Hg_d
dy(15) =  f_Tran_Hg_h_d - f_Tran_Hg_d_s  ;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Hg reservoir * d202hg calculations %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d202hg_a = - 1 ;
d202hg_s = y(17) / y(13) ;
d202hg_h = y(18) / y(14) ;
d202hg_d = y(19) / y(15) ;

%%%% fixed input
d202hg_veg = - 3 ;

%%%% Static atmospheric composition
dy(16) = 0 ;

%%%% Hg_s * d202hg_s
dy(17) = f_Hg_ocean_dep_s*d202hg_a - f_Hg_ocean_eva_s*d202hg_s + f_Tran_Hg_d_s*d202hg_d - f_Tran_Hg_s_h*d202hg_s + f_Hg_runoff*d202hg_veg - f_Hgb*d202hg_s;

%%%% Hg_h * d202hg_h
dy(18) = f_Hg_ocean_dep_h*d202hg_a - f_Hg_ocean_eva_h*d202hg_h + f_Tran_Hg_s_h*d202hg_s - f_Tran_Hg_h_d*d202hg_h  ;

%%%% Hg_d * d202hg_d
dy(19) =  f_Tran_Hg_h_d*d202hg_h - f_Tran_Hg_d_s*d202hg_d  ;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Save output as working   %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% record model dy states while working
workingstate.CO2_a(stepnumber,1) = CO2_a ;
workingstate.DIC_s(stepnumber,1) = DIC_s ;
workingstate.DIC_h(stepnumber,1) = DIC_h ;
workingstate.DIC_d(stepnumber,1) = DIC_d ;
workingstate.ALK_s(stepnumber,1) = ALK_s ;
workingstate.ALK_h(stepnumber,1) = ALK_h ;
workingstate.ALK_d(stepnumber,1) = ALK_d ;
workingstate.DIC_conc_s(stepnumber,1) = DIC_conc_s ;
workingstate.DIC_conc_h(stepnumber,1) = DIC_conc_h ;
workingstate.DIC_conc_d(stepnumber,1) = DIC_conc_d ;
workingstate.ALK_conc_s(stepnumber,1) = ALK_conc_s ;
workingstate.ALK_conc_h(stepnumber,1) = ALK_conc_h ;
workingstate.ALK_conc_d(stepnumber,1) = ALK_conc_d ;
workingstate.pH_s(stepnumber,1) = pH_s ;
workingstate.pH_h(stepnumber,1) = pH_h ;
workingstate.pH_d(stepnumber,1) = pH_d ;
workingstate.T_s(stepnumber,1) = T_s ;
workingstate.T_h(stepnumber,1) = T_h ;
workingstate.T_d(stepnumber,1) = T_d ;
workingstate.T_cont(stepnumber,1) = T_cont ;
workingstate.GAST(stepnumber,1) = GAST ;
workingstate.Atmospheric_CO2_ppm(stepnumber,1) = Atmospheric_CO2_ppm ;
workingstate.granw(stepnumber,1) = granw ;
workingstate.basw(stepnumber,1) = basw ;
workingstate.silw(stepnumber,1) = silw ;
workingstate.carbw(stepnumber,1) = carbw ;
workingstate.ccdeg(stepnumber,1) = ccdeg ;
workingstate.mccb(stepnumber,1) = mccb ;
workingstate.CO2_input_force(stepnumber,1) = CO2_input_force ;
workingstate.HCO3_s(stepnumber,1) = HCO3_s ;
workingstate.CO3_s(stepnumber,1) = CO3_s ;
workingstate.H_s(stepnumber,1) = H_s ;
workingstate.HCO3_h(stepnumber,1) = HCO3_h ;
workingstate.CO3_h(stepnumber,1) = CO3_h ;
workingstate.H_h(stepnumber,1) = H_h ;
workingstate.HCO3_d(stepnumber,1) = HCO3_d ;
workingstate.CO3_d(stepnumber,1) = CO3_d ;
workingstate.H_d(stepnumber,1) = H_d ;
workingstate.d13c_atm(stepnumber,1) = d13c_atm ;
workingstate.d13c_DIC_s(stepnumber,1) = d13c_DIC_s ;
workingstate.d13c_DIC_h(stepnumber,1) = d13c_DIC_h ;
workingstate.d13c_DIC_d(stepnumber,1) = d13c_DIC_d ;
workingstate.Hg_a(stepnumber,1) = Hg_a ;
workingstate.Hg_s(stepnumber,1) = Hg_s ;
workingstate.Hg_h(stepnumber,1) = Hg_h ;
workingstate.Hg_d(stepnumber,1) = Hg_d ;
workingstate.f_Hgb(stepnumber,1) = f_Hgb ;
workingstate.mocb(stepnumber,1) = mocb ;
workingstate.Hg_conc_s(stepnumber,1) = Hg_conc_s ;
workingstate.Hg_conc_h(stepnumber,1) = Hg_conc_h ;
workingstate.Hg_conc_d(stepnumber,1) = Hg_conc_d ;
workingstate.d202hg_a(stepnumber,1) = d202hg_a ;
workingstate.d202hg_s(stepnumber,1) = d202hg_s ;
workingstate.d202hg_h(stepnumber,1) = d202hg_h ;
workingstate.d202hg_d(stepnumber,1) = d202hg_d ;
workingstate.Hg_volc_force(stepnumber,1) = Hg_volc_force ;
workingstate.C_burial_force(stepnumber,1) = C_burial_force ;
workingstate.OXIDW_FORCE(stepnumber,1) = OXIDW_FORCE ;
workingstate.Hg_runoff_force(stepnumber,1) = Hg_runoff_force ;
workingstate.f_Hg_ocean_dep_s(stepnumber,1) = f_Hg_ocean_dep_s ;
workingstate.f_Hg_ocean_dep_h(stepnumber,1) = f_Hg_ocean_dep_h ;
workingstate.f_Hg_runoff(stepnumber,1) = f_Hg_runoff ;

%%%%%%%% record time
workingstate.time(stepnumber,1) = t ;
workingstate.time_myr(stepnumber,1) = t / 1e6 ;

%%%% final action: record current model step
stepnumber = stepnumber + 1 ;


end



