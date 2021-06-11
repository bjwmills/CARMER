%%%%%% CAR-MER
%%%%%% Coupled C-Hg 4-box ocean-atmosphere model
%%%%%% As used in Dal Corso et al. 2020 Nat. Comm.
%%%%%% Coded by Benjamin JW Mills // b.mills@leeds.ac.uk
%%%%%% Model frontend file: run this code to solve model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Define parameters   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% set up global structures
global stepnumber
global pars
global workingstate

%%%%%% water reservoir sizes in m3
pars.vol_s = 3.07e16 ;
pars.vol_h = 1.35e16 ;
pars.vol_d = 1.35e18 ;
pars.vol_ocean = pars.vol_s + pars.vol_h + pars.vol_d ;

%%%%%% inorganic carbon reservoirs in moles C
pars.CO2_a_0 = 5e16 ;
pars.DIC_s_0 = 6e16 ;
pars.DIC_h_0 = 3e16 ;
pars.DIC_d_0 = 3e18 ;
pars.ALK_s_0 = 6e16 ;
pars.ALK_h_0 = 3e16 ;
pars.ALK_d_0 = 3e18 ;

%%%%%% oxygen
pars.O0 = 3.7e19 ;

%%%%%% C isotope composition
pars.d13c_atm_0 = -7 ;
pars.d13c_DIC_s_0 = 0.1 ;
pars.d13c_DIC_h_0 = 0.1 ;
pars.d13c_DIC_d_0 = 0.1 ;

%%%%%% Mercury abundance
pars.total_ocean_Hg = 6.08e8 ; %%% Amos et al., 2013 BGC
pars.Hg_a_0 = 3.5e6 ;
pars.Hg_s_0 = pars.total_ocean_Hg * (pars.vol_s/(pars.vol_ocean)) ;
pars.Hg_h_0 = pars.total_ocean_Hg * (pars.vol_h/(pars.vol_ocean)) ;
pars.Hg_d_0 = pars.total_ocean_Hg * (pars.vol_d/(pars.vol_ocean)) ;

%%%%%% present day rates
pars.k_ccdeg = 8e12 ;
pars.k_carbw = 8e12 ;
pars.k_sfw = 0 ;
pars.k_mccb = pars.k_carbw + pars.k_ccdeg - pars.k_sfw ;
pars.k_silw = pars.k_mccb - pars.k_carbw ;
basfrac = 0.3 ;
pars.k_granw = pars.k_silw * (1-basfrac) ;
pars.k_basw = pars.k_silw * basfrac ;

%%%%%% organic C cycle
pars.k_ocdeg = 1.25e12 ;
pars.k_locb = 4.5e12 ;
pars.k_mocb = 4.5e12 ;
pars.k_oxidw = pars.k_mocb + pars.k_locb - pars.k_ocdeg ;

%%%% present day Hg fluxes
pars.k_Hgb = 2e6 ; %%%% set to equal volcanism + wildfire
pars.k_Hg_runoff = 2e6 ; %%%% set to equal burial
pars.k_Hg_volc = 1.5e6 ; %%%% from Amos 190 Mg/yr burial at present, x 1.5 for more degassing
pars.k_Hg_oceandep_total = 1.5e7 ; %%%% roughly Amos et al. 
pars.k_Hg_oceandep_s = pars.k_Hg_oceandep_total*0.85 ;
pars.k_Hg_oceandep_h = pars.k_Hg_oceandep_total*0.15 ;
pars.k_Hg_oceanevasion_total = 1.5e7 ;
pars.k_Hg_oceanevasion_s = pars.k_Hg_oceanevasion_total*0.85 ;
pars.k_Hg_oceanevasion_h = pars.k_Hg_oceanevasion_total*0.15 ;
pars.k_Hg_vegdep = 1e7 ; %%%% roughly Amos et lal. 
pars.k_Hg_vegevasion = 1e7 ;
pars.k_Hg_wildfire = 0.5e6 ;

%%% Hg isotopes
pars.d202hg_Hg_a_0 = -2 ;
pars.d202hg_Hg_s_0 = -2 ;
pars.d202hg_Hg_h_0 = -2 ;
pars.d202hg_Hg_d_0 = -2 ;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Initialise   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% set maximum step size for solver
options = odeset('maxstep',1e3) ;

%%%% run beginning
fprintf('Beginning run: \n')

%%%% set stepnumber to 1
stepnumber = 1 ;

%%%%%%% model timeframe in years (0 = present day)
pars.whenstart = -253e6 ;
pars.whenend = -251e6 ;

%%%% model start state
pars.startstate(1) = pars.CO2_a_0 ;
pars.startstate(2) = pars.DIC_s_0 ;
pars.startstate(3) = pars.DIC_h_0 ;
pars.startstate(4) = pars.DIC_d_0 ;
pars.startstate(5) = pars.ALK_s_0 ;
pars.startstate(6) = pars.ALK_h_0 ;
pars.startstate(7) = pars.ALK_d_0 ;
pars.startstate(8) = pars.CO2_a_0 * pars.d13c_atm_0 ;
pars.startstate(9) = pars.DIC_s_0 * pars.d13c_DIC_s_0 ;
pars.startstate(10) = pars.DIC_h_0 * pars.d13c_DIC_h_0 ;
pars.startstate(11) = pars.DIC_d_0 * pars.d13c_DIC_d_0 ;
pars.startstate(12) = pars.Hg_a_0 ;
pars.startstate(13) = pars.Hg_s_0 ;
pars.startstate(14) = pars.Hg_h_0 ;
pars.startstate(15) = pars.Hg_d_0 ;
pars.startstate(16) = pars.Hg_a_0 * pars.d202hg_Hg_a_0  ;
pars.startstate(17) = pars.Hg_s_0 * pars.d202hg_Hg_s_0 ;
pars.startstate(18) = pars.Hg_h_0 * pars.d202hg_Hg_h_0 ;
pars.startstate(19) = pars.Hg_d_0 * pars.d202hg_Hg_d_0 ;

%%%%%%% run the system 
[rawoutput.T,rawoutput.Y] = ode15s(@CARMER_equations,[pars.whenstart pars.whenend],pars.startstate,options);

%%%%% start time counter
tic



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Postprocessing   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% takes 'workingstate' from model and turns into 'state' %%%%%%%

%%%% size of output 
pars.output_length = length(rawoutput.T) ;
%%%%%%%%%% model finished output to screen
fprintf('Integration finished \t') ; fprintf('Total steps: %d \t' , stepnumber ) ; fprintf('Output steps: %d \n' , pars.output_length ) 
toc

%%%%%%%%% print final model states using final state for each timepoint
%%%%%%%%% during integration
fprintf('assembling state vectors... \t')
tic

%%%% trecords is index of shared values between ode15s output T vector and
%%%% model recorded workingstate t vector
[sharedvals,trecords] = intersect(workingstate.time,rawoutput.T,'stable') ;

%%%%%% assemble output state vectors
field_names = fieldnames(workingstate) ;
for numfields = 1:length(field_names)
    eval([' state.' char( field_names(numfields) ) ' = workingstate.' char( field_names(numfields) ) '(trecords) ; '])
end

%%%%%% done message
fprintf('Done: ')
endtime = toc ;
fprintf('time (s): %d \n', endtime )



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Plotting script   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% load data
load d13c_hg_data.mat

%%%% make one figure for everything, use subplots for each plot
figure

%%%% CO2 (ppm)
subplot(4,5,1)
hold on
box on
plot(state.time_myr,state.Atmospheric_CO2_ppm,'k')
xlabel('Time (Ma)')
ylabel('Atm. CO_{2} (ppm)')

%%%% Ocean DIC
subplot(4,5,2)
hold on
box on
plot(state.time_myr,state.DIC_conc_s,'r--')
plot(state.time_myr,state.DIC_conc_h,'c--')
plot(state.time_myr,state.DIC_conc_d,'b--')
xlabel('Time (Ma)')
ylabel('DIC (mol m^{-3})')

%%%% Ocean ALK
subplot(4,5,3)
hold on
box on
plot(state.time_myr,state.ALK_conc_s,'r--')
plot(state.time_myr,state.ALK_conc_h,'c--')
plot(state.time_myr,state.ALK_conc_d,'b--')
xlabel('Time (Ma)')
ylabel('ALK (mol m^{-3})')

%%%% Ocean HCO3
subplot(4,5,4)
hold on
box on
plot(state.time_myr,state.HCO3_s,'r--')
plot(state.time_myr,state.HCO3_h,'c--')
plot(state.time_myr,state.HCO3_d,'b--')
xlabel('Time (Ma)')
ylabel('HCO_{3} (mol m^{-3})')

%%%% Ocean CO3
subplot(4,5,5)
hold on
box on
plot(state.time_myr,state.CO3_s,'r--')
plot(state.time_myr,state.CO3_h,'c--')
plot(state.time_myr,state.CO3_d,'b--')
xlabel('Time (Ma)')
ylabel('CO_{3} (mol m^{-3})')

%%%% Ocean H+
subplot(4,5,6)
hold on
box on
plot(state.time_myr,state.H_s,'r--')
plot(state.time_myr,state.H_h,'c--')
plot(state.time_myr,state.H_d,'b--')
xlabel('Time (Ma)')
ylabel('H^{+} (mol m^{-3})')

%%%% Ocean pH
subplot(4,5,7)
hold on
box on
plot(state.time_myr,state.pH_s,'r--')
plot(state.time_myr,state.pH_h,'c--')
plot(state.time_myr,state.pH_d,'b--')
xlabel('Time (Ma)')
ylabel('pH')

%%%% Temperature
subplot(4,5,8)
hold on
box on
plot(state.time_myr,state.T_s,'r--')
plot(state.time_myr,state.T_h,'c--')
plot(state.time_myr,state.T_d,'b--')
plot(state.time_myr,state.T_cont,'m--')
plot(state.time_myr,state.GAST,'k--')
xlabel('Time (Ma)')
ylabel('T')

%%%% Degassing
subplot(4,5,9)
hold on
box on
plot(state.time_myr,state.ccdeg,'c')
xlabel('Time (Ma)')
ylabel('Degassing')

%%%% Continental weathering
subplot(4,5,10)
hold on
box on
plot(state.time_myr,state.basw,'k')
plot(state.time_myr,state.granw,'r')
plot(state.time_myr,state.silw,'b')
plot(state.time_myr,state.carbw,'c')
xlabel('Time (Ma)')
ylabel('Weathering')

%%%% Burial fluxes
subplot(4,5,11)
hold on
box on
plot(state.time_myr,state.mccb,'c')
xlabel('Time (Ma)')
ylabel('Burial')

%%%% d13C
subplot(4,5,12)
hold on
box on
plot(d13c_age,d13c_val,'g')
plot(state.time_myr,state.d13c_atm,'k')
plot(state.time_myr,state.d13c_DIC_s,'r')
plot(state.time_myr,state.d13c_DIC_h,'c')
plot(state.time_myr,state.d13c_DIC_d,'b')
xlabel('Time (Ma)')
ylabel('\delta^{13}C values')

%%%% input forcing
subplot(4,5,13)
hold on
box on
plot(state.time_myr,state.CO2_input_force,'k')
plot(state.time_myr,state.Hg_volc_force,'r')
xlabel('Time (Ma)')
ylabel('Input forcing')

%%%% biosphere forcing
subplot(4,5,14)
box on
semilogy(state.time_myr,state.C_burial_force,'k')
hold on
semilogy(state.time_myr,state.OXIDW_FORCE,'r')
semilogy(state.time_myr,state.Hg_runoff_force,'g')
xlabel('Time (Ma)')
ylabel('Input forcing')

%%%% Hg ocean
subplot(4,5,15)
hold on
box on
%%%% conveted to pM
plot(state.time_myr,state.Hg_conc_s*1e9,'r--')
plot(state.time_myr,state.Hg_conc_h*1e9,'c-')
plot(state.time_myr,state.Hg_conc_d*1e9,'b--')
xlabel('Time (Ma)')
ylabel('[Hg] (pM)')

%%%% Hg atmosphere
subplot(4,5,16)
hold on
box on
%%%% conveted to pM
plot(state.time_myr,state.Hg_a,'r')
xlabel('Time (Ma)')
ylabel('Atmospheric Hg (mol)')

%%%% Hg d202 isotopes
subplot(4,5,17)
hold on
box on
plot(d202hg_age,d202hg_val,'g')
plot(state.time_myr,state.d202hg_a,'k')
plot(state.time_myr,state.d202hg_s,'r--')
plot(state.time_myr,state.d202hg_h,'c-')
plot(state.time_myr,state.d202hg_d,'b--')
xlabel('Time (Ma)')
ylabel('d202hg Hg')

%%%% Hg burial versus data
subplot(4,5,18)
hold on
box on
plot(hg_age,hg_val,'g')
plot(state.time_myr,( state.f_Hgb./state.mocb ) * 8e7,'k')
xlabel('Time (Myr)')
ylabel('Hg/TOC vs model Hg/Corg molar')

%%%% Hg volc vs river input
subplot(4,5,19)
hold on
box on
plot(state.time_myr, ( state.f_Hg_ocean_dep_s + state.f_Hg_ocean_dep_h ) ./ state.f_Hg_runoff  ,'k')
xlabel('Time (Myr)')
ylabel('Hg atm input vs Hg runoff')


