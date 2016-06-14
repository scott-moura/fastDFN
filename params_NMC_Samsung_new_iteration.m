% %% Params for Electrochemical Model
% %   Created June 3, 2011 by Scott Moura
% %   Modified March 01, 2014 by Hector Perez 
%  
% %   (Feb 22, 2014) Adjusted for lengths to be in SI units [m]
% %                  These parameters are from Forman 2012, GA ParamID Study 
% %                  doi:10.1016/j.jpowsour.2012.03.009
%  
% %   (Mar 01, 2014) Added Parameter Sensitivity Nominal Parameters for rSPM 
%  
% %   (Jun 26, 2015) Modified to work with dfn_scott.m
%  
%  
% %   (Sep 23, 2015) Modified to use NMC parameters from [Fang DOI: 10.1002/er.1652  , Ji DOI:10.1149/2.047304jes  ,  Tanim DOI:10.1016/j.energy.2014.12.031 ]
%  
% %% Geometric Params
% % Thickness of each layer
% p.L_n = 12.3E-05;% <--FROM Samsung %4.00E-05; %2.885e-5;     % Thickness of negative electrode [m]
% p.L_s = 2E-5;% <--FROM Samsung %2.50E-05; %1.697e-5;     % Thickness of separator [m]
% p.L_p = 11.9E-5;% <--FROM Samsung %3.66E-05; %6.521e-5;     % Thickness of positive electrode [m]
%  
% % Particle Radii
% p.R_s_n = 5.00E-7;%5.00E-06; %3.596e-6;   % Radius of solid particles in negative electrode [m]
% p.R_s_p = 5.00E-7; %1.637e-7;   % Radius of solid particles in positive electrode [m]
%  
% % Volume fractions
% p.epsilon_s_n = 0.7215;%0.662; %0.3810;   % Volume fraction in solid for neg. electrode
% p.epsilon_s_p = 0.6516;%0.58; %0.4800;   % Volume fraction in solid for pos. electrode
%  
% p.epsilon_e_n = 0.3; %0.6190;   % Volume fraction in electrolyte for neg. electrode
% p.epsilon_e_s = 0.4; %0.3041;   % Volume fraction in electrolyte for separator
% p.epsilon_e_p = 0.3; %0.5200;   % Volume fraction in electrolyte for pos. electrode
%  
% % Specific interfacial surface area
% p.a_s_n = 3*p.epsilon_s_n / p.R_s_n;  % Negative electrode [m^2/m^3]
% p.a_s_p = 3*p.epsilon_s_p / p.R_s_p;  % Positive electrode [m^2/m^3]
%  
% %% Transport Params
% % Diffusion coefficient in solid
% p.D_s_n = 1.40E-14; %1.736e-14;  % Diffusion coeff for solid in neg. electrode, [m^2/s]
% p.D_s_p = 2.00E-14; %8.256e-14;  % Diffusion coeff for solid in pos. electrode, [m^2/s]
%  
% % Diffusion coefficient in electrolyte
% p.D_e = 1.50E-10; %6.911e-10;    % Diffusion coeff for electrolyte, [m^2/s]
%  
% p.brug = 1.5; %1.452;       % Bruggeman porosity
% p.D_e_n = p.D_e * p.epsilon_e_n^p.brug; % Effective diffusion coef. in neg. electrode, [m^2/s]
% p.D_e_s = p.D_e * p.epsilon_e_s^p.brug; % Effective diffusion coef. in neg. electrode, [m^2/s]
% p.D_e_p = p.D_e * p.epsilon_e_p^p.brug; % Effective diffusion coef. in neg. electrode, [m^2/s]
%  
% % Conductivity of solid
% p.sig_n = 1.00E+02; %100;    % Conductivity of solid in neg. electrode, [1/Ohms*m]
% p.sig_p = 1.00E+01; %100;    % Conductivity of solid in pos. electrode, [1/Ohms*m]
%  
% p.sig_eff_n = p.sig_n * p.epsilon_s_n;    % Eff. conductivity in neg. electrode, [1/Ohms*m]
% p.sig_eff_p = p.sig_p * p.epsilon_s_p;    % Eff. conductivity in pos. electrode, [1/Ohms*m]
%  
% % Conductivity of electrolyte
%  
% % Miscellaneous
% p.t_plus = 0.38; %0.2495;    % Transference number
% p.Faraday = 96487;    % Faraday's constant, [Coulumbs/mol]
% p.Area = 0.048965;%0.046;%<-- This is calculated from Samsung electrode geometry %1.02E-01; %0.3108;      % Electrode current collector area [m^2]
%  
% %% Kinetic Params
% p.R = 8.314472;       % Gas constant, [J/mol-K]
% p.alph = 0.5;         % Charge transfer coefficients
% p.R_SEI = 0;%3.391e-3;   % Resistivity of SEI layer, [Ohms*m^2]
% p.R_f_n = 15.00E-03; %p.R_SEI;    % Resistivity of SEI layer, [Ohms*m^2]
% p.R_f_p = 10.00E-03; %0;          % Resistivity of SEI layer, [Ohms*m^2]
%  
% % Reaction rates
% p.k_n = 1e-7; % Reaction rate in neg. electrode, [(A/cm^2)*(mol^3/mol)^(1+alpha)]
% p.k_p = 3e-7; % Reaction rate in pos. electrode, [(A/cm^2)*(mol^3/mol)^(1+alpha)]
%  
% %% Thermodynamic Params
% % Equilibrium potentials (ws/Uref_LiFEPO4.mat)
% load('ws/Uref_LiFePO4');
% p.theta_n = theta_n;
% p.theta_p = theta_p;
% p.Unref = Unref;
% p.Upref = Upref;
%  
% % Thermal dynamics
% p.C_p = 75;   % Heat capacity, [J/K]
% p.h = 12.4;   % Heat transfer coefficient, [W/K]
%  
% % Ambient Temperature
% p.T_amp = 298.15; % [K]
%  
% % Entropy coefficients
% p.dUref_dT = -0.4e-3; % [V/K] approx. from Al Hallaj et al 2000, JPS
%  
% % Lumped density [kg/m^2]
% p.rho_avg = 0.7934;
%  
% %% Concentrations
% % Maxima
% p.c_s_n_max = 3.1168e+04; % This is calculated from SAMSUNG OCP data, %2.948e4;    % Max concentration in anode, [mol/m^3]
% p.c_s_p_max = 4.2649e+04;% This is calculated from SAMSUNG OCP data, %5.18E+04; %1.035e4;    % Max concentration in cathode, [mol/m^3]
%  
% p.n_Li_s = 0.1391;%<-- This is calculated based on Samsung cell voltage limits and OCP maps %0.100687787424000;    % Total moles of lithium in solid phase [mol]
% p.c_e = 1e3;              % Fixed electrolyte concentration for SPM, [mol/m^3]
%  
% %% Cutoff voltages
% p.volt_max = 5;%3.8;
% p.volt_min = 2.0;%1.8;
%  
% %% Discretization parameters
% % Discrete time step
% p.delta_t = 1;
%  
% % Pade Order
% p.PadeOrder = 3;
%  
% % Finite difference points along r-coordinate
% % p.Nr = 10;
% % p.delta_r_n = p.R_s_n / p.Nr;
% % p.delta_r_p = p.R_s_p / p.Nr;
%  
% % Finite difference points along x-coordinate
% p.Nxn = 10;
% p.Nxs = 10;
% p.Nxp = 10;
% p.Nx = p.Nxn+p.Nxs+p.Nxp;
%  
% p.delta_x_n = 1 / p.Nxn;
% p.delta_x_s = 1 / p.Nxs;
% p.delta_x_p = 1 / p.Nxp;
%  
% %% Bessel Function Regression Params
% p.poly_I1 = [0.06994, -0.00641, 0.5016, -8.554e-5];
% p.poly_I2 = [0.02194, 0.1105, 0.003257, -0.0001637];
%  
% %% Parameter Sensitivity Nominal Parameters for rSPM Backstepping PDE Observer
% p.epsilon0 = 1;
% p.q0 = 1;
% p.alpha0 = 1;
% p.delta0=1;





















%% Params for Electrochemical Model
%   Created June 3, 2011 by Scott Moura
%   Modified March 01, 2014 by Hector Perez 
 
%   (Feb 22, 2014) Adjusted for lengths to be in SI units [m]
%                  These parameters are from Forman 2012, GA ParamID Study 
%                  doi:10.1016/j.jpowsour.2012.03.009
 
%   (Mar 01, 2014) Added Parameter Sensitivity Nominal Parameters for rSPM 
 
%   (Jun 26, 2015) Modified to work with dfn_scott.m
 
 
%   (Sep 23, 2015) Modified to use NMC parameters from [Fang DOI: 10.1002/er.1652  , Ji DOI:10.1149/2.047304jes  ,  Tanim DOI:10.1016/j.energy.2014.12.031 ]
 
%% Geometric Params
% Thickness of each layer
p.L_n = 12.3E-05;% <--FROM Samsung %4.00E-05; %2.885e-5;     % Thickness of negative electrode [m]
p.L_s = 2E-5;% <--FROM Samsung %2.50E-05; %1.697e-5;     % Thickness of separator [m]
p.L_p = 11.9E-5;% <--FROM Samsung %3.66E-05; %6.521e-5;     % Thickness of positive electrode [m]
 
% Particle Radii
p.R_s_n = 6.1708e-07;%5e-7;%1e-7;%5.00E-7;%5.00E-06; %3.596e-6;   % Radius of solid particles in negative electrode [m]
p.R_s_p = 5.0000e-07;%5e-7;%1e-7;%5.00E-7; %1.637e-7;   % Radius of solid particles in positive electrode [m]
 
% Volume fractions
p.epsilon_s_n = 0.7215;%0.662; %0.3810;   % Volume fraction in solid for neg. electrode
p.epsilon_s_p = 0.6516;%0.58; %0.4800;   % Volume fraction in solid for pos. electrode
 
p.epsilon_f_n = 0.0789;%0;%Filler
p.epsilon_f_p = 0.0998;%0;%Filler

p.epsilon_e_n = 1 - p.epsilon_s_n - p.epsilon_f_n;%0.3; %0.6190;   % Volume fraction in electrolyte for neg. electrode
p.epsilon_e_s = 0.4; %0.3041;   % Volume fraction in electrolyte for separator
p.epsilon_e_p = 1 - p.epsilon_s_p - p.epsilon_f_p; %0.5200;   % Volume fraction in electrolyte for pos. electrode
 
% Specific interfacial surface area
p.a_s_n = 3*p.epsilon_s_n/p.R_s_n;  % Negative electrode [m^2/m^3]
p.a_s_p = 3*p.epsilon_s_p/p.R_s_p;  % Positive electrode [m^2/m^3]
 
%% Transport Params
% Diffusion coefficient in solid
p.D_s_n = 8.8565e-13;%9.9945e-13;%1.40E-14; %1.736e-14;  % Diffusion coeff for solid in neg. electrode, [m^2/s]
p.D_s_p = 6.6027e-12;%1.0000e-14;%2.00E-14; %8.256e-14;  % Diffusion coeff for solid in pos. electrode, [m^2/s]
%p.D_s_n = 1.40E-14;%1.40E-14; %1.736e-14;  % Diffusion coeff for solid in neg. electrode, [m^2/s]
%p.D_s_p = 2.00E-14;%2.00E-14; %8.256e-14;  % Diffusion coeff for solid in pos. electrode, [m^2/s]
 
% Diffusion coefficient in electrolyte
p.D_e = 2.0864e-10;%1.50E-10; %6.911e-10;    % Diffusion coeff for electrolyte, [m^2/s]
 
% Diffusional conductivity in electrolyte
p.dactivity = 0;

p.brug = 1.5; %1.452;       % Bruggeman porosity
p.D_e_n = p.D_e * p.epsilon_e_n^p.brug; % Effective diffusion coef. in neg. electrode, [m^2/s]
p.D_e_s = p.D_e * p.epsilon_e_s^p.brug; % Effective diffusion coef. in neg. electrode, [m^2/s]
p.D_e_p = p.D_e * p.epsilon_e_p^p.brug; % Effective diffusion coef. in neg. electrode, [m^2/s]
 
% Conductivity of solid
p.sig_n = 1.00E+02; %100;    % Conductivity of solid in neg. electrode, [1/Ohms*m]
p.sig_p = 1.00E+01; %100;    % Conductivity of solid in pos. electrode, [1/Ohms*m]
 
p.sig_eff_n = p.sig_n * p.epsilon_s_n;    % Eff. conductivity in neg. electrode, [1/Ohms*m]
p.sig_eff_p = p.sig_p * p.epsilon_s_p;    % Eff. conductivity in pos. electrode, [1/Ohms*m]
 
% Conductivity of electrolyte
 
% Miscellaneous
 p.t_plus = 0.3662;%0.38; %0.2495;    % Transference number
p.Faraday = 96487;    % Faraday's constant, [Coulumbs/mol]
p.Area = 0.048965;%0.046;%<-- This is calculated from Samsung electrode geometry %1.02E-01; %0.3108;      % Electrode current collector area [m^2]
 
%% Kinetic Params
p.R = 8.314472;       % Gas constant, [J/mol-K]
p.alph = 0.5;         % Charge transfer coefficients
p.R_SEI = 0;%3.391e-3;   % Resistivity of SEI layer, [Ohms*m^2]
p.R_f_n = 15.00E-03; %p.R_SEI;    % Resistivity of SEI layer, [Ohms*m^2]
p.R_f_p = 10.00E-03; %0;          % Resistivity of SEI layer, [Ohms*m^2]
p.R_c = 5.1874e-05;%4.0000e-06;%7.00E-04; %0;          % Contact Resistance/Current Collector Resistance, [Ohms]
%p.R_c = 7.00E-04;
  
% Reaction rates
p.k_n = 1e-7; % Reaction rate in neg. electrode, [(A/cm^2)*(mol^3/mol)^(1+alpha)]
p.k_p = 3e-7; % Reaction rate in pos. electrode, [(A/cm^2)*(mol^3/mol)^(1+alpha)]
 
%% Thermodynamic Params
% Equilibrium potentials (ws/Uref_LiFEPO4.mat)
load('ws/Uref_LiFePO4');
p.theta_n = theta_n;
p.theta_p = theta_p;
p.Unref = Unref;
p.Upref = Upref;
 
% Thermal dynamics
p.C_p = 2000*0.2;%75;   % Heat capacity, [J/K]
p.h = 1;%12.4;   % Heat transfer coefficient, [W/K]
 
% Ambient Temperature
p.T_amp = 298.15; % [K]
 
% Entropy coefficients
p.dUref_dT = -0.4e-3; % [V/K] approx. from Al Hallaj et al 2000, JPS
 
% Lumped density [kg/m^2]
p.rho_avg = 0.7934;
 
%% Concentrations
% Maxima
p.c_s_n_max = 3.1168e+04; % This is calculated from SAMSUNG OCP data, %2.948e4;    % Max concentration in anode, [mol/m^3]
p.c_s_p_max = 4.2649e+04;% This is calculated from SAMSUNG OCP data, %5.18E+04; %1.035e4;    % Max concentration in cathode, [mol/m^3]
 
p.n_Li_s = 0.1391;%<-- This is calculated based on Samsung cell voltage limits and OCP maps %0.100687787424000;    % Total moles of lithium in solid phase [mol]
p.c_e = 1e3;              % Fixed electrolyte concentration for SPM, [mol/m^3]
 
%% Cutoff voltages
p.volt_max = 5;%3.8;
p.volt_min = 2.0;%1.8;
 
%% Discretization parameters
% Discrete time step
p.delta_t = 0.5;
 
% Pade Order
p.PadeOrder = 3;
 
% Finite difference points along r-coordinate
% p.Nr = 10;
% p.delta_r_n = p.R_s_n / p.Nr;
% p.delta_r_p = p.R_s_p / p.Nr;
 
% Finite difference points along x-coordinate
p.Nxn = 70;%10;
p.Nxs = 35;%10;
p.Nxp = 70;%10;
p.Nx = p.Nxn+p.Nxs+p.Nxp;
 
p.delta_x_n = 1 / p.Nxn;
p.delta_x_s = 1 / p.Nxs;
p.delta_x_p = 1 / p.Nxp;
 
%% Bessel Function Regression Params
p.poly_I1 = [0.06994, -0.00641, 0.5016, -8.554e-5];
p.poly_I2 = [0.02194, 0.1105, 0.003257, -0.0001637];
 
%% Parameter Sensitivity Nominal Parameters for rSPM Backstepping PDE Observer
p.epsilon0 = 1;
p.q0 = 1;
p.alpha0 = 1;
p.delta0=1;

