%% Params for Electrochemical Model
%   Created May 24, 2012 by Scott Moura
%
%   Most parameters are from DUALFOIL model
%   D_e is from Capiglia et al, for c_e = 1000 mol/m^3
%   Equilibrium potentials are from Bosch Klein TCST 2011
%   Electrode area chosen to correspond to 2.3 Ah cell

%% Geometric Params
% Thickness of each layer
p.L_n = 100e-6;     % Thickness of negative electrode [m]
p.L_s = 25e-6;     % Thickness of separator [m]
p.L_p = 100e-6;     % Thickness of positive electrode [m]

L_ccn = 25e-6;    % Thickness of negative current collector [m]
L_ccp = 25e-6;    % Thickness of negative current collector [m]

% Particle Radii
p.R_s_n = 10e-6;   % Radius of solid particles in negative electrode [m]
p.R_s_p = 10e-6;   % Radius of solid particles in positive electrode [m]

% Volume fractions
p.epsilon_s_n = 0.6;      % Volume fraction in solid for neg. electrode
p.epsilon_s_p = 0.5;      % Volume fraction in solid for pos. electrode

p.epsilon_e_n = 0.3;   % Volume fraction in electrolyte for neg. electrode
p.epsilon_e_s = 1.0;   % Volume fraction in electrolyte for separator
p.epsilon_e_p = 0.3;   % Volume fraction in electrolyte for pos. electrode

epsilon_f_n = 0.1;  % Volume fraction of filler in neg. electrode
epsilon_f_p = 0.2;  % Volume fraction of filler in pos. electrode

% Specific interfacial surface area
p.a_s_n = 3*p.epsilon_s_n / p.R_s_n;  % Negative electrode [m^2/m^3]
p.a_s_p = 3*p.epsilon_s_p / p.R_s_p;  % Positive electrode [m^2/m^3]

% Mass densities
rho_sn = 1800;    % Solid phase in negative electrode [kg/m^3]
rho_sp = 5010;    % Solid phase in positive electrode [kg/m^3]
rho_e =  1324;    % Electrolyte [kg/m^3]
rho_f = 1800;     % Filler [kg/m^3]
rho_ccn = 8954;   % Current collector in negative electrode
rho_ccp = 2707;   % Current collector in positive electrode

% Compute cell mass [kg/m^2]
m_n = p.L_n * (rho_e*p.epsilon_e_n + rho_sn*p.epsilon_s_n + rho_f*epsilon_f_n);
m_s = p.L_s * (rho_e*p.epsilon_e_n);
m_p = p.L_p * (rho_e*p.epsilon_e_p + rho_sp*p.epsilon_s_p + rho_f*epsilon_f_p);
m_cc = rho_ccn*L_ccn + rho_ccp*L_ccp;

% Lumped density [kg/m^2]
p.rho_avg = m_n + m_s + m_p + m_cc;

%% Transport Params
% Diffusion coefficient in solid
p.D_s_n = 3.9e-14;  % Diffusion coeff for solid in neg. electrode, [m^2/s]
p.D_s_p = 1e-13;  % Diffusion coeff for solid in pos. electrode, [m^2/s]

% Diffusion coefficient in electrolyte
% p.D_e = 2.7877e-10;    % Diffusion coeff for electrolyte, [m^2/s]

p.brug = 1.5;       % Bruggeman porosity

% Conductivity of solid
p.sig_n = 100;    % Conductivity of solid in neg. electrode, [1/Ohms*m]
p.sig_p = 10;    % Conductivity of solid in pos. electrode, [1/Ohms*m]

% p.sig_eff_n = p.sig_n * p.epsilon_s_n^p.brug;    % Eff. conductivity in neg. electrode, [1/Ohms*m]
% p.sig_eff_p = p.sig_p * p.epsilon_s_p^p.brug;    % Eff. conductivity in pos. electrode, [1/Ohms*m]

% Conductivity of electrolyte

% Miscellaneous
p.t_plus = 0.4;       % Transference number
p.Faraday = 96487;    % Faraday's constant, [Coulumbs/mol]
p.Area = 1;           % Electrode current collector area [m^2]

%% Kinetic Params
p.R = 8.314472;       % Gas constant, [J/mol-K]
p.alph = 0.5;         % Charge transfer coefficients
p.R_SEI = 1e-3;       % Resistivity of SEI layer, [Ohms*m^2]

% Reaction rates
p.k_n = 1e-5;  % Reaction rate in neg. electrode, [(A/m^2)*(mol^3/mol)^(1+alpha)]
p.k_p = 3e-7; % Reaction rate in pos. electrode, [(A/m^2)*(mol^3/mol)^(1+alpha)]

%% Thermodynamic Params
% Cubic Splines of Equilibrium potentials (ws/pp_LiFEPO4.mat)
% load('ws/BoschOCP');
% 
% ppn = spline(theta,OCPn);
% dppn = ppn;
% dppn.coefs = [3*ppn.coefs(:,1), 2*ppn.coefs(:,2), ppn.coefs(:,3)];
% dppn.order = 3;
% 
% ppp = spline(theta,OCPp);
% dppp = ppp;
% dppp.coefs = [3*ppp.coefs(:,1), 2*ppp.coefs(:,2), ppp.coefs(:,3)];
% dppp.order = 3;
% 
% p.Uppn = ppn;
% p.dUppn = dppn;
% p.Uppp = ppp;
% p.dUppp = dppp;

% Thermal dynamics
p.C_p = 2000;   % Heat capacity, [J/kg-K]
p.h = 0.36;   % Heat transfer coefficient, [W/K-m^2] 0

% Ambient Temperature
p.T_amp = 298; % [K]

% Entropy coefficients
p.dUref_dT = -0.4e-3; % [V/K] approx. from Al Hallaj et al 2000, JPS

%% Concentrations
% % Maxima based on 2.3Ah cell
% p.c_s_n_max = 1.4990e4;   % Max concentration in anode, [mol/m^3]
% p.c_s_p_max = 2.5609e4;    % Max concentration in cathode, [mol/m^3]

% Maxima based on DUALFOIL
p.c_s_n_max = 3.6e3 * 372 * 1800 / p.Faraday;   % Max concentration in anode, [mol/m^3]
p.c_s_p_max = 3.6e3 * 247 * 5010 / p.Faraday;    % Max concentration in cathode, [mol/m^3]

p.c_e = 1e3;              % Fixed electrolyte concentration for SPM, [mol/m^3]

%% Cutoff voltages
p.volt_max = 4.7;
p.volt_min = 2.0;

%% Discretization parameters
% Discrete time step
p.delta_t = 1;

% Pade Order
p.PadeOrder = 3;

% Finite difference points along r-coordinate
% p.Nr = 10;
% p.delta_r_n = p.R_s_n / p.Nr;
% p.delta_r_p = p.R_s_p / p.Nr;

% Finite difference points along x-coordinate
p.Nxn = 70;
p.Nxs = 35;
p.Nxp = 70;
p.Nx = p.Nxn+p.Nxs+p.Nxp;

p.delta_x_n = 1 / p.Nxn;
p.delta_x_s = 1 / p.Nxs;
p.delta_x_p = 1 / p.Nxp;


