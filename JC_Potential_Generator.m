% Jung-Cheatham model parameters adapted for three water models with
% Lorenz-Berthelot combining rules
% Use Watermodel = 'SPC/E', 'TIP3P', or 'TIP4PEW'
% Generates pair potential energy surfaces from 'Startpoint' up to a given
% length input 'Endpoint'. Plotswitch is a logical true or false to
% determine whether to plot the resulting potentials
% Recommended spacing is 0.005 angstroms or 0.0005 nm
% INPUT UNITS MUST BE ALL IN NANOMETERS. OUTPUTS ARE IN NANOMETERS AND
% kJ/mol
% DS is the dispersion scaling factor (should only affect the r6 term)
% ES is the epsilon scaling factor (increases well depth)
function [U_PM_out, U_PP_out, U_MM_out] = JC_Potential_Generator(Startpoint,Endpoint,Spacing,Salt,Watermodel,plotswitch,DS,ES)
%% Conversion factors
nm_per_Ang = 0.1; % nm per Angstrom
kJ_per_kcal = 4.184; % kj per kcal

%% Conversion factors and fundamental constants
kj_per_erg = 1e-10; % kJ per erg
Ang_per_cm = 1e+8; % Angstroms per cm
nm_per_cm = 1e+7; % nm per cm
Ang_per_m = 1e+10; % Angstroms per m
nm_per_m = 1e+9; % nm per m
NA = 6.0221409e23; % Molecules per mole
e_c = 1.60217662e-19; % Elementary charge in Coulombs
epsilon_0 = (8.854187817620e-12)*1000/(nm_per_m*NA); % Vacuum Permittivity C^2 mol kJ^-1 nm^-1
k_0 = 1/(4*pi*epsilon_0); % Coulomb constant in kJ nm C^-2 mol^-1
C_unit = 1e-60; % erg cm^6
D_unit = 1e-76; % erg cm^8

[Metal,Halide] = Separate_Metal_Halide(Salt);

%% JC Ion Parameters in SPC/E water
if strcmp(Watermodel,'SPC/E')
    Param.Li.sigma = (0.791*nm_per_Ang)*(2^(5/6)); % nm
    Param.Li.epsilon = (0.3367344)*kJ_per_kcal; % kJ/mol

    Param.Na.sigma = (1.212*nm_per_Ang)*(2^(5/6)); % nm
    Param.Na.epsilon = (0.3526418)*kJ_per_kcal; % kJ/mol

    Param.K.sigma = (1.593*nm_per_Ang)*(2^(5/6)); % nm
    Param.K.epsilon = (0.4297054)*kJ_per_kcal; % kJ/mol

    Param.Rb.sigma = (1.737*nm_per_Ang)*(2^(5/6)); % nm
    Param.Rb.epsilon = (0.4451036)*kJ_per_kcal; % kJ/mol

    Param.Cs.sigma = (2.021*nm_per_Ang)*(2^(5/6)); % nm
    Param.Cs.epsilon = (0.0898565)*kJ_per_kcal; % kJ/mol

    Param.F.sigma = (2.257*nm_per_Ang)*(2^(5/6)); % nm
    Param.F.epsilon = (0.0074005)*kJ_per_kcal; % kJ/mol

    Param.Cl.sigma = (2.711*nm_per_Ang)*(2^(5/6)); % nm
    Param.Cl.epsilon = (0.0127850)*kJ_per_kcal; % kJ/mol

    Param.Br.sigma = (2.751*nm_per_Ang)*(2^(5/6)); % nm
    Param.Br.epsilon = (0.0269586)*kJ_per_kcal; % kJ/mol

    Param.I.sigma = (2.919*nm_per_Ang)*(2^(5/6)); % nm
    Param.I.epsilon = (0.0427845)*kJ_per_kcal; % kJ/mol

elseif strcmp(Watermodel,'TIP3P')
    %% JC Ion Parameters in TIP3P water
    Param.Li.sigma = (1.025*nm_per_Ang)*(2^(5/6)); % nm
    Param.Li.epsilon = (0.0279896)*kJ_per_kcal; % kJ/mol

    Param.Na.sigma = (1.369*nm_per_Ang)*(2^(5/6)); % nm
    Param.Na.epsilon = (0.0874393)*kJ_per_kcal; % kJ/mol

    Param.K.sigma = (1.705*nm_per_Ang)*(2^(5/6)); % nm
    Param.K.epsilon = (0.1936829)*kJ_per_kcal; % kJ/mol

    Param.Rb.sigma = (1.813*nm_per_Ang)*(2^(5/6)); % nm
    Param.Rb.epsilon = (0.3278219)*kJ_per_kcal; % kJ/mol

    Param.Cs.sigma = (1.976*nm_per_Ang)*(2^(5/6)); % nm
    Param.Cs.epsilon = (0.4065394)*kJ_per_kcal; % kJ/mol

    Param.F.sigma = (2.303*nm_per_Ang)*(2^(5/6)); % nm
    Param.F.epsilon = (0.0033640)*kJ_per_kcal; % kJ/mol

    Param.Cl.sigma = (2.513*nm_per_Ang)*(2^(5/6)); % nm
    Param.Cl.epsilon = (0.0355910)*kJ_per_kcal; % kJ/mol

    Param.Br.sigma = (2.608*nm_per_Ang)*(2^(5/6)); % nm
    Param.Br.epsilon = (0.0586554)*kJ_per_kcal; % kJ/mol

    Param.I.sigma = (2.860*nm_per_Ang)*(2^(5/6)); % nm
    Param.I.epsilon = (0.0536816)*kJ_per_kcal; % kJ/mol

elseif strcmp(Watermodel,'TIP4PEW')
    %% JC Ion Parameters in TIP4P water
    Param.Li.sigma = (0.808*nm_per_Ang)*(2^(5/6)); % nm
    Param.Li.epsilon = (0.1039884)*kJ_per_kcal; % kJ/mol

    Param.Na.sigma = (1.226*nm_per_Ang)*(2^(5/6)); % nm
    Param.Na.epsilon = (0.1684375)*kJ_per_kcal; % kJ/mol

    Param.K.sigma = (1.590*nm_per_Ang)*(2^(5/6)); % nm
    Param.K.epsilon = (0.2794651)*kJ_per_kcal; % kJ/mol

    Param.Rb.sigma = (1.709*nm_per_Ang)*(2^(5/6)); % nm
    Param.Rb.epsilon = (0.4331494)*kJ_per_kcal; % kJ/mol

    Param.Cs.sigma = (1.888*nm_per_Ang)*(2^(5/6)); % nm
    Param.Cs.epsilon = (0.3944318)*kJ_per_kcal; % kJ/mol

    Param.F.sigma = (2.538*nm_per_Ang)*(2^(5/6)); % nm
    Param.F.epsilon = (0.0015752)*kJ_per_kcal; % kJ/mol

    Param.Cl.sigma = (2.760*nm_per_Ang)*(2^(5/6)); % nm
    Param.Cl.epsilon = (0.0116615)*kJ_per_kcal; % kJ/mol

    Param.Br.sigma = (2.768*nm_per_Ang)*(2^(5/6)); % nm
    Param.Br.epsilon = (0.0303773)*kJ_per_kcal; % kJ/mol

    Param.I.sigma = (2.952*nm_per_Ang)*(2^(5/6)); % nm
    Param.I.epsilon = (0.0417082)*kJ_per_kcal; % kJ/mol
else
    error(['Unknown Water Model: "' Watermodel ...
        '". Please choose one of "SPC/E", "TIP3P", or "TIP4PEW".'])
end

%% Parameter: q (charge)
q.Li =  1; % atomic
q.Na =  1; % atomic
q.K  =  1; % atomic
q.Rb =  1; % atomic
q.Cs =  1; % atomic

q.F  = -1; % atomic
q.Cl = -1; % atomic
q.Br = -1; % atomic
q.I  = -1; % atomic


%% Calculate parameters of interest for LJ potential
sigma_PP = Param.(Metal).sigma;
sigma_MM = Param.(Halide).sigma;
sigma_PM = ( sigma_PP + sigma_MM )/2;

epsilon_PP = Param.(Metal).epsilon;
epsilon_MM = Param.(Halide).epsilon;
epsilon_PM = sqrt(epsilon_PP*epsilon_MM);

% Change parameters into A/r12 - B/r6 format
A_PP = ES*4*epsilon_PP*(sigma_PP^12);
B_PP = DS*ES*4*epsilon_PP*(sigma_PP^6);

A_MM = ES*4*epsilon_MM*(sigma_MM^12);
B_MM = DS*ES*4*epsilon_MM*(sigma_MM^6);

A_PM = ES*4*epsilon_PM*(sigma_PM^12);
B_PM = DS*ES*4*epsilon_PM*(sigma_PM^6);

%% Generate range (r) in nm
r = Startpoint:Spacing:Endpoint;

%% Build PES: Plus-Minus

% Plus - Minus total potential
U_MetHal.Total = k_0*(e_c^2)*q.(Metal)*q.(Halide)./(r) ...
    + A_PM./(r.^12) ...
    - B_PM./(r.^6);

% components
U_MetHal.f = 1./r; % Electrostatics function f(r)
U_MetHal.g = -B_PM./(r.^6); % Dispersion g(r)
U_MetHal.h = A_PM./(r.^12);% Short range repulsion

% Plus - Minus total derivative
U_MetHal.dTotal = -k_0*(e_c^2)*q.(Metal)*q.(Halide)./(r.^2) ...
    - A_PM.*12./(r.^13) ...
    + B_PM.*6./(r.^7);

% components
U_MetHal.df = 1./(r.^2);% Electrostatics function
U_MetHal.dg = - B_PM.*6./(r.^7); % Dispersion
U_MetHal.dh = + A_PM.*12./(r.^13);% Short range repulsion

% remove infinities
U_MetHal = Remove_Infinities(U_MetHal);

%% Build PES: Plus - Plus

% Plus - Plus total potential
U_MetMet.Total = k_0*(e_c^2)*q.(Metal)*q.(Metal)./(r) ...
    + A_PP./(r.^12) ...
    - B_PP./(r.^6);

% components
U_MetMet.f = 1./r; % Electrostatics function f(r)
U_MetMet.g = -B_PP./(r.^6); % Dispersion g(r)
U_MetMet.h = A_PP./(r.^12);% Short range repulsion

% Plus - Plus total derivative
U_MetMet.dTotal = -k_0*(e_c^2)*q.(Metal)*q.(Metal)./(r.^2) ...
    - A_PP.*12./(r.^13) ...
    + B_PP.*6./(r.^7);

% components
U_MetMet.df = 1./(r.^2);% Electrostatics function
U_MetMet.dg = - B_PP.*6./(r.^7); % Dispersion
U_MetMet.dh = + A_PP.*12./(r.^13);% Short range repulsion

% remove infinities
U_MetMet = Remove_Infinities(U_MetMet);

%% Build PES: Minus - Minus

% Minus - Minus total potential
U_HalHal.Total = k_0*(e_c^2)*q.(Halide)*q.(Halide)./(r) ...
    + A_MM./(r.^12) ...
    - B_MM./(r.^6);

% components
U_HalHal.f = 1./r; % Electrostatics function f(r)
U_HalHal.g = -B_MM./(r.^6); % Dispersion g(r)
U_HalHal.h = A_MM./(r.^12);% Short range repulsion

% Minus - Minus total derivative
U_HalHal.dTotal = -k_0*(e_c^2)*q.(Halide)*q.(Halide)./(r.^2) ...
    - A_MM.*12./(r.^13) ...
    + B_MM.*6./(r.^7);

% components
U_HalHal.df = 1./(r.^2);% Electrostatics function
U_HalHal.dg = - B_MM.*6./(r.^7); % Dispersion
U_HalHal.dh = + A_MM.*12./(r.^13);% Short range repulsion

% remove infinities
U_HalHal = Remove_Infinities(U_HalHal);

%% Print
PM = [r ; U_MetHal.f ; U_MetHal.df ; U_MetHal.g ; U_MetHal.dg ; U_MetHal.h ; U_MetHal.dh];
U_PM_out = deblank( sprintf(['%16.10e   %16.10e %16.10e   %16.10e %16.10e   %16.10e %16.10e' newline],PM(:)) );

PP = [r ; U_MetMet.f ; U_MetMet.df ; U_MetMet.g ; U_MetMet.dg ; U_MetMet.h ; U_MetMet.dh];
U_PP_out = deblank( sprintf(['%16.10e   %16.10e %16.10e   %16.10e %16.10e   %16.10e %16.10e' newline],PP(:)) );

MM = [r ; U_HalHal.f ; U_HalHal.df ; U_HalHal.g ; U_HalHal.dg ; U_HalHal.h ; U_HalHal.dh];
U_MM_out = deblank( sprintf(['%16.10e   %16.10e %16.10e   %16.10e %16.10e   %16.10e %16.10e' newline],MM(:)) );


%% PLOT if plotswitch chosen
if plotswitch
    
    % Options
    lw=2;
    fs=25;

    h = cell(1,8);
    hold on

    h{1} = plot(r,U_MetHal.Total,'Color','k','LineWidth',lw,'LineStyle','-');
    h{2} = plot(r,U_MetMet.Total,'Color','r','LineWidth',lw,'LineStyle',':');
    h{3} = plot(r,U_HalHal.Total,'Color','g','LineWidth',lw,'LineStyle','-.');

    title(['Plot of JC Potentials for ' Salt],...
        'Interpreter','latex','fontsize',fs)

    set(gca,'box','on','TickLabelInterpreter','latex');
    set(gca,'XMinorTick','on','YMinorTick','on','FontSize',fs);
    xlabel('r (nm)','fontsize',fs,'Interpreter','latex');
    ylabel('Pair Potential (kJ mol$^{-1}$)','fontsize',fs,'Interpreter','latex');

    ylim([-900 1000]);
    xlim([Startpoint Endpoint]);

    % Blank line
    hline = refline([0 0]);
    hline.Color = 'k';
    hline.LineWidth = lw-1;
    hline.LineStyle = '--';
    leg1 = legend([h{:}],{['JC - ' Salt] ['JC - ' Metal Metal] ['JC - ' Halide Halide]});
    %leg1 = legend([h{:}],{['TF - ' Salt ' V(r)'] ['TF - ' Salt ' V''(r) Analytical'] ['TF - ' Salt ' V''(r) Numerical']});
%     leg1 = legend([h{:}],{['TF - ' Salt ' $V(r)$'] ['TF - ' Salt ' $f(r)$'] ...
%         ['TF - ' Salt ' $g(r)$'] ['TF - ' Salt ' $h(r)$'] ...
%         ['TF - ' Salt ' $\frac{d V (r)}{d r}$'] ['TF - ' Salt ' $\frac{d g (r)}{d r}$'] ...
%         ['TF - ' Salt ' $\frac{d h (r)}{d r}$'] ['TF - ' Salt ' $\frac{d h (r)}{d r}$']});
    leg1.Interpreter = 'latex';
end

% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% if plotswitch
%     %% Conversion Factors / Constants
%     m_per_nm = 1e-9; % nm per m
%     NA = 6.0221409e23; % Molecules per mole
%     e_c = 1.60217662e-19; % Elementary charge in Coulombs
%     J_per_kJ = 1000;
%     epsilon_0 = (8.854187817620e-12)*J_per_kJ*m_per_nm/NA; % Vacuum Permittivity C^2 mol kJ^-1 nm^-1
%     k_0 = 1/(4*pi*epsilon_0); % Coulomb constant in C^-2 nm kJ mol^-1
% 
%     % Use combining rules
%     epsilon_MM = Param.(Metal).epsilon; % kJ/mol
%     epsilon_HH = Param.(Halide).epsilon; % kJ/mol
%     epsilon_Salt = sqrt(Param.(Metal).epsilon*Param.(Halide).epsilon); % kJ/mol
%     
%     sigma_MM = Param.(Metal).sigma;
%     sigma_HH = Param.(Halide).sigma;
%     sigma_Salt = (1/2)*(Param.(Metal).sigma + Param.(Halide).sigma);
% 
%     q_M = 1*e_c; % charge of lithium in C
%     q_H = -1*e_c; % charge of fluoride in C
% 
%     %% Build PES for LiF
%     r = 0:0.001:1; % r vector
% 
%     U_Salt = k_0.*q_M.*q_H./(r) + 4*epsilon_Salt.*( ( sigma_Salt./r ).^12 - ( sigma_Salt./r ).^6 );
% 
%     U_MetMet= k_0.*q_M.*q_M./(r) + 4*epsilon_MM.*( ( sigma_MM./r ).^12 - ( sigma_MM./r ).^6 );
% 
%     U_HalHal  = k_0.*q_H.*q_H./(r) + 4*epsilon_HH.*( ( sigma_HH./r ).^12 - ( sigma_HH./r ).^6 );
% 
%     % Options
%     lw=2;
%     fs=25;
% 
%     hold on
%     h1 = plot(r,U_Salt,'Color','r','LineWidth',lw,'LineStyle','-');
%     h2 = plot(r,U_MetMet,'Color','b','LineWidth',lw,'LineStyle',':');
%     h3 = plot(r,U_HalHal,'Color','g','LineWidth',lw,'LineStyle','-.');
% 
%     title(['Plot of ' Metal Halide ' JC Potentials for ' Watermodel ...
%         ' Water Model'],'Interpreter','latex','fontsize',fs)
% 
%     set(gca,'box','on','TickLabelInterpreter','latex');
%     set(gca,'XMinorTick','on','YMinorTick','on','FontSize',fs);
%     xlabel('r (nm)','fontsize',fs,'Interpreter','latex');
%     ylabel('Pair Potential (kJ mol$^{-1}$)','fontsize',fs,'Interpreter','latex');
% 
%     ylim([-900 1000]); 
% 
%     % Blank line
%     hline = refline([0 0]);
%     hline.Color = 'k';
%     hline.LineWidth = 1;
%     hline.LineStyle = ':';
%     legend([h1 h2 h3],{['JC - ' Metal Halide] ['JC - ' Metal Metal] ['JC - ' Halide Halide]})
% end
end