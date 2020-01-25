% Generates Jung-Cheatham pair potential energy tables from 'Startpoint' up to a given
% length input 'Endpoint'. Plotswitch is a logical true or false to
% determine whether to plot the resulting potentials
% Recommended spacing is 0.005 angstroms or 0.0005 nm
% INPUT parameters units are: Energy = kJ/mol, length = nm
% OUTPUT units are: Energy = kJ/mol, length = nm
% These units are used in GROMACS by default.
% CRDamping is a boolean which adds a close range damping to the
% attractive R^-6 term
function [U_PM_out, U_PP_out, U_MM_out] = JC_Potential_Generator(Startpoint,...
    Endpoint,Spacing,Salt,Parameters,plotswitch,vdw_modifier,RVDW_Cutoff,CRDamping)

%% Conversion factors and fundamental constants
nm_per_m = 1e+9; % nm per m
NA = 6.0221409e23; % Molecules per mole
e_c = 1.60217662e-19; % Elementary charge in Coulombs
epsilon_0 = (8.854187817620e-12)*1000/(nm_per_m*NA); % Vacuum Permittivity C^2 mol kJ^-1 nm^-1
k_0 = 1/(4*pi*epsilon_0); % Coulomb constant in kJ nm C^-2 mol^-1
r_d = 0.10; % If close-range dispersion damping is on, damp at this distance
b = 150; % sigmoid "steepness" for damping

[Metal,Halide] = Separate_Metal_Halide(Salt);

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


%% Parameters from input for LJ potential
sigma_PP = Parameters(1,1);
sigma_MM = Parameters(1,2);
sigma_PM = Parameters(1,3);

epsilon_PP = Parameters(2,1);
epsilon_MM = Parameters(2,2);
epsilon_PM = Parameters(2,3);

% Change parameters into A/r12 - B/r6 format
A_PP = 4*epsilon_PP*(sigma_PP^12);
B_PP = 4*epsilon_PP*(sigma_PP^6);

A_MM = 4*epsilon_MM*(sigma_MM^12);
B_MM = 4*epsilon_MM*(sigma_MM^6);

A_PM = 4*epsilon_PM*(sigma_PM^12);
B_PM = 4*epsilon_PM*(sigma_PM^6);

%% Generate range (r) in nm
r = Startpoint:Spacing:Endpoint;
Nr = length(r);

%% Build PES: Plus-Minus

% If Damping at close range
if CRDamping
    
    f_r = 1./(1 + exp(-b.*(r - r_d))); % sigmoid damping function
    df_r = (b.*exp(-b.*(r - r_d)))./((1 + exp(-b.*(r - r_d))).^2); % sigmoid damping function derivative
    
    % Plus - Minus total potential (damp attractions at close range)
    U_MetHal.Total = f_r.*k_0*(e_c^2)*q.(Metal)*q.(Halide)./r ...
    + A_PM./(r.^12) ...
    - f_r.*B_PM./(r.^6);

    % components
    U_MetHal.f = 1./r; % Electrostatics function f(r)
    U_MetHal.g = -f_r.*B_PM./(r.^6); % Dispersion g(r)
    U_MetHal.h = A_PM./(r.^12) - k_0*(e_c^2)*q.(Metal)*q.(Halide)./(r) + ...
        f_r.*k_0*(e_c^2)*q.(Metal)*q.(Halide)./(r);% Short range repulsion
    
    % Plus - Minus total derivative
    U_MetHal.dTotal = -f_r.*k_0*(e_c^2)*q.(Metal)*q.(Halide)./(r.^2) ...
        + df_r.*k_0*(e_c^2)*q.(Metal)*q.(Halide)./(r)...
        - A_PM.*12./(r.^13) ...
        + f_r.*B_PM.*6./(r.^7) ...
        - df_r.*B_PM./(r.^6);

    % components
    U_MetHal.df = 1./(r.^2);% Electrostatics function
    U_MetHal.dg = - f_r.*B_PM.*6./(r.^7) + df_r.*B_PM./(r.^6); % Dispersion
    U_MetHal.dh = + A_PM.*12./(r.^13) - k_0*(e_c^2)*q.(Metal)*q.(Halide)./(r.^2) ...
        + f_r.*k_0*(e_c^2)*q.(Metal)*q.(Halide)./(r.^2) ...
        - df_r.*k_0*(e_c^2)*q.(Metal)*q.(Halide)./(r);% Short range repulsion
    
% No damping    
else
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
end

if contains(vdw_modifier,'potential-shift','IgnoreCase',true)
    EVDW_Cutoff = A_PM./(RVDW_Cutoff.^12) ...
    - B_PM./(RVDW_Cutoff.^6);
    
    % Shift by the dispersion energy at vdw cutoff radius. only affects one
    % energy component, not derivatives (i.e. forces)
    U_MetHal.Total = U_MetHal.Total - EVDW_Cutoff;
    U_MetHal.g = U_MetHal.g - EVDW_Cutoff;
end

% remove infinities
U_MetHal = Remove_Infinities(U_MetHal);

%% Build PES: Plus - Plus
% If Damping at close range
if CRDamping
    % Plus - Plus total potential
    U_MetMet.Total = k_0*(e_c^2)*q.(Metal)*q.(Metal)./(r) ...
        + A_PP./(r.^12) ...
        - f_r.*B_PP./(r.^6);

    % components
    U_MetMet.f = 1./r; % Electrostatics function f(r)
    U_MetMet.g = -f_r.*B_PP./(r.^6); % Dispersion g(r)
    U_MetMet.h = A_PP./(r.^12);% Short range repulsion

    % Plus - Plus total derivative
    U_MetMet.dTotal = -k_0*(e_c^2)*q.(Metal)*q.(Metal)./(r.^2) ...
        - A_PP.*12./(r.^13) ...
        + f_r.*B_PP.*6./(r.^7) ...
        - df_r.*B_PP./(r.^6);

    % components
    U_MetMet.df = 1./(r.^2);% Electrostatics function
    U_MetMet.dg = - f_r.*B_PP.*6./(r.^7) + df_r.*B_PP./(r.^6); % Dispersion
    U_MetMet.dh = + A_PP.*12./(r.^13);% Short range repulsion
else
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
end

if contains(vdw_modifier,'potential-shift','IgnoreCase',true)
    EVDW_Cutoff = A_PP./(RVDW_Cutoff.^12) ...
    - B_PP./(RVDW_Cutoff.^6);
    
    % Shift by the dispersion energy at vdw cutoff radius. only affects one
    % energy component, not derivatives (i.e. forces)
    U_MetMet.Total = U_MetMet.Total - EVDW_Cutoff;
    U_MetMet.g = U_MetMet.g - EVDW_Cutoff;
end

% remove infinities
U_MetMet = Remove_Infinities(U_MetMet);

%% Build PES: Minus - Minus

if CRDamping
    
    % Minus - Minus total potential
    U_HalHal.Total = k_0*(e_c^2)*q.(Halide)*q.(Halide)./(r) ...
        + A_MM./(r.^12) ...
        - f_r.*B_MM./(r.^6);

    % components
    U_HalHal.f = 1./r; % Electrostatics function f(r)
    U_HalHal.g = -f_r.*B_MM./(r.^6); % Dispersion g(r)
    U_HalHal.h = A_MM./(r.^12);% Short range repulsion

    % Minus - Minus total derivative
    U_HalHal.dTotal = -k_0*(e_c^2)*q.(Halide)*q.(Halide)./(r.^2) ...
        - A_MM.*12./(r.^13) ...
        + f_r.*B_MM.*6./(r.^7) ...
        - df_r.*B_MM./(r.^6);

    % components
    U_HalHal.df = 1./(r.^2);% Electrostatics function
    U_HalHal.dg = - f_r.*B_MM.*6./(r.^7) + df_r.*B_MM./(r.^6); % Dispersion
    U_HalHal.dh = + A_MM.*12./(r.^13);% Short range repulsion
else
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
end

if contains(vdw_modifier,'potential-shift','IgnoreCase',true)
    EVDW_Cutoff = A_MM./(RVDW_Cutoff.^12) ...
    - B_MM./(RVDW_Cutoff.^6);
    
    % Shift by the dispersion energy at vdw cutoff radius. only affects one
    % energy component, not derivatives (i.e. forces)
    U_HalHal.Total = U_HalHal.Total - EVDW_Cutoff;
    U_HalHal.g = U_HalHal.g - EVDW_Cutoff;
end

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


    dr = r(2) - r(1);

    h{1} = plot(r,U_MetHal.Total,'Color','k','LineWidth',lw,'LineStyle','-');
    h{2} = plot(r,U_MetMet.Total,'Color','r','LineWidth',lw,'LineStyle','-');
    h{3} = plot(r,U_HalHal.Total,'Color','g','LineWidth',lw,'LineStyle','-');

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