% Generates TF pair potential energy tables from 'Startpoint' up to a given
% length input 'Endpoint'. Plotswitch is a logical true or false to
% determine whether to plot the resulting potentials
% Recommended spacing is 0.005 angstroms or 0.0005 nm
% INPUT parameters units are: Energy = kJ/mol, length = nm
% OUTPUT units are: Energy = kJ/mol, length = nm
% These units are used in GROMACS by default.
% CRDamping is a boolean which adds a close range damping to the
% attractive R^-6 term
function [U_PM_out, U_PP_out, U_MM_out] = TF_Potential_Generator(Startpoint,...
    Endpoint,Spacing,Salt,Parameters,plotswitch,vdw_modifier,RVDW_Cutoff,CRDamping)

%% Conversion factors and fundamental constants
nm_per_m = 1e+9; % nm per m
NA = 6.0221409e23; % Molecules per mole
e_c = 1.60217662e-19; % Elementary charge in Coulombs
epsilon_0 = (8.854187817620e-12)*1000/(nm_per_m*NA); % Vacuum Permittivity C^2 mol kJ^-1 nm^-1
k_0 = 1/(4*pi*epsilon_0); % Coulomb constant in kJ nm C^-2 mol^-1

%% TF Parameter: q (charge)
q.Li =  1; % atomic
q.Na =  1; % atomic
q.K  =  1; % atomic
q.Rb =  1; % atomic
q.Cs =  1; % atomic

q.F  = -1; % atomic
q.Cl = -1; % atomic
q.Br = -1; % atomic
q.I  = -1; % atomic

%% Generate range (r) in nm
r = Startpoint:Spacing:Endpoint;

%% Split Salt Into Component Metal and Halide
[Metal,Halide] = Separate_Metal_Halide(Salt);

%% Parameters from input
alpha.PP = Parameters(1,1);
alpha.MM = Parameters(1,2);
alpha.PM = Parameters(1,3);

B.PP = Parameters(2,1);
B.MM = Parameters(2,2);
B.PM = Parameters(2,3);

C.PP = Parameters(3,1);
C.MM = Parameters(3,2);
C.PM = Parameters(3,3);

D.PP = Parameters(4,1);
D.MM = Parameters(4,2);
D.PM = Parameters(4,3);

%% Build Model Table: Plus-Minus (PM)

% Plus - Minus total potential
U_MetHal.Total = k_0*(e_c^2)*q.(Metal)*q.(Halide)./(r) ...
    + B.PM*exp(-alpha.PM.*r) ...
    - C.PM./(r.^6) ...
    - D.PM./(r.^8);

% Components of potential
U_MetHal.f = 1./r; % Electrostatics function f(r)
U_MetHal.g = - C.PM./(r.^6) - D.PM./(r.^8); % Dispersion g(r)
U_MetHal.h = B.PM*exp(-alpha.PM.*r);% Short range repulsion h(r)

% Plus - Minus total derivative
U_MetHal.dTotal = -k_0*(e_c^2)*q.(Metal)*q.(Halide)./(r.^2) ...
    - alpha.PM*B.PM*exp(-alpha.PM.*r) ...
    + 6.*C.PM./(r.^7) ...
    + 8.*D.PM./(r.^9);

% Components of derivative
U_MetHal.df = 1./(r.^2);% Electrostatics function -df(r)/dr
U_MetHal.dg = - 6.*C.PM./(r.^7) - 8.*D.PM./(r.^9) ; % Dispersion -dg(r)/dr
U_MetHal.dh = alpha.PM*B.PM*exp(-alpha.PM.*r);% Short range repulsion -dh(r)/dr

if contains(vdw_modifier,'potential-shift','IgnoreCase',true)
    EVDW_Cutoff = B.PM*exp(-alpha.PM.*RVDW_Cutoff) ...
    - C.PM./(RVDW_Cutoff.^6) ...
    - D.PM./(RVDW_Cutoff.^8);
    
    % Shift by the dispersion energy at vdw cutoff radius. only affects one
    % energy component, not derivatives (i.e. forces)
    U_MetHal.Total = U_MetHal.Total - EVDW_Cutoff;
    U_MetHal.g = U_MetHal.g - EVDW_Cutoff;
end

% Remove infinities (replace with zero)
U_MetHal = Remove_Infinities(U_MetHal);

%% Build Model Table: Plus - Plus (PP)

% Plus - Plus Total Potential
U_MetMet.Total = k_0*(e_c^2)*q.(Metal).*q.(Metal)./(r) ...
    + B.PP*exp(-alpha.PP.*r) ...
    - C.PP./(r.^6) ...
    - D.PP./(r.^8);

% Components of potential
U_MetMet.f = 1./r; % Electrostatics function f(r)
U_MetMet.g = - C.PP./(r.^6) - D.PP./(r.^8); % Dispersion g(r)
U_MetMet.h = B.PP*exp(-alpha.PP.*r); % Short range repulsion

% Plus - Plus Total Derivative
U_MetMet.dTotal = -k_0*(e_c^2)*q.(Metal).*q.(Metal)./(r.^2) ...
    - alpha.PP*B.PP*exp(-alpha.PP.*r)...
    + 6.*C.PP./(r.^7) ...
    + 8.*D.PP./(r.^9);

% Components of derivative
U_MetMet.df = 1./(r.^2); % Electrostatics function f(r)
U_MetMet.dg = - 6.*C.PP./(r.^7) - 8.*D.PP./(r.^9); % Dispersion g(r)
U_MetMet.dh = alpha.PP*B.PP*exp(-alpha.PP.*r); % Short range repulsion

if contains(vdw_modifier,'potential-shift','IgnoreCase',true)
    EVDW_Cutoff = B.PP*exp(-alpha.PP.*RVDW_Cutoff) ...
    - C.PP./(RVDW_Cutoff.^6) ...
    - D.PP./(RVDW_Cutoff.^8);
    
    % Shift by the dispersion energy at vdw cutoff radius. only affects one
    % energy component, not derivatives (i.e. forces)
    U_MetMet.Total = U_MetMet.Total - EVDW_Cutoff;
    U_MetMet.g = U_MetMet.g - EVDW_Cutoff;
end

% Remove infinities (replace with zero)
U_MetMet = Remove_Infinities(U_MetMet);

%% Build Model Table: Minus - Minus (MM)
% Minus - Minus Total Potential
U_HalHal.Total = k_0*(e_c^2)*q.(Halide).*q.(Halide)./(r)...
    + B.MM*exp(-alpha.MM.*r) ...
    - C.MM./(r.^6)...
    - D.MM./(r.^8);

% Components of potential
U_HalHal.f = 1./r; % Electrostatics function f(r)
U_HalHal.g = - C.MM./(r.^6) - D.MM./(r.^8); % Dispersion g(r)
U_HalHal.h = B.MM*exp(-alpha.MM.*r); % Short range repulsion

% Minus - Minus Total Derivative
U_HalHal.dTotal = -k_0*(e_c^2)*q.(Halide).*q.(Halide)./(r.^2)...
    - alpha.MM*B.MM*exp(-alpha.MM.*r)...
    + 6.*C.MM./(r.^7)...
    + 8.*D.MM./(r.^9);

% Components of derivative
U_HalHal.df = 1./(r.^2); % Electrostatics function f(r)
U_HalHal.dg = - 6.*C.MM./(r.^7) - 8.*D.MM./(r.^9); % Dispersion g(r)
U_HalHal.dh = alpha.MM*B.MM*exp(-alpha.MM.*r); % Short range repulsion

if contains(vdw_modifier,'potential-shift','IgnoreCase',true)
    EVDW_Cutoff = B.MM*exp(-alpha.MM.*RVDW_Cutoff) ...
    - C.MM./(RVDW_Cutoff.^6)...
    - D.MM./(RVDW_Cutoff.^8);
    
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

    h{1} = plot(r,U_MetHal.Total,'Color','k','LineWidth',lw,'LineStyle','-');
    h{2} = plot(r,U_MetMet.Total,'Color','r','LineWidth',lw,'LineStyle',':');
    h{3} = plot(r,U_HalHal.Total,'Color','g','LineWidth',lw,'LineStyle','-.');

    title(['Plot of TF Potentials for ' Salt],...
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
    leg1 = legend([h{:}],{['TF - ' Salt] ['TF - ' Metal Metal] ['TF - ' Halide Halide]});
    leg1.Interpreter = 'latex';
end

end