% Jung-Cheatham model parameters adapted for three water models with
% Lorenz-Berthelot combining rules
% Use Watermodel = 'SPC/E', 'TIP3P', or 'TIP4PEW'
% Output units: Sigma is nanometers, epsilon is kJ/mol
% DS is the dispersion scaling factor (should only affect the r6 term)
% ES is the epsilon scaling factor (increases well depth)
function [OutputMet,OutputHal] = JC_Potential_Parameters(Metal,Halide,Watermodel,plotswitch,DS,ES)
%% Conversion factors
nm_per_Ang = 0.1; % nm per Angstrom
kJ_per_kcal = 4.184; % kj per kcal

% Convert dispersion scaling factor in terms of epsilon and sigma scaling
Beta = 1/(DS^(1/6)); % Sigma scaling factor
Gamma = DS^2; % Epsilon scaling factor

%% JC Ion Parameters in SPC/E water
if strcmp(Watermodel,'SPC/E')
    Param.Li.sigma = Beta*(0.791*nm_per_Ang)*(2^(5/6)); % nm
    Param.Li.epsilon = ES*Gamma*(0.3367344)*kJ_per_kcal; % kJ/mol

    Param.Na.sigma = Beta*(1.212*nm_per_Ang)*(2^(5/6)); % nm
    Param.Na.epsilon = ES*Gamma*(0.3526418)*kJ_per_kcal; % kJ/mol

    Param.K.sigma = Beta*(1.593*nm_per_Ang)*(2^(5/6)); % nm
    Param.K.epsilon = ES*Gamma*(0.4297054)*kJ_per_kcal; % kJ/mol

    Param.Rb.sigma = Beta*(1.737*nm_per_Ang)*(2^(5/6)); % nm
    Param.Rb.epsilon = ES*Gamma*(0.4451036)*kJ_per_kcal; % kJ/mol

    Param.Cs.sigma = Beta*(2.021*nm_per_Ang)*(2^(5/6)); % nm
    Param.Cs.epsilon = ES*Gamma*(0.0898565)*kJ_per_kcal; % kJ/mol

    Param.F.sigma = Beta*(2.257*nm_per_Ang)*(2^(5/6)); % nm
    Param.F.epsilon = ES*Gamma*(0.0074005)*kJ_per_kcal; % kJ/mol

    Param.Cl.sigma = Beta*(2.711*nm_per_Ang)*(2^(5/6)); % nm
    Param.Cl.epsilon = ES*Gamma*(0.0127850)*kJ_per_kcal; % kJ/mol

    Param.Br.sigma = Beta*(2.751*nm_per_Ang)*(2^(5/6)); % nm
    Param.Br.epsilon = ES*Gamma*(0.0269586)*kJ_per_kcal; % kJ/mol

    Param.I.sigma = Beta*(2.919*nm_per_Ang)*(2^(5/6)); % nm
    Param.I.epsilon = ES*Gamma*(0.0427845)*kJ_per_kcal; % kJ/mol

elseif strcmp(Watermodel,'TIP3P')
    %% JC Ion Parameters in TIP3P water
    Param.Li.sigma = Beta*(1.025*nm_per_Ang)*(2^(5/6)); % nm
    Param.Li.epsilon = ES*Gamma*(0.0279896)*kJ_per_kcal; % kJ/mol

    Param.Na.sigma = Beta*(1.369*nm_per_Ang)*(2^(5/6)); % nm
    Param.Na.epsilon = ES*Gamma*(0.0874393)*kJ_per_kcal; % kJ/mol

    Param.K.sigma = Beta*(1.705*nm_per_Ang)*(2^(5/6)); % nm
    Param.K.epsilon = ES*Gamma*(0.1936829)*kJ_per_kcal; % kJ/mol

    Param.Rb.sigma = Beta*(1.813*nm_per_Ang)*(2^(5/6)); % nm
    Param.Rb.epsilon = ES*Gamma*(0.3278219)*kJ_per_kcal; % kJ/mol

    Param.Cs.sigma = Beta*(1.976*nm_per_Ang)*(2^(5/6)); % nm
    Param.Cs.epsilon = ES*Gamma*(0.4065394)*kJ_per_kcal; % kJ/mol

    Param.F.sigma = Beta*(2.303*nm_per_Ang)*(2^(5/6)); % nm
    Param.F.epsilon = ES*Gamma*(0.0033640)*kJ_per_kcal; % kJ/mol

    Param.Cl.sigma = Beta*(2.513*nm_per_Ang)*(2^(5/6)); % nm
    Param.Cl.epsilon = ES*Gamma*(0.0355910)*kJ_per_kcal; % kJ/mol

    Param.Br.sigma = Beta*(2.608*nm_per_Ang)*(2^(5/6)); % nm
    Param.Br.epsilon = ES*Gamma*(0.0586554)*kJ_per_kcal; % kJ/mol

    Param.I.sigma = Beta*(2.860*nm_per_Ang)*(2^(5/6)); % nm
    Param.I.epsilon = ES*Gamma*(0.0536816)*kJ_per_kcal; % kJ/mol

elseif strcmp(Watermodel,'TIP4PEW')
    %% JC Ion Parameters in TIP4P water
    Param.Li.sigma = Beta*(0.808*nm_per_Ang)*(2^(5/6)); % nm
    Param.Li.epsilon = ES*Gamma*(0.1039884)*kJ_per_kcal; % kJ/mol

    Param.Na.sigma = Beta*(1.226*nm_per_Ang)*(2^(5/6)); % nm
    Param.Na.epsilon = ES*Gamma*(0.1684375)*kJ_per_kcal; % kJ/mol

    Param.K.sigma = Beta*(1.590*nm_per_Ang)*(2^(5/6)); % nm
    Param.K.epsilon = ES*Gamma*(0.2794651)*kJ_per_kcal; % kJ/mol

    Param.Rb.sigma = Beta*(1.709*nm_per_Ang)*(2^(5/6)); % nm
    Param.Rb.epsilon = ES*Gamma*(0.4331494)*kJ_per_kcal; % kJ/mol

    Param.Cs.sigma = Beta*(1.888*nm_per_Ang)*(2^(5/6)); % nm
    Param.Cs.epsilon = ES*Gamma*(0.3944318)*kJ_per_kcal; % kJ/mol

    Param.F.sigma = Beta*(2.538*nm_per_Ang)*(2^(5/6)); % nm
    Param.F.epsilon = ES*Gamma*(0.0015752)*kJ_per_kcal; % kJ/mol

    Param.Cl.sigma = Beta*(2.760*nm_per_Ang)*(2^(5/6)); % nm
    Param.Cl.epsilon = ES*Gamma*(0.0116615)*kJ_per_kcal; % kJ/mol

    Param.Br.sigma = Beta*(2.768*nm_per_Ang)*(2^(5/6)); % nm
    Param.Br.epsilon = ES*Gamma*(0.0303773)*kJ_per_kcal; % kJ/mol

    Param.I.sigma = Beta*(2.952*nm_per_Ang)*(2^(5/6)); % nm
    Param.I.epsilon = ES*Gamma*(0.0417082)*kJ_per_kcal; % kJ/mol
else
    error(['Unknown Water Model: "' Watermodel ...
        '". Please choose one of "SPC/E", "TIP3P", or "TIP4PEW".'])
end

OutputMet = Param.(Metal);
OutputHal = Param.(Halide);

if plotswitch
    %% Conversion Factors / Constants
    m_per_nm = 1e-9; % nm per m
    NA = 6.0221409e23; % Molecules per mole
    e_c = 1.60217662e-19; % Elementary charge in Coulombs
    J_per_kJ = 1000;
    epsilon_0 = (8.854187817620e-12)*J_per_kJ*m_per_nm/NA; % Vacuum Permittivity C^2 mol kJ^-1 nm^-1
    k_0 = 1/(4*pi*epsilon_0); % Coulomb constant in C^-2 nm kJ mol^-1

    % Use combining rules
    epsilon_MM = Param.(Metal).epsilon; % kJ/mol
    epsilon_HH = Param.(Halide).epsilon; % kJ/mol
    epsilon_Salt = sqrt(Param.(Metal).epsilon*Param.(Halide).epsilon); % kJ/mol
    
    sigma_MM = Param.(Metal).sigma;
    sigma_HH = Param.(Halide).sigma;
    sigma_Salt = (1/2)*(Param.(Metal).sigma + Param.(Halide).sigma);

    q_M = 1*e_c; % charge of lithium in C
    q_H = -1*e_c; % charge of fluoride in C

    %% Build PES for LiF
    r = 0:0.001:1; % r vector

    U_Salt = k_0.*q_M.*q_H./(r) + 4*epsilon_Salt.*( ( sigma_Salt./r ).^12 - ( sigma_Salt./r ).^6 );

    U_MetMet= k_0.*q_M.*q_M./(r) + 4*epsilon_MM.*( ( sigma_MM./r ).^12 - ( sigma_MM./r ).^6 );

    U_HalHal  = k_0.*q_H.*q_H./(r) + 4*epsilon_HH.*( ( sigma_HH./r ).^12 - ( sigma_HH./r ).^6 );

    % Options
    lw=2;
    fs=25;

    hold on
    h1 = plot(r,U_Salt,'Color','r','LineWidth',lw,'LineStyle','-');
    h2 = plot(r,U_MetMet,'Color','b','LineWidth',lw,'LineStyle',':');
    h3 = plot(r,U_HalHal,'Color','g','LineWidth',lw,'LineStyle','-.');

    title(['Plot of ' Metal Halide ' JC Potentials for ' Watermodel ...
        ' Water Model'],'Interpreter','latex','fontsize',fs)

    set(gca,'box','on','TickLabelInterpreter','latex');
    set(gca,'XMinorTick','on','YMinorTick','on','FontSize',fs);
    xlabel('r (nm)','fontsize',fs,'Interpreter','latex');
    ylabel('Pair Potential (kJ mol$^{-1}$)','fontsize',fs,'Interpreter','latex');

    ylim([-900 1000]); 

    % Blank line
    hline = refline([0 0]);
    hline.Color = 'k';
    hline.LineWidth = 1;
    hline.LineStyle = ':';
    legend([h1 h2 h3],{['JC - ' Metal Halide] ['JC - ' Metal Metal] ['JC - ' Halide Halide]})
end
end