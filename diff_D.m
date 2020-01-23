
% clear
% clc

% basic parameter

lambda = linspace(400,1600,2000) * 1E-9;
theta = pi / 3;
n_core = 1.44;
n_clad = 1;
p_core = n_core^2;
p_clad = n_clad^2;
metal_t = 75E-9;
n_glass = 1.3;
p_glass = n_glass^2;
L = 50E-6;
D = 2E-6;


% metal

% permittivty
lambda_c = 8.9342E-6;
lambda_p = 1.6826E-7;
omega_p = 1.3388E16;
gamma = 7.07592E13;
eps_inf = 3.36174;

  
for d = linspace(2,10,6)*1E-6
    
%     n_core = d;
%     p_core = n_core^2;
    D = d;
    % layer 1: core

    q1 = sqrt(p_core - n_core^2 * sin(theta)) / p_core;
    % q1 = sqrt(1 / p_core) * cos(theta);
    
%     n_glass = d;
%     p_glass = n_glass^2;
    
    for i = 1:length(lambda)

    %     p_metal = 1 - (lambda(i)^2 * lambda_c) / (lambda_p^2 .* (lambda_c + 1i * lambda(i)));
        omega = 2 * pi * (2.99 * 1E8 / lambda(i));
        p_metal = eps_inf - omega_p^2 / (omega^2 + 1i * gamma * omega);

        % refractive
        n_metal = sqrt(p_metal);

        % layer 2: metal

        theta2 = asin(sin(theta) / (n_core / n_metal));
        q2 = sqrt(p_metal - n_core^2 * sin(theta)) / p_metal;
    %     q2 = sqrt(1 / p_metal) * cos(theta2);
        beta2 = 2 * pi * metal_t / lambda(i) * sqrt(p_metal - n_core^2 * sin(theta)^2);

        % layer 3: glass
        theta3 = asin(sin(theta2) / (n_metal / n_glass));
        q3 = sqrt(p_glass - n_core^2 * sin(theta)) / p_glass;
    %     q3 = sqrt(1 / p_glass) * cos(theta3);

        % transfer matrix

        m11 = cos(beta2);
        m12 = -1i * sin(beta2) / q2;
        m21 = -1i * q2 * sin(beta2);
        m22 = m11;

        % power

        p_trans = n_core^2 * sin(theta) * cos(theta) / (1 - n_core^2 * cos(theta)^2)^2;

        % reflection ratio

        r_p = ((m11 + m12 * q3) * q1 - (m21 + m22 * q3)) / ((m11 + m12 * q3) * q1 + (m21 + m22 * q3));
        R_p = abs(r_p)^2;
        Nef = L / (D * tan(theta));

        p(i) = R_p^Nef;
    end
    
    plot(lambda*1E9,p);
    hold on;
end
xlabel('Wavelength(nm)');
ylabel('Reflection');
axis([min(lambda)*1E9 max(lambda)*1E9 0 1]);

