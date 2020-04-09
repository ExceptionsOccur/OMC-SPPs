function [epsilon] = Cu_epsilon(lambda)

% Lorentz-Drude (LD) model parameters
omega_p = 10.83;
f_0 = 0.575;
gamma_0 = 0.030;

f_1 = 0.061;
gamma_1 = 0.378;
omega_1 = 0.291;

f_2 = 0.104;
gamma_2 = 1.056;
omega_2 = 2.957;

f_3 = 0.723;
gamma_3 = 3.213;
omega_3 = 5.300;

f_4 = 0.638;
gamma_4 = 4.305;
omega_4 = 11.18;

f_5 = 5.646;
gamma_5 = 2.419;
omega_5 = 20.29;

ev = 4.13566733e-1 .* 2.99792458 ./ (lambda .* 1e6);

Omega_p = sqrt(f_0) * omega_p;

epsilon = 1 - Omega_p^2 ./ (ev .* (ev + 1i * gamma_0)) + f_1 * omega_p^2 ./ ((omega_1^2 - ev.^2) - 1i .* ev .* gamma_1) + ...
    f_2 .* omega_p^2 ./ ((omega_2^2 - ev.^2) - 1i .* ev .* gamma_2) + f_3 .* omega_p^2 ./ ((omega_3^2 - ev.^2) - 1i .* ev .* gamma_3) + ...
    f_4 .* omega_p^2 ./ ((omega_4^2 - ev.^2) - 1i .* ev .* gamma_4) + f_5 .* omega_p^2 ./ ((omega_5^2 - ev.^2) - 1i .* ev .* gamma_5);
end

