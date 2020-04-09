function [epsilon] = Au_epsilon(lambda)

% Lorentz-Drude (LD) model parameters
omega_p = 9.03;
f_0 = 0.760;
gamma_0 = 0.053;

f_1 = 0.024;
gamma_1 = 0.241;
omega_1 = 0.415;

f_2 = 0.010;
gamma_2 = 0.345;
omega_2 = 0.830;

f_3 = 0.071;
gamma_3 = 0.870;
omega_3 = 2.969;

f_4 = 0.601;
gamma_4 = 2.494;
omega_4 = 4.304;

f_5 = 4.384;
gamma_5 = 2.214;
omega_5 = 13.32;

ev = 4.13566733e-1 .* 2.99792458 ./ (lambda .* 1e6);

Omega_p = sqrt(f_0) * omega_p;
epsilon = 1 - Omega_p^2 ./ (ev .* (ev + 1i * gamma_0)) + f_1 * omega_p^2 ./ ((omega_1^2 - ev.^2) - 1i .* ev .* gamma_1) + ...
    f_2 .* omega_p^2 ./ ((omega_2^2 - ev.^2) - 1i .* ev .* gamma_2) + f_3 .* omega_p^2 ./ ((omega_3^2 - ev.^2) - 1i .* ev .* gamma_3) + ...
    f_4 .* omega_p^2 ./ ((omega_4^2 - ev.^2) - 1i .* ev .* gamma_4) + f_5 .* omega_p^2 ./ ((omega_5^2 - ev.^2) - 1i .* ev .* gamma_5);
end

