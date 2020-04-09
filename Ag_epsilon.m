function [epsilon] = Ag_epsilon(lambda)

% Lorentz-Drude (LD) model parameters
omega_p = 9.01;
f_0 = 0.845;
gamma_0 = 0.048;

f_1 = 0.065;
gamma_1 = 3.886;
omega_1 = 0.816;

f_2 = 0.124;
gamma_2 = 0.452;
omega_2 = 4.481;

f_3 = 0.011;
gamma_3 = 0.065;
omega_3 = 8.185;

f_4 = 0.840;
gamma_4 = 0.916;
omega_4 = 9.083;

f_5 = 5.646;
gamma_5 = 2.419;
omega_5 = 20.29;

ev = 4.13566733e-1 .* 2.99792458 ./ (lambda .* 1e6);

Omega_p = sqrt(f_0) * omega_p;

epsilon = 1 - Omega_p^2 ./ (ev .* (ev + 1i * gamma_0)) + f_1 * omega_p^2 ./ ((omega_1^2 - ev.^2) - 1i .* ev .* gamma_1) + ...
    f_2 .* omega_p^2 ./ ((omega_2^2 - ev.^2) - 1i .* ev .* gamma_2) + f_3 .* omega_p^2 ./ ((omega_3^2 - ev.^2) - 1i .* ev .* gamma_3) + ...
    f_4 .* omega_p^2 ./ ((omega_4^2 - ev.^2) - 1i .* ev .* gamma_4) + f_5 .* omega_p^2 ./ ((omega_5^2 - ev.^2) - 1i .* ev .* gamma_5);
end

