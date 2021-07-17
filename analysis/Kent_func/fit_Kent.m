function [Parameter, gof] = fit_Kent(funcval)
ft = fittype( @(theta, phi, psi, kappa, beta, A, x, y) KentFunc(theta, phi, psi, kappa, beta, A, x, y), ...
    'independent', {'x', 'y'},'dependent',{'z'});
    %'problem', {'theta', 'phi', 'psi', 'kappa', 'beta', 'A'}, ...
% assume funcval is a 11*11 grid
ang_step = 18;
theta_arr = [-90:ang_step:90] / 180 * pi;
phi_arr = [-90:ang_step:90] / 180 * pi;
[phi_grid, theta_grid] = meshgrid(phi_arr, theta_arr); % return a grid of 
% degree of PC2 comes in the first index i, theta deg = (i-6) * 18
% degree of PC3 comes in the second index j, phi deg = (j-6) * 18
[Parameter, gof] = fit( [theta_grid(:), phi_grid(:)], [funcval(:)], ft, ...
                'StartPoint', [0, 0, pi/2, 0.1, 0.1, 0.1], ...
                'Lower', [-pi, -pi/2,  0, -Inf,   0,   0], ...
                'Upper', [ pi,  pi/2,  pi,  Inf, Inf, Inf]);%, ...
gof.coefname = string(coeffnames(Parameter)');
gof.coef = coeffvalues(Parameter);
gof.confint = confint(Parameter);
% "theta"    "phi"    "psi"    "kappa"    "beta"    "A" sequence
return