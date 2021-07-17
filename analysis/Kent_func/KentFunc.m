function fval = KentFunc(theta, phi, psi, kappa, beta, A, theta_z, phi_z)
% Assume theta_z, phi_z are column vectors ([0,2 pi])
% theta, phi, psi are scaler angular parameter ([0,2 pi])
% theta_z, phi_z
Z = [cos(theta_z).* cos(phi_z), sin(theta_z) .* cos(phi_z), sin(phi_z)]';
coord = SO3(theta, phi, psi);
mu1  = coord(:,1);
mu23 = coord(:,2:3);
fval = A * exp(kappa * mu1' * Z + beta * [1, -1] * (mu23' * Z).^2)';
end

function coord = SO3(theta, phi, psi)
    orig = [cos(theta)*cos(phi) sin(theta)*cos(phi)  sin(phi);
            -sin(theta)         cos(theta)                  0;
            cos(theta)*sin(phi) sin(theta)*sin(phi) -cos(phi)]';
    Rot23 = [1         0        0;
             0   cos(psi) sin(psi);
             0  -sin(psi) cos(psi)];
    coord = orig * Rot23;
end