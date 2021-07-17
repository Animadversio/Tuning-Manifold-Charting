function Kent_sphere_plot(ax, param_vec)
% theta, phi, psi, kappa, beta, A,
ang_step = 5;
L = (-90:ang_step:90) / 180 * pi;
T = (-90:ang_step:90) / 180 * pi;
[L, T] = meshgrid(L, T);
theta=param_vec(1); phi=param_vec(2); psi=param_vec(3); 
kappa=param_vec(4); beta=param_vec(5); A=param_vec(6);
fval = KentFunc(theta, phi, psi, kappa, beta, A, L(:), T(:));
% Visualize Kent distribution of 
fval = reshape(fval, size(L));
[X,Y,Z] = sph2cart(L,T,1);
if isempty(ax)
    figure;clf;
    surf(X,Y,Z,fval,'FaceAlpha',0.8,'EdgeColor','interp');
else
    surf(ax,X,Y,Z,fval,'FaceAlpha',0.8,'EdgeColor','interp');
end
colorbar; shading faceted; % interp flat
daspect([1 1 1]); axis tight; view([90,0]);
title('3D Plot of f on the sphere')
end