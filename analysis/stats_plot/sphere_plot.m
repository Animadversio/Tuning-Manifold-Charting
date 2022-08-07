function sphere_plot(ax,Theta,Phi,fval)
if max(max(Theta)) > pi || min(min(Theta)) < -pi
    Theta = Theta /180 *pi;
    Phi = Phi /180 *pi;
end
X = cos(Theta).* cos(Phi);
Y = sin(Theta) .* cos(Phi);
Z = sin(Phi);
% [X,Y,Z] = sph2cart(Theta,Phi,1);
if isempty(ax)
    figure;clf;
    surf(X,Y,Z,fval,'FaceAlpha',0.8,'FaceColor','interp');%,'EdgeColor','interp');%,'interp','EdgeColor','interp');
else
    surf(ax,X,Y,Z,fval,'FaceAlpha',0.8,'FaceColor','interp','EdgeColor','interp')%;,'EdgeColor','interp');%,[0.5,0.5,0.5]'interp','EdgeColor','interp');,[0.5,0.5,0.5]
end
colorbar; %shading faceted; % interp flat
daspect([1 1 1]); axis tight; view([90,0]);
title('3D Plot of f on the sphere')
end