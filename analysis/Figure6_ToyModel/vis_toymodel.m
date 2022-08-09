%% Manifold Paper Figure 6
%  2d toy model demonstrating the effect of tuning width and tuned dims on different stats.  

xmax = 20; ymax = 20;
Nx = 101; Ny = 101;
Cx = 15; Cy = 15;
FLATSTD = 1000;
xtick = linspace(-xmax,xmax,Nx);
ytick = linspace(-ymax,ymax,Ny);
gmap = gauss2d(xtick,ytick,0,0,5,5);

%% Demo of the evolution on the dummy function. 
centxy = zeros(1,2);
traj = [centxy];
for i = 1:10
xxyy = centxy + randn(5,2);
fvals = gauss2d_vec(xxyy,Cx,Cy,5,5);
[vsort, idxsort] = sort(fvals);
selxxyy = xxyy(idxsort(end-1:end),:);
centxy = mean(selxxyy,1);
traj = [traj;centxy];
end
%%
figure("pos", [680   322   640   650]);
tiledlayout(2,2,'padd','none','Tilesp','compact')
nexttile;
imagesc(xtick,ytick,gauss2d(xtick,ytick,Cx,Cy,5,5));
xticks([-20:10:20]);yticks([-20:10:20])
title("\sigma=5, D=2")
axis image; colormap hot;%colorbar;%copper
line([-20,20],[5,-5],'color','blue','linewidth',2)
hold on;traj = guided_randwalk([0,0],13,Cx,Cy,5,5);
plot(traj(:,1),traj(:,2),'-o','color','green','linewidth',1.5)
nexttile;
imagesc(xtick,ytick,gauss2d(xtick,ytick,Cx,Cy,10,10));
xticks([-20:10:20]);yticks([-20:10:20])
title("\sigma=10, D=2")
axis image; colormap hot;%colorbar;%copper
line([-20,20],[5,-5],'color','blue','linewidth',2)
hold on;traj = guided_randwalk([0,0],13,Cx,Cy,10,10);
plot(traj(:,1),traj(:,2),'-o','color','green','linewidth',1.5)
nexttile;
imagesc(xtick,ytick,gauss2d(xtick,ytick,Cx,Cy,5,FLATSTD));
xticks([-20:10:20]);yticks([-20:10:20])
title("\sigma=5, D=1")
axis image; colormap hot;%colorbar;%copper
line([-20,20],[5,-5],'color','blue','linewidth',2)
hold on;traj = guided_randwalk([0,0],13,Cx,Cy,5,FLATSTD);
plot(traj(:,1),traj(:,2),'-o','color','green','linewidth',1.5)
nexttile;
imagesc(xtick,ytick,gauss2d(xtick,ytick,Cx,Cy,10,FLATSTD));
xticks([-20:10:20]);yticks([-20:10:20])
title("\sigma=10, D=1")
axis image; colormap hot;%colorbar;%copper
line([-20,20],[5,-5],'color','blue','linewidth',2)
hold on;traj = guided_randwalk([0,0],13,Cx,Cy,10,FLATSTD);
plot(traj(:,1),traj(:,2),'-o','color','green','linewidth',1.5)
saveallform("E:\OneDrive - Harvard University\Manuscript_Manifold\Figure6Toy",...
            "toymodel_2dexample")

function G = gauss2d(xticks,yticks,mux,muy,sigx,sigy)
% 2d gaussian function
%  xticks, yticks are 1d vectors of the ticks 
[xx,yy]=meshgrid(xticks,yticks);
G=exp(-(((xx-mux)/sigx).^2+((yy-muy)/sigy).^2)/2);
end

function G = gauss2d_vec(xxyy,mux,muy,sigx,sigy)
% vectorized 2d Gaussian function. 
%  xxyy is N-by-2 array
xx = xxyy(:,1); yy = xxyy(:,2);
G=exp(-(((xx-mux)/sigx).^2+((yy-muy)/sigy).^2)/2);
end

function traj = guided_randwalk(centxy,steps,Cx,Cy,sigx,sigy)
% Random walk guided by the Gaussian tuning function 
%   specified by the parameter `Cx,Cy,sigx,sigy` 
%   Cx, Cy is the center coordinate
%   sigx, sigy is the std along x and y
traj = [centxy];
for i = 1:steps
xxyy = centxy + randn(5,2);
fvals = gauss2d_vec(xxyy,Cx,Cy,sigx,sigy);
[vsort, idxsort] = sort(fvals);
selxxyy = xxyy(idxsort(end-1:end),:);
centxy = mean(selxxyy,1);
traj = [traj;centxy];
end
end