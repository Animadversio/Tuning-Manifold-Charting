function [laplsEng_sp, TVEng_sp, dirEng_sp, laplsMap_sp, gradNormMap_sp, detJac_map] = sph_map_energy(actmap)
[PHI, THETA] = meshgrid(-90:18:90,-90:18:90); 
kerd2 = [[0,-1,0];[-1,4,-1];[0,-1,0]]/4;
% detJac = abs(cosd(PHI));
Ginv11 = 1./cosd(PHI).^2;Ginv11(:,1)=0;Ginv11(:,end)=0; % These components explode, set to 0. 
Ginv22 = ones(size(PHI));
% Sample points for the gradients are far away.
detJac_map = cosd([mean(PHI(:,[1,2]),2),PHI(:,2:end-1),mean(PHI(:,[end-1,end]),2)]);
actmap_sp = actmap;
actmap_sp(:,1)=mean(actmap(:,1));
actmap_sp(:,end)=mean(actmap(:,end));
laplsMap_sp = conv2(padarray(actmap_sp,[1,1],'replicate'), kerd2,'valid');%del2(actmap_sp);%
[actGx_sp,actGy_sp] = gradient(actmap_sp);
gradNormMap_sp = actGx_sp.^2 .* Ginv22 + actGy_sp.^2 .*Ginv11;
dirEng_sp = sum(gradNormMap_sp.*detJac_map,'all');
TVEng_sp = sum(sqrt(gradNormMap_sp).*detJac_map,'all');
laplsEng_sp = sum(abs(laplsMap_sp).*detJac_map,'all');
end