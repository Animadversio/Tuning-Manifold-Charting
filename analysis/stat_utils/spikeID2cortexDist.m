function [cortexDistmat, V1msk, V4msk, ITmsk] = spikeID2cortexDist(spikeID, unitID)
Extract_Channel_Organization;
IT_XY_Area = [IT_chan_XY, 3*ones(32,1)];
V1_XY_Area = [V1_chan_XY, 1*ones(16,1)];
V4_XY_Area = [V4_chan_XY, 2*ones(16,1)];
Ch2XY_Area = [IT_XY_Area;
              V1_XY_Area;
              V4_XY_Area];
XY_Area_arr = Ch2XY_Area(spikeID,:);
XY_Area_Unit_arr = [XY_Area_arr,unitID];
cortexDistmat = squareform(pdist(XY_Area_Unit_arr(:,1:2)));
sameArea = (XY_Area_arr(:,3)==XY_Area_arr(:,3)');
cortexDistmat(~sameArea) = nan;
ITmsk =  spikeID < 33;
V1msk = (spikeID < 49) & (spikeID > 32);
V4msk =  spikeID > 48;
end