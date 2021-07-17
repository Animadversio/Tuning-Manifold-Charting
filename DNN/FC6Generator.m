% FC6Generator under native matlab deeplearning framework 
classdef FC6Generator
   properties
      BGRMean = reshape([104.0, 117.0, 123.0],[1,1,3,1]);
      DeconvNet
      LinNet
   end
   methods
      function G = FC6Generator(loadpath)
         %D:\Github\Monkey_Visual_Experiment_Data_Processing
         if nargin == 0
             loadpath = "matlabGANfc6.mat";
         end
         %D:\Github\Monkey_Visual_Experiment_Data_Processing\DNN\matlabGANfc6.mat
         %E:\Github_Projects\Monkey_Visual_Experiment_Data_Processing\DNN\matlabGANfc6.mat
         data = load(loadpath,'DeconvNet','LinNet','BGRMean');
         G.DeconvNet = data.DeconvNet;
         G.LinNet = data.LinNet;
      end
      function acts = activations(G, codes, layername)
        if size(codes,2)==4096 
            codes = codes';
        end
        assert(size(codes,1)==4096, "Code shape error")%
        Z = dlarray(reshape(codes,1,1,4096,size(codes,2)), 'SSCB');
        hiddenout = G.LinNet.predict(Z);
        hiddenout_r = dlarray(hiddenout.reshape(4,4,256,[]).permute([2,1,3,4]),"SSCB"); % Note to flip your hiddenoutput
        acts = forward(G.DeconvNet, hiddenout_r, 'Outputs', layername);
        acts = extractdata(acts);
      end
      function acts = dlforward(G, codes)%, layername)
        % assume codes are dlarray with 4096 by n shape
        if size(codes,2)==4096 
            codes = codes';
        end
        assert(size(codes,1)==4096, "Code shape error")%
        Z = dlarray(reshape(codes,1,1,4096,size(codes,2)), 'SSCB');
        hiddenout = G.LinNet.predict(Z);
        hiddenout_r = dlarray(hiddenout.reshape(4,4,256,[]).permute([2,1,3,4]),"SSCB"); % Note to flip your hiddenoutput
%         acts = forward(G.DeconvNet, hiddenout_r, 'Outputs', layername);
        acts = G.DeconvNet.predict(hiddenout_r);
        acts = clip(acts + G.BGRMean, 0, 255); % extractdata(acts);
      end
      function imgs = visualize(G, codes)
        if size(codes,2)==4096 
            codes = codes'; 
        end 
        assert(size(codes,1)==4096, "Code shape error")% first dimension should have size 4096, second the batch
        Z = dlarray(reshape(codes,1,1,4096,size(codes,2)), 'SSCB');
        hiddenout = G.LinNet.predict(Z);
        hiddenout_r = dlarray(hiddenout.reshape(4,4,256,[]).permute([2,1,3,4]),"SSCB"); % Note to flip your hiddenoutput
        out = G.DeconvNet.predict(hiddenout_r);
        imgs = extractdata(out(:,:,:,:));
        imgs = uint8(min(max(imgs + G.BGRMean, 0), 255));
        imgs = imgs(:,:,[3,2,1],:);
      end
      function frame_cell = visualize_movie(G, latent_col)
       batchsize = 10; movieN = numel(latent_col); 
       frameN = cellfun(@(code)size(code, 1), latent_col);
       samplen = sum(frameN); 
       latent_arr = cat(1, latent_col{:});
       % Generate latent code using the current space
       matimgs = G.visualize(latent_arr);
       % sort the frames into cells corresponding to each movie. 
       frame_cell = arrayfun(@(iMv) matimgs(:,:,:,sum(frameN(1:iMv-1))+1:sum(frameN(1:iMv))), 1:movieN, "Uni", 0); 
     end
   end
end
