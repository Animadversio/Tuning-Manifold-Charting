classdef torchLPIPS
   % Compute dissimilarity between 2 images of the same size. Using
   % PerceptualSimilarity Metric from Pytorch
   % 
   % Setup Method : 
   % 
   % Install package `lpips` in a python environment. 
   %    `pip install lpips`. 
   % Make sure that python environment has a torch version ~ 1.3 or 1.1 so that torch and matlab can work together.
   % If there is issue loading `torch` or `numpy` check the `pyenv` and `path` variable.
   % For more info see `torchBigGAN` and `torchStyleGAN2` 
   % 
   % Usage: 
   % 
   % To compute a scalar difference of images
   %    D = torchImDist("squeeze")
   %    dsim = D.distance(img1, img2); % 0.525
   % 
   % To compute a spatial map of difference
   %    D = torchImDist("squeeze",true)
   %    dmap = D.distance(img1, img2); % 256 by 256 scalar mask
   %    figure;montage({img1, img2, dmap},'size',[1,3])
   % Change the distance metric use 
   %    D = D.select_metric("alex",1);
   properties
       D
       metric % a variable preset to specify which metric to use
       spatial % boolean for whether analyze spatial distribution of difference or not.
       colorspace % Measure distance in RGB or Lab space. 
   end
   methods
   function G = torchLPIPS(metric, spatial, colorspace)
       if nargin == 0
           metric = "squeeze";
           spatial = 0;
           colorspace = "RGB"; %'Lab'
       elseif nargin == 1
           spatial = 0;
           colorspace = "RGB"; %'Lab'
       elseif nargin == 2
           colorspace = "RGB"; %'Lab'
       end
       py.importlib.import_module("lpips");
       py.importlib.import_module("torch");
       py.importlib.import_module("numpy");
       G.colorspace = colorspace;
       G.spatial = spatial;
       G.metric = metric;
       switch metric
           case "SSIM"
               G.D = py.lpips.DSSIM(pyargs("colorspace",colorspace));
           case "L2"
               G.D = py.lpips.L2(pyargs("colorspace",colorspace));
           case "alex"
               G.D = py.lpips.LPIPS(pyargs("net", "alex", "spatial", spatial));
           case "vgg"
               G.D = py.lpips.LPIPS(pyargs("net", "vgg", "spatial", spatial));
           case "squeeze"
               G.D = py.lpips.LPIPS(pyargs("net", "squeeze", "spatial", spatial));
           otherwise
               fprintf("Metric unrecognized, use squeeze-net linear weighted image metric.\n")
               G.D = py.lpips.LPIPS(pyargs("net", "squeeze", "spatial", spatial));
               G.metric = "squeeze";
       end
       G.D.requires_grad_(false);
       G.D.eval();
       py.torch.set_grad_enabled(false);
   end
   function G = select_metric(G, metric, spatial, colorspace)
       if nargin <= 2, spatial = 0; end
       if nargin <= 3, colorspace = "RGB"; end
       G.colorspace = colorspace;
       G.spatial = spatial;
       G.metric = metric;
       switch metric
           case "SSIM"
               G.D = py.lpips.DSSIM(pyargs("colorspace",colorspace));
           case "L2"
               G.D = py.lpips.L2(pyargs("colorspace",colorspace));
           case "alex"
               G.D = py.lpips.LPIPS(pyargs("net", "alex", "spatial", spatial));
           case "vgg"
               G.D = py.lpips.LPIPS(pyargs("net", "vgg", "spatial", spatial));
           case "squeeze"
               G.D = py.lpips.LPIPS(pyargs("net", "squeeze", "spatial", spatial));
           otherwise
               fprintf("Metric unrecognized, use squeeze-net linear weighted image metric.\n")
               G.D = py.lpips.LPIPS(pyargs("net", "squeeze", "spatial", spatial));
               G.metric = "squeeze";
       end
       G.D.requires_grad_(false);
       G.D.eval();
   end
   
   function dists = distance(G, im1, im2)
       if max(im1,[],'all')>1.2, im1 = single(im1) / 255.0; end
       if max(im2,[],'all')>1.2, im2 = single(im2) / 255.0; end
       % interface with generate integrated code, cmp to FC6GAN
       im1_tsr = py.torch.tensor(py.numpy.array(permute(im1,[4,3,1,2]))).float();
       im2_tsr = py.torch.tensor(py.numpy.array(permute(im2,[4,3,1,2]))).float();
       dists = G.D.forward(im1_tsr, im2_tsr).squeeze().cpu().detach().numpy().double;
   end
   
   function distMat = distmat(G, imgs, B)
       if nargin==2, B = 100;end
       if max(imgs,[],'all')>1.2, imgs = single(imgs) / 255.0; end
       % interface with generate integrated code, cmp to FC6GAN
       imgn = size(imgs,4);
       distMat = zeros(imgn,imgn);
       if G.metric == "SSIM" % SSIM doesn't support batch distance computation, so do it one by one.
           for i=1:imgn
               for j=i+1:imgn
                   distMat(i,j)=G.distance(imgs(:,:,:,i), imgs(:,:,:,j));
                   distMat(j,i)=distMat(i,j);
               end
           end
       else % net-lin model support batch distance computation, so do it batch in row.
           for i=1:imgn
               csr=1;
               dist_row = [];
               while csr <= imgn
               csr_end = min(imgn, csr+B-1);
               dists = G.distance(imgs(:,:,:,i), imgs(:,:,:,csr:csr_end));
               dist_row = cat(2, dist_row, dists);
               csr = csr_end+1;
               end
               distMat(i, :) = dist_row;
           end
       end
   end
   function distMat = distmat2(G, imgs1, imgs2, B)
       if nargin==3, B = 50;end
       if max(imgs1,[],'all')>1.2, imgs1 = single(imgs1) / 255.0; end
       if max(imgs2,[],'all')>1.2, imgs2 = single(imgs2) / 255.0; end
       % interface with generate integrated code, cmp to FC6GAN
       imgn = size(imgs1,4);
       imgm = size(imgs2,4);
       distMat = zeros(imgn,imgm);
       if G.metric == "SSIM"
           for i=1:imgn
               for j=i+1:imgn
                   distMat(i,j)=G.distance(imgs1(:,:,:,i), imgs2(:,:,:,j));
               end
           end
       else
       for csr_i = 1:imgn
           csr_j=1;
           while csr_j <= imgm
               csr_end = min(imgm, csr_j+B-1);
               dists = G.distance(imgs1(:,:,:,csr_i), imgs2(:,:,:,csr_j:csr_end));
               distMat(csr_i, csr_j:csr_end) = dists';
               csr_j = csr_end+1;
           end
       end
       end
   end
   end
end