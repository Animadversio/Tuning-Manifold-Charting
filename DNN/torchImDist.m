classdef torchImDist
   % Usage: 
   % Compute dissimilarity between 2 images of the same size. Using
   % PerceptualSimilarity Metric from Pytorch
   % 
   % To compute a scalar difference of images
   %    D = torchImDist("squeeze")
   %    dsim = D.distance(img1, img2); % 0.525
   % 
   % To compute a spatial map of difference
   %    D = torchImDist("squeeze",1)
   %    dmap = D.distance(img1, img2); % 256 by 256 scalar mask
   %    figure;montage({img1, img2, dmap},'size',[1,3])
   % Change the distance metric use 
   %    D = D.select_metric("alex",1);
   properties
       D
       metric % a variable preset to specify which metric to use
       spatial % boolean for whether analyze spatial distribution of difference or not.
       net
   end
   methods
   function G = torchImDist(metric, spatial)
       if nargin == 0
           metric = "squeeze";
           spatial = 0;
       elseif nargin == 1
           spatial = 0;
       end
       switch getenv('COMPUTERNAME')
           case 'DESKTOP-9DDE2RH' % Office 3 Binxu's 
            repodir = "D:\Github\PerceptualSimilarity";
           case 'DESKTOP-MENSD6S' % Binxu's home work station
            repodir = "E:\Github_Projects\PerceptualSimilarity";
           case 'PONCELAB-ML2A' % MLa machine 
            repodir = "C:\Users\Poncelab-ML2a\Documents\Python\PerceptualSimilarity";
           case 'PONCELAB-ML2B' % MLb machine 
            repodir = "C:\Users\Ponce lab\Documents\Python\PerceptualSimilarity";
           otherwise
            repodir = "C:\Users\Poncelab-ML2a\Documents\Python\PerceptualSimilarity";
        end
       % install the torch 1.3.x and the biggan package like below.
        py.importlib.import_module("sys");
        syspath = py.sys.path(); % add the official stylegan2 repo. 
        syspath.append(repodir);
        py.importlib.import_module("models");
        py.importlib.import_module("torch");
        py.importlib.import_module("numpy");
        G.spatial = spatial;
        G.metric = metric;
       switch metric
           case "SSIM"
               G.D = py.models.PerceptualLoss(pyargs("model", "SSIM", "spatial", spatial,"colorspace","RGB"));
           case "L2"
               G.D = py.models.PerceptualLoss(pyargs("model", "L2", "spatial", spatial,"colorspace","RGB"));
           case "alex"
               G.D = py.models.PerceptualLoss(pyargs("model", "net-lin", "net", "alex", "spatial", spatial));
           case "vgg"
               G.D = py.models.PerceptualLoss(pyargs("model", "net-lin", "net", "vgg", "spatial", spatial));
           case "squeeze"
               G.D = py.models.PerceptualLoss(pyargs("model", "net-lin", "net", "squeeze", "spatial", spatial));
           otherwise
               fprintf("Metric unrecognized, use squeeze-net linear weighted image metric.\n")
               G.D = py.models.PerceptualLoss(pyargs("model", "net-lin", "net", "squeeze", "spatial", spatial));
               G.metric = "squeeze";
       end
       G.D.requires_grad_(false);
       G.D.eval(); % G.BGAN.to('cuda');
       py.torch.set_grad_enabled(false);
   end
   function G = select_metric(G, metric, spatial)
       if nargin <= 2
           spatial = 0;
       end
       G.spatial = spatial;
       G.metric = metric;
       switch metric
           case "SSIM"
               G.D = py.models.PerceptualLoss(pyargs("model", "SSIM", "spatial", spatial,"colorspace","RGB"));
           case "L2"
               G.D = py.models.PerceptualLoss(pyargs("model", "L2", "spatial", spatial,"colorspace","RGB"));
           case "alex"
               G.D = py.models.PerceptualLoss(pyargs("model", "net-lin", "net", "alex", "spatial", spatial));
           case "vgg"
               G.D = py.models.PerceptualLoss(pyargs("model", "net-lin", "net", "vgg", "spatial", spatial));
           case "squeeze"
               G.D = py.models.PerceptualLoss(pyargs("model", "net-lin", "net", "squeeze", "spatial", spatial));
           otherwise
               fprintf("Metric unrecognized, use squeeze-net linear weighted image metric.\n")
               G.D = py.models.PerceptualLoss(pyargs("model", "net-lin", "net", "squeeze"));
               G.metric = "squeeze";
               G.spatial = 0;
       end
   end
   function G = load_net(G)
       % this function is added to make a uniform interface for AlexNet
       % based distance metric. If not loaded D is more light weighted.
   G.net = alexnet;
   end
   function dists = distance(G, im1, im2)
       if max(im1,[],'all')>1.2, im1 = single(im1) / 255.0; end
       if max(im2,[],'all')>1.2, im2 = single(im2) / 255.0; end
       % interface with generate integrated code, cmp to FC6GAN
       im1_tsr = py.torch.tensor(py.numpy.array(permute(im1,[4,3,1,2])));
       im2_tsr = py.torch.tensor(py.numpy.array(permute(im2,[4,3,1,2])));
       dists = G.D.forward(im1_tsr, im2_tsr).squeeze().cpu().detach().numpy().double;
   end
   
   function dists = distance_op(G, im1, im2)
       if max(im1,[],'all')>1.2, im1 = single(im1) / 255.0; end
       if max(im2,[],'all')>1.2, im2 = single(im2) / 255.0; end
       % interface with generate integrated code, cmp to FC6GAN
       im1_tsr = py.torch.tensor(py.numpy.array(permute(im1,[4,3,1,2])));
       im2_tsr = py.torch.tensor(py.numpy.array(permute(im2,[4,3,1,2])));
       dists = G.D.forward(im1_tsr, im2_tsr, pyargs("distmat",true)).squeeze().cpu().detach().numpy().double;
   end
   
   function distMat = distmat(G, imgs, B)
       if nargin==2, B = 100;end
       if max(imgs,[],'all')>1.2, imgs = single(imgs) / 255.0; end
       % interface with generate integrated code, cmp to FC6GAN
       imgn = size(imgs,4);
       distMat = zeros(imgn,imgn);
       if G.metric == "SSIM"
       for i=1:imgn
           for j=i+1:imgn
               distMat(i,j)=G.distance(imgs(:,:,:,i), imgs(:,:,:,j));
               distMat(j,i)=distMat(i,j);
           end
       end
       elseif G.metric == "L2"
           L2dist = pdist(reshape(imgs,[],size(imgs,4))');
           distMat = squareform(L2dist);
       elseif G.metric == "FC6_corr"
           acts = G.net.activations(255*imgs,"fc6"); 
           FC6dist_corr = squareform(pdist(squeeze(acts)','correlation'));
           distMat = FC6dist_corr; 
       elseif G.metric == "FC6_L2"
           acts = G.net.activations(255*imgs,"fc6"); 
           FC6dist = squareform(pdist(squeeze(acts)','euclidean'));
           distMat = FC6dist; 
       else
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
   
   function distMat = distmat_B(G, imgs, B)
       if nargin==2, B = 70;end
       if max(imgs,[],'all')>1.2, imgs = single(imgs) / 255.0; end
       % interface with generate integrated code, cmp to FC6GAN
       imgn = size(imgs,4);
       distMat = zeros(imgn,imgn);
       csr_i = 1;
       while csr_i <= imgn
           csr_eni = min(imgn, csr_i+B-1);
           csr_j=csr_i;
           while csr_j <= imgn
               csr_end = min(imgn, csr_j+B-1);
               dists = G.distance_op(imgs(:,:,:,csr_i:csr_eni), imgs(:,:,:,csr_j:csr_end));
               distMat(csr_i:csr_eni, csr_j:csr_end) = dists';
               distMat(csr_j:csr_end, csr_i:csr_eni) = dists;
               csr_j = csr_end+1;
           end
           csr_i = csr_eni+1;
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
       csr_i = 1;
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