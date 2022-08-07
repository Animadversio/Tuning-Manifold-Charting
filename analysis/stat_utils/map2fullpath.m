function [imgfps, mapper] = map2fullpath(picnm_arr, imgdir)
% map the imgnm without suffix to the full path, for all images in a folder.
imgnm_wsfx = deblank([string(ls(imgdir+"\*.png"));...
                      string(ls(imgdir+"\*.jpg"));...
                      string(ls(imgdir+"\*.jpeg"));...
                      string(ls(imgdir+"\*.bmp"))]);
[~,imgnms,sfxs] = arrayfun(@fileparts, imgnm_wsfx); % older version matlab
mapper = containers.Map(imgnms,fullfile(imgdir, imgnm_wsfx));
imgfps = cellfun(@(I)mapper(I),picnm_arr,'uni',0); % image full path of tht input picnms 
end