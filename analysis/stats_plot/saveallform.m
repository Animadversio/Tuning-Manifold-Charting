function saveallform(figdirs,fignm,h,sfxlist)
% Save a (current) figure with all suffices in a figdir.
% signature:
%   saveallform(figdir,fignm,h,sfxlist)
% 
if nargin <=2, h=gcf; else, h = figure(h); end
if nargin <=3, sfxlist = ["fig","png","pdf"]; end
for figdir = reshape(figdirs,1,[])
for sfx = sfxlist
if strcmp(sfx, "fig")
   savefig(h,fullfile(figdir,fignm+"."+sfx))
elseif strcmp(sfx, "pdf")
%    set(h,'Units','Inches');
%    pos = get(h,'Position');
%    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%    print(h,fullfile(figdir,fignm+"."+sfx),'-dpdf','-bestfit')
   % function available in Matlab 2021, if older version use the commented part above. 
   exportgraphics(h,fullfile(figdir,fignm+"."+sfx),'ContentType','vector') 
else
   saveas(h,fullfile(figdir,fignm+"."+sfx))
end
end
end
end