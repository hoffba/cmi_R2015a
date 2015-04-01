function runTxformix(ifname,tfname,outfname)

bdir = uigetdir(pwd,'Select directory to process:');
cd(bdir)
elxdirs = dir(fullfile(bdir,'*_elxreg'));



[ifname,ipath] = uigetfile('*.mhd','Select image to transform');
[tfname,tpath] = uigetfile('*.txt','Select TransformParameter File');
[ofname,opath] = uiputfile('*.mhd','Saved transformed image as:');

outdir = fullfile(opath,'tempelx');
if ~exist(outdir,'dir')
    mkdir(outdir);
end
infn = fullfile(outdir,'temp-in.mhd');
img = ImageClass;
img.loadImg(0,{fullfile(ipath,ifname)});
img.saveImg(1,infn,1);

% Window Name:
namestr = ['Elastix Registration: ',ifname,' --> ',tfname];
% Call to Transformix:
tfxstr = ['/opt/elastix/bin/transformix',...
          ' -out ',outdir,...
          ' -in ',infn,...
          ' -tp ',fullfile(tpath,tfname),'; '];
% Cleanup string (removes temporary image files):
custr = '';
% custr = ['find ',outdir,' -name "elxtemp-*" -exec rm -f {} \; '];

stat = system(['xterm -geometry 170x50 -T "',namestr,'"',...
               ' -e ''',tfxstr,custr,';''']);
           
img.loadImg(0,fullfile(outdir,'result.mhd'));
img.saveImg(1,fullfile(opath,ofname),1);