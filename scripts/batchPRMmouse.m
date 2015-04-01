% CMI script
function batchPRMmouse(cmiObj)
%   Input: cmiObj = CMIclass object containing current settings

bdir = '/mnt/cmi/projects/CT_Lung/preclinical/LungCT-PRM/CTdata/Normals_C6/FLDs';
d = ' ';
istr = {'0','1','3','4','5','6','7','8','9'};
out = cell(2*length(istr) + 1,7);
out(1,:) = {'Exp','Ins','Emph','White','fSAD','Normal','Fibrosis'};

nids = length(istr);
for i = 1:nids % animal IDs
    inm = ['Norm',istr{i}];
    disp(inm)
    fexp = struct2cell(dir(fullfile(bdir,inm,'*Exp.fld')));
    if size(fexp,2)==2
        for j = 1:size(fexp,2)
            d(j) = fexp{1,j}(8);
        end
        for j = 1:2 % each animal has 2 time points
            nout = 2*(i-1) + j + 1;
            fins = dir(fullfile(bdir,inm,[inm,'*_d',d(j),'*_Ins*_R.fld']));
            if ~isempty(fins)
                fins = fins(end).name;
                % Load images
                cmiObj.loadImg(0,{fullfile(bdir,inm,fexp{1,j}),...
                                  fullfile(bdir,inm,fins)});
                out(nout,1:2) = {fexp{1,j} , fins};

                % Generate and save mask
                cmiObj.img.mask.merge('Replace',...
                    segLungMouse(cmiObj.img.mat(:,:,:,1)));
                cmiObj.img.saveMask(fullfile(bdir,inm,[fexp{1,j}(1:end-4),'_VOI.fld']));
                
                % Activate the PRM
                [~,v] = cmiObj.img.calcPRM(2);
                out(nout,3:end) = num2cell(v);
            end
        end
    end
end

assignin('base','results',out)
    

