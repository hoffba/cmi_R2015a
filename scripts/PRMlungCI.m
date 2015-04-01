% CMI script
function PRMlungCI(cmiObj)
%   Input: cmiObj = CMIclass object containing current settings

bdir = '/mnt/cmi/projects/CT_Lung/preclinical/LungCT-PRM/CTdata/Normals_C6/FLDs';
d = ' ';
istr = {'0','1','3','4','5','6','7','8','9'};
out = cell(2*length(istr) + 1,9);
out(1,:) = {'Exp','Ins','ExpMin','Exp-t','InsMin','InsMax','Ins-t','D','AUC'};

nids = length(istr);
for i = 1:nids % animal IDs
    inm = ['Norm',istr{i}];
    disp(inm)
    fexp = struct2cell(dir(fullfile(bdir,inm,'*Exp*NOVA.fld')));
    if size(fexp,2)==2
        for j = 1:size(fexp,2)
            d(j) = fexp{1,j}(8);
        end
        for j = 1:2 % each animal has 2 time points
            nout = 2*(i-1) + j + 1;
            fins = dir(fullfile(bdir,inm,[inm,'*_d',d(j),'*_Ins*_R*NOVA.fld']));
            if ~isempty(fins)
                % Load images
                cmiObj.loadImg(0,{fullfile(bdir,inm,fexp{1,j}),...
                                  fullfile(bdir,inm,fins(end).name)});
                out(nout,1:2) = {fexp{1,j} , fins(end).name};

                % Load mask
                fmsk = dir(fullfile(bdir,inm,[inm,'*_d',d(j),'*Exp_VOI.fld']));
                cmiObj.loadMask(fullfile(bdir,inm,fmsk(1).name));

                figure('Name',fexp{1,j}(1:8)); 
                
                % Expiration image (1)
                subplot(2,2,1)
                X = cmiObj.img.mat(:,:,:,1);
                X = X(cmiObj.img.mask.mat);
                [xmin,~,xt,n] = getDist(X);
                title(['Exp: ',num2str(n)])
                pause(0.1);

                % Inspiration image (2)
                subplot(2,2,2)
                Y = cmiObj.img.mat(:,:,:,2);
                Y = Y(cmiObj.img.mask.mat);
                [ymin,ymax,yt,n] = getDist(Y);
                title(['Ins: ',num2str(n)])
                pause(0.1);
                
                if any(isnan([ymin,ymax]))
                    disp('check')
                end
                
                subplot(2,2,3)
                [~,D,AUC] = PPplot(X,Y);

                out(nout,:) = { fexp{1,j} , fins(end).name , xmin-1000 , xt , ...
                                ymin-1000 , ymax-1000 , yt , D , AUC};
            end
        end
    end
end

assignin('base','results',out)


function [Dmin,Dmax,t,n] = getDist(X)
    % normalize CT lung data
    toff = 1; % offset to avoid 0's
    X = log(X+toff);
    max_logCT = max(X);
    X = X/max_logCT;
    % Mixed-Gaussian Fit
    [u,sig,t,n]=fit_mix_gaussian(X,1);
    i = find(t==max(t));
    % Unnormalize data
    u_plus=u(i)+1.6*sig(i);
    u_minus=u(i)-1.6*sig(i);
    Dmax=(exp(u_plus.*max_logCT))-toff;
    Dmin=(exp(u_minus.*max_logCT))-toff;
    t = t(i);
    

