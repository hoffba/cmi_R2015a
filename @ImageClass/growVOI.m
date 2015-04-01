% ImageClass function
% Grow VOI with thresholds
function growVOI(self,vec,opt)
if self.check && self.mask.check && nargin>=2
    if (nargin<3)
        opt = 0; % default growing option
    end
    step = inputdlg({'X-Y step:','Z step','Max # steps:','Display:'},...
        'Dilation step',1,{'3','1','500','on'});
    if ~isempty(step)
        dopt = strcmp(step{4},'on');
        step = str2double(step);
        maxct = step(3); step(3) = [];
        done = false;
        thmask = ((self.mat(:,:,:,vec)>=self.thresh(vec,1)) & (self.mat(:,:,:,vec)<=self.thresh(vec,2)));
        if opt==1 % Ends-in option
            tslc = find(squeeze(any(any(self.mask.mat)))); % To find ends
            if (length(tslc)<2)
                done = true;
            else
                thmask(:,:,[1:(min(tslc)-1),(max(tslc)+1):end]) = false;
            end
        end
        count = 0;
        warning('off', 'Images:initSize:adjustingMag');
        omask = self.mask.mat;
        if dopt
            td = round(size(omask)/2);
            h = figure; colormap(gray)
            subplot(1,3,1)
            h1 = imagesc(squeeze(self.mask.mat(:,:,td(3))));axis square
            subplot(1,3,2)
            h2 = imagesc(squeeze(self.mask.mat(:,td(2),:)));axis square
            subplot(1,3,3)
            h3 = imagesc(squeeze(self.mask.mat(td(1),:,:)));axis square
            pause(0.01)
        end
        tic
        while ~done
            %figure(h);imagesc(omask(:,:,round(size(omask,3)/2)));colormap gray;axis off;
            self.mask.morph('dilate',step);
            self.mask.merge('intersect',thmask);
            if all(omask(:)==self.mask.mat(:)) || (count==maxct)
                done = 1;
            end
            omask = self.mask.mat;
            if dopt
                set(h1,'CData',squeeze(omask(:,:,td(3))))
                set(h3,'CData',squeeze(omask(td(1),:,:)))
                set(h2,'CData',squeeze(omask(:,td(2),:)))
                title(['Completed step #: ' num2str(count)])
                pause(0.01)
            else
                disp(['Completed step #: ' num2str(count)])
            end
            count = count + 1;
        end
        disp(['Total Steps = ' num2str(count-1)])
        toc
        delete(h);
    end
end