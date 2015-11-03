% PRMclass function
% Calculate PRM from input values using current method
function [labels,prmpcts] = calcPRM(self,vals,vec,cvec,clabel,mask)
% Inputs: vals = [np x ndim] matrix of image values
%                       (first column is current (dynamic) image)
%         vec = [1 x ndim] vector of image #'s corresponding to vals
%         mask = 3D image mask

labels = self.prmmap(:,2);
prmpcts = zeros(1,length(labels));
if (nargin==6) % need all inputs
    [np,ndim] = size(vals);
    if (ndim==length(vec)) && (numel(find(mask))==np)
        tthresh = self.thresh; tthresh(tthresh(1:2*size(tthresh,1))==0) = cvec;
        tcut = self.cutoff; 
        if ~isempty(tcut)
            tcut(tcut(:,1)==0,1) = cvec;
        end
        [ind,C] = prm_seg(vals,vec,tthresh,tcut);
        if ~isempty(ind)
            % Check if # of thresholds matches current mapping scheme
            if size(C,2)==size(self.prmmap{1},2)
                prmvals = zeros(np,1);
                % Map prm C-index to desired PRM value
                % (note: any unmatched indexes are cropped)
                for i = 1:self.nprm
                    Li = find(ismember(C,self.prmmap{i,1},'rows'));
                    for j = 1:length(Li)
                        prmvals((ind==Li(j)) & (ind>0)) = i;
                    end
                end
            % Otherwise, change to arbitrary scheme
            else
                prmvals = ind;
                nC = size(C);
                self.nprm = nC(1);
                self.prmmap = [ mat2cell(C,ones(1,nC(1)),nC(2)) , ...
                                strcat('Region',cellstr(num2str((1:nC(1))')))];
                self.cmap = lines(nC(1));
            end
            if ~isempty(clabel) && iscell(clabel) && (length(clabel)==ndim)
                self.dlabels = clabel;
            else
                self.dlabels = strcat(repmat({'Dim'},ndim,1),num2str(vec(:)));
            end
            
            % Display scatterplot if desired
            %   2D only!
            if self.SPopts.show && (ndim>1)
                
%                 % Initialize PRM mask:
%                 self.mask.setDims(size(mask));
%                 self.mask.clear;
                
                if np > self.SPopts.Nmax
                    % Need to down-sample for large data sets
                    ind = round(1:(np/self.SPopts.Nmax):np);
                else
                    ind = 1:np;
                end
                % Plot the data
%                 txlim = []; tylim = [];
%                 if ~isempty(self.cutoff)
%                     ii = self.cutoff(:,1); ii(ii==0) = vec(1);
%                     ix = find(ii==vec(2),1);
%                     iy = find(ii==vec(1),1);
%                     if ix
%                         txlim = self.cutoff(ix,2:3);
%                     end
%                     if iy
%                         tylim = self.cutoff(iy,2:3);
%                     end
%                 end
                if isempty(self.hascatter) || ~(ishandle(self.hascatter) ...
                        && strcmp(get(self.hascatter,'Tag'),'PRMscatter'))
                    pos = get(0,'ScreenSize');
                    self.hfscatter = figure('Name','PRMscatter',...
                        'Position',[pos(3:4),pos(3:4)*2]/4);
                    self.hascatter = subplot(1,2,1);
                    set(self.hascatter,'Tag','PRMscatter');
                    hatext = subplot(1,2,2);
                    set(hatext,'YDir','reverse');
                    axis(hatext,'off');
                    title(hatext,'PRM Stats:','FontSize',12,'FontWeight','bold');
                    self.htscatter = text(0,0.5,'');
                    set(self.htscatter,'FontName','FixedWidth');
                end
                [ha,hs] = prm_plot(vals(ind,1),vals(ind,2),prmvals(ind),self.cmap,...
                                   self.prmmap(:,2),prmpcts,self.hascatter,...
                                   'xlabel',self.dlabels{1},'ylabel',self.dlabels{2},...
                                   'xlim',[self.SPopts.Xmin,self.SPopts.Xmax],...
                                   'ylim',[self.SPopts.Ymin,self.SPopts.Ymax]);
                
                               % Add threshold lines
                d = self.thresh(:,1:2);
                d(d==0) = cvec;
                m = self.thresh(:,3);
                b = self.thresh(:,4);
                ind = d(:,1)>d(:,2);
                m(ind) = 1./m(ind);
                plotlines(ha,m,b);
                
                % Update PRM stats:
                % Calculate statistics
                tot = zeros(self.nprm,1);
                for i = 1:self.nprm
                    tot(i) = sum(prmvals==i);
                end
                npcropped = nnz(prmvals);
                str = [ self.prmmap(:,2) , num2cell(tot/np*100) , ...
                        num2cell(tot/npcropped*100)]';
                set(self.htscatter,'String',...
                    ['                [% of VOI] ; [% of PRM]',char(10)...
                     sprintf('%16s = % 7.3f ; % 7.3f\n',str{:})]);
                
%   Commented out by CJG 20151103. Using check box for Wiener2 filter in ImageClass/calcPRM.                 
%                 % Calculate PRM stats to return:
%                 if self.normchk
%                     np = nnz(prmvals);
%                 end
                prmpcts = tot/np;
                
                self.hascatter = ha;
                self.hfscatter = get(ha,'Parent');
                self.hsscatter = hs;
            end
            
            % Update object values:
            if ~self.check
                self.mat = double(mask);
            end
            self.mat(mask) = prmvals;
            self.mat(~mask) = nan;
            [d(1),d(2),d(3),d(4)] = size(self.mat);
            self.dims = d;
            self.check = true;
            
        end
    end
end