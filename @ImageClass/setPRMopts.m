% ImageClass function
function stat = setPRMopts(self,cvec)
stat = false;
if (nargin==2) && isnumeric(cvec)
    % Load figure and retrieve figure handles
    hf = setRegions;
    h = guidata(hf); % grab figure object handles
    
    % Set GUIdata
    h.thresh = self.prm.thresh;
    h.cutoff = self.prm.cutoff;
    h.cmap = self.prm.cmap;
    h.prmmap = self.prm.prmmap;
    h.maxnp = self.prm.npmaxscat;
    h.prmscatter = self.prm.prmscatter;
    h.normchk = self.prm.normchk;
    guidata(hf,h);
    
    % Set GUI object properties
    if self.check
        str =  cellstr(strcat(num2str((1:self.dims(4))'),': ',self.labels(:)));
    else
        str = {''};
    end
    set(h.list_imgs,'String',str,'Value',cvec);
    set(h.list_regions,'String',self.prm.prmmap(:,2));
    set(h.axes_color,'Color',self.prm.cmap(h.i,:));
    set(h.edit_name,'String',self.prm.prmmap{h.i,2});
    set(h.table_thresh,'Data',num2cell(self.prm.thresh));
    set(h.table_C,'Data',num2cell(self.prm.prmmap{h.i,1})');
    set(h.table_crop,'Data',num2cell(self.prm.cutoff));
    set(h.edit_maxnp,'String',num2str(self.prm.npmaxscat));
    set(h.checkbox_showScatter,'Value',self.prm.prmscatter);
    set(h.checkbox_VolNorm,'Value',self.prm.normchk);
    
    uiwait(hf); % wait for options figure
    
    if ishandle(hf)
        % Grab PRM options from figure
        h = guidata(hf);
        close(hf);
        % Determine what labels to send in
        tval = h.thresh(:,1:2);
        tval = unique(tval(:));
        tval(tval==0) = [];
        % Set PRMclass properties
        if self.check
            labels = self.labels(tval(:));
        else
            labels = strcat('Dim',num2cell(num2str(tval')));
        end
        stat = self.prm.setOpts('thresh',h.thresh,'cutoff',h.cutoff,...
                                'cmap',h.cmap,'prmmap',h.prmmap,...
                                'labels',labels,...
                                'npmaxscat',h.maxnp,...
                                'prmscatter',h.prmscatter,...
                                'normchk',h.normchk);
    end
end