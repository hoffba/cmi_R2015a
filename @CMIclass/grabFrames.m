% CMIclass function
% Grab image frames from display (for montage / movie)
function F = grabFrames(self,opts)
% Save montage of current image(s) over defined dimension (typically 4, time, or 3, slice)
% Inputs: vdim = dimension to view
%         mdim = dimension to scroll over movie
%         mind = indexes of movie frames
oview = self.orient;
oslc = self.slc(oview);
ovec = self.vec;
F = struct('cdata',{},'colormap',{}); % initialize the frame structure
if nargin == 1 % Ask user to input required parameters
	opts = struct('vdim',{},'mdim',{},'mind',{});
    prompt = {'View dimension:',...
              'Increment dimension:',...
              'Increment indexes ("start":"increment":"end"):'};
    sstr = num2str(oview);
    ns = self.img.dims(oview);
    def = {num2str(oview),...       % View orientation
           sstr,...                 % Increment dimension
           ['1:1:' num2str(ns)]};          % Increment indexes
    answer = inputdlg(prompt,'Create Montage',1,def);
    if ~isempty(answer)
        vdim = str2double(answer{1}); % View orientation
        mdim = str2double(answer{2}); % Scroll dimension
        mind = eval(answer{3});       % Scroll indexes
        if ~((vdim > 3) || (vdim < 1) || (mdim > 4) || (mdim < 1) ||...
                (min(mind)<1) || (max(mind)>self.img.dims(mdim)) ||...
                ((mdim ~= vdim) && (mdim ~= 4)))
            opts.vdim = vdim;
            opts.mdim = mdim;
            opts.mind = mind;
        end
    end
end
if ~isempty(opts) && isstruct(opts) && all(isfield(opts,{'vdim','mdim','mind'}))
    v = [opts.vdim opts.mdim opts.mind];
    if (~any(~isnumeric(v) | isnan(v) | isempty(v)) || ~isempty(fname))
        if (opts.vdim ~= self.orient)
            self.setView(opts.vdim);
        end
    %     % Determine image location on figure:
    %     set(self.haxes,'Units','pixels');
    %     pos = get(self.haxes,'Position');
    %     set(self.haxes,'Units','normalized');
    %     % Accounting for plotting aspect ratio:
    %     aspn = pbaspect(self.haxes);
    %     aspn = aspn(1:2)/(max(aspn(1:2)));
    %     npos(3:4) = aspn.*pos(3:4);
    %     npos(1:2) = pos(1:2) + (pos(3:4) - npos(3:4))/2;
        currunits = get(self.haxes,'Units');
        set(self.haxes,'Units','pixels');
        npos = plotboxpos(self.haxes);
        set(self.haxes,'Units',currunits);
        npos(1:2) = npos(1:2)+1;
        npos(3:4) = npos(3:4)-2;
        for i = 1:length(opts.mind)
            % Change to next slice/vec
            if (opts.mdim == 4)
                self.setVec(opts.mind(i));
            else
                self.setSlice(opts.mind(i));
            end
            F(i) = getframe(self.hfig,npos);
        end
        % Set view back to original
        if self.orient ~= oview
            self.setView(oview);
        end
        if self.vec ~= ovec
            self.setVec(ovec);
        end
        if self.slc(oview) ~= oslc
            self.setSlice(oslc);
        end
    end
end