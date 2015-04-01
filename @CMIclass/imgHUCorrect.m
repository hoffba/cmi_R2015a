% CMIclass function
% HU Correct current vec image values using airHU and bloodHU
% (y = m*x + b)
% vec=0 means to process all vectors the same
% Revision history:
% 22 Oct 2013 jlb Sent to Ben for inclusion so we can have a menu item, but possibly not fully debugged 
%                 (although we were using in a previous version of the program)
function imgHUCorrect(self,vec,airHU,bloodHU)
if self.img.check
    

    % Retrieve current image scaling
    nvec = self.img.dims(4);
    oM = self.img.scaleM;
    oB = self.img.scaleB;
    
    do_correct = false;
    if (nargin==4) && isnumeric(airHU) && isnumeric(bloodHU) && ~any(vec<1) && ~any(vec>self.dims(4))
        tvec = vec;
        tairHU = airHU;
        tbloodHU = bloodHU;
        do_correct = true;
    else
        % Prompt for new airHU and bloodHU mean values for all dims
        str = {'Image vector:','Input Current Air HU (-995):','Input current blood HU (37):'};
        oairHU = -995; obloodHU = 37;
        oairHU = oairHU.*ones(1,nvec);
        obloodHU = obloodHU.*ones(1,nvec);
        def = {num2str(1:nvec),num2str(oairHU),num2str(obloodHU)};
        answer = inputdlg(str,'HU Correct:',1,def);
        if ~isempty(answer)
            tvec = str2num(answer{1});
            tairHU = str2num(answer{2});
            tbloodHU = str2num(answer{3});
            do_correct = true;
        end
    end
    if do_correct
        if (length(tvec)==length(tairHU)) && (length(tbloodHU)==length(tvec)) &&...
                ~any(isnan([tvec tairHU tbloodHU]))&& ~any(isempty([tvec tairHU tbloodHU]))
            if (tvec==0)
                tvec = 1:nvec;
                tairHU = tairHU.*ones(1,nvec);
                tbloodHU = tbloodHU.*ones(1,nvec);
            end
            %self.img.imgScale(tvec,tM,tB);
            self.img.imgHUCorrect(tvec,tairHU,tbloodHU);
            
            % update display scaling and PRM as needed
            tM = self.img.scaleM;
            tB = self.img.scaleB;
            oM = oM(tvec)'*[1 1];
            oB = oB(tvec)'*[1 1];
            tM = tM'*[1 1];
            tB = tB'*[1 1];
            self.clim(tvec,:) = tM./oM .* (self.clim(tvec,:) - oB) + tB;
            self.dispUDview;
            self.dispUDhist;
            if self.prmcheck
                self.img.calcPRM;
            end
            % update GUI objects:
            if self.guicheck
                [tmin,tmax] = self.getColorLimits;
                [imin,imax] = self.getColorMinMax;
                cpad = (imax - imin) / 2;
                set(self.h.slider_cmin,'Min',(imin-cpad),'Max',(imax+cpad),...
                                        'Value',tmin,'Enable','On');
                set(self.h.slider_cmax,'Min',(imin-cpad),'Max',(imax+cpad),...
                                        'Value',tmax,'Enable','On');
                set(self.h.edit_cmin,'String',num2str(tmin),'Enable','On');
                set(self.h.edit_cmax,'String',num2str(tmax),'Enable','On');
            end
        end
    end
end
