% CMIclass function
% Scale current image values (y = m*x + b)
function imgScale(self,~,~)
if self.img.check
    nvec = self.img.dims(4);
    oM = self.img.scaleM;
    oB = self.img.scaleB;
    str = {'Image vector:','Input scale factor (M):','Input scale offset (B):'};
    def = {num2str(1:nvec),num2str(oM),num2str(oB)};
    answer = inputdlg(str,'Scale Image (y=mx+b):',1,def);
    if ~isempty(answer)
        tvec = str2num(answer{1});
        tM = str2num(answer{2});
        tB = str2num(answer{3});
        if (length(tvec)==length(tM)) && (length(tB)==length(tvec)) &&...
                ~any(isnan([tvec tM tB]))&& ~any(isempty([tvec tM tB]))
            if (tvec==0)
                tvec = 1:nvec;
                tM = tM.*ones(1,nvec);
                tB = tB.*ones(1,nvec);
            end
            self.img.imgScale(tvec,tM,tB);
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