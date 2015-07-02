% CMIclass function
% Loads main figure and initialized handle structure
function stat = initMain(self)
stat = false;
% Check if figure needs to be loaded
if isempty(self.h) && self.guicheck
    self.h = guihandles(open('cmi.fig'));
end
% Update GUI callbacks
if ~isempty(self.h)
    % ADD LISTENERS/CALLBACKS FOR GUI CONTROLS
    set(self.h.mainFig,'CloseRequestFcn',@self.GUIexit);
    % ~~~~~ Sliders
    addlistener(self.h.slider_slice,'ContinuousValueChange', @self.setSlice);
    addlistener(self.h.slider_4D,   'ContinuousValueChange', @self.setVec);
    addlistener(self.h.slider_cmin, 'ContinuousValueChange', @self.setClim);
    addlistener(self.h.slider_cmax, 'ContinuousValueChange', @self.setClim);
    addlistener(self.h.slider_thMin, 'ContinuousValueChange', @self.imgThreshold);
    addlistener(self.h.slider_thMax, 'ContinuousValueChange', @self.imgThreshold);
    % ~~~~~ Edit Boxes
    set(self.h.edit_veclabel,'Callback',@self.setLabel);
    set(self.h.edit_vec,     'Callback',@self.setVec);
    set(self.h.edit_slc,     'Callback',@self.setSlice);
    set(self.h.edit_cmax,    'Callback',@self.setClim);
    set(self.h.edit_cmin,    'Callback',@self.setClim);
    % ~~~~~ Popups
    set(self.h.popup_voimode, 'Callback',@self.setVOImode);
    set(self.h.popup_colormap,'Callback',@self.setCMap);
    % ~~~~~ Buttons
    set(self.h.button_row,'Callback',@self.setView);
    set(self.h.button_col,'Callback',@self.setView);
    set(self.h.button_slc,'Callback',@self.setView);
    % ~~~~~ Checkboxes
    set(self.h.bgcheck,   'Callback',@self.setBGcheck);
    set(self.h.applytoall,'Callback',@self.setApplyAllCheck);
    % ~~~~~ UIMenu Items
    set(self.h.file_exit,        'Callback',@self.GUIexit);
    set(self.h.file_savemask,    'Callback',@self.saveMask);
    set(self.h.file_saveimg,     'Callback',@self.saveImg);
    set(self.h.file_mask2img,    'Callback',@self.mask2img);
    set(self.h.file_append,      'Callback',@self.loadImg);
    set(self.h.file_openimg,     'Callback',@self.loadImg);
    set(self.h.file_openmask,    'Callback',@self.loadMask);
    set(self.h.tools_eigentf,    'Callback',@self.eigOrient);
    set(self.h.tools_montage,    'Callback',@self.genMontage);
    set(self.h.tools_movie,      'Callback',@self.genMovie);
    set(self.h.tools_profile,    'Callback',@self.genProfile);
    set(self.h.tools_majorD,     'Callback',@self.getMaskMajorAxis);
    set(self.h.tools_drawVOI,    'Callback',@self.drawROI);
    set(self.h.tools_connToggle, 'Callback',@self.connToggle);
    set(self.h.tools_connThresh, 'Callback',@self.conn2voi);
    set(self.h.tools_connVOI,    'Callback',@self.conn2voi);
    set(self.h.tools_voi_fillholes,  'Callback',@self.fillHoles);
    set(self.h.tools_voi_segment,    'Callback',@self.growMask);
    set(self.h.tools_voi_erode,      'Callback',@self.erodeMask);
    set(self.h.tools_voi_dilate,     'Callback',@self.dilateMask);
    set(self.h.tools_voi_clearAll,   'Callback',@self.clearMask);
    set(self.h.tools_voi_clearSlice, 'Callback',@self.clearMask);
    set(self.h.tools_disp_bgSel,     'Callback',@self.setBGvec);
    set(self.h.tools_disp_overlay,   'Callback',@self.setOverCheck);
    set(self.h.tools_disp_histo,     'Callback',@self.setHistVis);
    set(self.h.tools_disp_overtransp,'Callback',@self.setOverTransp);
    set(self.h.tools_recon,          'Callback',@self.recon);
    set(self.h.tools_plot4D,         'Callback',@(varargin)self.img.plot4D);
    set(self.h.image_isosurf,    'Callback',@self.isosurf);
    set(self.h.image_showmip,    'Callback',@self.showMIP);
    set(self.h.image_flipdim,    'Callback',@self.imgFlip);
    set(self.h.image_rotate,     'Callback',@self.imgRotate);
    set(self.h.image_crop,       'Callback',@self.imgCrop);
    set(self.h.image_smooth,     'Callback',@self.imgSmooth);
    set(self.h.image_interp,     'Callback',@self.imgInterp);
    set(self.h.image_delete4D,   'Callback',@self.imgDelete);
    set(self.h.image_subsamp,    'Callback',@self.imgSubsample);
    set(self.h.image_scale,      'Callback',@self.imgScale);
    set(self.h.image_scaleHU,    'Callback',@(varargin)self.scale2HU);
    set(self.h.image_autoScale,  'Callback',@self.autoScale);
    set(self.h.image_threshmask, 'Callback',@self.thresh2mask);
    set(self.h.image_thresh,     'Callback',@self.imgThreshold);
    set(self.h.image_filter,     'Callback',@self.imgFilt);
    set(self.h.image_filtNOVA,   'Callback',@self.filtNOVA);
    set(self.h.image_surfcoilcorr,'Callback',@self.surfCoilCorrect);
    set(self.h.analysis_imgmath, 'Callback',@self.imgMath);
    set(self.h.analysis_prmopt,  'Callback',@self.setPRMopts);
    set(self.h.analysis_PRM,     'Callback',@self.activatePRM);
    set(self.h.analysis_calcCI,  'Callback',@self.calcPRMCI);
    set(self.h.analysis_PPplot,  'Callback',@(varargin)self.img.ppPlot(self));
    set(self.h.analysis_prmselect,'Callback',@self.selectPRMscatter);
    set(self.h.analysis_prmstats,'Callback',@(varargin)self.img.savePRM(varargin{:}));
    set(self.h.analysis_saveprm, 'Callback',@(varargin)self.img.savePRM(varargin{:}));
        strs = {'off','off','off'};
        strs{strcmp(class(self.img.model),{'FitClass','DCEclass','DiffClass'})} = 'on';
    set(self.h.analysis_genmod,  'Callback',@self.modCheck,...
                             'Checked',strs{1});
    set(self.h.analysis_perfmod, 'Callback',@self.modCheck,...
                             'Checked',strs{2});
    set(self.h.analysis_diffmod, 'Callback',@self.modCheck,...
                             'Checked',strs{3});
    set(self.h.analysis_modtype, 'Label',['Model Type: ',self.img.model.typeStr]);
    set(self.h.analysis_modsel,  'Callback',@self.modSelect,...
                                 'Label',['Model: ',self.img.model.getModName]);
    set(self.h.analysis_setx,    'Callback',@self.setFitX);
    set(self.h.analysis_curvefitopts,'Callback',@(varargin)self.img.model.setFitOptions);
    set(self.h.analysis_initguess,'Callback',@self.setPar0);
    set(self.h.analysis_dceopts, 'Callback',@self.setDCEopts);
    set(self.h.analysis_calcfit, 'Callback',@self.calcFit);
    set(self.h.analysis_lsf,     'Callback',@self.calcLinLSF);
    set(self.h.analysis_bmdcalc,'Callback',@self.CTdensityCal);
    set(self.h.analysis_MF,     'Callback',@self.calcMF);
    set(self.h.script_run,'Callback',@self.scriptRun);
    set(self.h.script_edit,'Callback',@self.scriptRun);
    set(self.h.script_create,'Callback',@self.scriptCreate);
    set(self.h.checkbox_thMin,'Callback',@self.imgThreshold);
    set(self.h.checkbox_thMax,'Callback',@self.imgThreshold);
    set(self.h.edit_thMin,'Callback',@self.imgThreshold);
    set(self.h.edit_thMax,'Callback',@self.imgThreshold);
    % ~~~~~ Context Menu Items
    set(self.h.tools_voi_stats,'Callback',@self.genStats);
    stat = true;
end
if self.img.check
    set(self.h.button_slc,'Value',self.orient==3);
    set(self.h.button_col,'Value',self.orient==2);
    set(self.h.button_row,'Value',self.orient==1);
    set(self.h.slider_4D,'Max',self.img.dims(4),...
                         'Value',self.vec,...
                         'SliderStep',[1,1]/(self.img.dims(4)-1))
    set(self.h.slider_slice,'Max',self.img.dims(3),...
                            'Value',self.slc(self.orient),...
                            'SliderStep',[1,1]/(self.img.dims(3)-1));
    set(self.h.edit_slc,'String',num2str(self.slc(self.orient)));
    set(self.h.edit_vec,'String',num2str(self.vec));
    set(self.h.edit_veclabel,'String',self.img.labels(self.vec));
    set(self.h.text_nv,'String',num2str(self.img.dims(4)));
    set(self.h.text_ns,'String',num2str(self.img.dims(3)));
    vext = self.getProp('ValLim');
    set([self.h.slider_thMin,self.h.slider_thMax],...
        'Min',vext(1),'Max',vext(2));
end