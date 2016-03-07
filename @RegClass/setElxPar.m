% RegClass function
% Function to manually set Elastix parameter and update GUI with new value
%   includes Parameter Name and Value validity checks
% Syntax:
%       RegClass.setElxPar(#,'Name',Value,...)
function setElxPar(self,i,varargin)
np = length(varargin)/2;
if (nargin>1) && isnumeric(i) && ~isempty(i) ...
        && ismember(i,1:length(self.elxObj.Schedule)) ...
        && (round(np)==np)
    
    % First set ElxClass properties
    if np~=0
        pname = varargin(1:2:end);
        vals = varargin(2:2:end);
        nres = self.elxObj.getPar(i,'NumberOfResolutions');
        C = {};
        for ipar = 1:np
            tval = vals{ipar};
            tC = {};
            switch pname{ipar}
                case 'jac'
                    tval = logical(tval);
                    self.jac = tval;
                    set(self.h.checkbox_jac,'Value',tval);
                case 'jacmat'
                    self.jacmat = logical(tval);
                    set(self.h.checkbox_jacmat,'Value',tval);
                case 'def'
                    self.def = logical(tval);
                    set(self.h.checkbox_def,'Value',tval);
                case 'FixedImagePyramidSchedule'
                    if isnumeric(tval) && (length(tval)==3*nres) && all(tval>0)
                        tC = [pname(ipar),{tval}];
                    end
                case 'MovingImagePyramidSchedule'
                    if isnumeric(tval) && (length(tval)==3*nres) && all(tval>0)
                        tC = [pname(ipar),{tval}];
                    end
                case 'MaximumNumberOfIterations'
                    if isnumeric(tval) && (length(tval)==nres) && all(tval>0)
                        tC = [pname(ipar),{tval}];
                    end
                case 'SP_A'
                    if isnumeric(tval) && (length(tval)==nres) && all(tval>0)
                        tC = [pname(ipar),{tval}];
                    end
                case 'SP_a'
                    if isnumeric(tval) && (length(tval)==nres) && all(tval>0)
                        tC = [pname(ipar),{tval}];
                    end
                case 'NumberOfSpatialSamples'
                    if isnumeric(tval) && (length(tval)==nres) && all(tval>0)
                        tC = [pname(ipar),{tval}];
                    end
                case 'GridSpacingSchedule'
                    if isnumeric(tval) && (length(tval)==3*nres) && all(tval>0)
                        tC = [pname(ipar),{tval}];
                    end
                case 'NumberOfResolutions'
                    % Shouldn't need more than 5 resolution steps
                    if isnumeric(tval) && ismember(tval,1:5) && (tval~=nres)
                        % Grab current values:
                        gcheck = strcmp('BSplineTransform',self.elxObj.getPar(i,'Transform'));
                        if gcheck
                            rdata = reshape(self.elxObj.getPar(i,'GridSpacingSchedule'),3,nres);
                        else
                            rdata = nan(3,nres);
                        end
                        rdata = [reshape(self.elxObj.getPar(i,'FixedImagePyramidSchedule'),3,[]);...
                                 reshape(self.elxObj.getPar(i,'MovingImagePyramidSchedule'),3,[]);...
                                 self.elxObj.getPar(i,'MaximumNumberOfIterations');...
                                 self.elxObj.getPar(i,'SP_A');...
                                 self.elxObj.getPar(i,'SP_a');...
                                 self.elxObj.getPar(i,'NumberOfSpatialSamples');...
                                 rdata];
                        % Add/Remove from front
                        if tval<nres
                            rdata(:,1:(nres-tval)) = [];
                        else
                            rdata = [repmat(rdata(:,1),1,tval-nres),rdata];
                        end
                        nres = tval;
                        % Pass parameters to ElxClass
                        tC = [pname(ipar),{tval},...
                              'FixedImagePyramidSchedule',reshape(rdata(1:3,:),1,[]),...
                              'MovingImagePyramidSchedule',reshape(rdata(4:6,:),1,[]),...
                              'MaximumNumberOfIterations',rdata(7,:),...
                              'SP_A',rdata(8,:),...
                              'SP_a',rdata(9,:),...
                              'NumberOfSpatialSamples',rdata(10,:)];
                        if gcheck
                            tC = [tC,'GridSpacingSchedule',reshape(rdata(11:13,:),1,[])];
                        end
                    end
                case 'ImageSampler'
                    if any(strcmp(tval,{'Full','Random','Grid','RandomCoordinate','RandomSparseMask'}))
                        tC = [pname(ipar),{tval}];
                    end
                case 'Metric'
                    if ischar(tval) && ismember(tval,{'AdvancedMeanSquares',...
                                                      'AdvancedNormalizedCorrelation',...
                                                      'AdvancedMattesMutualInformation',...
                                                      'NormalizedMutualInformation'})
                        if strncmp(self.elxObj.getPar(i,'Registration'),'MultiMetric',11)
                            tmetr = self.elxObj.getPar(i,'Metric');
                            tmetr{1} = tval;
                        else
                            tmetr = tval;
                        end
                        tC = [pname(ipar),{tmetr}];
                    end
                case 'TransformBendingEnergy'
                    % *** This is for the value Metric1Weight ***
                    % The Metric and other properties are automated
                    if isnumeric(tval) && ~isempty(tval) && ~isnan(tval)
                        tmetr = self.elxObj.getPar(i,'Metric');
                        if tval>0
                            tC = {'Metric1Weight',tval};
                            if ~strncmp(self.elxObj.getPar(i,'Registration'),'MultiMetric',11)
                                % Need to add a metric
                                tC = [tC,{'Metric',{self.elxObj.getPar(i,'Metric'),...
                                                    'TransformBendingEnergyPenalty'},...
                                          'Metric0Weight',1,...
                                          'UseRelativeWeights','false',...
                                          'Registration','MultiMetricMultiResolutionRegistration'}];
                            end
                        elseif iscell(tmetr) % Remove penalty from the Metrics
                            tC = {'Metric',tmetr{1},...
                                  'Metric0Weight',[],'Metric1Weight',[],...
                                  'UseRelativeWeights',[],...
                                  'Registration','MultiResolutionRegistration'};
                        end
                    end
                case 'FinalGridSpacingInVoxels'
                    if isnumeric(tval) && (length(tval)==3) && ~any(isnan(tval))
                        tC = [pname(ipar),{tval}];
                    end
                case 'DefaultPixelValue'
                    if isnumeric(tval) && ~isnan(tval)
                        tC = [pname(ipar),{tval}];
                    end
                case 'ResampleInterpolator'
                    if ischar(tval) && ismember(tval,{'FinalNearestNeighborInterpolator',...
                                                      'FinalLinearInterpolator',...
                                                      'FinalBSplineInterpolator'})
                        tC = [pname(ipar),{tval}];
                    end
                case {'WriteTransformParametersEachResolution','WriteResultImageAfterEachResolution'}
                    if ismember(tval,{'true','false'})
                        tC = [pname(ipar),{tval}];
                    end

            end
            C = [C,tC];
        end
        % Set all ElxClass parameters at once:
        if ~isempty(C)
            self.elxObj.setPar(i,C{:});
        end
    end
    
    % Update all GUI objects if updating selected Schedule step
    if i==self.ind
        self.selectTform(self.ind);
    end
end