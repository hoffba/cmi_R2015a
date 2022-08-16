function stat = check_info(self,info_in,tag)
% Check new info against what is currently loaded
% Inputs: info_in = new info to check
%         tag = tag of image to be loaded

stat = false;

% What info fields to match:
info_check = {'ImageSize','PixelDimensions'};

% Determine what geometry to match:
switch tag
    case {'ct_exp','seg_exp','scatnet','reg_ins','jac','prm','tprm'}
        eitag = 'ct_exp';
    case {'ct_ins','seg_ins'}
        eitag = 'ct_ins';
    case {'vessels','csa'}
        eitag = 'vessels';
    otherwise
        warning('Invalid input tag: %s',tag);
        return;
end

if isempty(self.dat.(eitag).info)
    if ismember(tag,{'ct_exp','ct_ins','vessels'})
        self.dat.(tag).info = info_in;
    else
        % Load info from CT file
        fname = self.getFileName(tag);
        if exist(fname,'file')
            stat = true;
            self.dat.(eitag).info = niftiinfo(fname);
        else
            self.getfile(eitag);
        end
    end
else
    stat = true;
    for i = 1:numel(info_check)
        stat = stat && all(self.dat.(eitag).info.(info_check{i}) == info_in.(info_check{i}));
    end
end