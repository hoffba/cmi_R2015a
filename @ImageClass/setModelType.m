% ImageClass function
function setModelType(self,val)
% Set type of model for self.model

if (nargin==2)
    % First, remember to delete old model
    delete(self.model)
    % Next, initialize new model
    switch val
        case 1 % General
            self.model = FitClass;
        case 2 % Perfusion
            self.model = DCEclass;
        case 3 % Diffusion
            self.model = DiffClass;
    end
end