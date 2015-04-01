% RegClass function
function setOptionsD(self, ~, ~)
if get(self.h.OptMenu, 'Value') == 3
%Enter in options for registration
        prompt = {'Enter p for monomodal, m for multimodal (autodetect default):', 'Choose Registration option (Rigid, Affine, NonRigid):', 'Enter Maximum number of grid refinement steps:',  'Display debug information (0, 1, 2):', 'Sigma fuid regularization constant:', 'Alpha constant (reduce edge influence):', 'Sigma Diff constant (diffusion regularization):', 'Interpolation method (Linear or Cubic):'};
        default = {'auto', 'Affine', 'auto', '2', '4', '4', '1', 'Linear'};
        answer = inputdlg(prompt, 'Registration Options', 1, default);
        modal = answer{1};
        if strcmp(answer{1}, 'auto')
            modal = [];
        end
        if strcmp(answer{3}, 'auto')
            mxref = [];
        else
            mxref = str2double(answer{3});
        end

        %create options structure
        self.demonOpt{1} = struct('Similarity',modal,'Registration',answer(2),'MaxRef',mxref,'Verbose',str2double(answer(4)),'SigmaFluid',str2double(answer(5)),'Alpha',str2double(answer(6)),'SigmaDiff',str2double(answer(7)),'Interpolation',answer(8)); 

        %affine optimizer options
        prompt = {'Set on if gradient is available, default is off:', 'Set linear search method, (0 or 1):', 'Set display to iter for every iteration, plot for line search results:', 'Max number of function iterations:', 'Max number of function evaluations:', 'Function Tolerance', 'Minimum xval tolerance', 'Minimum stepsize changes'};
        default = {'off', '1', 'off', '300', '1000', '1e-14', '1e-11', '1e-6'};
        answer = inputdlg(prompt, 'Affine/Rigid Optimizer options', 1, default);

        %run registration with the options given
        self.demonOpt{2} = struct('GradObj',answer(1),'GoalsExactAchieve',str2double(answer(2)), 'Display',answer(3),'MaxIter',str2double(answer(4)),'MaxFunEvals',str2double(answer(5)),'TolFun',str2double(answer(6)), 'TolX', str2double(answer(7)),'DiffMinChange',str2double(answer(8)));
end