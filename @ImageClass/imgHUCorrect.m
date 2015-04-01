% ImageClass function
% Scale current image HU values using blood & air values (y = m*x + b)
% Note: can correct multiple images in the set
function imgHUCorrect(self,vec,airHU,bloodHU)
if self.check
    if (nargin == 4)
        if isnumeric(airHU) && isnumeric(bloodHU) && ~any(vec<1) && ~any(vec>self.dims(4))
            %m = m(1); b = b(1); % only use one value
            vec = uint8(vec); % make sure they're ints
            
            % m = (CorrAir - CorrBlood)/(MeanAir - MeanBlood)
            % b = CorrBlood - m*MeanBlood
            nvec = self.dims(4);
            oM = self.scaleM;
            oB = self.scaleB;
            CorrAir = -995;
            CorrBlood = 37;
            CorrAir = CorrAir.*ones(1,nvec);
            CorrBlood = CorrBlood.*ones(1,nvec);
            m = (CorrAir - CorrBlood)./(airHU - bloodHU);
            b = CorrBlood - m.*bloodHU;
            
            for i = 1:length(vec)
                % Scale already xformed values by new scaling and offset
                % img = (I*oM + oB)*m + b;
                % new scale = oM*m and new offset = oB*m + b;
                tm = self.scaleM(vec(i));
                tb = self.scaleB(vec(i));
                %self.mat(:,:,:,vec(i)) = m(i)/tm * (self.mat(:,:,:,vec(i)) - tb) + b(i);
                self.mat(:,:,:,vec(i)) = m(i) * self.mat(:,:,:,vec(i)) + b(i);
                self.scaleM(vec(i)) = tm*m(i);
                self.scaleB(vec(i)) = tb*m(i)+b(i);
            end
        end
    end
end