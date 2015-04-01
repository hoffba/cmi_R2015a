% Mat3Dclass function
% Rotate 2D by multiples of 90 (CLOCKWISE)
function p = rotate90(self,dim,n)
p = 1:3;
if (nargin == 3) && ~isempty(n) && isnumeric(n) && ismember(dim,1:3)
    n = round(rem(n,4)); % minimize to 360degree rotation
    if n
        % Determine relative dimension order (determined by current view)
        switch dim
            case 1 % row view
                d = [2,3];
            case 2 % col view
                d = [1,3];
            case 3 % slice view
                d = [1,2];
        end
        
        % Determine direction
        if (abs(n)~=2) && (n ~= rem(n,2))
            n = -rem(n,2);
        end
        
        if (abs(n)==2) % 180deg
            self.flip(d(2));
        else
            p(p~=dim) = fliplr(p(p~=dim));
            self.permuteMat(p);
        end
        if (n==-1)
            td = d(2);
        else
            td = d(1);
        end
        self.flip(td);
    end
end