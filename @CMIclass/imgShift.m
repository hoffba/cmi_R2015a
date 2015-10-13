% CMIclass function
% Allows user to manually translate images
function imgShift(self,~,~)
if self.img.check
    
    % Determine what to shift:
    v = 1;
    if self.img.dims(4)>1
        answer = questdlg('Move all images or only the current?',...
            'Shift Image','All','Current','Cancel','Current');
        switch answer
            case 'All'
                v = 1:self.img.dims(4);
            case 'Current'
                v = self.img.vec;
            case 'Cancel'
                return;
        end
    end
    
    % Input points : (1) from, (2) to
    [x,y] = getpts(self.haxes);
    
    % Determine image dimensions to shift
    dx = zeros(1,3);
    dx(setxor(1:3,self.orient)) = [y(2)-y(1) , x(2)-x(1)];
    dx = dx./self.img.voxsz;
    
    % Shift the matrix:
    self.img.shiftMat(v,dx);
    
    % Update display:
    self.dispUDslice;
    
end