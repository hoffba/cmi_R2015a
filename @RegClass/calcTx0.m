% RegClass function
function calcTx0(self)
    if ~isempty(self.points{1}) && ~isempty(self.points{2})
        self.elObj.Tx0 = pts2T(self.points{2}. self.points{1});
    else
        disp('Points must be loaded to calculate Tx0');
    end
end