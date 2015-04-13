% RegClass function
function showTx0(self,M)
% Display initial transform in GUI

if nargin==1
    M = self.elxObj.Tx0.TransformParameters;
end
if ~isempty(M) && (numel(M)==12)
    M = M(:);
    set(self.h.edit_a11,'String',num2str(M(1)),'Enable','on');
    set(self.h.edit_a12,'String',num2str(M(4)),'Enable','on');
    set(self.h.edit_a13,'String',num2str(M(7)),'Enable','on');
    set(self.h.edit_a21,'String',num2str(M(2)),'Enable','on');
    set(self.h.edit_a22,'String',num2str(M(5)),'Enable','on');
    set(self.h.edit_a23,'String',num2str(M(8)),'Enable','on');
    set(self.h.edit_a31,'String',num2str(M(3)),'Enable','on');
    set(self.h.edit_a32,'String',num2str(M(6)),'Enable','on');
    set(self.h.edit_a33,'String',num2str(M(9)),'Enable','on');
    set(self.h.edit_t1,'String',num2str(M(10)),'Enable','on');
    set(self.h.edit_t2,'String',num2str(M(11)),'Enable','on');
    set(self.h.edit_t3,'String',num2str(M(12)),'Enable','on');
else
    set([self.h.edit_a11,self.h.edit_a12,self.h.edit_a13,...
        self.h.edit_a21,self.h.edit_a22,self.h.edit_a23,...
        self.h.edit_a31,self.h.edit_a32,self.h.edit_a33,...
        self.h.edit_t1,self.h.edit_t2,self.h.edit_t3],'Enable','off');
end