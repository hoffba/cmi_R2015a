% RegClass function
function showTx0(self,M)
% Display initial transform in GUI

if ~isempty(M)
    set(self.h.edit_a11,'String',num2str(M(1,1)),'Enable','on');
    set(self.h.edit_a12,'String',num2str(M(1,2)),'Enable','on');
    set(self.h.edit_a13,'String',num2str(M(1,3)),'Enable','on');
    set(self.h.edit_a21,'String',num2str(M(2,1)),'Enable','on');
    set(self.h.edit_a22,'String',num2str(M(2,2)),'Enable','on');
    set(self.h.edit_a23,'String',num2str(M(2,3)),'Enable','on');
    set(self.h.edit_a31,'String',num2str(M(3,1)),'Enable','on');
    set(self.h.edit_a32,'String',num2str(M(3,2)),'Enable','on');
    set(self.h.edit_a33,'String',num2str(M(3,3)),'Enable','on');
    set(self.h.edit_t1,'String',num2str(M(1,4)),'Enable','on');
    set(self.h.edit_t2,'String',num2str(M(2,4)),'Enable','on');
    set(self.h.edit_t3,'String',num2str(M(3,4)),'Enable','on');
end