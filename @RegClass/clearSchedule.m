function clearSchedule(self)

stat = self.elxObj.rmStep;
if stat && self.guicheck
    set(self.h.listbox_Tforms,'String',{});
    self.selectTform(0);
end
