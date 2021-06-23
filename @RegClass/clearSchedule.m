function clearSchedule(self)

stat = self.elxObj.rmStep;
if stat
    set(self.h.listbox_Tforms,'String',{});
    self.selectTform(0);
end
