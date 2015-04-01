% CMIclass function
% Set current view
function setView(self,x,~)
if (nargin>1) % Determine inputs
   if (nargin==3) && ishandle(x)
       val = find(strcmp(get(x,'Tag'),{'button_row','button_col','button_slc'}));
   else
       val = round(x);
   end
    if ~isempty(val) && any(val==(1:3))
        if ~isempty(self.h)
            % Must make sure exactly one of the buttons is pressed at any time
            switch val
                case 1
                    set([self.h.button_col,self.h.button_slc],'Value',0)
                    set(self.h.button_row,'Value',1)
                case 2
                    set([self.h.button_row,self.h.button_slc],'Value',0)
                    set(self.h.button_col,'Value',1)
                case 3
                    set([self.h.button_col,self.h.button_row],'Value',0)
                    set(self.h.button_slc,'Value',1)
            end
        end
        self.orient = val;
        
        % Update the display
        self.dispUDview;
        if self.guicheck
            % Update GUI objects
            ns = self.img.dims(val);
            enslc = 'Off';
            if ns>1
                enslc = 'On';
            end
            val = self.slc(val);
            set(self.h.edit_slc,'String',num2str(val),...
                                'Enable',enslc);
            set(self.h.slider_slice,'Value',val,'Max',ns,...
                                    'SliderStep',[1 1]/max(ns-1,1),...
                                    'Enable',enslc);
            set(self.h.text_ns,'String',num2str(ns));
        end
    end
end