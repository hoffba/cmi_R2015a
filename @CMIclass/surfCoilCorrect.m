% CMIclass function
function surfCoilCorrect(self,~,~)
% Corrects signal dropoff in MRI images

self.img.surfCoilCorrect(self.vec);
self.dispUDimg;
