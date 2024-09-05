% Create morphological structuring element
% Function 'strel' might not meet the requirement

function index = calcSe(px,r)
  xdim = floor(r/px(1))*2 + 1;
  ydim = floor(r/px(2))*2 + 1;
  zdim = floor(r/px(3))*2 + 1;
  o = [xdim+1 ydim+1 zdim+1]/2;
  index = zeros(zdim*ydim*xdim,3);
  num = 1;
  for x = 1:xdim
    for y = 1:ydim
      for z = 1:zdim
        if ( ((o(1)-x)*px(1))^2 +((o(2)-y)*px(2))^2 + ((o(3)-z)*px(3))^2 ) <= r^2
          index(num,:) = [x-o(1), y-o(2), z-o(3)];
          num = num+1;
        end
      end
    end
  end
  index = index(1:num-1,:);
    