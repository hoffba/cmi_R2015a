function P = calcProgression(exp,ins,J,seg)

% Preprocessing:
seg = logical(seg);
exp = imgaussfilt3(exp,1)+1000;
ins = imgaussfilt3(ins,1)+1000;

P = zeros(size(exp));
P(seg) = exp(seg) - ins(seg) .* J(seg);

