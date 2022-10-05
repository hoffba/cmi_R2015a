function P = calcProgression(exp,ins,J,seg)

% Preprocessing:
seg = logical(seg);
exp = imgaussfilt(exp);
ins = imgaussfilt(ins);

P = zeros(size(exp));
P(seg) = exp(seg) - ins(seg) ./ J(seg);

