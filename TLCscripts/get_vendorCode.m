function vndr = get_vendorCode(mnfr_str)
% Auxiliary function to genrate three-letter vendor-code, depending on "Manufactureer" input string.
% 
% MD 20160523: extracted from "build_acringroup"
%

% Initialize possibel string patterns:
PMS1 = 'Phi'; PMS2 = 'PHI';
GMS1 = 'GE';  GMS2 = 'Ge';
SMS1 = 'SIE'; SMS2 = 'Sie';

% check manufacturer:
p1 = strfind(mnfr_str, PMS1); p2 = strfind(mnfr_str, PMS2);
g1 = strfind(mnfr_str, GMS1); g2 = strfind(mnfr_str, GMS2);
s1 = strfind(mnfr_str, SMS1); s2 = strfind(mnfr_str, SMS2);

if ( ~isempty(p1) || ~isempty(p2) ) % True for Philips
    vndr = 'PMS';
elseif ( ~isempty(g1) || ~isempty(g2) ) % True for GE
    vndr = 'GMS';
elseif ( ~isempty(s1) || ~isempty(s2) ) % True for Siemens
    vndr = 'SMS';
else
    vndr = 'UNK'
end % if ~isempty

return