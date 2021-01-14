function stvar = make_safe_fieldname(tagStr)
% this function checks potential fieldnames and makes sure they are valid
% for MATLAB syntax, e.g. 2DInterpolation -> x2DInterpolation (must begin
% with a letter)

tagStr = strtrim(tagStr);

if (isletter(tagStr(1)))
    stvar = tagStr;
else
    stvar = strcat('x', tagStr);
end

if strfind(stvar, ';'), stvar = strrep(stvar, ';', '_'); end
if strfind(stvar, '@'), stvar = strrep(stvar, '@', '_'); end % VD13
if strfind(stvar, '-'), stvar = strrep(stvar, '-', '_'); end