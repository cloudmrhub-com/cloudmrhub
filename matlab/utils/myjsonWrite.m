function [] =myjsonWrite(J,fn)
%before you have to jsonencode(the struct)
%myjsonWrite(jsonencode(O),filenameof the json file);
fid = fopen(fn, 'w');
if fid == -1, error('Cannot create JSON file'); end
fwrite(fid, J, 'char');
fclose(fid);
end