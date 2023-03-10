function J=readanddecodejson(x)


fid = fopen(x);
raw = fread(fid,inf);
str = char(raw');
fclose(fid);
J = jsondecode(str);