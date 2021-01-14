function [tempArr,mrprot,count] = walkdown(mrprot,structpos,tempArr,fillin,count)


fieldnames = eval(['fieldnames(' structpos ');']);   

for i = 1:length(fieldnames)
    if (eval(['isstruct(' structpos '.' fieldnames{i} ')']))
        ind = strfind(tempArr,'{');
        tempArr = tempArr((ind(1)+1):end);
        [tempArr,mrprot,count] = walkdown(mrprot,[structpos '.' fieldnames{i}],tempArr,fillin,count);
        ind2 = strfind(tempArr,'}');
        if (isempty(ind2))
            tempArr = '';
        else
            tempArr = tempArr((ind2(1)+1):end);
        end;
        count = count +1;
    else
        ind = strfind(tempArr,'{');
        ind2 = strfind(tempArr,'}');
        if (fillin == 1)
            try
                eval([structpos '.' fieldnames{i} ' = ' ( tempArr((ind(1)+1):(ind2(1)-1)) ) ';']);
            catch
                eval([structpos '.' fieldnames{i} ' = 0;']);
            end;
        end;
        tempArr = tempArr((ind2(1)+1):end);
        count = count +1;
    end;
end;