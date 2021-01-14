function value = get_xprot_value(stubStr, tagType)
% stubStr contains the values of interest, but also might include
% modifiers. we will just ignore the modifiers, which seem to always be
% terminated by CR.

[tagStr, newtagType, remStr] = findNextTag(stubStr);
while (~isempty(remStr)) % found a tag
    if (~isempty(newtagType)) % not a modifier tag???
        error('parse_xprot()::get_xprot_value(): ERROR: found unknown tag!')
    else
        % these are the modifier tags we know about
        if ((strncmpi(tagStr, 'Precision', 9))      || ...
            (strncmpi(tagStr, 'LimitRange', 10))    || ...
            (strncmpi(tagStr, 'MinSize', 7))        || ...
            (strncmpi(tagStr, 'MaxSize', 7))        || ...
            (strncmpi(tagStr, 'Limit', 5))          || ...
            (strncmpi(tagStr, 'Default', 7))        || ...
            (strncmpi(tagStr, 'InFile', 6))         || ...
            (strncmpi(tagStr, 'Context', 7))        || ...
            (strncmpi(tagStr, 'Dll', 3))            || ...
            (strncmpi(tagStr, 'Class', 5))          || ...
            (strncmpi(tagStr, 'Comment', 7))        || ...
            (strncmpi(tagStr, 'Label', 5))        || ...
            (strncmpi(tagStr, 'Tooltip', 7))        || ...
            (strncmpi(tagStr, 'Visible', 7))        || ...
            (strncmpi(tagStr, 'Unit', 4)))
            % acknowledge that we know about these tags
        else
%             fprintf('parse_xprot(): WARNING: found unknown modifier %s\n',tagStr);
        end
        
        % remove the line containing the tag
        lstart = strfind(stubStr, ['<' tagStr '>']);
        %lend = strfind(stubStr(lstart+1:end), char(10));
        lend = lstart + length(tagStr) + 1;
        if (lstart > 1)
            stubStr = [stubStr(1:(lstart-1)) stubStr((lend+1):end)];
        else
            stubStr = stubStr((lend+1):end);
        end
    end

    [tagStr, newtagType, remStr] = findNextTag(stubStr);
end

if (strncmpi(tagType, 'ParamChoice', 11))
    value = '';%getQuotString(stubStr);
    ok = true;
end;
if (strncmpi(tagType, 'ParamString', 11))
    value = getQuotString(stubStr);
    ok = true;
else
    stubStr = strrep(stubStr, char(10), ' '); % remove newlines
    if (strncmpi(tagType, 'ParamBool', 9))
        stubStr = strrep(stubStr, '"true"', '1');
        stubStr = strrep(stubStr, '"false"', '0');
        [value,ok] = str2num(stubStr); %#ok<ST2NM>
        if (length(stubStr) > 1)
            value = false;
        else
            value = (value ~= 0);
        end;
    elseif (strncmpi(tagType, 'ParamLong', 9))
        [value,ok] = str2num(stubStr); %#ok<ST2NM>
    elseif (strncmpi(tagType, 'ParamDouble', 11))
        [value,ok] = str2num(stubStr); %#ok<ST2NM>
        ind = strfind(stubStr,'.');
        if (length(ind) > 0)
            if (length(ind) == length(value) -1) % there is an indication of bitnumber
                value = value(2:end);
            end;
        end;
    end
end

if (~ok)
    %fprintf('WARNING: get_xprot_value failed: %s\n', stubStr);
    %disp(value);
end