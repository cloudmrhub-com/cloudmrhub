function mrprot = parse_loop(mrprot, workarr, tagName, level)
% this function is called recursively to parse the xprotocol

[tagStr, tagType, workarr] = findNextTag(workarr);
while (~isempty(workarr))
    % for parammap, add the name as another level of the current array name
    % and spawn off another copy of this function to deal with them
    if (strncmpi(tagType, 'ParamMap', 8))
        if (~isempty(tagStr))
            level = level + 1;
            tagName{level} = make_safe_fieldname(tagStr);
        else
            level = level + 1;
            if (level == 1)
                tagName{level} = 'x';
            else
                tagName{level} = tagName{level-1};
            end;
        end;
        stubArr = extractBraceString(workarr);
        mrprot = parse_loop(mrprot, stubArr, tagName, level);
        workarr = workarr(length(stubArr)+1:end);
        %if (~isempty(tagStr))
            level = level - 1;
        %end

    % for paramarray, special processing is required
    % for now, though, just skip them since the most important ones are
    % still in the text mrprot and it would be a lot more work to process
    % them here
    elseif (strncmpi(tagType, 'ParamArray', 10))
        if (~isempty(tagStr))
            level = level + 1;
            tagName{level} = make_safe_fieldname(tagStr);
        end
        %fprintf('found array %s\n',build_tag_name(tagName,level+1));
        %fprintf('found array %s\n',tagStr);
        stubArr = extractBraceString(workarr);
%         workarr = workarr(length(stubArr)+1:end);

%         %skip the <default> tag, which always seems to be there
        [tagStr2, tagType, stubArr] = findNextTag(stubArr);
%         if (strcmpi(tagStr2,'MinSize') || strcmpi(tagStr2,'MaxSize'))
%             stubArr = extractBraceString(workarr);
%             temparr = workarr(length(stubArr)+1:end);
%             ind = findstr(temparr,'}');
%             workarr = temparr((ind(1)+1):end);
%         else
            if (~strcmpi(tagStr2,'Default'))
                %fprintf('parse_xprot(): unknown ParamArray format!')
                stubArr = extractBraceString(workarr);
            end
            % now map the fields
    %         [tagStr, tagType, stubArr] = findNextTag(stubArr);

            mrprot = parse_loop(mrprot, stubArr, tagName, level);

            temparr = workarr;
            % walk down the string along the structure
            structpos = ['mrprot'];
            for i = 1:level
                structpos = [structpos '.' tagName{i}];
            end;
            [temparr,~,count] = walkdown(mrprot,structpos,temparr,0,0);

            if (length(strfind(temparr,'{')) >= count )
                % walk down the structure again, along the values of the stucture, and fill the values in
                [temparr,mrprot] = walkdown(mrprot,structpos,temparr,1,0);
            end;


            workarr = workarr(length(stubArr)+1:end);
%         end;
        
%         remArr = extractBraceString(stubArr);
%         stubArr = stubArr(length(remArr)+1:end);

        if (~isempty(tagStr))
            level = level - 1;
        end

    elseif (strncmpi(tagType, 'PipeService', 11))
        if (~isempty(tagStr))
            level = level + 1;
            tagName{level} = make_safe_fieldname(tagStr);
        end
        stubArr = extractBraceString(workarr);
        
        %skip the <Class> tag, which always seems to be there
        [tagStr, tagType, stubArr] = findNextTag(stubArr);
        if (~strcmpi(tagStr,'Class'))
            fprintf('parse_xprot(): unknown PipeService format!')
            stubArr = extractBraceString(workarr);
        end
        
        mrprot = parse_loop(mrprot, stubArr, tagName, level);
        workarr = workarr(length(stubArr)+1:end);
        if (~isempty(tagStr))
            level = level - 1;
        end
        
    elseif (strncmpi(tagType, 'Pipe', 4))
        if (~isempty(tagStr))
            level = level + 1;
            tagName{level} = make_safe_fieldname(tagStr);
        end
        stubArr = extractBraceString(workarr);
        mrprot = parse_loop(mrprot, stubArr, tagName, level);
        workarr = workarr(length(stubArr)+1:end);
        if (~isempty(tagStr))
            level = level - 1;
        end
        
    elseif (strncmpi(tagType, 'ParamFunctor', 12))
        if (~isempty(tagStr))
            level = level + 1;
            tagName{level} = make_safe_fieldname(tagStr);
        end
        stubArr = extractBraceString(workarr);
        
        %skip the <Class> tag, which always seems to be there
        [tagStr, tagType, stubArr] = findNextTag(stubArr);
        if (~strcmpi(tagStr,'Class'))
            fprintf('parse_xprot(): unknown PipeService format!')
            stubArr = extractBraceString(workarr);
        end
        
        mrprot = parse_loop(mrprot, stubArr, tagName, level);
        workarr = workarr(length(stubArr)+1:end);
        if (~isempty(tagStr))
            level = level - 1;
        end
        
    % these are values that we know how to process, so do so
    elseif ((strncmpi(tagType, 'ParamBool', 9))     || ...
            (strncmpi(tagType, 'ParamLong', 9))     || ...
            (strncmpi(tagType, 'ParamDouble', 11))  || ...
            (strncmpi(tagType, 'ParamString', 11))  || ...
            (strncmpi(tagType, 'ParamChoice', 11)))
        if (isempty(tagStr))
            tagStr = 'x';
        end;
        tagName{level+1} = make_safe_fieldname(tagStr);
        ind = strfind(workarr,'}');
        ind2 = strfind(strtrim(workarr((ind(1)+1):end)),'{');
        if (length(ind2) == 0)
            ind2(1) = 0;
        end;
        if (ind2(1) == 1)
            value = [];
            ind = strfind(workarr,'{');
            l = 1;
            for i = 1:length(ind)
                stubStr = extractBraceString(workarr(ind(i):end));
                %fprintf('%s = %s\n',tagName{level+1},stubStr);
                temp = get_xprot_value(stubStr, tagType);  
                if (~isempty(temp) && length(temp) == 1)
                    value(l) = temp;
                    l = l+1;
                end;
            end;
        else
            stubStr = extractBraceString(workarr);
            %fprintf('%s = %s\n',tagName{level+1},stubStr);
            value = get_xprot_value(stubStr, tagType);
        end;
        fields = {tagName{1:level+1}, value};
        mrprot = setfield(mrprot, fields{:});
        workarr = workarr(length(stubStr)+1:end);

    % we don't care about the things below, but acknowledge that we know
    % about them
    elseif (...%(strncmpi(tagType, 'ParamChoice', 11))      || ...
            (strncmpi(tagType, 'Class', 5))             || ...
            (strncmpi(tagType, 'Connection', 10))       || ...
            (strncmpi(tagType, 'Event', 5))             || ...
            (strncmpi(tagType, 'ParamFunctor', 12))     || ...
            (strncmpi(tagType, 'Method', 6))            || ...
            (strncmpi(tagType, 'ProtocolComposer', 16)) || ...
            (strncmpi(tagType, 'Dependency', 10))       || ...
            (strncmpi(tagType, 'ParamCardLayout', 15))  || ...
            (strncmpi(tagType, 'EVACardLayout', 13)))
        stubStr = extractBraceString(workarr);
        workarr = workarr(length(stubStr)+1:end);

        % ignore these also
    elseif ((strncmpi(tagStr, 'Name', 4))             || ...
            (strncmpi(tagStr, 'ID', 2))               || ...
            (strncmpi(tagStr, 'Comment', 7))          || ...
            (strncmpi(tagStr, 'Label', 5))            || ...
            (strncmpi(tagStr, 'Visible', 7))          || ...
            (strncmpi(tagStr, 'Userversion', 11)))
        % skip to end of line for these
        lend = strfind(workarr, char(10));
        workarr = workarr(lend+1:end);
    elseif (strncmpi(tagStr, 'EVAStringTable', 14))
        % skip brace string for this one
        stubStr = extractBraceString(workarr);
        workarr = workarr(length(stubStr)+1:end);
        
    % something unknown has happened if we reach this point, so throw a warning
    else
        if (~isempty(tagType))
%             fprintf('parse_xprot(): WARNING: found unknown tag %s (%s)\n',tagStr,tagType);
        end;
    end

    [tagStr, tagType, workarr] = findNextTag(workarr);
end