function [tagStr, tagType, remStr] = findNextTag(inStr)
% returns <tag> name and the remainder of the string following the tag.
% for e.g. <Tag>, returns tagStr='Tag', tagType=''
% for e.g. <ParamLong."Tag">, returns tagStr='Tag', tagType='ParamLong'
% if no tag is found, returns null strings

tagStr = [];
tagType = [];
remStr = [];

startPos = strfind(inStr,'<'); % look for start of tag
if (startPos)
    endPos = strfind(inStr,'>'); % look for end of tag
    if (endPos)
        % found complete tag
        fullTag = inStr(startPos+1:endPos-1);
        
        % now check for name/type
        dotPos = strfind(fullTag,'."');
        if (dotPos)
            tagStr = getQuotString(fullTag(dotPos+1:end));
            tagType = fullTag(1:dotPos-1);
        else
            tagStr = fullTag;
        end
        
        % return remainder
        remStr = inStr(endPos+1:end);
    end
end
