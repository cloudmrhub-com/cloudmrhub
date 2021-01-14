function stvar = extractBraceString(text)
% extracts string from within curly braces, including nested braces

% tlen = length(text);
% tstart = strfind(text,'{');
% 
% stvar = [];
% if (tstart)
%     tnests = 1;
%     for x=tstart+1:tlen
%         if (text(x) == '{')
%             tnests = tnests + 1;
%         elseif (text(x) == '}')
%             tnests = tnests - 1;
%         end
%         if (tnests == 0)
%             stvar = text(tstart+1:x-1);
%             break;
%         end
%     end
% end

tstart = strfind(text,'{');
tend = strfind(text,'}');

stvar = [];
if (~isempty(tstart) & ~isempty(tend))
    [brackind,ind] = sort([tstart,tend]);
    pos = [ones(1,length(tstart)),-1*ones(1,length(tend))];
    endind = find(cumsum(pos(ind))==0,1);
    if (~isempty(endind))
        endpos = brackind(endind);
        stvar = text(tstart(1)+1:endpos-1);
    end;
end;