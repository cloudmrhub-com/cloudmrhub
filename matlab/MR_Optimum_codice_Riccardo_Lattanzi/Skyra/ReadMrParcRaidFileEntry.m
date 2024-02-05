function MrParcRaidFileEntry=ReadMrParcRaidFileEntry(fid)

        meas.ulmeasID = fread(fid,1,'uint32');     % number or meas in file
        meas.ulfileID = fread(fid,1,'uint32') ;    % number or meas in file
        meas.uloff = fread(fid,1,'uint64')     ;   % number or meas in file
        meas.ullen = fread(fid,1,'uint64')      ;  % number or meas in file
%        meas.ch1 = (fread(fid,[1,64],'*uchar'));        % number or meas in file
        meas.ch1 = char(fread(fid,[1,64],'uchar'));        % number or meas in file
        meas.ch2 = char(fread(fid,[1,64],'uchar')) ;       % number or meas in file
        MrParcRaidFileEntry = meas;
end