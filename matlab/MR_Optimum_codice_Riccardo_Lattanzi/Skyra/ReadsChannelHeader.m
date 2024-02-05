function sChannelHeader = ReadsChannelHeader(fid)

    ulTypeandChannelLength=fread(fid,1,'uint32');     
    % 0 : (1)  8bit type (0x02 => ChannelHeader)
                                                           
    %     (3) 24bit channel length (header+data) in byte
    meas.ulType=bitand(ulTypeandChannelLength,hex2dec('000000FF'));
    meas.ulChannelLength=bitand(ulTypeandChannelLength,hex2dec('FFFFFF00'));                                                     
    meas.lMeasUID=fread(fid,1,'int32');            % 4 line index                   
    meas.ulScanCounter=fread(fid,1,'uint32');      % 8  scan counter [1...] 
    meas.ulReserved1=fread(fid,1,'uint32');        % 12
    meas.ulUnused1=fread(fid,1,'uint32');          % 16
    meas.ulUnused2=fread(fid,1,'uint32');          % 20
    meas.ulChannelId=fread(fid,1,'uint16');        % 24
    meas.ulUnused3=fread(fid,1,'uint16');          % 26
    meas.ulCRC=fread(fid,1,'uint32');              % 28 CRC32 Checksum

    sChannelHeader = meas; 
end