function [sScanHeader,sLoopCounter] = ReadsScanHeader(fid)
    ulDMALength=fread(fid,1,'uint32');                     % 0-3   bit  0..24: DMA length [bytes]
%       bit     25: pack bit
%       bit 26..31:
%       pci_rx enable flags                   4 byte
    mask=bin2dec('00000000111111111111111111111111');           % New with VB13
    meas.ulDMALength=bitand(ulDMALength,mask);

    mask=bin2dec('00000001000000000000000000000000');           % New with VB13
    meas.ulPackBit=bitand(ulDMALength,mask);
    
    mask=bin2dec('11111110000000000000000000000000');           % New with VB13
    meas.ulPCI_RXenable=bitand(ulDMALength,mask);

    
    meas.lMeasUID=fread(fid,1,'int32');                         % 4   measurement user ID
    meas.ulScanCounter=fread(fid,1,'uint32');                   % 8  scan counter [1...]
    meas.ulTimeStamp=fread(fid,1,'uint32');                     % 12 time stamp [2.5 ms ticks since 00:00]
    meas.ulPMUTimeStamp=fread(fid,1,'uint32');                  % 16 PMU time stamp [2.5 ms ticks since last trigger]
    meas.ushSystemType=fread(fid,1,'uint16');                   % 20 System Type
    meas.ushPTabPosDelay=fread(fid,1,'uint16');                 % 22
    meas.lPTABPosX=fread(fid,1,'int32');                        % 24-27   absolute Position in 0.1mm
    meas.lPTABPosY=fread(fid,1,'int32');                        % 28   absolute Position in 0.1mm
    meas.lPTABPosZ=fread(fid,1,'int32');                        % 32   absolute Position in 0.1mm
    meas.ulReserved1=fread(fid,1,'uint32');                     % 36
    
    meas.ulEvalInfoMask1=fread(fid,1,'uint32');                % 40 evaluation info mask field
    meas.ulEvalInfoMask2=fread(fid,1,'uint32');                %  evaluation info mask field

    meas.ushSamplesInScan=fread(fid,1,'uint16');                % 48 # of samples acquired in scan
    meas.ushUsedChannels=fread(fid,1,'uint16');                 % 50 # of channels used in scan
%   LoopCounter                                       % 52 Loop Counter
        data.ushLine = fread(fid,1,'uint16');
        data.ushAcquisition = fread(fid,1,'uint16');
        data.ushSlice = fread(fid,1,'uint16');
        data.ushPartition = fread(fid,1,'uint16');
        data.ushEcho = fread(fid,1,'uint16');
        data.ushPhase = fread(fid,1,'uint16');
        data.ushRepetition = fread(fid,1,'uint16');
        data.ushSet = fread(fid,1,'uint16');
        data.ushSeg = fread(fid,1,'uint16');
        data.ushIda = fread(fid,1,'uint16');
        data.ushIdb = fread(fid,1,'uint16');
        data.ushIdc = fread(fid,1,'uint16');
        data.ushIdd = fread(fid,1,'uint16');
        data.ushIde = fread(fid,1,'uint16');
    fread(fid,4,'uchar');                                        % 80 Cut-off values
    meas.ushKSpaceCentreColumn=fread(fid,1,'uint16');           % 84 centre of echo
    meas.ushCoilSelect=fread(fid,1,'uint16');                   % 86 centre of echo
    meas.fReadOutOffcentre=fread(fid,1,'float32');              % 88 ReadOut offcenter value
    meas.ulTimeSinceLastRF=fread(fid,1,'uint32');               % 92 Sequence time stamp since last RF pulse
    meas.ushKSpaceCentreLineNo=fread(fid,1,'uint16');           % 96 number of K-space centre line
    meas.ushKSpaceCentrePartitionNo=fread(fid,1,'uint16');      % 98 number of K-space centre partition
    fread(fid,28,'char');                                       % 100 Slice Data
    fread(fid,4,'uint16');                                      % 128 free parameter for IceProgram
    fread(fid,4,'uint16');                                      % 136 free parameter
    fread(fid,44,'uint8');                                      % 144 Application parameter
    meas.ulCRC=fread(fid,1,'uint32');                           % 188 CRC32 Checksum

    sScanHeader = meas;
    sLoopCounter = data;
end
