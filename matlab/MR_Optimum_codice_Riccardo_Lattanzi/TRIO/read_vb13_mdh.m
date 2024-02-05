% READ_VB13_MDH - Read a meas.asc measurement data header from an open file
%                 Works for VB13 software version
%
% [header, ok, nBytesRead] = read_vb13_mdh(fid)
%
% header: struct with entries
% - lMeasUID
% - ulScanCounter
% - ulTimeStamp
% - ulPMUTimeStamp
% - aulEvalInfoMask
% - ushSamplesInScan
% - ushUsedChannels
% - sLC
% - sCutOff
% - ushKSpaceCentreColumn
% - ushCoilSelect
% - fReadOutOffcentre
% - ulTimeSinceLastRF
% - ushKSpaceCentreLineNo
% - ushKSpaceCentrePartitionNo
% - aushIceProgramPara
% - aushFreePara
% - sSD
% - ushChannelId
% - ushPTABPosNeg
%
% ok: true if header read was successful, false if not
%
% nBytesRead: number of bytes actually read; full mdh length if successful, potentially less if not
%
% fid: file ID for reading
%
% Stephan.Kannengiesser@siemens.com
%
% -----------------------------------------------------------------------------
%   Copyright (C) Siemens AG 2003  All Rights Reserved. Confidential.
% -----------------------------------------------------------------------------

function [header, ok] = read_vb13_mdh(fid, nBytesRead)

b                    = 0;

ok                   = true;

% typedef struct
% {
%  PACKED_MEMBER( uint32_t,     ulFlagsAndDMALength           );    // bit  0..24: DMA length [bytes]
%                                                                   // bit     25: pack bit
%                                                                   // bit 26..31: pci_rx enable flags                   4 byte
[header.ulFlagsAndDMALength, count] = fread(fid, 1,      'ulong');  
if count == 1
    %disp('* read ok')
    b = b + 4;
else
    %disp('* read NOT ok')
    ok = false;
    return
end

%  PACKED_MEMBER( int32_t,      lMeasUID                      );    // measurement user ID                               4
[header.lMeasUID, count]            = fread(fid, 1,      'long');
if count == 1
    b = b + 4;
else
    ok = false;
    return
end

%  PACKED_MEMBER( uint32_t,     ulScanCounter                 );    // scan counter [1...]                               4
[header.ulScanCounter, count]       = fread(fid, 1,      'ulong');
if count == 1
    b = b + 4;
else
    ok = false;
    return
end

%  PACKED_MEMBER( uint32_t,     ulTimeStamp                   );    // time stamp [2.5 ms ticks since 00:00]             4
[header.ulTimeStamp, count]         = fread(fid, 1,      'ulong');
if count == 1
    b = b + 4;
else
    ok = false;
    return
end

%  PACKED_MEMBER( uint32_t,     ulPMUTimeStamp                );    // PMU time stamp [2.5 ms ticks since last trigger]  4
[header.ulPMUTimeStamp, count]      = fread(fid, 1,      'ulong');
if count == 1
    b = b + 4;
else
    ok = false;
    return
end

%  PACKED_MEMBER( uint32_t,     aulEvalInfoMask[MDH_NUMBEROFEVALINFOMASK] ); // evaluation info mask field           8
[header.aulEvalInfoMask, count]     = fread(fid, [1 2],  'ulong');
if count == 2
    b = b + 8;
else
    ok = false;
    return
end

%  PACKED_MEMBER( uint16_t,     ushSamplesInScan              );    // # of samples acquired in scan                     2
[header.ushSamplesInScan, count]    = fread(fid, 1,      'ushort');
if count == 1
    b = b + 2;
else
    ok = false;
    return
end

%  PACKED_MEMBER( uint16_t,     ushUsedChannels               );    // # of channels used in scan                        2   =32
[header.ushUsedChannels, count]     = fread(fid, 1,      'ushort');
if count == 1
    b = b + 2;
else
    ok = false;
    return
end

%  PACKED_MEMBER( sLoopCounter, sLC                           );    // loop counters                                    28   =60
%typedef struct {
%  PACKED_MEMBER( uint16_t,  ushLine         ); /* line index                   */
%  PACKED_MEMBER( uint16_t,  ushAcquisition  ); /* acquisition index            */
%  PACKED_MEMBER( uint16_t,  ushSlice        ); /* slice index                  */
%  PACKED_MEMBER( uint16_t,  ushPartition    ); /* partition index              */
%  PACKED_MEMBER( uint16_t,  ushEcho         ); /* echo index                   */
%  PACKED_MEMBER( uint16_t,  ushPhase        ); /* phase index                  */
%  PACKED_MEMBER( uint16_t,  ushRepetition   ); /* measurement repeat index     */
%  PACKED_MEMBER( uint16_t,  ushSet          ); /* set index                    */
%  PACKED_MEMBER( uint16_t,  ushSeg          ); /* segment index  (for TSE)     */
%  PACKED_MEMBER( uint16_t,  ushIda          ); /* IceDimension a index         */
%  PACKED_MEMBER( uint16_t,  ushIdb          ); /* IceDimension b index         */
%  PACKED_MEMBER( uint16_t,  ushIdc          ); /* IceDimension c index         */
%  PACKED_MEMBER( uint16_t,  ushIdd          ); /* IceDimension d index         */
%  PACKED_MEMBER( uint16_t,  ushIde          ); /* IceDimension e index         */
%} sLoopCounter;                                /* sizeof : 28 byte             */
[header.sLC, count]                        = fread(fid, [1 14], 'ushort');
if count == 14
    b = b + 28;
else
    ok = false;
    return
end

%  PACKED_MEMBER( sCutOffData,  sCutOff                       );    // cut-off values                                    4
%typedef struct{
%  PACKED_MEMBER( uint16_t,  ushPre          );    /* write ushPre zeros at line start */
%  PACKED_MEMBER( uint16_t,  ushPost         );    /* write ushPost zeros at line end  */
% } sCutOffData;
[header.sCutOff, count]                    = fread(fid, [1 2],  'ushort');
if count == 2
    b = b + 4;
else
    ok = false;
    return
end

%  PACKED_MEMBER( uint16_t,     ushKSpaceCentreColumn         );    // centre of echo                                    2
[header.ushKSpaceCentreColumn, count]      = fread(fid, 1,      'ushort');
if count == 1
    b = b + 2;
else
    ok = false;
    return
end

%  PACKED_MEMBER( uint16_t,     ushCoilSelect                 );    // Bit 0..3: CoilSelect                              2
[header.ushCoilSelect, count]              = fread(fid, 1,      'ushort');
if count == 1
    b = b + 2;
else
    ok = false;
    return
end

%  PACKED_MEMBER( float,        fReadOutOffcentre             );    // ReadOut offcenter value                           4
[header.fReadOutOffcentre, count]          = fread(fid, 1,      'float');
if count == 1
    b = b + 4;
else
    ok = false;
    return
end

%  PACKED_MEMBER( uint32_t,     ulTimeSinceLastRF             );    // Sequence time stamp since last RF pulse           4
[header.ulTimeSinceLastRF, count]          = fread(fid, 1,      'ulong');
if count == 1
    b = b + 4;
else
    ok = false;
    return
end

%  PACKED_MEMBER( uint16_t,     ushKSpaceCentreLineNo         );    // number of K-space centre line                     2
[header.ushKSpaceCentreLineNo, count]      = fread(fid, 1,      'ushort');
if count == 1
    b = b + 2;
else
    ok = false;
    return
end

%  PACKED_MEMBER( uint16_t,     ushKSpaceCentrePartitionNo    );    // number of K-space centre partition                2
[header.ushKSpaceCentrePartitionNo, count] = fread(fid, 1,      'ushort');
if count == 1
    b = b + 2;
else
    ok = false;
    return
end

%  PACKED_MEMBER( uint16_t,     aushIceProgramPara[MDH_NUMBEROFICEPROGRAMPARA] ); // free parameter for IceProgram   8   =88
[header.aushIceProgramPara, count]         = fread(fid, [1 4],  'ushort');
if count == 4
    b = b + 8;
else
    ok = false;
    return
end

%  PACKED_MEMBER( uint16_t,     aushFreePara[MDH_FREEHDRPARA] );    // free parameter                          4 * 2 =   8
[header.aushFreePara, count]               = fread(fid, [1 4],  'ushort');
if count == 4
    b = b + 8;
else
    ok = false;
    return
end

%  PACKED_MEMBER( sSliceData,   sSD                           );    // Slice Data                                       28   =124
%typedef struct {
%  PACKED_MEMBER( float,  flSag          );
%  PACKED_MEMBER( float,  flCor          );
%  PACKED_MEMBER( float,  flTra          );
%} sVector;
%typedef struct{
%  PACKED_MEMBER( sVector,         sSlicePosVec     ); /* slice position vector        */
%  PACKED_MEMBER( float,           aflQuaternion[4] ); /* rotation matrix as quaternion*/
%} sSliceData;                                         /* sizeof : 28 byte             */
[header.sSD, count]                        = fread(fid, [1 7],  'float');
if count == 7
    b = b + 28;
else
    ok = false;
    return
end

%  PACKED_MEMBER( uint16_t,	   ushChannelId                  );    // channel Id must be the last parameter             2
[header.ushChannelId, count]               = fread(fid, 1,      'ushort');
if count == 1
    b = b + 2;
else
    ok = false;
    return
end

%  PACKED_MEMBER( uint16_t,	   ushPTABPosNeg                 );    // negative, absolute PTAB position in [0.1 mm]      2
[header.ushPTABPosNeg, count]              = fread(fid, 1,      'ushort');
if count == 1
    b = b + 2;
else
    ok = false;
    return
end

% } sMDH;                                        // total length: 32 * 32 Bit (128 Byte)            128
% 
% #endif   /* MDH_H */

%disp(['* rest bytes ' num2str(128-b)])

dummy                         = fread(fid, [1 128-b], 'char');
