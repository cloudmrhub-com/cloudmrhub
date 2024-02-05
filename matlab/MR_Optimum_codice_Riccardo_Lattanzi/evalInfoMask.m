% EVALINFOMASK - Analyze the evalInfoMask entry in a meas.out measurement data 
%                header. Works for VB11 (and possibly for VA21) software
%
% [mask, taglist] = evalInfoMask(value)
%
% Stephan.Kannengiesser@siemens.com
%
% -----------------------------------------------------------------------------
%   Copyright (C) Siemens AG 2003  All Rights Reserved. Confidential.
% -----------------------------------------------------------------------------

function [mask, taglist] = evalInfoMask(value)

% const MdhBitField MDH_ACQEND            ((unsigned long)0);
mask.MDH_ACQEND            = bitand(value(1), 2^0);

% const MdhBitField MDH_RTFEEDBACK        (1);
mask.MDH_RTFEEDBACK        = bitand(value(1), 2^1);

% const MdhBitField MDH_HPFEEDBACK        (2);
mask.MDH_HPFEEDBACK        = bitand(value(1), 2^2);

% const MdhBitField MDH_ONLINE            (3);
mask.MDH_ONLINE            = bitand(value(1), 2^3);

% const MdhBitField MDH_OFFLINE           (4);
mask.MDH_OFFLINE           = bitand(value(1), 2^4);

% const MdhBitField MDH_SYNCDATA          (5);       // readout contains synchroneous data
mask.MDH_SYNCDATA          = bitand(value(1), 2^5);

% GAP

% const MdhBitField MDH_LASTSCANINCONCAT  (8);       // Flag for last scan in concatination
mask.MDH_LASTSCANINCONCAT  = bitand(value(1), 2^8);

% GAP

% const MdhBitField MDH_RAWDATACORRECTION (10);      // Correct the rawadata with the rawdata correction factor
mask.MDH_RAWDATACORRECTION = bitand(value(1), 2^10);

% const MdhBitField MDH_LASTSCANINMEAS    (11);      // Flag for last scan in measurement
mask.MDH_LASTSCANINMEAS    = bitand(value(1), 2^11);

% const MdhBitField MDH_SCANSCALEFACTOR   (12);      // Flag for scan specific additional scale factor
mask.MDH_SCANSCALEFACTOR   = bitand(value(1), 2^12);

% const MdhBitField MDH_2NDHADAMARPULSE   (13);      // 2nd RF exitation of HADAMAR
mask.MDH_2NDHADAMARPULSE   = bitand(value(1), 2^13);

% const MdhBitField MDH_REFPHASESTABSCAN  (14);      // reference phase stabilization scan         
mask.MDH_REFPHASESTABSCAN  = bitand(value(1), 2^14);

% const MdhBitField MDH_PHASESTABSCAN     (15);      // phase stabilization scan
mask.MDH_PHASESTABSCAN     = bitand(value(1), 2^15);

% const MdhBitField MDH_D3FFT             (16);      // execute 3D FFT         
mask.MDH_D3FFT             = bitand(value(1), 2^16);

% const MdhBitField MDH_SIGNREV           (17);      // sign reversal
mask.MDH_SIGNREV           = bitand(value(1), 2^17);

% const MdhBitField MDH_PHASEFFT          (18);      // execute phase fft     
mask.MDH_PHASEFFT          = bitand(value(1), 2^18);

% const MdhBitField MDH_SWAPPED           (19);      // swapped phase/readout direction
mask.MDH_SWAPPED           = bitand(value(1), 2^19);

% const MdhBitField MDH_POSTSHAREDLINE    (20);      // shared line               
mask.MDH_POSTSHAREDLINE    = bitand(value(1), 2^20);

% const MdhBitField MDH_PHASCOR           (21);      // phase correction data    
mask.MDH_PHASCOR           = bitand(value(1), 2^21);

% const MdhBitField MDH_PATREFSCAN        (22);      // additonal scan for PAT reference line/partition
mask.MDH_PATREFSCAN        = bitand(value(1), 2^22);

% const MdhBitField MDH_PATREFANDIMASCAN  (23);      // additonal scan for PAT reference line/partition that is also used as image scan
mask.MDH_PATREFANDIMASCAN  = bitand(value(1), 2^23);

% const MdhBitField MDH_REFLECT           (24);      // reflect line              
mask.MDH_REFLECT           = bitand(value(1), 2^24);

% const MdhBitField MDH_NOISEADJSCAN      (25);      // noise adjust scan --> Not used in NUM4        
mask.MDH_NOISEADJSCAN      = bitand(value(1), 2^25);

% const MdhBitField MDH_SHARENOW          (26);      // all lines are acquired from the actual and previous e.g. phases
mask.MDH_SHARENOW          = bitand(value(1), 2^26);

% const MdhBitField MDH_LASTMEASUREDLINE  (27);      // indicates that the current line is the last measured line of all succeeding e.g. phases
mask.MDH_LASTMEASUREDLINE  = bitand(value(1), 2^27);

% const MdhBitField MDH_FIRSTSCANINSLICE  (28);      // indicates first scan in slice (needed for time stamps)
mask.MDH_FIRSTSCANINSLICE  = bitand(value(1), 2^28);

% const MdhBitField MDH_LASTSCANINSLICE   (29);      // indicates  last scan in slice (needed for time stamps)
mask.MDH_LASTSCANINSLICE   = bitand(value(1), 2^29);

% const MdhBitField MDH_TREFFECTIVEBEGIN  (30);      // indicates the begin time stamp for TReff (triggered measurement)
mask.MDH_TREFFECTIVEBEGIN  = bitand(value(1), 2^30);

% const MdhBitField MDH_TREFFECTIVEEND    (31);      // indicates the   end time stamp for TReff (triggered measurement)
mask.MDH_TREFFECTIVEEND    = bitand(value(1), 2^31);

% Second EvalInfoMask Value...
if length(value) > 1

% Gap

%const MdhBitField MDH_FIRST_SCAN_IN_BLADE       (40);  // Marks the first line of a blade
mask.MDH_FIRST_SCAN_IN_BLADE = bitand(value(2), 2^(40-32));

%const MdhBitField MDH_LAST_SCAN_IN_BLADE        (41);  // Marks the last line of a blade
mask.MDH_LAST_SCAN_IN_BLADE  = bitand(value(2), 2^(41-32));

%const MdhBitField MDH_LAST_BLADE_IN_TR          (42);  // Set for all lines of the last BLADE in each TR interval
mask.MDH_LAST_BLADE_IN_TR    = bitand(value(2), 2^(42-32));

% Gap

%const MdhBitField MDH_RETRO_LASTPHASE           (45);  // Marks the last phase in a heartbeat
mask.MDH_RETRO_LASTPHASE           = bitand(value(2), 2^(45-32));

%const MdhBitField MDH_RETRO_ENDOFMEAS           (46);  // Marks an ADC at the end of the measurement
mask.MDH_RETRO_ENDOFMEAS           = bitand(value(2), 2^(46-32));

%const MdhBitField MDH_RETRO_REPEATTHISHEARTBEAT (47);  // Repeat the current heartbeat when this bit is found
mask.MDH_RETRO_REPEATTHISHEARTBEAT = bitand(value(2), 2^(47-32));

%const MdhBitField MDH_RETRO_REPEATPREVHEARTBEAT (48);  // Repeat the previous heartbeat when this bit is found
mask.MDH_RETRO_REPEATPREVHEARTBEAT = bitand(value(2), 2^(48-32));

%const MdhBitField MDH_RETRO_ABORTSCANNOW        (49);  // Just abort everything
mask.MDH_RETRO_ABORTSCANNOW        = bitand(value(2), 2^(49-32));

%const MdhBitField MDH_RETRO_LASTHEARTBEAT       (50);  // This adc is from the last heartbeat (a dummy)
mask.MDH_RETRO_LASTHEARTBEAT       = bitand(value(2), 2^(50-32));

%const MdhBitField MDH_RETRO_DUMMYSCAN           (51);  // This adc is just a dummy scan, throw it away
mask.MDH_RETRO_DUMMYSCAN           = bitand(value(2), 2^(51-32));

%const MdhBitField MDH_RETRO_ARRDETDISABLED      (52);  // Disable all arrhythmia detection when this bit is found
mask.MDH_RETRO_ARRDETDISABLED      = bitand(value(2), 2^(52-32));

end % if length(value) > 1


% Make a list of tags

%taglist = '';
taglist = [];

if mask.MDH_ACQEND            
  taglist = [taglist; ' MDH_ACQEND                   '];
end
if mask.MDH_RTFEEDBACK        
  taglist = [taglist; ' MDH_RTFEEDBACK               '];
end
if mask.MDH_HPFEEDBACK        
  taglist = [taglist; ' MDH_HPFEEDBACK               '];
end
if mask.MDH_ONLINE            
  taglist = [taglist; ' MDH_ONLINE                   '];
end
if mask.MDH_OFFLINE           
  taglist = [taglist; ' MDH_OFFLINE                  '];
end
if mask.MDH_SYNCDATA
  taglist = [taglist; ' MDH_SYNCDATA                 '];
end
if mask.MDH_LASTSCANINCONCAT  
  taglist = [taglist; ' MDH_LASTSCANINCONCAT         '];
end
if mask.MDH_RAWDATACORRECTION 
  taglist = [taglist; ' MDH_RAWDATACORRECTION        '];
end
if mask.MDH_LASTSCANINMEAS    
  taglist = [taglist; ' MDH_LASTSCANINMEAS           '];
end
if mask.MDH_SCANSCALEFACTOR   
  taglist = [taglist; ' MDH_SCANSCALEFACTOR          '];
end
if mask.MDH_2NDHADAMARPULSE   
  taglist = [taglist; ' MDH_2NDHADAMARPULSE          '];
end
if mask.MDH_REFPHASESTABSCAN  
  taglist = [taglist; ' MDH_REFPHASESTABSCAN         '];
end
if mask.MDH_PHASESTABSCAN     
  taglist = [taglist; ' MDH_PHASESTABSCAN            '];
end
if mask.MDH_D3FFT             
  taglist = [taglist; ' MDH_D3FFT                    '];
end
if mask.MDH_SIGNREV           
  taglist = [taglist; ' MDH_SIGNREV                  '];
end
if mask.MDH_PHASEFFT          
  taglist = [taglist; ' MDH_PHASEFFT                 '];
end
if mask.MDH_SWAPPED           
  taglist = [taglist; ' MDH_SWAPPED                  '];
end
if mask.MDH_POSTSHAREDLINE    
  taglist = [taglist; ' MDH_POSTSHAREDLINE           '];
end
if mask.MDH_PHASCOR           
  taglist = [taglist; ' MDH_PHASCOR                  '];
end
if mask.MDH_PATREFSCAN        
  taglist = [taglist; ' MDH_PATREFSCAN               '];
end
if mask.MDH_PATREFANDIMASCAN  
  taglist = [taglist; ' MDH_PATREFANDIMASCAN         '];
end
if mask.MDH_REFLECT           
  taglist = [taglist; ' MDH_REFLECT                  '];
end
if mask.MDH_NOISEADJSCAN      
  taglist = [taglist; ' MDH_NOISEADJSCAN             '];
end
if mask.MDH_SHARENOW          
  taglist = [taglist; ' MDH_SHARENOW                 '];
end
if mask.MDH_LASTMEASUREDLINE  
  taglist = [taglist; ' MDH_LASTMEASUREDLINE         '];
end
if mask.MDH_FIRSTSCANINSLICE  
  taglist = [taglist; ' MDH_FIRSTSCANINSLICE         '];
end
if mask.MDH_LASTSCANINSLICE   
  taglist = [taglist; ' MDH_LASTSCANINSLICE          '];
end
if mask.MDH_TREFFECTIVEBEGIN  
  taglist = [taglist; ' MDH_TREFFECTIVEBEGIN         '];
end
if mask.MDH_TREFFECTIVEEND    
  taglist = [taglist; ' MDH_TREFFECTIVEEND           '];
end


if length(value) > 1


if mask.MDH_FIRST_SCAN_IN_BLADE
  taglist = [taglist; ' MDH_FIRST_SCAN_IN_BLADE      '];
end

if mask.MDH_LAST_SCAN_IN_BLADE
  taglist = [taglist; ' MDH_LAST_SCAN_IN_BLADE       '];
end

if mask.MDH_LAST_BLADE_IN_TR
  taglist = [taglist; ' MDH_LAST_BLADE_IN_TR         '];
end


if mask.MDH_RETRO_LASTPHASE
  taglist = [taglist; ' MDH_RETRO_LASTPHASE          '];
end

if mask.MDH_RETRO_ENDOFMEAS
  taglist = [taglist; ' MDH_RETRO_ENDOFMEAS          '];
end

if mask.MDH_RETRO_REPEATTHISHEARTBEAT
  taglist = [taglist; ' MDH_RETRO_REPEATTHISHEARTBEAT'];
end

if mask.MDH_RETRO_REPEATPREVHEARTBEAT
  taglist = [taglist; ' MDH_RETRO_REPEATPREVHEARTBEAT'];
end

if mask.MDH_RETRO_ABORTSCANNOW
  taglist = [taglist; ' MDH_RETRO_ABORTSCANNOW       '];
end

if mask.MDH_RETRO_LASTHEARTBEAT
  taglist = [taglist; ' MDH_RETRO_LASTHEARTBEAT      '];
end

if mask.MDH_RETRO_DUMMYSCAN
  taglist = [taglist; ' MDH_RETRO_DUMMYSCAN          '];
end

if mask.MDH_RETRO_ARRDETDISABLED
  taglist = [taglist; ' MDH_RETRO_ARRDETDISABLED     '];
end


end % if length(value) > 1
