% READ_MEAS_VB13 - read meas.out file along all dimensions for software version VB13
%
% function [raw noise ref phasecor pc_order centerlines lastheader rest] = ...
%           read_meas_vb13(filename [,mode] [,undoOS] [,repWise] [,average] [, extRefOnly])
%
%   raw             : rawdata
%
%   noise           : noise adjust scan (for iPAT measurements)
%
%   ref             : iPAT reference lines (*currently not filled*)
%
%   phasecor        : TSE / EPI phase correction scans (*not thoroughly
%                       tested, use with care!*)
%
%   pc_order        : experimental output
%
%   centerlines     : line index / indices of k-space center lines
%                       (may contain an additional 1 for TSE due to refScan)
%
%   lastheader      : last measurement data header as read by this function
%                       (for writing / duplicating meas.out)
%
%   rest            : binary data for the rest of the meas.out file
%                       (for writing / duplicating meas.out)
%
%   modes
%   - 'STD'         : all loop counters are used
%   - 'TSE'         : segments ignored since redundant with lines
%   - 'DRY'         : only print max loop counters
%
%   undoOS          : 0 = preserve ADC oversampling (2x); 1 = remove
%
%   repWise         : 1 = put multiple measurements (repetitions) into
%                       separate files (for large data sets)
%
%   average         : 0 = preserve averages in their own matrix dimension
%                     1 = accumulate raw data across AVE dimension
%
%   extRefOnly      : 1 = extract only "EXTRA" reference lines for iPAT
%
%
% Stephan.Kannengiesser@siemens.com
%
% -----------------------------------------------------------------------------
%   Copyright (C) Siemens AG 2006  All Rights Reserved. Confidential.
% -----------------------------------------------------------------------------
%
function [raw, noise, ref, phasecor, pc_order, centerlines, lastheader, rest, phasestab, ch_id_track] = ...
    read_meas_vb13(filename, mode, undoOS, repWise, average, extRefOnly)

% mode as default
if nargin < 2
  mode = 'TSE';
end
if nargin < 3
%   undoOS = 1;
  undoOS = 0; %RL Jan/8/2010
end
if nargin < 4
  repWise = 0;
end
if nargin < 5
  average = 1;
end
if nargin < 6
  extRefOnly = 0;
end

% first pass: open file and find out about dimensions
%frontofFile = 32;

maxLC  = zeros([1 14]);
maxCOL = 0;
maxCHA = 0;

centerlines = zeros([1 0]);
centerpars  = zeros([1 0]);

fid         = fopen(filename, 'rb');
% First long contains the size of the file header
frontofFile = fread(fid, [1 1], 'int32'); % NOT long ?
disp(['* Front of file reads length ' num2str(frontofFile)])
fHeadRest   = fread(fid, [1 frontofFile-4], 'char');
readMore    = 1;
nScansRead  = 0;
while readMore
  nScansRead = nScansRead + 1;
  header   = read_vb13_mdh(fid);
  mask     = evalInfoMask(header.aulEvalInfoMask(1));
  readMore = ~mask.MDH_ACQEND;
  % RL Jan/6/2010: need to keep track of the used channels to trim the raw data
  if maxCHA == 0
      used_ch = header.ushUsedChannels;
      ch_id_track = zeros(1,used_ch);
      track_index = 1;
  else
      track_index = track_index + 1;
  end
  
  if readMore

    % just dummy read the data!
    a      = fread(fid, [2 header.ushSamplesInScan], 'float');
    if any(size(a) ~= [2 header.ushSamplesInScan])
      disp('*** ERROR reading from file')
      readMore = 0;
    end % if

    if extRefOnly & ~mask.MDH_PATREFSCAN
      continue
    end
    
    % Update max. indices
    maxLC  = max([header.sLC ; maxLC]);
    maxCOL = max(maxCOL, header.ushSamplesInScan);
    maxCHA = max(maxCHA, header.ushChannelId);
    
    % RL Jan/6/2010: need to keep track of the used channels to trim the raw data
    if track_index <= used_ch
        ch_id_track(track_index) = header.ushChannelId + 1; %RL this was causing the channel shift
%         ch_id_track(track_index) = header.ushChannelId;
    end
    
    %centerline = header.ushKSpaceCentreLineNo;
    %disp(centerline)
    
    if (length(centerlines) == 0) | ...
          ~any(centerlines == repmat(header.ushKSpaceCentreLineNo+1, ...
                                     size(centerlines)))
      centerlines = [centerlines header.ushKSpaceCentreLineNo+1];
    end
    
    if (length(centerpars) == 0) | ...
          ~any(centerpars == repmat(header.ushKSpaceCentrePartitionNo+1, ...
                                     size(centerpars)))
      centerpars = [centerpars header.ushKSpaceCentrePartitionNo+1];
    end
    
  end % readMore
end % while

disp([num2str(nScansRead) ' scans read'])

lastheader = header;

% remember dimensions
dimNames = ['LIN'; 'ACQ'; 'SLC'; 'PAR'; 'ECO'; 'PHS'; 'REP'; 'SET'; 'SEG'; ...
            'IDA'; 'IDB'; 'IDC'; 'IDD'; 'IDE'];
dimNames = ['COL'; 'CHA'; dimNames];

%size(dimNames)
%size(maxLC)

disp('* sizes:')
disp([dimNames repmat(': ',[16 1]) num2str([maxCOL maxCHA+1 maxLC+1].')])
centerlines
centerpars

if mode == 'DRY'
  
  % Try to read the rest of the file
  [rest count] = fread(fid, [1 inf], 'char');

  disp(['* read ' num2str(count) ' extra bytes after last header'])
  
  fclose(fid);
  
  return
  
else
  fclose(fid);
  rest = NaN; % RL Jan/05/2010
end % if DRY

if mode == 'EPI'
  disp('- EPI: ignoring segments and sets')
  maxLC(9) = 0;
  maxLC(8) = 0;
end
if mode == 'TSE'
  disp('- TSE: ignoring segments')
  maxLC(9) = 0;
end
if average
  disp('- average: return only averaged data')
  maxLC(2) = 0;
end
if repWise
  disp(['- repWise: returning only last repetition, writing 1 to ' ...
        num2str(repWise) ' to disk!'])
  maxLC(7) = 0;
end
currRepCtr = 0;

% Now read for real (flags???)
if undoOS
  raw_size = [maxCOL/2 maxCHA+1 maxLC+1];
else
  raw_size = [maxCOL maxCHA+1 maxLC+1];
end
nCmplx = prod(raw_size);
disp(['* trying to allocate MBytes: ' ...
      num2str(16*nCmplx/1024/1024)])

% NOTE: some loop counters are not distinct labels, but
%       provide additional information without requiring
%       more storage, e.g., segments!
try
  raw      = zeros(raw_size);
  noise    = zeros([maxCOL 0]);
  ref      = zeros([maxCOL 0]);
  phasecor = zeros([maxCOL 0]);
  phasestab = zeros([maxCOL 0]);
  pc_order = zeros([1      0]); % slc/par
catch
  disp('*** Error allocating memory - aborting')
  return
end

fid      = fopen(filename, 'rb');
fHead    = fread(fid, [1 frontofFile], 'char');
readMore = 1;
while readMore
  header         = read_vb13_mdh(fid);
  [mask taglist] = evalInfoMask(header.aulEvalInfoMask(1));
  readMore       = ~mask.MDH_ACQEND;
  
  % check for new repetition
  if repWise & (header.sLC(7) ~= currRepCtr)
    disp(['- dumping repetition ' num2str(currRepCtr+1)])
    save(sprintf('raw_rep%05d',currRepCtr+1), 'raw');
    currRepCtr = header.sLC(7);

    % External references usually come first
    if extRefOnly
      break
    end
    
    % New: we might only want the first repWise measurements
    if currRepCtr >= repWise
      break
    end

    raw = zeros(raw_size);
    
  end % if repWise & (header.sLC(7) ~= currRepCtr)
  
  if readMore
    a      = fread(fid, [2 header.ushSamplesInScan], 'float');
    if any(size(a) ~= [2 header.ushSamplesInScan])
      disp('*** ERROR reading from file')
      readMore = 0;
      break
    end
    b      = a(1,:) + i.*a(2,:);
    % reverse if necessary
    if mask.MDH_REFLECT
      b = b(1,end:-1:1);
    end % reverse
    % distinguish between attributes
    if mask.MDH_NOISEADJSCAN
      % very few, and first: don't worry about copying
      noise = [noise b.'];
    elseif mask.MDH_PHASCOR
      % process the same as noise
      phasecor = [phasecor b.'];
      % - but remember the slices
      pc_order = [pc_order header.sLC(3)+1];
    elseif mask.MDH_REFPHASESTABSCAN | mask.MDH_PHASESTABSCAN
      phasestab = [phasestab b.'];
    else % if 0
      % assume this is an imaging scan
      if undoOS
        lb = length(b);
        b  = fftshift(ifft(fftshift(b(:))));
        %b  = b(lb/4+1:3*lb/4);
        b  = fftshift(fft(fftshift(b(lb/4+1:3*lb/4))));
      end
      lb = length(b);
      f = (raw_size(1)-lb)/2+1;
      t = (raw_size(1)+lb)/2;
      if mode == 'TSE'
        header.sLC(9) = 0;
      end
      if mode == 'EPI'
        header.sLC(9) = 0;
        header.sLC(8) = 0;
      end
      if repWise
        header.sLC(7) = 0;
      end
      if extRefOnly 
        if (~mask.MDH_PATREFSCAN)
          continue
        end
      else
        %disp([mask.MDH_PATREFSCAN mask.MDH_PATREFANDIMASCAN])
        if mask.MDH_PATREFSCAN & (~mask.MDH_PATREFANDIMASCAN)
          continue
        end
      end
      if average
        header.sLC(2) = 0;
      end
      %if any([header.ushSamplesInScan header.ushChannelId+1 header.sLC] ...
      if any([lb header.ushChannelId+1 header.sLC+1] ...
             > raw_size)
        % index out of range - ignoring
        disp('** ignoring index')
        continue
      end
      try
        if average & header.sLC(2) > 0
          raw( f:t, header.ushChannelId+1, ...
               header.sLC(1)+1, ...
               header.sLC(2)+1, ...
               header.sLC(3)+1, ...
               header.sLC(4)+1, ...
               header.sLC(5)+1, ...
               header.sLC(6)+1, ...
               header.sLC(7)+1, ...
               header.sLC(8)+1, ...
               header.sLC(9)+1, ...
               header.sLC(10)+1, ...
               header.sLC(11)+1, ...
               header.sLC(12)+1, ...
               header.sLC(13)+1, ...
               header.sLC(14)+1 ) = ...
          raw( f:t, header.ushChannelId+1, ...
               header.sLC(1)+1, ...
               header.sLC(2)+1, ...
               header.sLC(3)+1, ...
               header.sLC(4)+1, ...
               header.sLC(5)+1, ...
               header.sLC(6)+1, ...
               header.sLC(7)+1, ...
               header.sLC(8)+1, ...
               header.sLC(9)+1, ...
               header.sLC(10)+1, ...
               header.sLC(11)+1, ...
               header.sLC(12)+1, ...
               header.sLC(13)+1, ...
               header.sLC(14)+1 ) + b(:);
        else
          raw( f:t, header.ushChannelId+1, ...
               header.sLC(1)+1, ...
               header.sLC(2)+1, ...
               header.sLC(3)+1, ...
               header.sLC(4)+1, ...
               header.sLC(5)+1, ...
               header.sLC(6)+1, ...
               header.sLC(7)+1, ...
               header.sLC(8)+1, ...
               header.sLC(9)+1, ...
               header.sLC(10)+1, ...
               header.sLC(11)+1, ...
               header.sLC(12)+1, ...
               header.sLC(13)+1, ...
               header.sLC(14)+1 ) = b(:);
        end
      catch
        disp('*** SIZE too large, aborting')
        readMore = 0;
        f
        t
        header.ushChannelId+1
        header.sLC.'+1
        break
      end % try
    % else
    end % if mask....
    
  end % if readMore
  
end % while readMore

fclose(fid);
