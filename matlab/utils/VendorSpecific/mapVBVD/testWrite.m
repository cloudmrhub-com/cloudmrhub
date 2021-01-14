



f='meas_MID00024_FID188178_Multislice.dat';
% 
 F=CLOUDMRRD(f);
% 
 K=F.getNoiseKSpace();
 
 F.write2DCartesianKspacedatainISMRMRDv1(K,'prova.H5');






% % It is very slow to append one acquisition at a time, so we're going
% % to append a block of acquisitions at a time.
% % In this case, we'll do it one repetition at a time to show off this
% % feature.  Each block has nYsamp aquisitions
% 
% 
% 
% 
% 
% 
% f='meas_MID00024_FID188178_Multislice.dat';
% 
% F=CLOUDMRRD(f);
% 
% 
% K=F.getNoiseKSpace();
% 
% 
% 
% size(K)
% % 1     1     1    96    96     5    16
% %   O={'1: Average',
% %       '2: contrast',
% %       '3: repetition',
% %       '4: Frequency Encode',
% %       '5: Phase Encode',
% %       '6: Slice',
% %       '7: Coils'};
% 
% 
% 
% 
% 
% filename ='test.H5';
% dset = ismrmrd.Dataset(filename);
% 
% 
% nX=size(K,4);
% nCoils=size(K,7);
% nYsamp=prod(size(K))/(nX*nCoils);
% nAvg=size(K,1);
% nCnt=size(K,2);
% nRep=size(K,3);
% nPh=size(K,5);
% nSl=size(K,6);
% 
% %nYsamp is number of actually
% acqblock = ismrmrd.Acquisition(nYsamp);
% 
% 
% % Set the header elements that don't change
% acqblock.head.version(:) = 1;
% acqblock.head.number_of_samples(:) = nX;
% acqblock.head.center_sample(:) = floor(nX/2);
% acqblock.head.active_channels(:) = nCoils;
% acqblock.head.available_channels(:) =nCoils;
% acqblock.head.read_dir  = repmat([1 0 0]',[1 nYsamp]);
% acqblock.head.phase_dir = repmat([0 1 0]',[1 nYsamp]);
% acqblock.head.slice_dir = repmat([0 0 1]',[1 nYsamp]);
% 
% 
% 
% 
% 
% counter=0;
% 
% for avg=1:nAvg
%     
%     
%     
%     
%     for cnt=1:nCnt
%         for rep = 1:nRep
%             
%             for sl = 1:nSl
%                 for p = 1:nPh
%                     counter=counter+1;
%                     
%                     if (p==1)
%                         acqblock.head.flagClearAll(counter);
%                         acqblock.head.flagSet('ACQ_FIRST_IN_ENCODE_STEP1',counter);
%                         
%                         if avg == 1
%                             acqblock.head.flagSet('ACQ_FIRST_IN_AVERAGE',counter);
%                         end
%                         if rep == 1
%                             acqblock.head.flagSet('ACQ_FIRST_IN_REPETITION',counter);
%                         end
%                         if cnt == 1
%                             acqblock.head.flagSet('ACQ_FIRST_IN_CONTRAST',counter);
%                         end
%                         if sl == 1
%                             acqblock.head.flagSet('ACQ_FIRST_IN_SLICE',counter);
%                         end
%                         
%                         
%                         
%                     end
%                     
%                     % Set the header elements that change from acquisition to the next
%                     % c-style counting
%                     acqblock.head.scan_counter(counter) = counter;
%                     % Note next entry is k-space encoded line number (not acqno which
%                     % is just the sequential acquisition number)
%                     acqblock.head.idx.kspace_encode_step_1(counter) = p-1;
%                     acqblock.head.idx.repetition(counter) = rep-1;
%                     acqblock.head.idx.average(counter) = avg-1;
%                     
%                     acqblock.head.idx.slice(counter)=sl-1;
%                     acqblock.head.idx.contrast(counter)=cnt-1;
%                     acqblock.head.idx.kspace_encode_step_2(counter)=0;
%                     acqblock.head.idx.segment(counter)=0;
%                     acqblock.head.idx.set(counter)=0;
%                     
%                     
%                     
%                     
%                     
%                     
%                     
%                     
%                     
%                     % Set the flags
%                     %                     acqblock.head.flagClearAll(counter);
%                     
%                     
%                     %                     if rep == 1
%                     %                         acqblock.head.flagSet('ACQ_FIRST_IN_REPETITION',counter);
%                     %                     end
%                     %
%                     %                     if rep ==nRep
%                     %                              acqblock.head.flagSet('ACQ_LAST_IN_REPETITION',counter);
%                     %
%                     %                     end
%                     %
%                     %
%                     %
%                     %
%                     %
%                     %                             acqblock.head.flagSet('ACQ_FIRST_IN_ENCODE_STEP1', counter);
%                     %                         acqblock.head.flagSet('ACQ_FIRST_IN_SLICE', counter);
%                     %
%                     %
%                     %
%                     %
%                     %                     elseif counter==nPh
%                     %                         acqblock.head.flagSet('ACQ_LAST_IN_ENCODE_STEP1', counter);
%                     %                         acqblock.head.flagSet('ACQ_LAST_IN_SLICE', counter);
%                     %                         acqblock.head.flagSet('ACQ_LAST_IN_REPETITION',counter);
%                     %                     end
%                     
%                     
%                     
%                     % fill the data
%                     acqblock.data{counter} = squeeze(K(avg,cnt,rep,:,p,sl,:));
%                     % Append the acquisition block
%                     
%                     
%                     if counter==nPh
%                         %                         acqblock.head.flagSet('ACQ_LAST_IN_ENCODE_STEP1', counter);
%                         %                         acqblock.head.flagSet('ACQ_LAST_IN_SLICE', counter);
%                         %                         acqblock.head.flagSet('ACQ_LAST_IN_REPETITION',counter);
%                         %                     end
%                         
%                         
%                         
%                         %increase counter;
%                         
%                         
%                         
%                         acqblock.head.flagSet('ACQ_LAST_IN_ENCODE_STEP1',counter);
%                         
%                         if avg == nAvg
%                             acqblock.head.flagSet('ACQ_LAST_IN_AVERAGE',counter);
%                         end
%                         if rep == nRep
%                             acqblock.head.flagSet('ACQ_LAST_IN_REPETITION',counter);
%                         end
%                         if cnt == nCnt
%                             acqblock.head.flagSet('ACQ_LAST_IN_CONTRAST',counter);
%                         end
%                         
%                         if sl == nSl
%                             acqblock.head.flagSet('ACQ_LAST_IN_SLICE',counter);
%                         end
%                         
%                         
%                         
%                         
%                     end
%                     
%                     
%                     
%                     %                     if avg ==nAvg
%                     %                         acqblock.head.flagSet('ACQ_LAST_IN_AVERAGE',counter-1);
%                     %
%                     %                     end
%                     %
%                     %
%                     %                     if rep == 1
%                     %                         acqblock.head.flagSet('ACQ_FIRST_IN_REPETITION',counter);
%                     %                     end
%                     %
%                     %                     if rep ==nRep
%                     %                         acqblock.head.flagSet('ACQ_LAST_IN_REPETITION',counter);
%                     %
%                     %                     end
%                     
%                     
%                     
%                     
%                     
%                 end
%                 
%                 if sl ==nSl
%                     acqblock.head.flagSet('ACQ_LAST_IN_SLICE',counter);
%                 end
%             end
%         end
%     end
% end
% acqblock.head.flagSet('ACQ_LAST_IN_MEASUREMENT',counter);
% 
% dset.appendAcquisition(acqblock);
% 
% 
% %get the resolution ffrom somewhere
% 
% r=[1 1 1];
% 
% %%%%%%%%%%%%%%%%%%%%%%%%
% %% Fill the xml header %
% %%%%%%%%%%%%%%%%%%%%%%%%
% % We create a matlab struct and then serialize it to xml.
% % Look at the xml schema to see what the field names should be
% 
% header = [];
% 
% % Experimental Conditions (Required)
% header.experimentalConditions.H1resonanceFrequency_Hz = 128000000; % 3T
% 
% % Acquisition System Information (Optional)
% header.acquisitionSystemInformation.systemVendor = 'CLOUDMR www.cloudmrhub.com';
% header.acquisitionSystemInformation.systemModel = 'CM scanner v01';
% header.acquisitionSystemInformation.receiverChannels = nCoils;
% header.acquisitionSystemInformation.systemFieldStrength_T=2.893620;
% header.acquisitionSystemInformation.relativeReceiverNoiseBandwidth=0793;
% 
% % The Encoding (Required)
% header.encoding.trajectory = 'cartesian';
% header.encoding.encodedSpace.fieldOfView_mm.x = nX*r(1);
% header.encoding.encodedSpace.fieldOfView_mm.y = nPh*r(2);
% header.encoding.encodedSpace.fieldOfView_mm.z = nSl*r(3);
% header.encoding.encodedSpace.matrixSize.x = nX;
% header.encoding.encodedSpace.matrixSize.y = nPh;
% header.encoding.encodedSpace.matrixSize.z = 1;
% % Recon Space
% % (in this case same as encoding space)
% header.encoding.reconSpace = header.encoding.encodedSpace;
% % Encoding Limits
% header.encoding.encodingLimits.kspace_encoding_step_0.minimum = 0;
% header.encoding.encodingLimits.kspace_encoding_step_0.maximum = nX-1;
% header.encoding.encodingLimits.kspace_encoding_step_0.center = floor(nX/2);
% header.encoding.encodingLimits.kspace_encoding_step_1.minimum = 0;
% header.encoding.encodingLimits.kspace_encoding_step_1.maximum = nPh-1;
% header.encoding.encodingLimits.kspace_encoding_step_1.center = floor(nPh/2);
% 
% header.encoding.encodingLimits.kspace_encoding_step_2.minimum = 0;
% header.encoding.encodingLimits.kspace_encoding_step_2.maximum = 0;
% header.encoding.encodingLimits.kspace_encoding_step_2.center = 0;
% 
% header.encoding.encodingLimits.repetition.minimum = 0;
% header.encoding.encodingLimits.repetition.maximum = nRep-1;
% header.encoding.encodingLimits.repetition.center = 0;
% 
% header.encoding.encodingLimits.average.minimum = 0;
% header.encoding.encodingLimits.average.maximum = nAvg-1;
% header.encoding.encodingLimits.average.center = 0;
% 
% 
% header.encoding.encodingLimits.phase.minimum = 0;
% header.encoding.encodingLimits.phase.maximum = nPh-1;
% header.encoding.encodingLimits.phase.center = 0;
% 
% 
% header.encoding.encodingLimits.contrast.minimum = 0;
% header.encoding.encodingLimits.contrast.maximum = nCnt-1;
% header.encoding.encodingLimits.contrast.center = 0;
% 
% header.encoding.encodingLimits.slice.minimum = 0;
% header.encoding.encodingLimits.slice.maximum = nSl-1;
% header.encoding.encodingLimits.slice.center = 0;
% 
% header.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_1 = 1 ;
% header.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_2 = 1 ;
% % header.encoding.parallelImaging.calibrationMode = 'embedded' ;
% 
% % Commented code below appears not necessary - saw this parameter after converting
% % a scanner file using siemens_to_ismrmrd
% % header.userParameters.userParameterLong.name = 'EmbeddedRefLinesE1' ;
% % header.userParameters.userParameterLong.value = ACShw *2  ;
% 
% %% Serialize and write to the data set
% xmlstring = ismrmrd.xml.serialize(header);
% dset.writexml(xmlstring);
% 
% %% Write the dataset
% dset.close();



