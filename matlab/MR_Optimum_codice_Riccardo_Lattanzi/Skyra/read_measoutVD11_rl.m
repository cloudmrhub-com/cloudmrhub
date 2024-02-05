function [Raw4,nChannels,nslices,nphases,nacq,ADCchannels] = read_measoutVD11(datatextname)
%
% Funktion zum Einlesen von Rohdaten aus Measout Files
% Die Rohdaten werden in Bilddaten transformiert
% Die Bilddaten werden als n-dimensionales Array an die rufende Funktion
% zurückgegeben.
% 

oversampling=1;  % Oversampling always on
if oversampling==1
    fprintf(' Oversampling On \n');
end

%------------------------------------------------


fid=fopen(datatextname,'r');
data=fread(fid,1,'uint32');                % First four bytes contains length
if data<32
    % VD11 mdh
    hhr = waitbar(0,['Reading VD11 data. ','Please wait!' '...']);
    MrParcRaidFileHeader.ulhdSize = data;
    MrParcRaidFileHeader.ulcount = fread(fid,1,'uint32');     % number or meas in file
    disp(sprintf('Anzahl an Messungen im Meas.out File_: %d', MrParcRaidFileHeader.ulcount));
    
    for ii=1:(MrParcRaidFileHeader.ulcount)
 %       disp(sprintf('Lese bei Position: %d',ftell(fid)));
        MrParcRaidFileEntry(ii)=ReadMrParcRaidFileEntry(fid);     % number or meas in file
    end


    % Number of measurements
    for ii=1:(MrParcRaidFileHeader.ulcount)
        fseek(fid,MrParcRaidFileEntry(ii).uloff,'bof');                   
        %  read "header_lengths" bytes before 1. MDH
        % disp(sprintf('Lese bei Position: %d',ftell(fid)));
        header_length=fread(fid,1,'uint32');                 % First four bytes contains length
        fseek(fid,header_length-4,'cof');                   % 4 byte für header Länge abziehen

        % 
        % disp(sprintf('Lese bei Position: %d',ftell(fid)));
        [sScanHeader, sLC]=ReadsScanHeader(fid);    %sLC = sLoopCounter
        % Number of channels
        nChannels=sScanHeader.ushUsedChannels;
        [mask,tag]=evalInfoMask(sScanHeader.ulEvalInfoMask1);disp(tag);
        disp(sprintf('Anzahl Kanäle: %d', sScanHeader.ushUsedChannels));
        kk=0;
        while ~bitand(sScanHeader.ulEvalInfoMask1, 2^0);
            kk=kk+1;
            % Iterate all number of channels
            for jj=1:sScanHeader.ushUsedChannels
                sChannelHeader=ReadsChannelHeader(fid);

                % *********************************************
                % Hier werden die DAten gelesen
                % * 2 fuer re/im; 
                 raw_block = fread(fid,sScanHeader.ushSamplesInScan*2,'float32');
                % *********************************************
                ADCchannels(jj)=sChannelHeader.ulChannelId; 
                if mask.MDH_REFLECT==0
                    raw_new(:,jj,sLC.ushSlice+1,sLC.ushLine+1,sLC.ushAcquisition+1)=raw_block;
     
                else                    % MDH_REFLECT==1 : Flip line (for EPI)
                    raw_new(:,jj,sLC.ushSlice+1,sLC.ushLine+1,sLC.ushAcquisition+1)=fliplr(raw_block);
                end
            end   % Iterate all number of channels
            % *********************************************
            % *********************************************
            % disp(sprintf('Lese bei Position: %d',ftell(fid)));
            [sScanHeader, sLC]=ReadsScanHeader(fid);
            [mask,tag]=evalInfoMask(sScanHeader.ulEvalInfoMask1);
            %disp(tag);
        end
        disp(sprintf('Anzahl Scans gelesen: %d', kk));
        disp('Gelesene ADC Kanäle auf RX Kassette:'); disp(strcat('RX',num2str(ADCchannels')));
                
    end

else
    % VA11 mdh
    header_length=data;
end

close(hhr);
fclose(fid);

%frewind(fid)                                        % Rewind and restart
%fread(fid,header_length,'uchar');                   % read "header_lengths" bytes before 1. MDH




%% *********************************************************************
% ab hier alter Code
% 2D oder 3D Messung ??
% blocks=nch*np*sl*acq;   % 2D ?
% blocks=nch*np           % 3D ??


[nlines nChannels nslices nphases nacq]=size(raw_new);
nlines=nlines/2;
raw_new2=reshape(raw_new,[2 nlines nChannels nslices nphases nacq]);
clear raw_new
raw_new3=raw_new2(1,:,:,:,:,:)+j*raw_new2(2,:,:,:,:,:);
clear raw_new2
raw_new3=reshape(raw_new3,[nlines nChannels nslices  nphases nacq ]);
raw_new4=zeros(nlines, nChannels, nslices, nphases);
for i=1:nacq
    raw_new4=raw_new4+raw_new3(:,:,:,:,i);
end
clear raw_new3
Raw4=shiftdim(raw_new4,3); 
clear raw_new4


% 
% %------------------------------------------------
% % MDH_REFLECT=1;        
% % if MDH_REFLECT==1
% %     for i1=1:nch
% %         for i2=1:sl
% %             [Raw4(:,:,i1,i2)]=MDH_REFLECT_Correct(Raw4(:,:,i1,i2));     % Flip every second line (for EPI)
% %         end
% %     end
% %     
% % end
% 
% for i1=1:nChannels
%     for i2=1:nslices
%         ima2(:,:,i1,i2)=fftshift(fft2(fftshift(Raw4(:,:,i1,i2)))); % Calculate images
%     end
% end
% 
% clear Raw4
% if oversampling==1
%     for i1=1:nChannels
%         for i2=1:nslices
%             ima(:,:,i1,i2)=ima2(:,nlines/4+1:3*nlines/4,i1,i2)*1000;
%         end
%     end
% end
% clear ima2;
% 
% %whos
% % *********************************************************************
% debug=1;
% if debug==1
%     [nlines nChannels nslices nphases nacq]=size(ima);
% 
%     disp(strcat('nChannels         : ', num2str(nChannels)));
%     disp(strcat('nlines            : ', num2str(nlines)));
%     disp(strcat('nphases           : ', num2str(nphases)));
%     disp(strcat('nslices           : ', num2str(nslices)));
%     disp(strcat('nacq              : ', num2str(nacq)));
% end

end





% *********************************************************************
% 
% 
% 
% 
% 
