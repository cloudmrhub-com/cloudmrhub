function [origin,cosines,info]=mapVBVDutils(X)
%X can be a string with  name of a file or a 

switch(class(X))
    case 'char'
            DATA=mapVBVD(X);
    case 'struct'
        if(isfield(X,'hdr'))
            DATA=X;
        end
    case 'twix_map_obj'
            DATA.image=X;
    case 'cell'
              DATA=DATA{2};
   
        
end


if(~isstruct(DATA))
    DATA=DATA{2};
end

%number of acquisition in the image
NACQ=DATA.image.NAcq;


cosines=zeros(NACQ,9);
origin=zeros(NACQ,3);


% for acq=1:NACQ
% 
% %     from SIEMENS TO ISMRMRD     // Convert Siemens quaternions to direction cosines. In
% % the Siemens convention the quaternion corresponds to a rotation matrix with columns P R S.
% %Siemens stores the quaternion as (W,X,Y,Z)
%     
% % slicepos[Ox,Oy,Oz,QW,QX,QY,QZ];
% %tested against ismrtmrd converter
%     Q=quaternion(DATA.image.slicePos([5:7 4],acq)');
%    
%     cosines(acq,:)= reshape(quat2rotm(Q),1,[]);
%     origin(acq,:)=DATA.image.slicePos(1:3,acq);
% 
% end

%FOV=[DATA.hdr.Config.ReadFoV DATA.hdr.Config.PeFOV  ]

%size=[ DATA.hdr.Config.NImageCols  DATA.hdr.Config.NImageLins  DATA.hdr.Config.NSlc];


%GlobalImageScaleFactor
D=DATA.hdr.Dicom;
C=DATA.hdr.Config;
MY=DATA.hdr.MeasYaps;


info.SiemensVersion       = DATA.image.softwareVersion;

info.seq.siz=[DATA.hdr.Config.RoFOV DATA.hdr.Config.PeFOV ];

info.BaseResolution=DATA.hdr.Config.BaseResolution;
try;info.GlobalImageScaleFactor=DATA.hdr.Config.GlobalImageScaleFactor;end;
try;info.ImaScaleFactor=DATA.hdr.Config.ImaScaleFactor;end;
info.SIZE=[DATA.hdr.Config.ImageColumns DATA.hdr.Config.ImageLines];
info.NFEncoding=DATA.hdr.Config.NColMeas;
try;info.ReadoutOversampling=DATA.hdr.Config.ReadoutOversamplingFactor;end;

info.NAve=DATA.hdr.Config.NAve;
info.NSlc=DATA.hdr.Config.NSlc;
info.NAve=DATA.hdr.Config.NAve;
info.NRep=DATA.hdr.Config.NRep;
info.NSeg=DATA.hdr.Config.NSeg;


info.NoOfFourierLines=DATA.hdr.Config.NoOfFourierLines;
try;info.NoiseScaleFactor=DATA.hdr.Config.NoiseScaleFactor;end;
try;info.ReceiverNoiseEquivalentBandwidth=DATA.hdr.Config.ReceiverNoiseEquivalentBandwidth;end;

info.seq.TI=DATA.hdr.Config.TI*1e-3;
info.seq.TR_=MY.alTR{1}*1e-3;
info.seq.TE=MY.alTE{1}*1e-3;
info.seq.FA=DATA.hdr.Dicom.adFlipAngleDegree;
info.seq.name=D.tScanningSequence;
info.seq.owner=D.tSequenceOwner;
info.seq.variant=D.tSequenceVariant;
try;info.seq.description=C.SequenceDescription;end;
try;info.seq.string=C.SequenceString;end;
info.seq.resolution =[     DATA.hdr.Dicom.dSliceResolution];
try;info.seq.ETL=D.EchoTrainLength;end;
try;info.seq.GRETL=D.GradientEchoTrainLength;;end;
info.seq.acceleration = MY.sFastImaging.lTurboFactor;
if(isempty(DATA.hdr.Config.Is3D));info.seq.is3D=false;else;info.seq.is3D=true;end


try;info.sPatPosition=DATA.hdr.Config.sPatPosition;end;

info.slice.Order=DATA.hdr.Config.relSliceNumber(find(DATA.hdr.Config.relSliceNumber>-1));

for t=1:numel(info.NSlc)
    info.slice.Thickness(t)=MY.sSliceArray.asSlice{t}.dThickness;
    info.slice.normal(t)=MY.sSliceArray.asSlice{t}.sNormal
    info.slice.Order(t)=C.relSliceNumber(t);
    info.slice.Position=MY.sSliceArray.asSlice{t}.sPosition;
    
end
    
    
try ;info.spacingBetweenslices=C.SpacingBetweenSlices;end;


info.scanner.B0=DATA.hdr.Dicom.flMagneticFieldStrength;


nf=DATA.hdr.Config.RoFOV;
removeOS=[1:nf/4, 1+nf*3/4:nf];
totake=setdiff([1:nf],removeOS);
info.removeOSIndexes=totake;














