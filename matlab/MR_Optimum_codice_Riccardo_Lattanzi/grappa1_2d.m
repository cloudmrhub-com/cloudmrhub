function [krecon,ws] = grappa1_2d(sig,acs,af,ksize,flag_acs)
%-------------------------------------------------------------------------
% GRAPPA reconstruction of a 2D image with uniform acceleration along the 
% y dimension. Based on Griswold MA et al. Generalized autocalibrating
% partially parallel acquisitions (GRAPPA). Magn Reson Med. 2002
% Jun;47(6):1202-10. The borders of k-space are zero-filled to accomodate 
% the GRAPPA convolution operator. 
%-------------------------------------------------------------------------
%	Input:
%	sig: accelerated acquisition [ny/af,nx,nc]. 
%	acs: autocalibration signal  [nyacs,nxacs,nc]. 
%	af: acceleration factor (integers only).
%   ksize: grappa kernel size [ksy,ksx]
%   flag_acs: 1, include the acs data into the final reconstruction
%
%	Output:
%	krecon: coil-by-coil k-space recon [ny,nx,nc].
%   ws: grappa weights [af,nc*ksx*ksy,nc].
%--------------------------------------------------------------------------
% Ricardo Otazo, NYU-CBI
% January 17, 2008
%
% Riccardo Lattanzi
% March 10, 2010 - added the weights to the output
%--------------------------------------------------------------------------
[ny,nx,nc]=size(sig);           
[nyacs,nxacs,nc]=size(acs);
% acs data into the final reconstruction?
if nargin<5,flag_acs=0;end
% grappa kernelsize
if nargin<4,ksy=4;ksx=3;else ksy=ksize(1);ksx=ksize(2);end

fprintf(strcat('\n 1D-GRAPPA reconstruction (',int2str(ksy),'x',int2str(ksx),' kernel) ........................................'));
fprintf(strcat('\n Gradient-encoding data: ',int2str(ny),'x',int2str(nx),'(Ry=',int2str(af),')'));
fprintf(strcat('\n Autocalibration data: ',int2str(nyacs),'x',int2str(nxacs)));
fprintf(strcat('\n Number of channels: ',int2str(nc),'\n'));

src=zeros((nyacs-(ksy-1)*af)*(nxacs-(ksx-1)),nc*ksx*ksy);
targ=zeros(af-1,(nyacs-(ksy-1)*af)*(nxacs-(ksx-1)),nc);

src_idx=0;                         
% reordering data
fprintf(' blocks:');
for xi=(ksx-1)/2+1:nxacs-(ksx-1)/2,
for yi=1:nyacs-(ksy-1)*af,
    src_idx=src_idx+1;
    % collects a ksyxksx block of source points around the target point
    block=[];
    for bxi=-(ksx-1)/2:(ksx-1)/2,
    for byi=0:(ksy-1),
        block=cat(1,block,squeeze(acs(yi+byi*af,xi+bxi,:)));
    end
    end
    src(src_idx,:)=block;
    % target point for these source points
    for afi=1:af-1,
        targ(afi,src_idx,:)=squeeze(acs(yi+(ksy-2)/2*af+afi,xi,:)); 
    end                                                         
end
end

fprintf('done! - weights:');
% weights using pseudoinverse
ws=zeros(af-1,nc*ksx*ksy,nc);
for afi=1:af-1,
    ws(afi,:,:)=inv(src'*src)*src'*squeeze(targ(afi,:,:));      
end
% reconstructed k-space matrix
krecon=zeros(ny*af+ksy*af,nx+ksx-1,nc);
% krecon=zeros(ny*af+(ksy-1)*af,nx+ksx-1,nc);

% source data in the reconstructed k-space matrix
krecon(ksy/2*af+1:af:end-ksy/2*af,(ksx-1)/2+1:nx+(ksx-1)/2,:)=sig;                         

% %************************OLD RICARDO'S CODE******************************
% fprintf('done! - recon:');
% for xi=(ksx-1)/2+1:nx+(ksx-1)/2,
% for yi=1:af:ny*af,
%     src=[];
%     for bxi=-(ksx-1)/2:(ksx-1)/2,
%     for byi=0:(ksy-1),
%         src=cat(1,src,squeeze(krecon(yi+byi*af,xi+bxi,:)));
%     end
%     end
%     % recon=source*weights
%     for afi=1:af-1,
%         krecon(yi+(ksy-2)/2*af+afi,xi,:)=transpose(src)*squeeze(ws(afi,:,:));
%     end                                                         
% end
% end
% % extracting the reconstructed data
% if mod(af,2)==1,
%     krecon=krecon(ksy/2*af+1:ny*af+ksy/2*af-1,(ksx-1)/2+1:nx+(ksx-1)/2,:);                         
% else
%     krecon=krecon(ksy/2*af+1:ny*af+ksy/2*af,(ksx-1)/2+1:nx+(ksx-1)/2,:);
% end
% %*************************************************************************

%************************RICCARDO LATTANZI March 23,2010********************
fprintf('done! - recon:');
for xi=(ksx-1)/2+1:nx+(ksx-1)/2,
for yi=1:af:ny*af,
    src=[];
    for bxi=-(ksx-1)/2:(ksx-1)/2,
    for byi=1:ksy,
        src=cat(1,src,squeeze(krecon(yi+byi*af,xi+bxi,:)));
    end
    end
    % recon=source*weights
    for afi=1:af-1,
        krecon(ksy/2*af+yi+afi,xi,:)=transpose(src)*squeeze(ws(afi,:,:));
    end                                                         
end
end
% extracting the reconstructed data
if mod(af,2)==1,
    krecon=krecon(ksy/2*af+1:ny*af+ksy/2*af-1,(ksx-1)/2+1:nx+(ksx-1)/2,:);                         
else
    krecon=krecon(ksy/2*af+1:ny*af+ksy/2*af,(ksx-1)/2+1:nx+(ksx-1)/2,:);
end
%*************************************************************************

if flag_acs,
    [ny,nx,nc]=size(krecon);
    krecon(ny/2-nyacs/2+1:ny/2+nyacs/2,nx/2-nxacs/2+1:nx/2+nxacs/2,:)=acs;
end

fprintf('done!\n');