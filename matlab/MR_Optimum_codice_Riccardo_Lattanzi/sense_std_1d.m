function [s_r,g]=sense_std_1d(y,c,method,msvp,psi)
%
%   1D-SENSE: reconstruction for acceleration along 1 dimension (2D image)
%             acceleration: uniform subsampling along 1st dimension of y 
%-------------------------------------------------------------------------
%	Input:
%	y: accelerated acquisition [ny_acc, nx, nc]. 
%	c: coil sensitivity maps [ny, nx, nc].
%   method: 'pinv', 'ssvd','tsvd'
%   msvp: maximum singular value percentage (for svd method) 
%
%	Output:
%	s_r: reconstructed image [ny, nx].
%	g: g-factor map [ny, nx].
%--------------------------------------------------------------------------
% Ricardo Otazo
%--------------------------------------------------------------------------
%
ref=c;
if strcmp(method,'pinv')==0 && strcmp(method,'ssvd')==0 && strcmp(method,'tsvd')==0,
    fprintf('\n Error! Method name is incorrect ...\n\n');
    return;
end
% checking input parameters
if nargin<5,psi=eye(size(y,3));end
if nargin<4,msvp=1;end
msvp=msvp/100;

% data preparation
y=permute(y,[3 1 2]);
c=permute(c,[3 1 2]);
[nc,ny,nx]=size(c);
[nc,nyfold,nxfold]=size(y);
af=round(ny/nyfold);
if mod(af,2)==0,
    y=ifftshift(y,2);
end

% whitening matrix
[u,s,v]=svd(psi);
w=pinv(sqrt(s))*v';	

step_y=size(y,2);
% iteration for each point in the aliased image
for xi=1:nx,
   for yi=1:step_y,
       % prewhitening
       y_w=w*squeeze(y(:,yi,xi));           
       E=w*squeeze(c(:,yi:step_y:ny,xi));   
       %solution
       if strcmp(method,'pinv'),
           %s_r(yi:step_y:ny,xi)=pinv(E)*y_w; %#ok<AGROW>
           s_r(yi:step_y:ny,xi)=inv(E'*E)*E'*y_w; %#ok<AGROW>
           % g-factor
           g(yi:step_y:ny,xi)=sqrt(diag(pinv(E'*E)).*diag(E'*E)); %#ok<AGROW>
       else
            % svd-based solution
            [U,S,V] = svd(E);E_sv=diag(S);E_sv0=E_sv;
            if strcmp(method,'ssvd'),
                E_sv=E_sv+max(E_sv)*msvp;
            elseif strcmp(method,'tsvd'),
                E_sv(find(E_sv<max(E_sv)*msvp))=0;
            end
            Ei=V*pinv(diag(E_sv))*U(:,1:size(V,1))';
            s_r(yi:step_y:ny,xi)=Ei*y_w; %#ok<AGROW>
            % g-factor
           g(yi:step_y:ny,xi)=sqrt(diag(Ei*Ei').*diag(E'*E));  %#ok<AGROW>
       end
   end
end
s_r=sqrt(af)*s_r;