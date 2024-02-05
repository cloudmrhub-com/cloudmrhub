function new_csensset=coil_im2csensset(csensset,sens_type,reference_dataset)

% This function will take a raw coil sensitivity set and produce a more processed version
%
% function new_csensset=coil_im2csensset(csensset,sens_type,reference_dataset)
%
% Input:
%	new_csensset:		raw sensitivity set
%	sens_type:			type of smoothing algorithm to be used-
%								'sense':	SENSE type sensitivity (csensset/reference image)
%								'sum of squares': divide csensset by Sum of Squares of csensset
%								'Wang and Reykowski': Not implemented yet
%								'standard': no filtering applied but csensset is multiplied by mask function
%	reference_image:	body coil reference image
%
% Output: 
%	new_csensset:	processed sensitivity set
%
% Charles A. McKenzie  8/9/99  

num_coil=size(csensset,1);
new_csensset=zeros(size(csensset));

%switch lower(sens_type)
switch sens_type

case 1 %'standard'
   new_csensset=csensset;
   
case 2 %'sense'
   %divide csensset by reference_dataset (body coil image)
   new_csensset=csensset./repmat(reference_dataset(:).',num_coil,1);
   
case 3 %'sum of squares'
   %calculate Sum of Squares (SoS) sensitivity reference
   SoS_sens_ref=sqrt(sum(abs(csensset).^2,1));
   
   %divide each coil sensitivity by the sum of squares sensitivity reference
   new_csensset=csensset./repmat(SoS_sens_ref,num_coil,1);
   
case 4 %'Wang and Reykowski'
   
   new_csensset=csensset;
   
otherwise
   new_csensset=csensset;
   disp('Unknown Sensitivity Construction Type: No construction done')
   
end % switch

