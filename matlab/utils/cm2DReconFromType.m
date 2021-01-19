function [L,type]=cm2DReconFromType(type)

if(isstruct(type))
TYPE=type.Type; 
else
    TYPE=type;
end

if nargin>0
    %    switch lower(m)
    switch lower(type)
        case 'rss'
            L=cm2DRecon();
        case {'sense','msense'}
            L=cm2DReconSENSE();
        case 'b1'
            L=cm2DReconB1();
        case 'grappa'
            L=cm2DReconGRAPPA();
        case 'adapt'
            L=cm2DReconAdapt();
        otherwise
            L=cm2DRecon();
            
    end
    
    
else
    METHODS={'rss','sense','grappa','b1','msense','adapt'};
    fprintf(1,'available methods ')
    for m=1:numel(METHODS)
        fprintf(1,[METHODS{m} ', ']);
    end
    fprintf(1,'\b\b  \neros.montin@gmail.com\n');
    
    
end