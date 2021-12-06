
clear all

A= linspace(1,128,128).';

IM=[]

SL=10

for s=1:SL
for a=1:256
    IM(:,a,s)=complex(a*A,a*A);   
end
end


A= linspace(1,20,20).';
for s=1:SL
for a=1:20
    IM2(:,a,s)=complex(A,A);
end
end


A=cmOutput()for s=1:SL

A.add2DImagetoExport(IM,'gradient 1-128')
A.add2DImagetoExport(IM2,'gradient 1-20')

