slicenum=45;

nx=640;
ny=460;
nc=12;



if slicenum<10,
    filenameweights=['kspaceweights_00' int2str(slicenum) '.dat'];
    filenamekrecon=['recondkspace_00' int2str(slicenum) '.dat'];
else
    filenameweights=['kspaceweights_0' int2str(slicenum) '.dat'];
    filenamekrecon=['recondkspace_0' int2str(slicenum) '.dat'];
end;

cd I:\computer_backup_23apr2011\Pippas_files\subrawpaper\from_Stephan\Rep0\kspaceweights

kweightsimages=readkweightsimagedomain(filenameweights,nx,ny,nc);%order: (y,coilin,coilout,x). 

cd I:\computer_backup_23apr2011\Pippas_files\subrawpaper\from_Stephan\Rep0\recondkspace

image3d=readkreconfunction(filenamekrecon,nx,ny,nc);%order: nx,ny,nc
gmap=zeros(nx,ny);
sosimage=zeros(nx,ny);
for countx=1:nx,
    for county=1:ny,
        mW=zeros(nc,nc);
        mPT=zeros(1,nc);
        for coilout=1:nc,
            for coilin=1:nc,
                mW(coilout,coilin)=kweightsimages(county,coilin,coilout,countx);
            end;
            mPT(coilout)=conj(image3d(countx,county,coilout));
        end;
        mPTmW=mPT*mW;
        numerator=sqrt(abs(mPTmW*mPTmW'));
        denominator=sqrt(abs(mPT*mPT'));
        gmap(countx,county)=numerator/denominator;
        sosimage(countx,county)=denominator;
    end;
end;

figf(1)
imshow(gmap,[0 max(gmap(:))])
figf(2)
imshow(sosimage,[0 max(sosimage(:))])


