classdef cmOutput<handle
    %main class of cloudmr
    %Main reconstrucion class
    %author:eros.montin@gmail.com
    %v012021
    %46&2 just ahead of me
    
    
    properties(Access=private)
        Exporter % a cell array of name and alues to be exported
        Log
        Type
        SubType
        OUTPUTLOGFILENAME
        OUTPUTFILENAME
    end
    
    
    
    methods(Access=protected)
        
        function addToExporter(this,type,name,value)
            this.Exporter=[this.Exporter;[{type},{name},{double(value)}]];
        end
        
    end
    methods
        function this = cmOutput()
            this.logIt(['class ' class(this) ' instantiated'],'start');
            try
                TurnOffWarnings();
                
            catch
                
            end
            
            
            
            
            
        end
        
        
        function setOutputLogFileName(this,L)
            this.OUTPUTLOGFILENAME=L;
        end
        
        function o=getOutputLogFileName(this)
            o=this.OUTPUTLOGFILENAME;
        end
        
        
        function setOutputFileName(this,L)
            this.OUTPUTFILENAME=L;
        end
        
        function o=getOutputFileName(this)
            o=this.OUTPUTFILENAME;
        end
        
        
        function setLog(this,L)
            this.Log=L;
        end
        
        function o=getLog(this)
            o=this.Log;
        end
        
        
        
             function setTypeOutput(this,L)
            this.Type=L;
        end
        
        function o=getTypeOutput(this)
            o=this.Type;
        end
        
        
        
        function appendLog(this,L)
            this.setLog(cat(2,this.getLog(),L));
        end
        
        function logIt(this,W,t)
            if nargin<3
                t='null';
            end
            j.time=datestr(clock);
            j.text=W;
            j.type=t;
            this.Log=[this.Log j];
            
        end
        
        
        
        
        
        
        function O=getResults(this)
            O=this.exportResults();
        end
        
        
        function O=getResultsAsImages(this)
            O=this.getImagesFromResuls(this.getResults());
        end
        
        
        
        
        
        function O=exportResults(this,fn)
            
            O.version='20201223';
            O.author='eros.montin@gmail.com';
            
            if isempty(this.Type)
                
                O.type='DATA';
            else
                O.type=this.Type;
            end
            
            
            
            if isempty(this.SubType)
                
                O.subtype='';
            else
                O.subtype=this.SubType;
            end
            
            
            
            
            for t=1:size(this.Exporter,1)
                if(strcmp(this.Exporter{t,1},'image2D'))
                    im.slice=this.image2DToJson(this.Exporter{t,3});
                    im.imageName=this.Exporter{t,2};
                    O.images(t)=im;
                    clear im;
                    
                end
                
                
            end
            if (nargin>1)
                myjsonWrite(jsonencode(O),fn);
            else
                try
                    myjsonWrite(jsonencode(O),this.getOutputFileName());
                end
            end
        end
        
        
        
        
        
        
        function exportLog(this,fn)
            if(nargin>1)
                myjsonWrite(jsonencode(this.Log),fn);
            else
                try
                    myjsonWrite(jsonencode(this.Log),this.getOutputLogFileName());
                catch
                    display('error');
                end
                
            end
        end
        
        function whatHappened(this)
            %where's bugo?
            for t=1:numel(this.Log)
                display(this.Log(t));
            end
        end
        
        
        
        function errorMessage(this)
            this.logIt('ERROR','error');
        end
        
        function outputError(this,fn)
            this.errorMessage();
            this.exportLog(fn);
        end
        
        
        
        
        
        
        function add2DImagetoExport(this,im,name)
            %adds to results the image and its name
            % a.add2DImagetoExport(randn(20),'imname')
            if (isnumeric(im))
                this.addToExporter('image2D',name,im);
            else
                display('nothing added')
            end
        end
    end
    
    methods (Static)
        
        function Oout=image2DToJson(im)
            for sl=1:size(im,3)
                o=im(:,:,sl).';
                O.w=size(o,2);
                O.h=size(o,1);
                O.type= 'double complex';
                O.Vr=real(reshape(o,1,[]));
                O.Vi=imag(reshape(o,1,[]));
                Oout(sl)=O;
            end
        end
        
     
               function noise_coeff=calculateNoiseCoefficientsMatrix(noisecov)
            noise_coeff=zeros(size(noisecov));
            for itemp = 1:size(noisecov,1)
                for jtemp = 1:size(noisecov,1)
                    noise_coeff(itemp,jtemp) = noisecov(itemp,jtemp)/sqrt(noisecov(itemp,itemp)*noisecov(jtemp,jtemp));
                end
            end
        end
        
        function o=rescale01(bla)
            MI=min(bla(:));
            MA=max(bla(:));
            o= (bla(:) - MI)./ ( MA - MI );
        end
        
        
        
        
        function [O]=getJsonResultFromJsonFile(file)
            %from jsonfilename you get the struct of data
            O=readanddecodejson(file);
        end
        
        
        
        function [O]=getImagsesFromJsonResultFile(file)
            %from the json results file get the entire image
            %set in a struct variable
            data=CLOUDMROutput.getJsonResultFromJsonFile(file);
            O=CLOUDMROutput.getImagesFromResuls(data);
        end
        
        function [O] =getImagesFromResuls(data)
            %from the struct derived by the json results file got the entire image
            %set
            
                for imnumber=1:numel(data.images)
                h=data.images(imnumber).slice.h;
                w= data.images(imnumber).slice.w;
                NSL=numel(data.images(imnumber).slice);
                clear image_;
                %we switch i know!:)
                image_=NaN(w,h,NSL);
                O(imnumber).ImageName=data.images(imnumber).imageName;
                for slnumber=1:NSL
                    a=data.images(imnumber).slice(slnumber).Vr(:)+data.images(imnumber).slice(slnumber).Vi(:)*1i;
                    
                    image_(:,:,slnumber)=reshape(a,h,w).';
                end
                O(imnumber).image=image_;
                
                
            end
        end
        
        
        
        function TEST= test()
            %instantiate the class
            TEST=cmOutput();
            TEST.logIt('the test start now','start')
            T0=randn(5)+i*randn(5);
            T1=randn(5);
            TEST.logIt('created 2 results','ok')
            %add two results
            TEST.logIt('add first result','ok')
            try
                TEST.add2DImagetoExport(T0,'test0');
            catch
                TEST.errorMessage()
            end
            
            TEST.logIt('add second result','ok')
            try
                TEST.add2DImagetoExport(T1,'test0');
            catch
                TEST.errorMessage()
            end
            %get the results
            TEST.logIt('retrieve the results json','ok')
            try
                R= TEST.getResults();
            catch
                TEST.errorMessage();
            end
            TEST.logIt('results retrieved','ok')
            
            %get the results
            TEST.logIt('retrieve the results json','ok')
            
            TEST.logIt('retrieve the results as images','ok')
            try
                O=TEST.getResultsAsImages();
                TEST.logIt('retrieve the results json','ok')
            catch
                TEST.errorMessage();
            end
            %check the results
            TEST.logIt('tested the results ','?')
            try
                f=sum(abs(O(1).image(:)-T0(:)))+sum(abs(O(2).image(:)-T1(:)));
                TEST.logIt(['tested the results  and result error is :' num2str(f) ],'ok')
            catch
                TEST.errorMessage();
            end
            
            
            if(f==0)
                TEST.logIt(['tested the results  and result error is :' num2str(f) ],'ok')
            else
                TEST.errorMessage();
            end
            
            
            
           
            
            
            
            
        end
        
        
      
        
        function [KOUT]=plotImageAfterTest(IM,tit)
            
            figure()
            subplot(121)
            imshow(IM,[]);
            title(tit);
            colorbar();
            subplot(122)
            imshow(abs(IM),[]);
            title([tit  ' abs']);
            colorbar();
            
            ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
            text(0.5, 0.98,[num2str(size(IM,1)) 'x' num2str(size(IM,2))])
            
            
        end
        
        
        
        function [KOUT]=plotTwoImagesAfterTest(IM1,IM2,tit1,tit2)
            
            figure()
            subplot(121)
            imshow(IM1,[]);
            title(tit1);
            colorbar();
            subplot(122)
            imshow(IM2,[]);
            title([tit2]);
            colorbar();
            
            ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
            text(0.5, 0.98,['IM1: ' num2str(size(IM1,1)) 'x' num2str(size(IM1,2))])
            text(0.5, 0.88,['IM2: ' num2str(size(IM2,1)) 'x' num2str(size(IM2,2))])
            
       
            RD=abs(IM1(:)-IM2(:))./IM1(:)*100;
                 G=find(~isinf(RD(:)) & ~isnan(RD(:)));
            text(0.5, 0.95,['%RD mean: ' num2str(nanmean(abs(RD(G)))) ', std: ' num2str(nanstd(abs(RD(G))))])
            
        end
        
        
        
      
        
                function NK =resize2DMatrixbyInterpolationVoxelSpace(theK,x1,x2)
    
        [nX, nY]=size(theK);
       [newxg,newyg]=ndgrid(x1,x2);
       %old grid
       [oldxg,oldyg]=ndgrid([1:nX],[1:nY]);
       %scattered interpolation
        if(isreal(theK))
         F=scatteredInterpolant([oldxg(:) oldyg(:)],double(theK(:)));
        %composition of the channels
        NK=F(newxg,newyg);
        else
        FR=scatteredInterpolant([oldxg(:) oldyg(:)],double(real(theK(:))));
        FI=scatteredInterpolant([oldxg(:) oldyg(:)],double(imag(theK(:))));
        %composition of the channels
        NK=FR(newxg,newyg) +1i*FI(newxg,newyg);
        end    
        
                end
        
    end
end

