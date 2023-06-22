function wien()


A=[linspace(1,9,8);zeros(1,8);linspace(11,19,8);zeros(1,8);linspace(21,29,8);zeros(1,8);linspace(31,39,8);zeros(1,8)]
A=A+1j*A;

wienc(A,2)

function NOISE=wienc(noiseonly,box)
            %noiseonly is the reconstructed noiseonly image 2D 1 slice, box the size
            %of the box
            
            [NC,NR]=size(noiseonly);
            
            kx=box;
            ky=box; 
            NOISE=nan(NC,NR);
             PADDEDNOISE = padarray(noiseonly,[kx ky],NaN);
                for ic=1:NC
                    for ir=1:NR
                        pic=kx+ic+[-kx:kx];
                        pir=ky+ir+[-ky:ky];
                        try
                            NOISE(ic,ir)=nanstd(reshape(PADDEDNOISE(pic,pir),[],1));
                        catch
                        end
                    end
                end
                
        end
    end