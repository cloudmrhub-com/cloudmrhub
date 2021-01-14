


fsignal='/data/MYDATA/TestSNR_15Apr2019_multislice/RAWDATA/meas_MID00024_FID188178_Multislice.dat';


fs='/data/MYDATA/TestSNR_15Apr2019_multislice/RAWDATA/meas_MID00027_FID188181_Multislice_no_RF.dat';

fm='/data/MYDATA/TestSNR_15Apr2019_multislice/RAWDATA/multi_RAID/Signal_Data_Multislice_MULTIRAID.dat';

% 
 FS=CLOUDMRRD(fs);
% 
 
 
  KSN=CLOUDMRSpreadthenoiseinto2DKSpace(FS,1);

  FM=CLOUDMRRD(fm);

 
  KMN=CLOUDMRSpreadthenoiseinto2DKSpace(FM,0);
 
 ACMm=CLOUDMRgetclassfromOptions(CLOUDMRgetOptions('rssbart'));

 
 
 
 ACMs=CLOUDMRgetclassfromOptions(CLOUDMRgetOptions('rssbart'));
 
 
 ACMs.setNoiseKSpace(KSN);
 
 ACMm.setNoiseKSpace(KMN);

 
 
 ACMm2=ACMm;
 
covs=ACMs.getNoiseCovariance() ;


covm=ACMm.getNoiseCovariance();
 


  

 
 
 
 
 DIFF=abs((covm(:)-covs(:))./covs(:))
 mean(DIFF*100)
 
display('done!');






%calculate SN




SIGNAL=CLOUDMRRD(fsignal);

thesignal=SIGNAL.getKSpaceImageSlice(1,1,1,1);


thesignalm=FM.getKSpaceImageSlice(1,1,1,1);


% subplot(131); imshow(singleSNR,[]); colorbar();subplot(132); imshow(multiSNR,[]); colorbar();subplot(132); imshow(multiSNR2,[]); colorbar();



clear o;
o=CLOUDMRgetOptions('rssbart');

SNRs=testSNR___(o,thesignal,KSN);

SNRm=testSNR___(o,thesignalm,KMN);


T=nanmean(100*abs((SNRs(:)-SNRm(:))./SNRs(:)))

subplot(121); imshow(SNRs,[]); colorbar();title('singlei RAID RSS');
subplot(122); imshow(SNRm,[]); colorbar();title(['multi RAID RSS, ME: ' num2str(T) '% ' ]);






clear o;
o=CLOUDMRgetOptions('b1espirit');

SNRs=testSNR___(o,thesignal,KSN,thesignal);

SNRm=testSNR___(o,thesignalm,KMN,thesignalm);


T=nanmean(100*abs((SNRs(:)-SNRm(:))./SNRs(:)))

subplot(121); imshow(SNRs,[]); colorbar();title('singlei RAID B1');
subplot(122); imshow(SNRm,[]); colorbar();title(['multi RAID B1, ME: ' num2str(T) '% ' ]);

 



clear o;
o=CLOUDMRgetOptions('b1espirit');

SNRs=testSNR___(o,thesignal,KSN,thesignal);

SNRm=testSNR___(o,thesignalm,KMN,thesignalm);


T=nanmean(100*abs((SNRs(:)-SNRm(:))./SNRs(:)))

subplot(121); imshow(SNRs,[]); colorbar();title('singlei RAID B1');
subplot(122); imshow(SNRm,[]); colorbar();title(['multi RAID B1, ME: ' num2str(T) '% ' ]);
 
 


clear o;
o=CLOUDMRgetOptions('msensesimplesense');

o.AccelerationF=1;
o.AccelerationP=4;


o.Autocalibration=0;


thesignald= undersamplemSense2D(thesignal,o.AccelerationF,o.AccelerationP,o.Autocalibration);

SNRs=testSNR___(o,thesignald,KSN,thesignal);

thesignalmd= undersamplemSense2D(thesignalm,o.AccelerationF,o.AccelerationP,o.Autocalibration);

SNRm=testSNR___(o,thesignalmd,KMN,thesignalm);


T=nanmean(100*abs((SNRs(:)-SNRm(:))./SNRs(:)))

subplot(121); imshow(abs(SNRs),[]); colorbar();title('singlei RAID sense');
subplot(122); imshow(abs(SNRm),[]); colorbar();title(['multi RAID sense, ME: ' num2str(T) '% ' ]);
 




 
 
 
 
 