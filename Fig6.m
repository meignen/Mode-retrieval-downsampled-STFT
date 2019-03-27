close all
 %bat signal  
 load -ascii batsig.txt
 s = batsig;
 s = s(145:end)';
 s = hilbert(s);
 s = s(:);
 N = length(s);
 
 if (N == 2^(floor(log2(N))))
  Nfft = N;
 else
  Nfft = 2^(floor(log2(N))+1);
 end
 sigma_opt = 0.13;   
 Lg = 65;
 clwin = 10;
  
  
 %Hamming window
 
 hlength=floor(Lg);%optimal window determined by Rényi entropy
 hlength=hlength+1-rem(hlength,2);%the length of the filter has to be odd
 h = tftb_window(hlength,'Hamming');
 [hrow,hcol]=size(h); 
 Lh=(hrow-1)/2;
  
 [tfr,norm2h] = tfrstft_three_case_down(s,Nfft,2,h,Lh,1,0); 
 Abstfr = abs(tfr);
 %ridge extraction
 nr = 4;
 [Cs] = exridge_mult(tfr,nr,0,0,clwin);
 Cs = Cs';
 
 figure
 imagesc(Abstfr(1:150,1:230))
 set(gca,'ydir','normal');
 hold on;
 for k=1:nr
  plot(Cs(1:230,k));
 end
 hold off;
 
 
 [tfr,norm2h] = tfrstft_three_case_down(s,Nfft,2,h,Lh,4,0); 
 Abstfr = abs(tfr);
 
 %ridge extraction
 nr = 4;
 [Cs] = exridge_mult(tfr,nr,0,0,clwin);
 Cs = Cs';
 figure
 imagesc(Abstfr(1:150,1:floor(230/4)))
 set(gca,'ydir','normal');
 hold on;
 for k=1:nr
  plot(Cs(1:floor(230/4),k));
 end
 hold off;
 
 [tfr,norm2h] = tfrstft_three_case_down(s,Nfft,2,h,Lh,8,0); 
 Abstfr = abs(tfr);
 %ridge extraction
 nr = 4;
 [Cs] = exridge_mult(tfr,nr,0,0,clwin);
 Cs = Cs';

 figure
 imagesc(Abstfr(1:150,1:floor(230/8)))
 set(gca,'ydir','normal');
 hold on;
 for k=1:nr
  plot(Cs(1:floor(230/8),k));
 end
 hold off;
 
 [tfr,norm2h] = tfrstft_three_case_down(s,Nfft,2,h,Lh,16,0); 
 Abstfr = abs(tfr);
 %ridge extraction
 nr = 4;
 [Cs] = exridge_mult(tfr,nr,0,0,clwin);
 Cs = Cs';
 figure
 imagesc(Abstfr(1:150,1:floor(230/16)))
 set(gca,'ydir','normal');
 hold on;
 for k=1:nr
   plot(Cs(1:floor(230/16),k));
 end
 hold off;
  


 [X,XX,coeffX,Y,YY,coeffY,Z,ZZ,coeffZ]=recons_bat(3);
 [X1,XX1,coeffX1,Y1,YY1,coeffY1,Z1,ZZ1,coeffZ1]=recons_bat(4);

 d=0:10;
 %SNR 0
 
[snr_out,nbcoeff]= snrout_down(3,'Hamming',16,0);
[snr_out1,nbcoeff1]= snrout_down(3,'Hamming',8,0);
[snr_out2,nbcoeff2]= snrout_down(3,'Hamming',4,0);
coeff0 = [nbcoeff nbcoeff1 nbcoeff2];
snrout0 = [snr_out snr_out1 snr_out2];

%SNR 10
[snr_out,nbcoeff]= snrout_down(3,'Hamming',16,10);
[snr_out1,nbcoeff1]= snrout_down(3,'Hamming',8,10);
[snr_out2,nbcoeff2]= snrout_down(3,'Hamming',4,10);
coeff = [nbcoeff nbcoeff1 nbcoeff2];
snrout = [snr_out snr_out1 snr_out2];

figure
plot(coeffX,XX,coeffX1,XX1,'--',coeff0,snrout0,'-.',coeffZ,ZZ,':',coeffZ1,ZZ1,'->',coeff,snrout,'-<');

figure
plot(3*(2*d+1),X,4*(2*d+1),X1,'--',3*(2*d+1),Z,'-.',4*(2*d+1),Z1,':');
