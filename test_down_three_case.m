function test_down_three_case(cas,window,downsamp,shift)
 
 %cas    : 1, first STFT technique
 %         2, second STFT technique (filter with unit energy on its support)
 %         3, third STFT technique (filter L1 normed on its support)  
 %window : choice for COLA(R) window
 %downsamp : downsampling factor between 1 and Lh
 %shift  : we use the sample m Lh + shift to reconstruct (shift < Lh) 
 
 N = 4096;
 t =(0:N-1)/N;
 a  = 2;
 s1 = a.*exp(2*pi*1i*(1000*t+60*cos(3*pi*t)));
 s2 = a.*exp(2*pi*1i*(400*t+30*cos(3*pi*t)));
 s  = s1+s2;
 s = s(:);
 Nfft = 512;
 s    = s(:);
 
 %we build the filter h
 if strcmp(window,'hamming')
  hlength=floor(161);
  hlength=hlength+1-rem(hlength,2);%the length of the filter has to be odd
  h = tftb_window(hlength,window);
  [hrow,~ ]=size(h); 
  Lh=(hrow-1)/2;
 else
  prec = 10^(-3);
  sigma_w = 0.03;
  L =  N*sigma_w;
  Lh = floor(L*sqrt(-log(prec)/pi))+1;
  h = amgauss(2*Lh+1,Lh+1,L);    
 end
 
 [tfr,norm2h] = tfrstft_three_case_down(s,Nfft,cas,h,Lh,downsamp,shift);
 close all
 figure()
 B = size(tfr);
 imagesc((0:B(2)-1)/B(2),0:B(1)/2-1,abs(tfr(1:B(1)/2,:)));
 set(gca,'ydir','normal');
 axis square
 [x] = itfrstft_three_case_down(tfr,cas,N,h,shift);
 figure()
 subplot(2,1,1)
 plot(t,imag(s),'--',t,imag(x))
 subplot(2,1,2)
 plot(t,real(s),'--',t,real(x))
 %test if the reconstruction is correct
 Y = max(abs(imag(s)-imag(x)))
 Z = max(abs(real(s)-real(x)))
end 