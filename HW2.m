% Connor Schleicher AMATH 582 HW 2

% initialize workspace
clear all; close all; clc;

% initialize recording from HW prompt
load handel
v = y'/2;
figure(1)
subplot(2,1,1)
plot((1:length(v))/Fs,v);
xlabel('Time [sec]');
ylabel('Amplitude');
title('Signal of Interest, v(n)');

% initialize fourier modes 
L = length(v)/Fs;
n = length(v);
t2=linspace(0,L,n+1); t=t2(1:n); 
k=(2*pi/L)*[0:n/2 -n/2:-1];  
ks=fftshift(k);


subplot(2,1,2)
vt = fft(v);
plot(ks,abs(fftshift(vt))/max(abs(vt)));
ylabel('FFT(v)'), xlabel('Frequency (\omega)')
title('FFT of Signal');

% p8 = audioplayer(v,Fs);
% playblocking(p8);

% initialize fourier modes 
L = length(v)/Fs;
n = length(v);
t2=linspace(0,L,n+1); t=t2(1:n); 
k=(2*pi/L)*[0:n/2 -n/2:-1];  
ks=fftshift(k);


%% Spectogram of Handel with Gabor Transform Varying the Window Width

slidew = [5, 30, 70];
slidet = 0:1/L:L;
Specto = [];
%Specto = zeros(length(slidew)*length(slidet),length(v));
for i = 1:length(slidew)
    figure(2)
    for j = 1:length(slidet)
       g = exp(-slidew(i)*(t-slidet(j)).^2); % Gabor transform
       vg = g.*v; % apply the Gabor transform
       vgt = fft(vg); % take the Fourier transform 
       Specto = [Specto; abs(fftshift(vgt))]; %storing data for plotting
       subplot(3,1,1), plot(t,v,'k',t,g,'r'), title('Gabor Filtering and signal'), legend('v','Gabor filter')
       xlabel('Time [sec]'), ylabel('Amplitude')
       subplot(3,1,2), plot(t,vg,'k'), title('Gabor Filter of Signal')
       xlabel('Time [sec]'), ylabel('Amplitude')
       subplot(3,1,3), plot(ks, abs(fftshift(vgt))/max(abs(vgt))), title('Transformation of Signal')
       xlabel('Frequency (\omega)'), ylabel('FFT')
       drawnow
       %pause(0.05)  
    end
figure(3)
subplot(1,3,i)
pcolor(slidet,ks,Specto((i-1)*length(slidet) + 1:i*length(slidet),:).'),
titletext = 'Window Width: ' + string(slidew(i));
title(titletext), xlabel('Time [sec]'), ylabel('Frequency [Hz]')
%set(gca,'Ylim',[0 1500],'Fontsize',[14]) 
shading interp
colormap(hot)
end

%% Spectogram of Handel with Gabor Transform Varying the Sampling Rate
slidew = 70;
sampling_rate = [30, 80, 150];
Specto = [];
%Specto = zeros(length(slidew)*length(slidet),length(v));
for i = 1:length(sampling_rate)
    slidet = linspace(0,L,sampling_rate(i));
    figure(3)
    Specto = [];
    for j = 1:length(slidet)
       g = exp(-slidew*(t-slidet(j)).^2); % Gabor transform
       vg = g.*v; % apply the Gabor transform
       vgt = fft(vg); % take the Fourier transform 
       Specto = [Specto; abs(fftshift(vgt))]; %storing data for plotting
       subplot(3,1,1), plot(t,v,'k',t,g,'r'), title('Gabor Filtering and signal'), legend('v','Gabor filter')
       xlabel('Time [sec]'), ylabel('Amplitude')
       subplot(3,1,2), plot(t,vg,'k'), title('Gabor Filter of Signal')
       xlabel('Time [sec]'), ylabel('Amplitude')
       subplot(3,1,3), plot(ks, abs(fftshift(vgt))/max(abs(vgt))), title('Transformation of Signal')
       xlabel('Frequency (\omega)'), ylabel('FFT')
       drawnow
       %pause(0.05)  
    end
    
    figure(4)
    subplot(1,3,i)
    pcolor(slidet,ks,Specto.'),
    titletext = 'Sampling Rate: ' + string(sampling_rate(i));
    title(titletext), xlabel('Time [sec]'), ylabel('Frequency [Hz]')
    %set(gca,'Ylim',[0 1500],'Fontsize',[14]) 
    shading interp
    colormap(hot)
end

%% Spectogram of Handel with mexican hat wavelet

width = 70;
slidet = 0:1/L:L;
SpectoMH = zeros(length(slidet),length(v));
%Specto = zeros(length(slidew)*length(slidet),length(v));
figure(5)
for j = 1:length(slidet)
   mh = (1-width*t.^2).*exp(-width*(t-slidet(j)).^2); % mexican hat wavelet
   mh = mh/max(abs(mh)); % scale to 1
   vmh = mh.*v; % apply the Gabor transform
   vmht = fft(vmh); % take the Fourier transform 
   SpectoMH(j,:) = abs(fftshift(vmht)); %storing data for plotting
   subplot(3,1,1), plot(t,v,'k',t,mh,'r'), title('Gabor Filtering and signal'), legend('v','Gabor filter')
   xlabel('Time [sec]'), ylabel('Amplitude'), ylim([-1 1])
   subplot(3,1,2), plot(t,vmh,'k'), title('Gabor Filter of Signal')
   xlabel('Time [sec]'), ylabel('Amplitude')
   subplot(3,1,3), plot(ks, abs(fftshift(vmht))/max(abs(vmht))), title('Transformation of Signal')
   xlabel('Frequency (\omega)'), ylabel('FFT')
   drawnow
   %pause(0.05)  
end

%% Spectogram of Handel with Shannon wavelet

width = 700;
slidet = 0:1/L:L;
SpectoS = zeros(length(slidet),length(v));
figure(6)
for j = 1:length(slidet)
   k = j-1; 
   num_steps = floor(length(v)/length(slidet));
   s = zeros(1,length(t));
   s((k*num_steps)+1:k*num_steps+width) = 1;
   vs = s.*v; % apply the Gabor transform
   vst = fft(vs); % take the Fourier transform 
   SpectoS(j,:) = abs(fftshift(vst)); %storing data for plotting
   subplot(3,1,1), plot(t,v,'k',t,s,'r'), title('Gabor Filtering and signal'), legend('v','Gabor filter')
   xlabel('Time [sec]'), ylabel('Amplitude')
   subplot(3,1,2), plot(t,vs,'k'), title('Gabor Filter of Signal')
   xlabel('Time [sec]'), ylabel('Amplitude')
   subplot(3,1,3), plot(ks, abs(fftshift(vst))/max(abs(vst))), title('Transformation of Signal')
   xlabel('Frequency (\omega)'), ylabel('FFT')
   drawnow
   %pause(0.05)  
end

%%
figure(6)
subplot(1,2,1)
pcolor(slidet,ks,SpectoMH.'),
title('Spectogram with Mexican Hat Wavelet')
xlabel('Time [sec]'), ylabel('Frequency [Hz]')
shading interp
colormap(hot)

subplot(1,2,2)
pcolor(slidet,ks,SpectoS.'),
title('Spectogram with Shannon filter')
xlabel('Time [sec]'), ylabel('Frequency [Hz]')
shading interp
colormap(hot)

%% Frequency graphs of piano and recorder
figure()
subplot(2,2,1)
tr_piano=16; % record time in seconds
y_p=audioread('music1.wav'); Fs_p=length(y_p)/tr_piano;
plot((1:length(y_p))/Fs_p,y_p);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Mary had a little lamb (piano)'); drawnow
%p8 = audioplayer(y,Fs); playblocking(p8);

subplot(2,2,3)
tr_rec=14; % record time in seconds
y_r=audioread('music2.wav'); Fs_r=length(y_r)/tr_rec;
plot((1:length(y_r))/Fs_r,y_r);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Mary had a little lamb (recorder)');
%p8 = audioplayer(y,Fs); playblocking(p8);

S_p = y_p'/2;
L_p = length(S_p)/Fs_p;
n_p = length(S_p);
t2_p=linspace(0,L_p,n_p+1); t_p=t2_p(1:n_p); 
k_p=(2*pi/L_p)*[0:n_p/2-1 -n_p/2:-1];  
ks_p=fftshift(k_p);
St_p = fft(S_p);

S_r = y_r'/2;
L_r = length(S_r)/Fs_r;
n_r = length(S_r);
t2_r = linspace(0,L_r,n_r+1); t_r = t2_r(1:n_r);
k_r = (2*pi/L_r)*[0:n_r/2-1 -n_r/2:-1];
ks_r = fftshift(k_r);
St_r = fft(S_r);

subplot(2,2,2)
plot(ks_p,fftshift(St_p)/max(abs(St_p)))
ylabel('FFT(v)'), xlabel('Frequency (\omega)')
title('FFT of Piano')%, xlim([-0.5 0.5])

subplot(2,2,4)
plot(ks_r,fftshift(St_r)/max(abs(St_r)))
ylabel('FFT(v)'), xlabel('Frequency (\omega)')
title('FFT of Recorder')%, xlim([-0.5 0.5])

%% Spectorgram of Piano

numstep = 100;
width = 150;

slidet_p = linspace(0,tr_piano,numstep);
Specto_p = zeros(length(slidet_p),length(y_p));
figure()
for p = 1:length(slidet_p)
   g_p = exp(-width*(t_p-slidet_p(p)).^2); % Gabor transform
   Sg_p = g_p.*S_p; % apply the Gabor transform
   Sgt_p = fft(Sg_p); % take the Fourier transform 
   Specto_p(p,:) = abs(fftshift(Sgt_p)); %storing data for plotting
   subplot(3,1,1), plot(t_p,S_p,'k',t_p,g_p,'r'), title('Gabor Filtering and signal'), legend('v','Gabor filter')
   xlabel('Time [sec]'), ylabel('Amplitude')
   subplot(3,1,2), plot(t_p,Sg_p,'k'), title('Gabor Filter of Signal')
   xlabel('Time [sec]'), ylabel('Amplitude')
   subplot(3,1,3), plot(ks_p, abs(fftshift(Sgt_p))/max(abs(Sgt_p))), title('Transformation of Signal')
   xlabel('Frequency (\omega)'), ylabel('FFT')
   drawnow
end

figure()
pcolor(slidet_p,ks_p,log(Specto_p.'+1)), shading interp
xlabel('Time [sec]'), ylabel('Frequency [Hz]'),title('Spectogram of piano')
set(gca,'Ylim',[1000 2500],'Fontsize',[14]) 
colormap(hot)

%% Spectogram of Recorder

numstep = 100;
width = 150;

slidet_r = linspace(0,tr_rec,numstep);
Specto_r = zeros(length(slidet_r),length(y_r));
figure()
for r = 1:length(slidet_r)
   g_r = exp(-width*(t_r-slidet_r(r)).^2); % Gabor transform
   Sg_r = g_r.*S_r; % apply the Gabor transform
   Sgt_r = fft(Sg_r); % take the Fourier transform 
   Specto_r(r,:) = abs(fftshift(Sgt_r)); %storing data for plotting
   subplot(3,1,1), plot(t_r,S_r,'k',t_r,g_r,'r'), title('Gabor Filtering and signal'), legend('v','Gabor filter')
   xlabel('Time [sec]'), ylabel('Amplitude')
   subplot(3,1,2), plot(t_r,Sg_r,'k'), title('Gabor Filter of Signal')
   xlabel('Time [sec]'), ylabel('Amplitude')
   subplot(3,1,3), plot(ks_r, abs(fftshift(Sgt_r))/max(abs(Sgt_r))), title('Transformation of Signal')
   xlabel('Frequency (\omega)'), ylabel('FFT')
end

figure ()
pcolor(slidet_r,ks_r,log(Specto_r.'+1)), shading interp
xlabel('Time [sec]'), ylabel('Frequency [Hz]'),title('Spectogram of recorder')
set(gca,'Ylim',[4000 8000],'Fontsize',[14]) 
colormap(hot)

%% Comparison of full spectogram for piano and recorder

figure ()
subplot(2,1,1)
pcolor(slidet_p,ks_p,log(Specto_p.'+1)), shading interp
xlabel('Time [sec]'), ylabel('Frequency [Hz]'),title('Spectogram of piano')
set(gca,'Ylim', [0 (10^5)/2],'Fontsize',[14]) 
colormap(hot)

subplot(2,1,2)
pcolor(slidet_r,ks_r,log(Specto_r.'+1)), shading interp
xlabel('Time [sec]'), ylabel('Frequency [Hz]'),title('Spectogram of recorder')
set(gca,'Ylim', [0 (10^5)/2],'Fontsize',[14]) 
colormap(hot)
