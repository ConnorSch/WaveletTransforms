% Connor Schleicher AMATH 582 HW 2

% initialize workspace
clear all; close all; clc;

% initialize recording from HW prompt
load handel
v = y'/2;
plot((1:length(v))/Fs,v);
xlabel('Time [sec]');
ylabel('Amplitude');
title('Signal of Interest, v(n)');

% p8 = audioplayer(v,Fs);
% playblocking(p8);

% initialize fourier modes 
L = length(v)/Fs;
n = length(v);
t2=linspace(0,L,n+1); t=t2(1:n); 
k=(2*pi/L)*[0:n/2 -n/2:-1];  
ks=fftshift(k);


%%
slidew = [2, 5, 7];
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
       subplot(3,1,1), plot(t,v,'k',t,g,'r')
       subplot(3,1,2), plot(t,vg,'k')
       subplot(3,1,3), plot(ks, abs(fftshift(vgt))/max(abs(vgt)))
       drawnow
       %pause(0.05)  
    end
figure(i+2)
pcolor(slidet,ks,Specto((i-1)*length(slidet) + 1:i*length(slidet),:).'),
titletext = 'Spectogram at filter width ' + string(slidew(i));
title(titletext)
shading interp
colormap(hot)
end

%%
figure()
tr_piano=16; % record time in seconds
y_p=audioread('music1.wav'); Fs_p=length(y_p)/tr_piano;
plot((1:length(y_p))/Fs_p,y_p);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Mary had a little lamb (piano)'); drawnow
%p8 = audioplayer(y,Fs); playblocking(p8);

figure()
tr_rec=14; % record time in seconds
y_r=audioread('music2.wav'); Fs_r=length(y_r)/tr_rec;
plot((1:length(y_r))/Fs_r,y_r);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Mary had a little lamb (recorder)');
%p8 = audioplayer(y,Fs); playblocking(p8);

%%
S_p = y_p'/2;
S_r = y_r'/2;

L_p = length(S_p)/Fs_p;
n_p = length(S_p);
t2_p=linspace(0,L_p,n_p+1); t_p=t2_p(1:n_p); 
k_p=(2*pi/L_p)*[0:n_p/2-1 -n_p/2:-1];  
ks_p=fftshift(k_p);

L_r = length(S_r)/Fs_r;
n_r = length(S_r);
t2_r = linspace(0,L_r,n_r+1); t_r = t2_r(1:n_r);
k_r = (2*pi/L_r)*[0:n_r/2-1 -n_r/2:-1];
ks_r = fftshift(k_r);

slidet_p = 0:1/tr_piano:tr_piano;
Specto_p = [];
figure()
for p = 1:length(slidet_p)
   g_p = exp(-0.5*(t_p-slidet_p(p)).^2); % Gabor transform
   Sg_p = g_p.*S_p; % apply the Gabor transform
   Sgt_p = fft(Sg_p); % take the Fourier transform 
   Specto_p = [Specto_p; abs(fftshift(Sgt_p))]; %storing data for plotting
   subplot(3,1,1), plot(t_p,S_p,'k',t_p,g_p,'r')
   subplot(3,1,2), plot(t_p,Sg_p,'k')
   subplot(3,1,3), plot(ks_p, abs(fftshift(Sgt_p))/max(abs(Sgt_p)))
   drawnow
end
figure()
pcolor(slidet_p,ks_p,Specto_p.'),
title('Spectogram of piano')
shading interp
colormap(hot)

slidet_r = 0:1/tr_rec:tr_rec;
Specto_r = [];
figure()
for r = 1:length(slidet_r)
   g_r = exp(-0.5*(t_r-slidet_r(r)).^2); % Gabor transform
   Sg_r = g_r.*S_r; % apply the Gabor transform
   Sgt_r = fft(Sg_r); % take the Fourier transform 
   Specto_r = [Specto_r; abs(fftshift(Sgt_r))]; %storing data for plotting
   subplot(3,1,1), plot(t_p,S_p,'k',t_p,g_p,'r')
   subplot(3,1,2), plot(t_p,Sg_p,'k')
   subplot(3,1,3), plot(ks_p, abs(fftshift(Sgt_p))/max(abs(Sgt_p)))
end
figure ()
pcolor(slidet_r,ks_r,Specto_r.'),
title('Spectogram of recorder')
shading interp
colormap(hot)