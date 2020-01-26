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

p8 = audioplayer(v,Fs);
playblocking(p8);

