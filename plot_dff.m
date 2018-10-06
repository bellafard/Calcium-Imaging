clear; clc; close all;

%% Adding paths
addpath(genpath('/media/s2v2/Arash/caim/Sources/CNMF_E/ca_source_extraction'));
load('/media/s2v2/Arash/C142/2p/2018_09_17_am/20180917_13_13_40_c142_s01/20180917_13_13_40_c142_s01_XYT_rig_source_extraction/frames_1_24000/LOGS_02-Oct_12_24_07/02-Oct_14_23_21.mat');

%% 
Tr = neuron.C_raw;
Tp = neuron.C;
N  = neuron.A;

%%
n = 50;

[yr,yp,c] = plot_activity(n,Tr,Tp,N);
x = [0:1/30:(length(yr)-1)/30];

figure(1);
plot(x,yr);
xlabel('Time [Sec]')
hold on;
plot(x,yp);
hold off;

figure(2);
imagesc(c);

%%
function [a,b,c] = plot_activity(n,Tr,Tp,N)
    a = Tr(n,:);
    b = Tp(n,:);
    c = reshape(N(:,n),[512,512]);
end


