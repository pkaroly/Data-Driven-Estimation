%% generateData
% plots simulated data from the neural mass model at different input values


%%
% Dean Freestone, Philippa Karoly 2016
% This code is licensed under the MIT License 2018

%%
clear
clc
close all

addpath(genpath('../src/'));

time = 60;
Fs = 1e3;
x = 1/Fs:1/Fs:time;
sigma_R = 0;
for input = 0:10:320
    [A,B,C,N_states,N_syn,N_inputs,N_samples,xi, ...
        v0,varsigma,Q,R,H,y] = set_params(input,[],time,Fs, sigma_R);
    plot(x,y,'k');
    set(gca,'box','off','xtick',[0 time]);
    xlabel('Time (s)');
    ylabel('ECoG (mV)');
    title(sprintf('input / drive: %d',input));
    drawnow;
    pause(0.5);
end

%%
function plotPotentials(xi)
figure
figure('name','parameter estimates' ,'units','normalized','position',[0 0 1 1] )
subplot(411),plot(xi(1,:))
title('Inhibitory -> Pyramidal');
subplot(412),plot(xi(3,:))
title('Pyramidal -> Inhibitory');
subplot(413),plot(xi(5,:))
title('Pyramidal -> Excitatory');
subplot(414),plot(xi(7,:))
title('Excitatory -> Pyramidal');
end

function plotAlpha(xi)
    figure('name','parameter estimates' ,'units','normalized','position',[0 0 1 1] )
    subplot(511),plot(xi(9,:))
    title('Input');
    subplot(512),plot(xi(10,:))
    title('Inhibitory -> Pyramidal');
    subplot(513),plot(xi(11,:))
    title('Pyramidal -> Inhibitory');
    subplot(514),plot(xi(12,:))
     title('Pyramidal -> Excitatory');
    subplot(515),plot(xi(13,:))
    title('Excitatory -> Pyramidal');
end