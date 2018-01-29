clear all
clc
close all

time = 5;
Fs = 0.4e3;
x = 1/Fs:1/Fs:time;
for input = 0:50:1000
    [A,B,C,N_states,N_syn,N_inputs,N_samples,xi, ...
        v0,varsigma,Q,R,H,y] = set_params(input,[],time,Fs);
    plot(x,y,'k');
    set(gca,'box','off','xtick',[0 time]);
    xlabel('Time (s)');
    ylabel('ECoG (mV)');
    title(sprintf('input: %d',input));
    pause(0.5);
    drawnow;
end