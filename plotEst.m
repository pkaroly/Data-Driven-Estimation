clear
close all
clc

color_pk;

ParamName = {'input-Py','In-Py','Py-In','Py-Ex','Ex-Py'};

% Patients
Patient{1} = '23_002';
Patient{2} = '23_003';
Patient{3} = '23_004';
Patient{4} = '23_005';
Patient{5} = '23_006';
Patient{6} = '23_007';

Patient{7} = '24_001';
Patient{8} = '24_002';
Patient{9} = '24_004';
Patient{10} = '24_005';

Patient{11} = '25_001';
Patient{12} = '25_002';
Patient{13} = '25_003';
Patient{14} = '25_004';
Patient{15} = '25_005';

Fs = 399.6097561;
Tbefore = 30;

upperThresh = 25;
lowerThresh = 1;

%% plot stuff
col2 = [204, 0, 51]/255;
col1 = [0, 102, 204]/255;
font = 'arial';
fsize = 12;

for iSz = 1:13
    
    savename = sprintf('Sz%d.png',iSz);
    
    if exist(savename, 'file')
        continue;
    end
    
    %% plot
    count = 0;
    
    try
        for iParam = [1,2,5,4,3]
            count = count + 1;
            for iCh = 1:16
                
                load(sprintf('status seizures/estimates/Seizure%04d/Ch%d',iSz,iCh));
                y = xi_hat(8+iParam,:);
                
                y = log(y.^2);
                background = y(1:round(Tbefore*Fs));
                segment = y(round(Tbefore*Fs):end);
                
                if iCh == 1
                    Len = size(segment,2);
                    chMat = zeros(16,Len);
                    t = 1/Fs:1/Fs:Len/Fs;
                end
                
                % thresh
                %         segment(segment < lowerThresh & segment ~= -inf) = lowerThresh;
                %         background(background < lowerThresh & background ~= -inf) = lowerThresh;
                %         segment(segment > upperThresh) = upperThresh;
                %         background(background  > upperThresh) = upperThresh;
                
                % background normalise
                mu = mean(background);
                segment = (segment - mu) / mu;
                
                chMat(iCh,:) = segment;
                
            end
            chMat(isnan(chMat)) = 0;
            
            subplot(5,1,count);
            imagesc(chMat);
            colormap(cmap);
            caxis([-0.5 0.5]);
            
            set(gca,'box','off','xtick',[1, Len],'xticklabel',[],...
                'ytick',[1 16],'yticklabel',[],...
                'fontname',font,'fontsize',fsize)
            title(ParamName{count})
            if count == 5
                if t(end) > 600
                    set(gca,'xticklabel',round([0 t(end)]/60),'yticklabel',[1 16])
                    xlabel('Duration (m)');
                else
                    set(gca,'xticklabel',round([0 t(end)]),'yticklabel',[1 16])
                    xlabel('Duration (s)');
                end
                ylabel('Channel');
            end
            drawnow
        end
        %%
        set(gcf,'paperunits','centimeters','paperposition',[0 0 20 20])
        print(gcf,'-dpng',savename);
        
    catch
        continue;
    end
    
end