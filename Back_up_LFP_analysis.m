%% LFP analysis
clear all
clc
close all
%... Load LFP.....
LFP = struct;
data = load('ampSortLFPMs.mat');
LFP.LCstim = data.ampSortLFPMs{1, 12}{1, 3}(:,:);

trial_num = size(LFP.LCstim,1);
Time =  1:size(LFP.LCstim,2);
%% Power Spectogeram Analysis
%.. Put all the trials together
power= [];
for trial = 1:trial_num
    trial
    %........................Scalogeram (Wavelet transform power)..........
    fs = 1000;% Sampling frequency
    y = LFP.LCstim(trial,:);
    [cfs,frq] = cwt(y,fs);
    power(:,:,trial) = (cfs);
end
%% Single average trial

figure
plot(mean(LFP.LCstim),'-k','LineWidth',1.5)



%% Single trial Example
figure
for itter = [1,2,3,4,5,95,96,97,98,99,100]
    % plot(LFP.LCstim(itter,:)+0.2*itter)
    subplot(5,2,(mod(itter,5)+1)*2-(itter>50))
    eqlzr1 = min(min(abs(nanmean(power(:,:,itter)))));
    eqlzr2 = max(max(abs(nanmean(power(:,:,itter)))));
    eqlzr3 = mean(mean(abs(nanmean(power(:,1:200,itter)))));
    
    avg_pwr = abs((power(:,:,itter)))';
    baseline = mean(avg_pwr(1:100,:));
    avg_pwr = avg_pwr - ones(1401,1)*baseline;
    fs = 1000;
    [~,start_freq] = min(abs(frq-300)); %find the the element of 50Hz
    A=avg_pwr(:,start_freq:end)';
    % K = 4^-1*ones(2,2);
    % A = conv2(A,K,'same');
    pcolor(Time,frq(start_freq:end),A);
    shading interp
    
    AX = gca;
    set(gca,'YScale','log');
    set(AX,'Ydir','normal');
    hold on
    freq = 2.^(round(log2(2)):round(log2(300)));
    AX.YTickLabelMode = 'auto';
    AX.YTick = freq;
    AX.XTick = []
    AX.YTick = [];
    % TimePoints = [-1:0.5:2];
    % ylabel('Frequency')
    % xlabel('Time')
    colormap(gca,brewermap(256, '*RdYlBu'));
    % colormap(gca,'jet');
    caxis([-max(eqlzr1,eqlzr2),max(eqlzr1,eqlzr2)])
    
    % axis square
end
% close all
%%
% figure
% counter = 1;
% for itter = [1,2,3,4,5,95,96,97,98,99,100]
%     counter = counter+1;
%     plot(LFP.LCstim(itter,:)+0.2*counter)
%     hold on
% end
% 
% figure
% subplot(1,2,1);
%     plot(Time,mean(LFP.LCstim(:,:)),'-o')
%     hold on
%     xlim([-800,200])
%     subplot(1,2,2);
%     plot(Time,mean(LFP.LCstim(:,:)),'-o')
%     hold on
%     xlim([500,1500])
    
%     errorbar(Time,mean(LFP.LCstim(:,:)),std(LFP.LCstim(:,:))/sqrt(size(LFP.LCstim,1)))
%% 
% figure
% trial = 93;
% subplot(1,2,1);
%     plot(Time,(LFP.LCstim(trial,:)),'-k','LineWidth',1.5)
%     hold on
%     xlim([-800,200])
%     ylim([-0.2 0.2])
%     subplot(1,2,2);
%     
%     plot(Time,(LFP.LCstim(trial,:)),'-k',...
%         'LineWidth',1.5)
%     hold on
%     xlim([500,1500])
%        ylim([-0.2 0.2])
% close all
%%
%  plot(Time,(LFP.LCstim(93,:)),'-k','LineWidth',1.5)
%  hold on 
% yl = ylim;
% upb = [yl(2), yl(2)];%Vector of random data
% lowb = [yl(1),yl(1)];%2nd vector of data points;
% jbfill([200,300],upb,lowb,[0.5 0.5 0.5],[0 0 0],0,1)
% ylim(yl)
% 
% hold on 
% 
% upb = [0.1, 0.1];%Vector of random data
% lowb = [-0.2,-0.2];%2nd vector of data points;
% jbfill([400,1400],upb,lowb,[0 0 1],[1 1 1],0,0.1)
% 
% upb = [0.05, 0.05];%Vector of random data
% lowb = [-0.1,-0.1];%2nd vector of data points;
% jbfill([0,200],upb,lowb,[1 0 0],[1 1 1],0,0.1)
% xlim([0 1400])
% ylabel('Power (v)')
% xlabel('Time (Stimulus onset aligned)')
% save2pdf(['Single trial example'],gcf,800)
%% Plot Power Spectrum

h= figure(2);
scrsz = get(0,'ScreenSize');
scrsz(4) = scrsz(4)/1.5;
scrsz(3) = scrsz(3)/2;
set(h, 'Position',scrsz);

ax_main = axes('Position', [0.1 0.2 0.55 0.55]);
% ax_right = axes('Position', [0.7 0.2 0.1 0.55]);
ax_top = axes('Position', [0.1 0.8 0.55 0.1]);

axes(ax_main)
eqlzr1 = min(min(abs(nanmean(power(:,:,:),3))));
eqlzr2 = max(max(abs(nanmean(power(:,:,:),3))));
eqlzr3 = mean(mean(abs(nanmean(power(:,1:200,:),3))));

avg_pwr = abs(nanmean(power,3))';
baseline = mean(avg_pwr(1:100,:));
avg_pwr = avg_pwr - ones(1401,1)*baseline;
fs = 1000;
[~,start_freq] = min(abs(frq-300)); %find the the element of 50Hz
A=avg_pwr(:,start_freq:end)';
K = 16^-1*ones(10,2);
A = conv2(A,K,'same');
pcolor(Time,frq(start_freq:end),A);
shading interp

AX = gca;
set(gca,'YScale','log');
set(AX,'Ydir','normal');
hold on
freq = 2.^(round(log2(2)):round(log2(300)));
AX.YTickLabelMode = 'auto';
AX.YTick = freq;
% TimePoints = [-1:0.5:2];
ylabel('Frequency')
xlabel('Time from LC stimulation')
colormap(gca,brewermap(256, '*RdBu'));
% colormap(gca,'jet');
caxis([-max(eqlzr1,eqlzr2),max(eqlzr1,eqlzr2)])
% caxis([eqlzr1,eqlzr2])

% ......Colorbar.........
h = colorbar(ax_main, 'manual', 'Position', [0.85 0.2 0.02 0.55]);
ylabel(h,'Power');
axes(ax_main);
% hold on;
% [main_x,main_y] = max(x);
% xl = xlim;
% plot([xl(1) xl(2)],[freq_axis(right_y) freq_axis(right_y)],...
%     'Color', [1 0 0 0.2], 'LineWidth', 1.5);
hold on
yl = ylim;
Time_Onset = 200;
plot([Time(Time_Onset(1)),Time(Time_Onset(1))],[yl(1) yl(2)],'--',...
    'Color', [1 0 0 0.5], 'LineWidth', 1.5);

%........... Top Axis.........
axes(ax_top);
area(Time, mean(A),...
    'EdgeColor', 'none', 'FaceColor', [0.5 0.5 0.5]);
box off;
ax_top.XTickLabel = [];
ylabel('Power (dB)');
hold on;
ylm =ylim
plot([200,200],ylim,'--','Color', [1 0 0 0.5], 'LineWidth', 1.5)
xlim([1,1400])

% save2pdf(['LFP spectogram post stimulation'],gcf,600)
% print(,'-dpdf',pdfFileName,sprintf('-r%d',dpi))
%% Compare early and late trials

power_early= [];
early_num = 50;
counter = 0;
for trial = 1:early_num
    counter = 1+counter;
    %........................Scalogeram (Wavelet transform power)..........
    fs = 1000;% Sampling frequency
    y = LFP.LCstim(trial,:);
    [cfs,frq] = cwt(y,fs);
    power_early(:,:,counter) = (cfs);
end

power_late= [];
late_num = 50;
counter = 0;
for trial = late_num:trial_num
    counter = 1+counter;
    %........................Scalogeram (Wavelet transform power)..........
    fs = 1000;% Sampling frequency
    y = LFP.LCstim(trial,:);
    [cfs,frq] = cwt(y,fs);
    power_late(:,:,counter) = (cfs);
end


figure
for type =1:2 % early and late trials
    
    if type==1
        pwr = power_early;
    else
        pwr = power_late;
    end
    
    subplot(1,3,type)
    eqlzr1 = min(min(abs(nanmean(pwr,3))));
    eqlzr2 = max(max(abs(nanmean(pwr,3))));
    eqlzr3 = mean(mean(abs(nanmean(pwr(:,1:200,:)))));
    
    avg_pwr = abs(nanmean(pwr,3))';
    baseline = mean(avg_pwr(1:100,:));
    avg_pwr = avg_pwr - ones(1401,1)*baseline;
    fs = 1000;
    [~,start_freq] = min(abs(frq-300)); %find the the element of 50Hz
    A=avg_pwr(:,start_freq:end)';
    % K = 4^-1*ones(2,2);
    % A = conv2(A,K,'same');
    pcolor(Time,frq(start_freq:end),A);
    shading interp
    
    AX = gca;
    set(gca,'YScale','log');
    set(AX,'Ydir','normal');
    hold on
    freq = 2.^(round(log2(2)):round(log2(300)));
    AX.YTickLabelMode = 'auto';
    AX.YTick = freq;
    AX.XTick = []
    AX.YTick = [];
    % TimePoints = [-1:0.5:2];
    % ylabel('Frequency')
    % xlabel('Time')
    % colormap(gca,brewermap(256, '*RdYlBu'));
    colormap(gca,'jet');
    caxis([-max(eqlzr1,eqlzr2),max(eqlzr1,eqlzr2)])
    
    % axis square
end
close


figure(2)
pwr = mean((power_early),3)' - mean((power_late),3)';
eqlzr1 = min(min(abs(pwr)));
eqlzr2 = max(max(abs(pwr)));

avg_pwr = abs(pwr)';
% baseline = mean(avg_pwr(:,1:100));
% avg_pwr = avg_pwr - ones(1401,1)*baseline;

fs = 1000;
[~,start_freq] = min(abs(frq-300)); %find the the element of 50Hz
A=avg_pwr(start_freq:end,:);
% K = 4^-1*ones(2,2);
% A = conv2(A,K,'same');
pcolor(Time,frq(start_freq:end),A);
shading interp

AX = gca;
set(gca,'YScale','log');
set(AX,'Ydir','normal');
hold on
freq = 2.^(round(log2(2)):round(log2(300)));
AX.YTickLabelMode = 'auto';
AX.YTick = freq;
AX.XTick = []
AX.YTick = [];
% TimePoints = [-1:0.5:2];
% ylabel('Frequency')
% xlabel('Time')
colormap(gca,brewermap(256, '*RdBu'));
% colormap(gca,'jet');
caxis([-max(eqlzr1,eqlzr2),max(eqlzr1,eqlzr2)])
colorbar
%  axis square
%% theta oscillation time modulation 
figure(3)

subplot(1,3,[1:2])
%................For Gamma Oscillation.......
% start_freq = 50;
% end_freq = 80;
%.............................


%............. For theta Oscillation...........
start_freq = 1;
end_freq = 4;
band = intersect(find(frq<end_freq),find(frq>start_freq));

start_time = 1;
end_time = 1400;
time_step = 1;

power_band_early = squeeze(mean(abs(power_early(band,start_time:time_step:end_time,:)),1));
power_band_late = squeeze(mean(abs(power_late(band,start_time:time_step:end_time,:)),1));

mean_early = mean(power_band_early,2);
std_early = std(power_band_early')'/sqrt(size(power_band_early,2));

p1 = plot(Time(start_time:time_step:end_time),mean_early,'Color',[0.5 0.5 0.5],'LineWidth',1.5)
hold on
a= mean_early-std_early;
b = mean_early+std_early;
[ph,msg]=jbfill(Time(start_time:time_step:end_time),a',b',[0.5 0.5 0.5],[0 0 0],0,0.2);




mean_late = mean(power_band_late,2);
std_late = std(power_band_late')'/sqrt(size(power_band_late,2));


p2 = plot(Time(start_time:time_step:end_time),mean_late,'Color',[0 0.5 0],'LineWidth',1.5)
a= mean_late-std_late;
b = mean_late+std_late;
[ph,msg]=jbfill(Time(start_time:time_step:end_time),a',b',[0 0.5 0],[0 0 0],0,0.2);

ylm = ylim;
plot([200,200],ylm,'--','Color',[0.5 0 0 0.5],'LineWidth',1.5)
ylim(ylm)
legend([p1,p2],{'Early (Low NE level)','Late (High NE level)'},'Location','southwest')
ylabel('Theta band average power')
xlabel('Time from LC stimulation (ms)')
xlim([1 1400])


%.. Show the individual theta differeneces are correlated with temporal differences
subplot(1,3,[3])
%..... For theta
start_freq = 1;
end_freq = 4;
band = intersect(find(frq<end_freq),find(frq>start_freq));
theta_power = []
% power_early
trial_time = [1:100].*1.4*13/60;

for itter = 1:100
    theta_power(itter) = mean(mean(abs(power(band,400:1400,itter))));
%     if itter ==13
%         theta_power(13) = nan;
% trial_time(13) = nan;
%     end
    hold on 
    p2 = plot(trial_time(itter),theta_power(itter),'Ok','MarkerFaceColor',[0 itter/100 0])
end


hold on
[brob,stats] = robustfit(trial_time,theta_power);
x = trial_time;
p1 = plot(x,brob(1)+brob(2)*x,'-','LineWidth',1.5,'Color',[0.9 0 0]);
xlm = xlim;
plot(xlm,[mean(theta_power),mean(theta_power)],':','LineWidth',1,'Color',[0.5 0.5 0.5]);
grid on 
xlim([0,30])
xlabel('Trial Latency(mintue)')
legend([p1,p2],{['r =',num2str(floor(100*stats.coeffcorr(1,2))/100)],'trial'})
% save2pdf(['Ne level effect'],gcf,800)