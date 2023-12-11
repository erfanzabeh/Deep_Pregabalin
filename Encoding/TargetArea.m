%% Comapre channels

clear all
clc
close all

%... Load iEEG.....
lfp = struct;
data = load('Complt_sig.mat');
lfp.pre = data.pre;
lfp.are = data.are;
lfp.gab = data.gab;
clear data1 data2 data3

%........... data properties
trl_num = size(lfp.pre,1);
cnhl_num = size(lfp.pre,2);
Fs = 500; %hz each 2ms
Time =  [1:size(lfp.pre,3)]/Fs;

% CLR{3} = [0.8,0.5,0.5];
% CLR{2} = [0.5,0.5,0.5];
CLR{1} = [0.5,0.8,0.8];
CLR{2} = [0.5,0.8,0.8]*0.5;
CLR{3} = [0.5,0.8,0.8]*0.25;
CLR{4} = [0.5,0.8,0.8]*0.1;

drug_strn{1} = 'Pregabalin';
drug_strn{2} = 'Gabapentin';
drug_strn{3} = 'Arecoline';
%% demo
trl = 901;
figure
for chnl = 1:2
    x = squeeze(lfp.pre(trl,chnl,:));
%     clr = [0, 1-1/(chnl), 1-1/(chnl)];
    clr = [0 0 0]
    plot(2*chnl+x,'Color',clr,'LineWidth',1.5);
    hold on
%     grid on
%     grid minor
end
 save2pdf(['waveform sample'],gcf,800)

%% Calculate Power Spectogeram
%.. Put all the trials together
power= [];
for chan = 1:4
            signal = lfp.pre;
    for trl = 1:trl_num/2
        trl
        %........................Scalogeram (Wavelet transform power)..........
        y = squeeze(signal(trl,chan,:));
        [cfs,frq] = cwt(y,Fs);
        power{chan}(:,:,trl) = (cfs);
    end
end

%% Plot power spectrum
h= figure(3);
scrsz = get(0,'ScreenSize');
scrsz(4) = scrsz(4)/2;
scrsz(3) = scrsz(3)/3;
set(h, 'Position',scrsz);
Ax = gcf;
Ax.Color = [1 1 1];

for chan = 1:4
    subplot(1,2,1)
    avg_pwr = abs(nanmean(power{chan},3))';
    [~,start_freq] = min(abs(frq-300)); %find the the element of 50Hz
    A=avg_pwr(:,start_freq:end)';
    
    clr = [0, 1-1/(chan), 1-1/(chan)];
    CLR{chan} = clr;
    plot(log2(frq(start_freq:end)), smooth(mean(A'),20),'-.','Color',clr,'LineWidth',2);
    % area(log(frq(start_freq:end)), mean(A'),...
    %     'EdgeColor', 'none', 'FaceColor', [0.5 0.5 0.5]);
    AX = gca;
    % set(gca,'XScale','log'); set(AX,'Xdir','normal');
    AX.XTickLabel = freq;
    freq = 2.^(round(log2(2)):round(log2(300)));
     AX.XTick =[1:8];
    AX.XTickLabel = freq;
    xlim(log2([3,256]))
    ylabel('Power (dB)');
    hold on;
    xlabel('Frequency')
    legend('Visual cortex',' Somatosensory cortex','Motor cortex','Olfactory bulbs',...
        'Location','northoutside')
    box off
    subplot(1,2,2)
    
    
    A=avg_pwr;
    [~,lw_f1]  = min(abs(frq-10));
    [~,frq50]  = min(abs(frq-50));
    low_pwr = mean(A(:,lw_f1:end),2);
    high_pwr = mean(A(:,frq50:lw_f1),2);
    
%     x = (log(1000*low_pwr));
%     y = (log(1000*high_pwr));
    
    LH(chan,:) = high_pwr./low_pwr;
    plot(chan,mean(LH(chan,:)),'O','LineWidth',3,'color',clr)
    hold on
    errorbar(chan,mean(LH(chan,:)),std(LH(chan,:)),'Color',[0.5 0.5 0.5],'LineWidth',1)
    box off
end

b = bar(mean(LH'),'FaceColor','flat');
b.CData(1,:) = CLR{1};
b.CData(2,:) = CLR{2};
b.CData(3,:) = CLR{3};
b.CData(4,:) = CLR{4};
b.EdgeColor = [0 0 0]
b.FaceAlpha = 0.6;

Ax = gca;
Ax.XTick = [1:4];
Ax.XTickLabel = {'Visual cortex';' Somatosensory cortex';'Motor cortex';'Olfactory bulbs'};
Ax.XTickLabelRotation = 45;
ylabel('High to Low frequency power ratio')
save2pdf(['electrode comparison'],gcf,800)

% plot([200,200],ylim,'--','Color', [1 0 0 0.5], 'LineWidth', 1.5)
% xlim([1,1400])
