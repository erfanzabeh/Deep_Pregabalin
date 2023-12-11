%% General analysis

clear all
clc
close all
%... Load iEEG.....
lfp = struct;
data = load('Complt_sig.mat');
lfp.pre = data.pre;
% lfp.phe = data.phe;
lfp.are = data.are;
lfp.gab = data.gab;

clear data1 data2 data3

trl_num = size(lfp.pre,1);
cnhl_num = size(lfp.pre,2);
Fs = 500; %hz each 2ms
Time =  [1:size(lfp.pre,3)]/Fs;

%% demo
trl = 901;
figure
for chnl = 1:4
    x = squeeze(lfp.pre(trl,chnl,:));
    clr = [0, 1-1/(chnl), 1-1/(chnl)];
    plot(chnl+x,'Color',clr);
    hold on
end
%% Drug comparison

h= figure(2);
scrsz = get(0,'ScreenSize');
scrsz(4) = scrsz(4)/1.5;
scrsz(3) = scrsz(3)/1;
set(h, 'Position',scrsz);
Ax = gcf;
Ax.Color = [1 1 1];

CLR{3} = [0.8,0.5,0.5];
CLR{2} = [0.5,0.5,0.5];
CLR{1} = [0.5,0.8,0.8];

drug_strn{1} = 'Pregabalin';
drug_strn{2} = 'Gabapentin';
drug_strn{3} = 'Arecoline';

trl = 90;
chnl = 2;

subplot(3,3,[1,2]) %Pregabalin
x = squeeze(lfp.pre(trl,chnl,:));
clr = CLR{1};

plot(Time,x,'Color',clr,'LineWidth',2);
title(drug_strn{1})
box off

subplot(3,3,[4,5]) %Phenazepam
x = squeeze(lfp.gab(trl,chnl,:));
clr = CLR{2} ;
plot(Time,x,'Color',clr,'LineWidth',2);
title(['Similar Drug (',drug_strn{2},')'])
box off

subplot(3,3,[7,8]) %Arecoline
x = squeeze(lfp.are(trl,chnl,:));
clr = CLR{3};
plot(Time,x,'Color',clr,'LineWidth',2);
title('Other Drug (Arecoline)')
box off

%% Calculate Power Spectogeram
%.. Put all the trials together
power= [];
for drug = 1:3
    
    switch drug
        case 1
            signal = lfp.pre;
        case 2
            signal = lfp.gab;
        case 3
            signal = lfp.are;
    end
    
    for trl = 1:trl_num/5
        trl
        %........................Scalogeram (Wavelet transform power)..........
        y = squeeze(signal(trl,chnl,:));
        [cfs,frq] = cwt(y,Fs);
        power{drug}(:,:,trl) = (cfs);
    end
end

%% Plot Power Spectrum

figure(2);
EQLZR1 = [];
EQLZR2 = [];

for drug = 1:3
    subplot(3,3,3*drug)
    
    avg_pwr = abs(nanmean(power{drug},3))';
    [~,start_freq] = min(abs(frq-300)); %find the the element of 50Hz
    A=avg_pwr(:,start_freq:end)';
    % K = 16^-1*ones(10,2); A = conv2(A,K,'same');
    pcolor(Time,frq(start_freq:end),A);
    
    eqlzr1 = max(max(A));
    eqlzr2 = min(min(A));
    
    shading interp
    AX = gca; set(gca,'YScale','log'); set(AX,'Ydir','normal');
    freq = 2.^(round(log2(2)):round(log2(300)));
    AX.YTickLabelMode = 'auto'; AX.YTick = freq;
    ylabel('Frequency');xlabel('Trial Time (s)');
    colormap(gca,brewermap(256, '*RdBu'));
    colorbar
    %......... Save color range .....
    EQLZR1(drug) = eqlzr1;
    EQLZR2(drug) = eqlzr2;
end

% save2pdf(['iEEG signal comparison'],gcf,800)
% for drug = 1:3
%     subplot(3,3,3*drug)
%     caxis([min(EQLZR2),max(EQLZR1)])
% end
