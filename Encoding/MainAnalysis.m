%%  iEEG analysis
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

CLR{3} = [0.8,0.5,0.5];
CLR{2} = [0.5,0.5,0.5];
CLR{1} = [0.5,0.8,0.8];

drug_strn{1} = 'Pregabalin';
drug_strn{2} = 'Gabapentin';
drug_strn{3} = 'Arecoline';
%% demo
trl = 901;
figure
for chnl = 1:4
    x = squeeze(lfp.pre(trl,chnl,:));
    clr = [0, 1-1/(chnl), 1-1/(chnl)];
    plot(chnl+x,'Color',clr);
    hold on
end
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

%% Plot power spectrum
figure

for drug = 1:3
    avg_pwr = abs(nanmean(power{drug},3))';
    [~,start_freq] = min(abs(frq-300)); %find the the element of 50Hz
    A=avg_pwr(:,start_freq:end)';
    
    
    plot(log2(frq(start_freq:end)), (mean(A')),'-.','Color',CLR{drug},'LineWidth',3);
    % area(log(frq(start_freq:end)), mean(A'),...
    %     'EdgeColor', 'none', 'FaceColor', [0.5 0.5 0.5]);
    AX = gca;
    % set(gca,'XScale','log'); set(AX,'Xdir','normal');
    % AX.XTickLabelMode = 'auto'; AX.XTickLabel = freq;
    AX.XTickLabel = freq;
    xlim(log2([2,256]))
    ylabel('Power (dB)');
    hold on;
    xlabel('Frequency')
end
legend(drug_strn{1},drug_strn{2},drug_strn{3})
save2pdf(['spectrum comparison'],gcf,800)

% plot([200,200],ylim,'--','Color', [1 0 0 0.5], 'LineWidth', 1.5)
% xlim([1,1400])

%% plot TSNE

figure
for drug = [3,2,1]
    avg_pwr = squeeze(abs(nanmean(power{drug},2)))';
    [~,start_freq] = min(abs(frq-300)); %find the the element of 50Hz
    A=avg_pwr(:,start_freq:end);
    
    
    [~,lw_f1]  = min(abs(frq-10));
    [~,frq50]  = min(abs(frq-50));
    low_pwr = mean(A(:,lw_f1:end),2);
    high_pwr = mean(A(:,frq50:lw_f1),2);
    
    x = (log(1000*low_pwr));
    y = (log(1000*high_pwr));
    
    p{drug} = plot(x,y,'.','MarkerSize',20,...
        'color',[CLR{drug}],...
        'LineWidth',2)
    hold on
end

X = [];Y = [];
for drug = [3,2,1]
    avg_pwr = squeeze(abs(nanmean(power{drug},2)))';
    A=avg_pwr;
    [~,lw_f1]  = min(abs(frq-10));
    [~,frq50]  = min(abs(frq-50));
    low_pwr = mean(A(:,lw_f1:end),2);
    high_pwr = mean(A(:,frq50:lw_f1),2);
    
    x = (log(1000*low_pwr));
    y = (log(1000*high_pwr));
    
    X = [X,x];%save for later analysis
    Y = [Y,y];
    hold on
    plot(mean(x),mean(y),'+','Color',0.8*CLR{drug},'MarkerSize',150*sqrt(std(x).^2+std(y).^2),'LineWidth',1)
    plot(mean(x),mean(y),'O','Color',0.8*CLR{drug},'MarkerSize',150*sqrt(std(x).^2+std(y).^2),'LineWidth',1)
    grid on
    
    
    [x_d,y_d] = ginput(1)
    plot(x_d,y_d,'*k')
    [~,select_trl{drug}]= min(abs(x_d-log(1000*low_pwr)))
end
legend([p{1},p{2},p{3}],{drug_strn{1},drug_strn{2},drug_strn{3}},'Location','northwest')
xlabel('Low frequncy power')
ylabel('High frequncy power')


% save2pdf(['Drug similarity space'],gcf,800)

%% I selected examples here :D
% save(['selected_trials.mat'],'select_trl')
load('selected_trials.mat')

%%  select individual trial examples
h= figure(3);
scrsz = get(0,'ScreenSize');
scrsz(4) = scrsz(4)/1.5;
scrsz(3) = scrsz(3)/2;
set(h, 'Position',scrsz);
Ax = gcf;
Ax.Color = [1 1 1];

CLR{3} = [0.8,0.5,0.5];
CLR{2} = [0.5,0.5,0.5];
CLR{1} = [0.5,0.8,0.8];

drug_strn{1} = 'Pregabalin';
drug_strn{2} = 'Gabapentin';
drug_strn{3} = 'Arecoline';

chnl = 2;
t1 = 1;
t2 = find(Time==2);

subplot(3,3,[1,2]) %Pregabalin
x = squeeze(lfp.pre(select_trl{1},chnl,t1:t2));
clr = CLR{1};

plot(Time(t1:t2),x,'Color',clr,'LineWidth',2);
title(drug_strn{1})
box off
ylim([-2.5 2.5])

subplot(3,3,[4,5]) %Phenazepam
x = squeeze(lfp.gab(select_trl{2},chnl,t1:t2));
clr = CLR{2} ;
plot(Time(t1:t2),x,'Color',clr,'LineWidth',2);
title(['Similar Drug (',drug_strn{2},')'])
box off
ylim([-2.5 2.5])

subplot(3,3,[7,8]) %Arecoline
x = squeeze(lfp.are(select_trl{3},chnl,t1:t2));
clr = CLR{3};
plot(Time(t1:t2),x,'Color',clr,'LineWidth',2);
title('Other Drug (Arecoline)')
box off
ylim([-2.5 2.5])

pwr= [];
for drug = 1:3
    
    switch drug
        case 1
            signal = lfp.pre;
        case 2
            signal = lfp.gab;
        case 3
            signal = lfp.are;
    end
    
    trl = select_trl{drug};
    %........................Scalogeram (Wavelet transform power)..........
    y = squeeze(signal(trl,chnl,:));
    [cfs,frq] = cwt(y,Fs);
    pwr{drug}(:,:) = (cfs);
end


EQLZR1 = [];
EQLZR2 = [];

for drug = 1:3
    subplot(3,3,3*drug)
    
    avg_pwr = abs(pwr{drug})';
    [~,start_freq] = min(abs(frq-300)); %find the the element of 50Hz
    A=avg_pwr(t1:t2,start_freq:end)';
    % K = 16^-1*ones(10,2); A = conv2(A,K,'same');
    pcolor(Time(t1:t2),frq(start_freq:end),A);
    
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
    caxis([min(EQLZR2),max(EQLZR1)])
end
save2pdf(['Single trial example'],gcf,800)

%% Coeffcient of similarity
h= figure(3);
scrsz = get(0,'ScreenSize');
scrsz(4) = scrsz(4)/1.5;
scrsz(3) = scrsz(3)/2;
set(h, 'Position',scrsz);
Ax = gcf;
Ax.Color = [1 1 1];

A = [X(:),Y(:)]; 
clust = kmeans(A,3);
[s,h] = silhouette(A,clust,'Euclidean')


% D12 = pdist([1,2,3; 5,6,7]')

l3 = [X(:,1),Y(:,1)];
l2 = [X(:,2),Y(:,2)];
l1 = [X(:,3),Y(:,3)];

d1 = mean(pdist([l1]));
d2 = mean(pdist([l2]));
d3 = mean(pdist([l3]));

D12 =mean(pdist([l1-l2]));
D13 =mean(pdist([l1-l3]));
D23 =mean(pdist([l2-l3]));


sc1(1) = D12;
sc1(2) = D13;

sc2(1) = D12;
sc2(2) = D23;

sc3(1) = D13;
sc3(2) = D23;


c = categorical({drug_strn{1},drug_strn{2},drug_strn{3}});
b = bar(c,[D12;D13;D23],'FaceColor','flat');
b.CData(1,:) = CLR{1};
b.CData(2,:) = CLR{2};
b.CData(3,:) = CLR{3};
