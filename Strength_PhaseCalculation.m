fs        = 512;

sub = 63:81;
sub([5,15]) = [];

load('Chan128.mat');
chanlocs = cell(1,128);
for i=1:length(chanlocs)
    chanlocs{i} = chan(i).labels;
end

%% data preparation
SigS001 = cell(length(sub),1);
for s=1:length(sub)   
    address = ['\N',num2str(sub(s)),'\PreFFTData.mat'];
    load(address);
    
    EEG = pop_eegfiltnew(EEG,7,0);
    EEG = pop_eegfiltnew(EEG,0,12);
    
    LS001 = 0;
    for i=1:length(EEG.event)        
        if strcmp(EEG.event(i).type,'S001')
            LS001 = LS001 + 1;
        end
    end
    
    sigS001 = zeros(EEG.nbchan,round(38*EEG.srate),LS001);
    k1 = 0;
    k2 = 0;
    for i=1:length(EEG.event)
        if strcmp(EEG.event(i).type,'S001')
            k2 = k2 + 1;
            sigS001(:,:,k2) = EEG.data(:,EEG.event(i).latency:...
                EEG.event(i).latency + round(38*EEG.srate)-1);
        end
    end
    
    tmp = cell(length(EEG.chanlocs),1);
    for i=1:length(EEG.chanlocs)
        tmp{i} = EEG.chanlocs(i).labels;
    end
    temp = ismember(chanlocs,tmp);
    tmp1 = find(temp==0);
    tmp2 = find(temp==1);
    NewSig(tmp2,:,:) = sigS001;
    NewSig(tmp1,:,:) = NaN;
    sigS001 = NewSig;    
    clear NewSig
    
    SigS001{s} = sigS001; 
end
%% TFR calculation
freq = cell(length(sub),1);
cfg = [];
cfg.method = 'wavelet';
cfg.output = 'pow';
cfg.foi = 7:.25:12;
cfg.pad = 'nextpow2';
cfg.width = 3;
cfg.keeptrials = 'no';
cfg.baseline = 'no';

[b, a, fCenterList, nTap, fEdges] = filterBank_prop(3, 3, 2.65, 2, 512, 'fir1', .25);
b = b{1};
a = a{1};
for k=1:length(sub)     
    for j=1:size(SigS001{k},3)        
        sig = SigS001{k}(:,:,j);
    
        sig = pop_importdata('data',sig,'srate',fs);
        sig.chanlocs = chan;
        sig = eeglab2fieldtrip(sig,'preprocessing');   
        for i=1:length(sig.time)
            sig.time{i} = linspace(0,38,size(sig.trial{i},2));
        end
        cfg.toi      = sig.time{1}(sig.time{1}>=0 & sig.time{1}<=38);
        pow          = ft_freqanalysis(cfg,sig);
        pow          = pow.powspctrm;
        
        pow = squeeze(mean(pow,2,'omitnan'));
        ind = isnan(mean(pow,2,'omitnan'));
        pow(ind,:) = [];
        pow(:,isnan(mean(pow,1))) = 0;
        pow = pow';        
        pow = filtfilt(b,a,pow);
        pow = pow';
        Pow(find(~ind),:) = pow;
        Pow(find(ind),:) = NaN;
        pow = Pow;
        clear Pow
        pow = pow(:,round(2*fs)+1:end-round(2*fs));
            
        freq{k}(:,:,j) = pow;
    end
end
%% phase difference calculation
t = linspace(0,38,length(sig.data));
stimulus = .6*sin(2*pi*3*t+pi/2);
stimulus = angle(hilbert(stimulus));
stimulus = stimulus(:,round(2*fs)+1:end-round(2*fs));

ang = cell(length(sub),1);
for k=1:length(sub)
    for j=1:size(freq{k},3)
        temp = hilbert(freq{k}(:,:,j)');
        temp = angle(temp)';
        
        for jj=1:size(temp,1)
            temp(jj,:) = circ_dist(stimulus,temp(jj,:));
        end
        temp = mean(exp(1i*temp),2);
        
        ang{k}(:,j) = angle(temp);
    end
end

%% strength and phase plot
sub = 1:17;
clear tmp
for i=1:17
    tmp(:,i) = mean(exp(1i*ang{i}),2,'omitnan');
end
aa = abs(tmp);
bb = angle(tmp);
aa = mean(aa,2,'omitnan');
figure
topoplot(aa,chan,'electrodes','on','style','map')
colormap jet

a = angle(mean(exp(1i*bb),2,'omitnan'));
figure
topoplot(a,chan,'electrodes','on','style','map')
clim([-pi,pi])
phasemap;
phasebar;

%% phase and strength plot
ch = [10,11,12,16,18,19,20,22,23,24]; % cluster

sub = 1:17;
clear tmp
for i=1:17
    tmp(:,i) = mean(exp(1i*ang{i}),2,'omitnan');
end

cc = mean(tmp(ch,:),1,'omitnan');
cc1 = abs(cc);
cc2 = angle(cc);

figure
for j=1:length(sub)
    circ_plotWithR(cc2(j),cc1(j),'pretty',[],colorMap,[min(cc1),max(cc1)], 20,true,true,'linewidth',6,'iter',j)
    hold on
end

figure
colorMap = [0,0,0;.25,.25,.25;.45,.45,.45;.65,.65,.65;.85,.85,.85];
colormap(colorMap)
t = linspace(-.333/2,.333/2,340);
z = sin(2*pi*3*t+pi/2);
t = linspace(-333/2,333/2,340);
plot(t,z,'LineWidth',6)
set(gca,'FontSize',40,'FontWeight','Bold')
hold on
tpi = linspace(-pi,pi,340);
for i=1:length(cc2)
    [~,x] = min(abs(tpi-cc2(i)));
    val = cc1(i);
    scatter(t(x),z(x),3000,[.5,.5,.5],'filled');
end
clim([min(cc1), max(cc1)])

%% Surrogate
t = linspace(0,38,19438);
stimulus = .6*sin(2*pi*3*t+pi/2);
stimulus = angle(hilbert(stimulus));
stimulus = stimulus(:,round(2*fs)+1:end-round(2*fs));

Stimulus = stimulus;
q = cell(length(sub),1);
for k=1:length(sub)
    temp = hilbert(permute(freq{k},[2,1,3]));
    q{k} = angle(temp);    
end
ang = cell(length(sub),1000);
for i=1:1000
    tic
    temp = Stimulus;
    temp = temp(:,randperm(length(temp),length(temp)));
    
    stimulus = temp;
    for k=1:length(sub)
        temp = repmat(stimulus',1,128,size(q{k},3));
        
        temp = exp(1i*q{k})./exp(1i*temp);
        temp = squeeze(mean(temp,1));

        ang{k,i} = angle(temp);        
    end
    toc
end