fs        = 512;

sub = 63:81;
sub([5,15]) = [];

freq = cell(length(sub),1);
cfg = [];
cfg.method = 'wavelet';
cfg.output = 'pow';
cfg.foi = 2:.5:30;
cfg.pad = 'nextpow2';
% cfg.width = linspace(2.6,10,length(cfg.foi));
cfg.width = 3;
cfg.keeptrials = 'no';
cfg.baseline = 'no';


load('H:\Amiens 2021\Edalati\New data\Chan128.mat');
chanlocs = cell(1,128);
for i=1:length(chanlocs)
    chanlocs{i} = chan(i).labels;
end

%% Calculating the TFR of silence period of each subject and averaging over time to find baseline
Mean = zeros(length(chan),length(cfg.foi),length(sub));
for s=1:length(sub)
    address = ['\N',num2str(sub(s)),'\PreFFTOfSilencePeriod.mat'];
    load(address);
    
    sig = EEG.data;
    
    temp = floor(size(sig,2)/38/fs);    
    sig = reshape(sig,size(sig,1),[],temp);
    sig = sig(:,1:floor(38*fs)-1,:);    
    
    tmp = cell(length(EEG.chanlocs),1);
    for i=1:length(EEG.chanlocs)
        tmp{i} = EEG.chanlocs(i).labels;
    end
    temp = ismember(chanlocs,tmp);
    tmp1 = find(temp==0);
    tmp2 = find(temp==1);
    NewSig(tmp2,:,:) = sig;
    NewSig(tmp1,:,:) = NaN;
    sig = NewSig;    
    clear NewSig

    sig = sig(:,:,randperm(size(sig,3),10));
    
    sig = pop_importdata('data',sig,'srate',fs);
    sig.chanlocs = chan;
    sig = eeglab2fieldtrip(sig,'preprocessing');   
    for i=1:length(sig.time)
        sig.time{i} = linspace(0,38,size(sig.trial{1},2));
    end
    cfg.toi      = sig.time{1}(sig.time{1}>=0 & sig.time{1}<=38);
    pow          = ft_freqanalysis(cfg,sig);

    pow.powspctrm = pow.powspctrm(:,:,2*fs+1:end-2*fs);
    Mean(:,:,s) = mean(pow.powspctrm,3,'omitnan');
end
%% TFR of the main data
SigS001 = cell(length(sub),1);
Baseline = cell(length(sub),1);
for s=1:length(sub)   
    address = ['\N',num2str(sub(s)),'\PreFFTdata.mat'];
    load(address);
    
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

for k=1:length(sub)  
    tic
    sig = SigS001{k};   
    
    sig(:,:,sum(isnan(squeeze(sum(sig,2))))==128) = [];
    
    sig = pop_importdata('data',sig,'srate',fs);
    sig.chanlocs = chan;
    sig = eeglab2fieldtrip(sig,'preprocessing');   
    for i=1:length(sig.time)
        sig.time{i} = linspace(0,38,size(sig.trial{i},2));
    end
    
    cfg.toi      = sig.time{1}(sig.time{1}>=0 & sig.time{1}<=38);

    power = zeros(length(sig.label),length(cfg.foi),length(cfg.toi)-2*round(2*fs)-1);
    for i=1:length(sig.trial)
        temp       = sig;
        temp.trial = temp.trial(i);
        temp.time  = temp.time(i);
        pow        = ft_freqanalysis(cfg,temp);
        temp       = pow.powspctrm;
        temp       = temp(:,:,round(2*fs)+1:end-round(2*fs)-1);
        
        power      = power + temp;
    end
    power = power./length(sig.trial);
        
    freq{k,1} = power;
toc
end
%% Normalizing the achieved TFR
Freq = freq;

for s=1:length(sub)
    Freq{s} = (Freq{s} -  repmat(Mean(:,:,s),1,1,size(Freq{s},3))) ./ ...
        repmat(Mean(:,:,s),1,1,size(Freq{s},3));%./Std(:,:,s);   
end
%%
cfg             = [];
cfg.showlabels  = 'yes';
cfg.showoutline = 'yes';
cfg.masknans    = 'no';

temp = zeros(size(Freq{1}));
for s=1:length(sub)
    temp(:,:,:,s) = Freq{s};
end
temp = mean(temp,4,'omitnan');
temp = reshape(temp,size(temp,1),size(temp,2),1023,[]);
temp = mean(temp,4,'omitnan');
pow.powspctrm = temp;

temp = pow.elec.chanpos(:,1);
pow.elec.chanpos(:,1) = -pow.elec.chanpos(:,2);
pow.elec.chanpos(:,2) = temp;    
figure
pow.time = linspace(0,2,1023);
pow.freq = 2:.5:30;
ft_multiplotTFR(cfg,pow);