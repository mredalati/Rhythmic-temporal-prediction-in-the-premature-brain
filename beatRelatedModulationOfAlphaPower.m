%%
Freq = freq;

for s=1:length(sub)
    Freq{s} = (Freq{s} -  repmat(Mean(:,:,s),1,1,size(Freq{s},3))) ./ ...
        repmat(Mean(:,:,s),1,1,size(Freq{s},3));%./Std(:,:,s);   
end

for k=1:17
    Freq{k} = Freq{k}(:,:,round(2*fs)+1:end-1*round(2*fs));
end

f = 2:.5:30;
ff = f>=7 & f<=12;
for k=1:17
    Freq{k} = squeeze(mean(Freq{k}(:,ff,:),2,'omitnan'));
    temp = Freq{k};
    Freq{k} = temp - repmat(mean(temp,2),1,size(temp,2));
end
%% Time course of alpha power
clear sig
t1 = round(.333*fs);
t2 = round(.333*fs);
for s=1:length(sub)
    ind = [1:6:16*6+1,5:6:16*6+5];
%     ind = [6:6:16*6,2:6:16*6+2];
%     ind = 1:16*6+5;
    
    temp = zeros(size(Freq{s},1),t1+t2,size(Freq{s},3)*length(ind));
    for i=1:length(ind)        
        temp(:,:,(i-1)*size(Freq{s},3)+1:i*size(Freq{s},3)) = ...
            Freq{s}(:,round(ind(i)*.333*fs)-t1+1:...
            round(ind(i)*.333*fs)+t2,:);
    end
    sig(:,:,s) = mean(temp,3,'omitnan');
    sig(:,:,s) = sig(:,:,s) - repmat(mean(sig(:,:,s),2),1,size(sig,2),1);    
end
t = linspace(-(t1/fs)*1000,(t2/fs)*1000,size(sig,2));
ch = [10,11,12,16,18,19,20,22,23,24];

waveforms = squeeze(mean(sig(ch,:,:),'omitnan'))';
ci = bootci(1000,@mean,waveforms);
CI = ci;
figure
plot(t, mean(waveforms,'omitnan'), 'k', 'LineWidth', 6)
hold on
fill([t, fliplr(t)], [CI(1,:), fliplr(CI(2,:))], [.5,.5,.5], 'FaceAlpha', 0.3, 'EdgeAlpha', 0)
set(gca,'FontSize',60,'FontWeight','Bold')

%% Statistical analysis
q = squeeze(mean(sig(ch,:,:),'omitnan'));
q = q';
tt = nan(length(sub),1);
tt2 = nan(length(sub),1);
tprim = t(t<=333/2 & t>=-333/2);


for i=1:length(sub)
    if ~isempty(findpeaks(q(i,t<=333/2 & t>=-333/2)))
        [a,b] = findpeaks(q(i,t<=333/2 & t>=-333/2));        
        [~,idx] = max(a);
        tt(i) = tprim(b(idx));
    end
end

for i=1:length(sub)
    tprim2 = t(t<=333 & t>=tt(i));
    if ~isempty(findpeaks(-q(i,t<=333 & t>=tt(i))))
        [a,b] = findpeaks(-q(i,t<=333 & t>=tt(i)));        
        tt2(i) = tprim2(b(1));
    end
end
tt = tt2;

figure
ax = boxplot(tt,'Orientation', 'horizontal','color','k');
bx = findobj('Tag','boxplot');
bx.Children(2).LineStyle = ':';
set(gca,'FontSize',60,'FontWeight','Bold');
set(ax,'LineWidth', 6);
axis([-400,400,.5,1.5])
h = lillietest(tt);
if h==1
    [p,h,stat] = signrank(tt,0,'Tail','Right');
else
    [h,p,~,stat] = ttest(tt,0,'Tail','Right');
end