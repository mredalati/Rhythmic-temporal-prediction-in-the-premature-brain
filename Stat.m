load('ChannelsNeighbouring.mat')
load('ForStrength&Phase.mat')

p = nan(128,1);
for i=1:128
    a = mean(sig{1}(i,:),'omitnan');
    b = mean(sig{2}(i,:,:),2,'omitnan');
    b = squeeze(b);    

    s = abs(a-squeeze(mean(b)))./squeeze(std(b,[]));
    p(i) = erfc(s./sqrt(2));

end
p(isnan(mean(sig{1},2,'omitnan'))) = NaN;
cluster = ft_findcluster(p<.05,connmat,3);
cluster(isnan(p)) = [];

%% Rayleigh test
load('EEG_StimPhaseDifference.mat')
for i=1:128
    [p(i),z(i)] = circ_rtest(bb(i,:));
end