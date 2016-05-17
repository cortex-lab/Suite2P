

%lam0 = lamF-repmat(mean(lamF')',1,size(lamF,2));
%[ul sl vl] = svd(lam0*lam0');

[~,iu] = sort(ul(:,1),'ascend');

bin = 30;
dd  = lamF(:,1:bin:end);
%dd=squeeze(sum(reshape(lamF(:,1:floor(size(lamF,2)/bin)*bin),size(lamF,1),bin,floor(size(lamF,2)/bin)),2));
dd  = zscore(dd',1,1)';

dp=dd(iu(2:25:end),:);
imagesc([0:1:(size(dp,2)-1)*1],[1:size(dp,1)],dp,[-3 3]);
set(gca,'fontsize',16)
xlabel('time (min)');
ylabel('neuron id');
%title('binned data (1 min bins)');
title('sampled every minute');
