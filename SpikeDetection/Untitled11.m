% figure; 
clf
iNN = ceil(rand * size(Winit,2));
plot(zscore(W(:,iNN)))
hold all
plot(zscore(V(:,1)))
