function [cc,mfactor] = regSVDcorr(A,B)

[nX nY numBlocks] = size(A);

eps0 = single(1e-20);
m1 = fft(fft(A,[],1),[],2);
m1 = m1./(abs(m1)+eps0);

m2 = fft(fft(B,[],1),[],2);
m2 = m2./(abs(m2)+eps0);

% make m1 and m2 into 4D matrices for comparison across 10x10
m1 = repmat(m1,1,1,numBlocks);
m1 = reshape(m1,[nX nY numBlocks numBlocks]);
m2 = repmat(m2,1,1,numBlocks);
m2 = reshape(m2,[nX nY numBlocks numBlocks]);
m2 = permute(m2,[1 2 4 3]);

% only compute for lower triangular of matrix
indL = find(tril(ones(numBlocks,numBlocks)));
indS = find(eye(numBlocks,numBlocks));
for i=1:length(indS)
    indSL = find(indL==indS(i));
end
m1 = m1(:,:,indL);
m2 = m2(:,:,indL);
% multiply by this factor to get back to 10x10
mfactor = 2*ones(length(indL),1);
mfactor(indSL) = 1;

% compute correlation matrix
cc = real(ifft(ifft(m1 .* conj(m2),[],1),[],2));
