function goodcell = remove_doubles2(F, med)

goodF = 1 - mean(diff(F,1).^2, 1)./(2*var(F, [], 1));

CC = corrcoef(my_conv2(F, 2, 1));
%
[NT NN] = size(F);
Ls = zeros(NN, NN, 3);
for i = 1:size(med,2)
    xs = med(:,i) * ones(1,NN);
    Ls(:,:,i) = (xs - xs').^2;
end
ll = Ls(:,:,3);
ll(ll>64.^2) = Inf;
ll(ll<64.^2) = 0;
Ls(:,:,3) = ll;
Ls = sum(Ls,3).^.5;
Ls = Ls + diag(Inf * ones(NN,1));
%
goodcell = true(NN,1);
ix = Ls<5;
for i = 1:NN
    ifs = find(ix(i,:) & CC(i, :)>.6, 1);
    if ~isempty(ifs)
       if goodF(i)>goodF(ifs)
          goodcell(ifs) = 0; 
       else
           goodcell(i) = 0;
       end
    end
end
%

% F = F(:, goodcell);
% med = med(goodcell, :);
% %
% [NT NN] = size(F);
% 
% clear Ls
% for i = 1:size(med,2)
%     xs = med(:,i) * ones(1,NN);
%     Ls(:,:,i) = (xs - xs').^2;
% end
% Ls = sum(Ls,3).^.5;
% Ls = Ls + diag(Inf * ones(NN,1));
% 
