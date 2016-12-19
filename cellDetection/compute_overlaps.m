function LtM = compute_overlaps(ipix1, ilam1, xy1, ipix2, ilam2, xy2, d)

n1 = size(ipix1,2);
n2 = size(ipix2,2);

LtM = zeros(n1, n2);

d2 = zeros(n1,n2);
for j = 1:2
    d2 = d2 + bsxfun(@minus, xy1(:,j), xy2(:,j)').^2;
end
hasOverlap = d2 < d^2;


for j = 1:n2
    ineigh = find(hasOverlap(:,j));
    for k = 1:length(ineigh)
        tile = zeros(d/2, d/2);
        
        ipix = ipix1(:,ineigh(k));
        
        % relative position?
        ipix = ipix + 
        
    end
end