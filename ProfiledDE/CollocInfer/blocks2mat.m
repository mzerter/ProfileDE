function Hmat = blocks2mat(H)  
% Cell array of matrices -> large matrix

[l,k] = size(H);
[m,n] = size(H{1,1});
Hmat = zeros(l*m,k*n);
for i=1:l
    for j=1:k
        indi = ((i-1)*m+1):(i*m);
        indj = ((j-1)*n+1):(j*n);
        Hmat(indi,indj) = H{i,j};
    end
end

end

