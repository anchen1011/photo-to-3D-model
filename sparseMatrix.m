function A = sparseMatrix(i, j, Aij, imsize)
% Get a set of 2D linear neiborhood operators and build a sparse matrix
% encoding all the operators.
%
% Input:
%   Aij = [ni, nj, nc] nc = number of neirbourhoods with constraints
%   i = row index
%   j = column index
%   imsize = [nrows ncols]
%
% Output:
%   A = sparse matrix. Each row contains one 2D linear operator.

nrows = imsize(1);

[ni, nj, nc] = size(Aij);
nij = ni*nj;

a = zeros([nc*nij 1]);
m = zeros([nc*nij 1]);
n = zeros([nc*nij 1]);

% Center pixel:
[jj, ii] = meshgrid(-(ni-1)/2:(ni-1)/2);
ii = ii(:); 
jj = jj(:);
%  jj = column index
%  ii = row index


k = 0;
for c = 1:nc    
    % Translate pixels into matrix indices:
    x = (i(c)-1+ii) + (j(c)-1+jj)*nrows + 1;
    
    % Insert data:
    a(k+1:k+nij) = Aij(:,:,c);
    m(k+1:k+nij) = c;
    n(k+1:k+nij) = x;
    
    k = k + nij;
end


% Create sparse matrices
A = sparse(m(1:k),n(1:k),a(1:k));


