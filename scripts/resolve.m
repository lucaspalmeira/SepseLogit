function [alpha, x] = resolve(A,b)
% RESOLVE
% obtem alpha = (I + A'A)\(A'b) sem computar a matriz A'A
% marcos@dcc.ufmg.br

    [m,n] = size(A);
    Im = speye(m); % alterado
    In = speye(n); % alterado
    M = sparse(m+n,m+n);
    M = [Im,-A; A',In];
    nb = zeros(m+n,1);
    nb(1:m) = -b;
    x = M\nb;
    alpha = x(m+1:end);
end
