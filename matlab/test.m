% matrix test
n = 500;
A = randn(n,n);
b = randn(n,1);
tic,
x = A\b;
toc
tic,
xx = lu_copt(A,b);
toc
tic,
xxx = qr_copt(A,b);
toc