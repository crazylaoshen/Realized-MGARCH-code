A = A_true;
B = B_true;
x = tril(A,-1);
y = triu(A,1);

z =setdiff(A, diag(A));
mean(z)
mean(diag(A))
mean(diag(B))


L = EliminationM(10);
D = DuplicationM(10);
Av_true =  L*kron(A_true,A_true)*D;
mean(diag(Av_true))
Bv_true =  L*kron(B_true,B_true)*D;
mean(diag(Bv_true))