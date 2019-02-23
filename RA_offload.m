function Q_RA=RA_offload(m)
global W_sum D_sum;
D=rand(1,m);
D=D_sum*D/sum(D);

W=rand(1,m);
W=W_sum*W/sum(W);

Q_RA=Q_need(W,D);