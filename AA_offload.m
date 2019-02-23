function Q_AA=AA_offload(m)
global W_sum D_sum;
W=W_sum/m*ones(1,m);
D=D_sum/m*ones(1,m);
Q_AA=Q_need(W,D);