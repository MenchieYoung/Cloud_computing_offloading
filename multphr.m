%function [x,mu,lam,ouput]=multphr(fun,hf,gf,dfun,dhf,dgf,x0) 
%fun dfun分别是目标函数及其梯度
%hf dhf分别是等式约束函数及其Jacobi矩阵的转置
%gf dgf分别是不等式约束及其Jacobi矩阵的转置
%x0是初始点
%x是近似最优解
%mu lam分别是相应于等式约束和不等式约束的乘子向量（lamda*）
%output是结构变量，输出近似极小值f，迭代次数，内迭代次数等

% function multphr(m)

global alpha P f SNR beta  gama0  W_sum D_sum rho ant;
    m=8;    %从节点数 目前参数设置下应大于3,否则不收敛；应小于10，否则目前迭代次数不够
    alpha=0.5;
    a1=50;     %mW
    b1=100;    %mW
    a3=4;  %8000kB/s  % 8000khz  %8Mhz
    b3=16;  %Mhz
    a4=20;
    b4=35;
    gama0=20;
    W_sum=80;       %总子载波数 %先考虑等式，学姐论文开始也一直认为是=   %拉格朗日（回去找之前论文，没找到）  罚函数(万方论文)      80    %总带宽吧
    D_sum=40;      
    %%%%%%%%% m=2应该要考虑一下，8*m
    %%%%%%%%%%%%%%%%%%%%这一块参数设置emmm%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Q0=zeros(1,m);  %初始化
    beta=1-alpha ;
    rho=2^(-1/2);
    P = a1 + (b1-a1).*rand([1 m]);  %在区间（a，b）内产生均匀分布的随机数   %m行 n列
% %    g = a2 + (b2-a2).*rand([1 m]);
    f = a3 + (b3-a3).*rand([1 m]);
    SNR = a4 + (b4-a4).*rand([1 m]);

maxk=80;       %最大迭代次数 %外迭代 %外迭代太多反而会导致目标函数值变大？20就不好  % 80
%Mk=2.0;          %罚因子     %原先较sigma
Mk=2.0;     
%epsilon=1e-5;    %终止误差值
epsilon=0.0001;
epsilon2=0.1;
epsilon3=1;
eta=10;          %放大罚因子值，一般取10
theta=0.8;       %PHR算法中的实参数
k=0;             %外迭代次数
%ink=0;           %内迭代次数

 ant=200;   % 蚂蚁数量
 %he=zeros(2,ant);
 he=zeros(1,2);
 %gi=zeros(1,ant);
% X0{1,:}=zeros(1,m); 
% X0{2,:}=zeros(1,m);
 X0=cell(2,ant);
%%%%%%%%%%%%%%%%%%%给定初始1000个点%%%%%%%%其实只设定一个也行，因为最后解出的值也都是最值处的x，而不是1000只蚂蚁的，设一个的话，感觉应该先运行一下ACO得到最优值。。所以没必要给he初值%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%不行！！！！还是要以1000个点作为初始值以及变量的！！！！！！但是把1000个作为一个整体%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 for i=1:ant 
   X0{1,i}=(0+(W_sum-0)*rand(1,m)); %随机设置每只蚂蚁的初值位置  %横 W  %就该用{}而不是（），（）拿出的是cell，而不是给这个cell赋值
   X0{2,i}=(0+(D_sum-0)*rand(1,m)); %纵  D
  % ta(i)= Q(X0{1,i},X0{2,i},Mk,lamdak);     %是Q没关系，要的就是index算he
 end  
   %  [ta_best,best_index]=max(ta);
   %  he_X={X0{1,best_index};X0{2,best_index}};             %W   D   
%    % % he_X={X{1,i};X{2,i}};
%    % % he(:,ant)=feval(hf,he_X);   %等式约束  %试了一下 X{：,i}取不出来一列？？？   %实验了可以横向量可以给一列赋值  he第一行：W的差值  he第二行：D的差值
   %  he=feval(hf,he_X); 
%     % %gi(i)=feval(gf,X{2,i});    %不等式约束 
% he_X=cell(2,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%妈的绕了好久感觉白绕了%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%n=length(x);
% l=length(he);       %现在是2
% m=length(gi);    
% mu=0.1*ones(l,1); 
% lam=0.1*ones(m,1);
%lamdak=0.1*ones(2,ant); 
 lamdak=0.1*ones(1,2);
% lam=0.1*ones(ant,1);
betak=10;     
betaold=10;     %用来检验终止条件的两个值
h=zeros(1,maxk);
min_help=zeros(1,maxk);
Qobj_min=zeros(1,maxk);
maxxyt=cell(2,maxk); %上行x（W），下行y（D）
times=100;
while betak>epsilon && k<maxk 
   %调用BFGS程序求解无约束子问题  [ik,x,val]=bfgs('mpsi','dmpsi',x0,fun,hf,gfun,dfun,dhf,dgf,mu,lam,Mk);   
    % ink=ink+ik;  
    [lastW_distri,lastD_distri,maxx,maxy,minvalue]=antcolony_offload(m,X0,Mk,lamdak,k,times) ;    %改成返回值最终一次分布的1000个点也觉得不妥，因为有的点在局部最优，需要的还是maxx，maxy  
     
%      %he_X={maxx;maxy};
%      %he=feval(hf,he_X);       %y
       he=hf(maxx,maxy);
     %he=hf(lastW_distri,lastD_distri);
     %he=cell2mat(he);
%     % gi=feval(gf,x);
     betak=sqrt(norm(he,2)^2);            %+norm(min(gi,lam/Mk),2)^2);  %norm(he,2)返回矩阵的‘2’范数  %betak相当于h（x^k）,betaold相当于h（x^(k-1）)
    
     h(k+1)=betak;   %下标不能从0开始 %h(1)对应k=0
     min_help(k+1)=minvalue;   %辅助函数
     Qobj_min(k+1)=Q_need(maxx,maxy); %真正目标函数
     maxxyt{1,k+1}=maxx;
     maxxyt{2,k+1}=maxy;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
     if betak>epsilon              %没达到hf趋于0，更新乘子向量
         lamdak=lamdak-Mk*he;      %mu改成lamdak了  %主要步骤
     %    lam=max(0.0,lam-Mk*gi);   %不等式的
            if k>=2 && betak>theta*betaold     %判断收敛快慢
                Mk=eta*Mk;
            end
     end 
     if betak<epsilon3    %在可行域里适当多迭代一些，找最值
         times=200;    
     end
%    h(k)=betak;
     k=k+1;    %我理解的迭代次数是对Mk和lamdak而言的
     %k=0，对Mk等迭代了0次，但是对h min obj_min等已经发生了第一次
     %也就是说%%%%%是分界线
     %最后一次迭代的Mk和Lmadak并没有用到
     betaold=betak;
     X0=[lastW_distri;lastD_distri];  %擦擦不是大括号，里面已经是cell了  %不行还是要1000个点，不然没法迭代，但是可以再输出1000个点专门为迭代用，maxx，maxy为he所用
                                  %竟然还有中括号...
     
end
 k_plot=(1:1:k);
figure,plot(k_plot-1,h(k_plot)),title('h(k)每次迭代变化历程');grid on;
figure,plot(k_plot-1,Qobj_min(k_plot)),title('真正目标函数每次迭代最小值变化历程');grid on;hold on;
Q_AA=AA_offload(m);
Q_RA=RA_offload(m);
plot([0,k],[Q_AA,Q_AA],'--r');
figure,plot(k_plot-1,min_help(k_plot)),title('辅助函数每次外迭代最小值变化历程');grid on;
%  k_plot2=(1:1:k-10);
% figure,plot(k_plot2-1,min_help(k_plot2)),title('辅助函数每次外迭代最小值变化历程');grid on;

j_index=0;
j=zeros(1,k);
for jj=1:1:k
  if h(jj)<epsilon2
      j(j_index+1)=jj;
      j_index=j_index+1;
  end
end
temp=Qobj_min(j(1:j_index));
[Qobj_minvalue,Q_obj_index]=min(temp);  %b  = [a(1:10),a(20:25),a(51:60)];
  W_need=maxxyt{1,j(Q_obj_index)};
  D_need=maxxyt{2,j(Q_obj_index)};
% % f=feval(fun,x);
% % output.fval=f;
% output.iter=k;    %外迭代次数
% % %output.inner_iter=ink;
% output.beta=betak;
% output.Q_minvalue=minvalue;
% %output.Qobj_minvalue=Q_need(maxx,maxy);
output.ACO_ALMM=Qobj_minvalue;
output.AA=Q_AA;
output.RA=Q_RA;

