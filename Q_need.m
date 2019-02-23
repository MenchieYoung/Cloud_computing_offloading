function [Q1_obj]= Q_need(W,D) 
    global alpha P f SNR beta gama0 rho ;     
    Q0=(D.*(alpha+beta.*P))./(rho*(W.*log2(1+SNR)))+(alpha*gama0.*D)./f;    
%        global f SNR gama0 rho ; 
%        Q0=D./(rho*(W.*log2(1+SNR)))+(gama0*D)./f; 
    Q1_obj=sum(Q0);