

tic;
clear;
Iex=0;%外部电流，单位 pA
M_area=100;%定义膜面积，单位是um^2
Cm=0.01*M_area;%单位是电容 pF
gL=0.003*M_area;%单位是电导率 nS
gK=0.02;%单位是电导率 nS
VL=-54.4;%单位是 mV
V_K=-77;%单位是 mV


Nch=18*M_area;

delta_t=10^-6;   % dimension  ms
T_step=2^26;
% 
ge=zeros(1,Nch+1);
for i_ch=1:Nch+1
    ge(i_ch)=gL+(i_ch-1)*gK;
end
Ve=zeros(1,Nch+1);
for i_ch=1:Nch+1
    Ve(i_ch)=(gL*VL+(i_ch-1)*gK*V_K+Iex)/ge(i_ch);
end
ge=ge/Cm;


nfft=T_step;
%df_V=zeros(1,T_step);

length_PS=5000;
AV_df_V=zeros(1,length_PS);

Num_aver=10000;
parpool(10)

parfor i_aver=1:Num_aver

    V=zeros(1,T_step);
    S_k=zeros(4,T_step);

%
    V(1)=-66;
    S_k(:,1)=[round(0.4*Nch);round(0.2*Nch);round(0.05*Nch);0];
    n_i=1;


while n_i<T_step-1
                  % 开始，计算dwell time，停留时间。
                  U1=rand;
                  F=1;
                  integ_lambda=0;        %the integrel of lambda(s)
                  while F>U1   &&   n_i<T_step-1
                       V(n_i+1)=V(n_i)-ge(S_k(4,n_i)+1)*(V(n_i)-Ve(S_k(4,n_i)+1))*delta_t;
                       S_k(:,n_i+1)=S_k(:,n_i);
                       alfa_n=0.01*(V(n_i)+55)/(1-exp(-(V(n_i)+55)/10));
                       bata_n=0.125*exp(-(V(n_i)+65)/80);
                       lambda1=(Nch-sum(S_k(:,n_i)))*4*alfa_n;
                       lambda2=S_k(1,n_i)*bata_n;
                       lambda3=S_k(1,n_i)*3*alfa_n;
                       lambda4=S_k(2,n_i)*2*bata_n;
                       lambda5=S_k(2,n_i)*2*alfa_n;
                       lambda6=S_k(3,n_i)*3*bata_n;
                       lambda7=S_k(3,n_i)*alfa_n;
                       lambda8=S_k(4,n_i)*4*bata_n;
                       lambda_T=lambda1+lambda2+lambda3+lambda4+lambda5+lambda6+lambda7+lambda8;
                       integ_lambda=integ_lambda+lambda_T*delta_t;
                       F=exp(-integ_lambda);
                       n_i=n_i+1;
                   end
                  
                  %%%%%%%%%%%%%%%%%%%%%%%%%Jump
                  
                  U2=rand;
                  interv1=        lambda1;%/lambda_T;
                  interv2=interv1+lambda2;%/lambda_T;
                  interv3=interv2+lambda3;%/lambda_T;
                  interv4=interv3+lambda4;%/lambda_T;
                  interv5=interv4+lambda5;%/lambda_T;
                  interv6=interv5+lambda6;%/lambda_T;
                  interv7=interv6+lambda7;%/lambda_T;
                  
                  
                  Rates=[0 interv1 interv2 interv3 interv4 interv5 interv6 interv7 lambda_T]/lambda_T;
                  Jumpto=FindLoca2(U2,Rates);
                  
                  
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%Jump to the next channel state
                  switch    Jumpto
                      case   1
                          S_k(:,n_i+1)=S_k(:,n_i)+[1;0;0;0];
                      case   2
                          S_k(:,n_i+1)=S_k(:,n_i)+[-1;0;0;0];
                      case   3
                          S_k(:,n_i+1)=S_k(:,n_i)+[-1;1;0;0];
                      case   4
                          S_k(:,n_i+1)=S_k(:,n_i)+[1;-1;0;0];
                      case   5
                          S_k(:,n_i+1)=S_k(:,n_i)+[0;-1;1;0];
                      case   6
                          S_k(:,n_i+1)=S_k(:,n_i)+[0;1;-1;0];
                      case   7
                          S_k(:,n_i+1)=S_k(:,n_i)+[0;0;-1;1];
                      case   8
                          S_k(:,n_i+1)=S_k(:,n_i)+[0;0;1;-1];
                  end
                  V(n_i+1)=V(n_i);
                  
                  n_i=n_i+1;
                  
end          

    df_V=fft(V,nfft)/nfft;   %note  there are T_step+1 numbers in V
    df_V1=df_V(1:length_PS).*conj(df_V(1:length_PS));
    df_V1(1)=0;
    AV_df_V=AV_df_V+df_V1;

end

delete(gcp('nocreate'))
AV_df_V=AV_df_V/Num_aver;

Freq=0:1:length_PS-1;
Freq=Freq/(nfft*delta_t); 


save PowerSpectVM100deltt6order AV_df_V   Freq T_step nfft Num_aver delta_t M_area  length_PS Cm gL gK V_K VL;


% subplot(2,1,1);
% plot([1:length(V)]*delta_t,S_k(4,:));
% set(gca,'ylim',[-0.5,18])
% subplot(2,1,2);                                                                                 
 %plot([1:length(V)]*delta_t,V);
 %set(gca,'ylim',[-77,-50])