

tic;
clear;
Iex=0;% external current  pA per  squre micrometer，
M_area=150;%定义膜面积，单位是um^2
Cmd=0.01; % 0.01 pF per squre micrometer
%Cm=Cmd*M_area;%单位是电容 pF
gLd=0.003; %0.003 nS per squre micrometer
%gL=gLd*M_area;%单位是电导率 nS
gK=0.02;%单位是电导率 the conductivity of per Potassium channel nS
Dens=18;  % Distribution density of Potassium Channels 18 per squre micrometer
gKm=gK*Dens; %the maxumum conductivity of the Potassium channel per  squre micrometer
N_Ch=Dens*M_area;
ScalingF=1/sqrt(N_Ch);

VL=-54.4;%单位是 mV
V_K=-77;%单位是 mV





gL=gLd/Cmd;
gK=gKm/Cmd;

delta_t=10^-3;   % dimension  ms
T_step=2^20;


x_vec=zeros(4,T_step);
V=zeros(1,T_step);
Kesha_vec=zeros(4,T_step);

% 中间量，conduct，difference of Voltage
    gK_t=0;
    ge=zeros(1,T_step);
    Ve=0;
    Dif_v=0;
    
    x_vec(:,1)=[0.2677; 0.3741; 0.2323;  0.0541];
    V(1)=-66;
    Kesha_vec(:,1)=zeros(4,1);

for i_t=2:T_step
        j_t=i_t-1;
        alfa_n=0.01*(V(j_t)+55)/(1-exp(-(V(j_t)+55)/10));
        bata_n=0.125*exp(-(V(j_t)+65)/80);
        
        x_0=1-sum(x_vec(:,j_t));     %x_vec(1,i_t-1)(i_t-1)-x_vec(2,i_t-1)(i_t-1)-x_vec(3,i_t-1)(i_t-1)-x_vec(4,i_t-1)(i_t-1);
        
        A=[7*alfa_n+bata_n   4*alfa_n-2*bata_n   4*alfa_n           4*alfa_n
           -3*alfa_n         2*alfa_n+2*bata_n   -3*bata_n           0
            0                -2*alfa_n           alfa_n+3*bata_n    -4*bata_n
            0                  0                 -alfa_n            4*bata_n];
        
        
            b13=sqrt(3*alfa_n*x_vec(1,j_t));  b14=sqrt(2*bata_n*x_vec(2,j_t));
            b25=sqrt(2*alfa_n*x_vec(2,j_t));  b26=sqrt(3*bata_n*x_vec(3,j_t));
            b37=sqrt(alfa_n*x_vec(3,j_t));    b38=sqrt(4*bata_n*x_vec(4,j_t));
     Gama=[sqrt(4*alfa_n*x_0)     -sqrt(bata_n*x_vec(1,j_t))    -b13       b14   0      0     0       0
           0                      0                              b13      -b14   -b25   b26   0       0                            
           0                      0                              0         0     b25    -b26  -b37    b38 
           0                      0                              0         0     0      0     b37    -b38];
        
        
        
        gK_t=gK*(x_vec(4,j_t)+ScalingF*Kesha_vec(4,j_t));
        ge(j_t)=gL+gK_t;
        Ve=(gL*VL+gK_t*V_K)/ge(j_t);
        Dif_v=Ve-V(j_t);
        
        
        x_vec(:,i_t)=x_vec(:,j_t)+delta_t*(-A*x_vec(:,j_t)+[4*alfa_n;0;0;0]);
        V(i_t)=V(j_t)+delta_t*ge(j_t)*Dif_v;      %(-gL*(V(i_t-1)-VL)-gK*(x_vec(4,i_t-1)+1/sqrt(N_Ch)*Kesha_vec(4,i_t-1))*(V(i_t-1)-V_K)+Iex);
        Kesha_vec(:,i_t)=Kesha_vec(:,j_t)-A*Kesha_vec(:,j_t)*delta_t+Gama*normrnd(0,1,8,1)*sqrt(delta_t);
        

end
    gK_t=gK*(x_vec(4,T_step)+ScalingF*Kesha_vec(4,T_step));
    ge(T_step)=gL+gK_t;   %ge中，没有包含Cmd
   


% subplot(2,1,1);
% plot([1:length(V)]*delta_t,S_k(4,:));
% set(gca,'ylim',[-0.5,18])
% subplot(2,1,2);                                                                                 
% plot([1:length(V)]*delta_t,V);
% set(gca,'ylim',[-77,-50])