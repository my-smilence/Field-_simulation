clc
clear
close all;
%%
%初始化
field_init(-1);

%%    %%%%%%参数设定%%%%%%%%%
f0                      =   4e6;              % 中心频率 [Hz]
c                       =   1540;             % Speed of sound [m/s]
fs                      =   100e6;            % 采样频率 [Hz]
Ts                      =   1/fs;             % 采样时间间隔 [s]
lambda                  =   c/f0;             % 波长 [m]
height                  =   6e-3;             % 阵元宽度6mm [m]
pitch                   =   lambda/2;         % 阵元中心间距半波长[m]    
kerf                    =   pitch/10;         % 阵元空隙间距 m
width                   =   pitch-kerf;       % 阵元长度[m]    
N_elements              =   96;               % 阵元数

%%  %%%%%%%%%%%%%%%%%%%%%% 初始化阵列设定 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
focus = [0 0 100];

Th = xdc_linear_array (N_elements, width, height, kerf, 1, 1,focus);
Th2 = xdc_linear_array (N_elements, width, height, kerf, 1, 1,focus);
%发射和接收阵列

BW=1*f0;
impulse_response=sin(2*pi*f0*(0:Ts:2/BW));
impulse_response=impulse_response.*hanning(max(size(impulse_response)))';
%脉冲响应函数

    xdc_impulse (Th, impulse_response);
    xdc_impulse (Th2, impulse_response);
    %发射接收阵列相同

excitation=sin(2*pi*f0*(0:Ts:2/BW));
excitation=excitation.*hanning(max(size(excitation)))';
xdc_excitation (Th, excitation);
%激励信号  激励信号和滤波器相同即为匹配滤波，匹配滤波频域再时间上平移t1

LengthConv=(length(excitation) + 2*length(impulse_response)-2)*0.5;
%这是半个脉冲的长度   卷积：M+L-1=75+75-1

%% %%%%%%%%%%%%% 设定显示区域 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lx            = lambda/10;  %波长1/10步进
% imaging_start = 28/1000;
% imaging_depth = 4/1000; % 30/1000; 
% image_width   = 10e-3;   
% Lz             = c*Ts/2; % L是z轴步长，为半波长
% z_start       = imaging_start;
% z_stop        = imaging_start+imaging_depth;  %Z=（28：L：28+4）;
% X             = -image_width/2: Lx: image_width/2;  %中心左右各2.5mm宽度范围
% 
% Z             = z_start: Lz : z_stop;  %仅显示29-31mm深度；此时聚焦点100mm

%% %%%%%%%%%%%%%%%%%%%%%%%%%% 设定散射点 %%%%%%%%%%%%%%%%%%%%%%%%%%
dz=10/1000;        %  散射点之间轴间距[m]
dx=[repmat(-1/1000,1,5),repmat(1/1000,1,5)];        %  x轴方向投影为两个点
z_start=10/1000;  %  开始深度 [m]
Npoints=5;       %  共计10*2个点

positions = [dx; zeros(1,Npoints*2); repmat((1:Npoints)*dz+z_start,1,2)]';
amp=ones(Npoints*2,1);
%[-1 0 20] [1 0 20] ------- [-1 0 60][1 0 60]

%% %%%%%%%%%%%%%%%%% 
td        = LengthConv*2; %单个脉冲的长度/点数
tdd       = (td-1)/2;   %

Lp    = 4;       % Lp 子阵列长度，一半N
P     = N_elements-Lp+1;    % M-L+1= 96/2+1=48+1=49  空间平滑系数
e     = ones(Lp,1);         % 单位矢量，一般为1

yy    = zeros(25000,10000); %初始化地址矩阵

W     = (N_elements*pitch); % -kerf

Tapo  = hamming(N_elements)'; 
xdc_apodization (Th, 0, Tapo); %发射变迹

%% %%%%%%%%%%%%%%%%%%% 聚焦设定延时 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 1:N_elements;
sector=60*pi/180;   %20度扇扫;-10~10
% theta= -sector/2;
i=1;

for theta = -sector/2:pi/180:sector/2
         xdc_center_focus (Th,[0 0 0]);
         xdc_focus (Th, 0, [100*sin(theta)/1000 0 100*cos(theta)]/1000);
         xdc_center_focus (Th2,[0 0 0]);
         xdc_focus (Th2, 0, [100*sin(theta)/1000 0 100*cos(theta)]/1000);
        [v1,times]=calc_scat_multi(Th, Th2, positions, amp);    %(0:N-1)/fs+t  times=第一个接收到的散射点的信号时刻
        
        v1_zeropadded= [zeros(round(times*fs) , N_elements); v1];%补齐起始时间
        
        v(1:length(v1_zeropadded),:)=v1_zeropadded;
       
        %Pass the signal through an AWGN channel having a 60 dB signal-to-noise ratio (SNR)
        %adding white Gaussian noise
        RF_data=awgn( v, 60, 'measured'); %添加高斯白噪声；60=snr指定了每一个采样点信号与噪声的比率，单位为dB
                  
       
%          x_length=(N_elements-1)*pitch;
%          x=-x_length/2-(pitch/2):pitch:x_length/2+(pitch/2);     %阵列位置坐标xy,单位m
%          
%          x_position=50*sin(theta)/1000;       %此时焦点坐标 单位mm
%          y_position=50*cos(theta)/1000;
%          
%          delay(:,i)=1000000000*sqrt((x-x_position).^2+y_position.^2)/1540;    %焦点和每个阵元的距离 ns延时 
%          delay(:,i)=delay(:,i)-delay(48,i);                                  %延时差 ns
%          delay(:,i)=round(delay(:,i)/(1e9/fs));          %移位点数，个数  10ns=dt
%          delay_point(1:96,i)=[delay(1:47,i);delay(49:97,i)]; %拼接96阵元
%           for j=1:N_elements              
%               RF_data(:,j)=circshift(RF_data(:,j),delay_point(j,i));
%           end
        
                
          RF_data=RF_data/max(max(RF_data));%查看回波信号是否对齐
          [N,M]=size(RF_data);
          
          if (theta == -sector/2 && theta == 0)
          figure;  
          for k=1:N_elements
         plot((0:N-1)/fs+times,RF_data(:,k)+k), hold on
          end
         hold off
          end
         
        RF_data=sum(RF_data,2);     %延时叠加
        display(:,i)=RF_data(1:7952);       %扫描线存储
           
        N=0:7932-1;
        p_display=N*(1e9/fs)*c/1000000;     %n*dt*c  单位mm
        x_display(:,i)=p_display*sin(theta);     %扫面线各个点坐标
        y_display(:,i)=p_display*cos(theta);

          i=i+1;
 
          
end
        H_display=abs(hilbert(display)); 
        H_display=20*log10(H_display)/N_elements;
         figure;       
        for k=1:61
         plot((0:N-1)/fs+times,H_display(:,k)+k), hold on
         end
         hold off
       







