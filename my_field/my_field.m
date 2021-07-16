clc
clear
close all;
%%
%��ʼ��
field_init(-1);

%%    %%%%%%�����趨%%%%%%%%%
f0                      =   4e6;              % ����Ƶ�� [Hz]
c                       =   1540;             % Speed of sound [m/s]
fs                      =   100e6;            % ����Ƶ�� [Hz]
Ts                      =   1/fs;             % ����ʱ���� [s]
lambda                  =   c/f0;             % ���� [m]
height                  =   6e-3;             % ��Ԫ���6mm [m]
pitch                   =   lambda/2;         % ��Ԫ���ļ��벨��[m]    
kerf                    =   pitch/10;         % ��Ԫ��϶��� m
width                   =   pitch-kerf;       % ��Ԫ����[m]    
N_elements              =   96;               % ��Ԫ��

%%  %%%%%%%%%%%%%%%%%%%%%% ��ʼ�������趨 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
focus = [0 0 100];

Th = xdc_linear_array (N_elements, width, height, kerf, 1, 1,focus);
Th2 = xdc_linear_array (N_elements, width, height, kerf, 1, 1,focus);
%����ͽ�������

BW=1*f0;
impulse_response=sin(2*pi*f0*(0:Ts:2/BW));
impulse_response=impulse_response.*hanning(max(size(impulse_response)))';
%������Ӧ����

    xdc_impulse (Th, impulse_response);
    xdc_impulse (Th2, impulse_response);
    %�������������ͬ

excitation=sin(2*pi*f0*(0:Ts:2/BW));
excitation=excitation.*hanning(max(size(excitation)))';
xdc_excitation (Th, excitation);
%�����ź�  �����źź��˲�����ͬ��Ϊƥ���˲���ƥ���˲�Ƶ����ʱ����ƽ��t1

LengthConv=(length(excitation) + 2*length(impulse_response)-2)*0.5;
%���ǰ������ĳ���   �����M+L-1=75+75-1

%% %%%%%%%%%%%%% �趨��ʾ���� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lx            = lambda/10;  %����1/10����
% imaging_start = 28/1000;
% imaging_depth = 4/1000; % 30/1000; 
% image_width   = 10e-3;   
% Lz             = c*Ts/2; % L��z�Ჽ����Ϊ�벨��
% z_start       = imaging_start;
% z_stop        = imaging_start+imaging_depth;  %Z=��28��L��28+4��;
% X             = -image_width/2: Lx: image_width/2;  %�������Ҹ�2.5mm��ȷ�Χ
% 
% Z             = z_start: Lz : z_stop;  %����ʾ29-31mm��ȣ���ʱ�۽���100mm

%% %%%%%%%%%%%%%%%%%%%%%%%%%% �趨ɢ��� %%%%%%%%%%%%%%%%%%%%%%%%%%
dz=10/1000;        %  ɢ���֮������[m]
dx=[repmat(-1/1000,1,5),repmat(1/1000,1,5)];        %  x�᷽��ͶӰΪ������
z_start=10/1000;  %  ��ʼ��� [m]
Npoints=5;       %  ����10*2����

positions = [dx; zeros(1,Npoints*2); repmat((1:Npoints)*dz+z_start,1,2)]';
amp=ones(Npoints*2,1);
%[-1 0 20] [1 0 20] ------- [-1 0 60][1 0 60]

%% %%%%%%%%%%%%%%%%% 
td        = LengthConv*2; %��������ĳ���/����
tdd       = (td-1)/2;   %

Lp    = 4;       % Lp �����г��ȣ�һ��N
P     = N_elements-Lp+1;    % M-L+1= 96/2+1=48+1=49  �ռ�ƽ��ϵ��
e     = ones(Lp,1);         % ��λʸ����һ��Ϊ1

yy    = zeros(25000,10000); %��ʼ����ַ����

W     = (N_elements*pitch); % -kerf

Tapo  = hamming(N_elements)'; 
xdc_apodization (Th, 0, Tapo); %����伣

%% %%%%%%%%%%%%%%%%%%% �۽��趨��ʱ %%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 1:N_elements;
sector=60*pi/180;   %20����ɨ;-10~10
% theta= -sector/2;
i=1;

for theta = -sector/2:pi/180:sector/2
         xdc_center_focus (Th,[0 0 0]);
         xdc_focus (Th, 0, [100*sin(theta)/1000 0 100*cos(theta)]/1000);
         xdc_center_focus (Th2,[0 0 0]);
         xdc_focus (Th2, 0, [100*sin(theta)/1000 0 100*cos(theta)]/1000);
        [v1,times]=calc_scat_multi(Th, Th2, positions, amp);    %(0:N-1)/fs+t  times=��һ�����յ���ɢ�����ź�ʱ��
        
        v1_zeropadded= [zeros(round(times*fs) , N_elements); v1];%������ʼʱ��
        
        v(1:length(v1_zeropadded),:)=v1_zeropadded;
       
        %Pass the signal through an AWGN channel having a 60 dB signal-to-noise ratio (SNR)
        %adding white Gaussian noise
        RF_data=awgn( v, 60, 'measured'); %��Ӹ�˹��������60=snrָ����ÿһ���������ź��������ı��ʣ���λΪdB
                  
       
%          x_length=(N_elements-1)*pitch;
%          x=-x_length/2-(pitch/2):pitch:x_length/2+(pitch/2);     %����λ������xy,��λm
%          
%          x_position=50*sin(theta)/1000;       %��ʱ�������� ��λmm
%          y_position=50*cos(theta)/1000;
%          
%          delay(:,i)=1000000000*sqrt((x-x_position).^2+y_position.^2)/1540;    %�����ÿ����Ԫ�ľ��� ns��ʱ 
%          delay(:,i)=delay(:,i)-delay(48,i);                                  %��ʱ�� ns
%          delay(:,i)=round(delay(:,i)/(1e9/fs));          %��λ����������  10ns=dt
%          delay_point(1:96,i)=[delay(1:47,i);delay(49:97,i)]; %ƴ��96��Ԫ
%           for j=1:N_elements              
%               RF_data(:,j)=circshift(RF_data(:,j),delay_point(j,i));
%           end
        
                
          RF_data=RF_data/max(max(RF_data));%�鿴�ز��ź��Ƿ����
          [N,M]=size(RF_data);
          
          if (theta == -sector/2 && theta == 0)
          figure;  
          for k=1:N_elements
         plot((0:N-1)/fs+times,RF_data(:,k)+k), hold on
          end
         hold off
          end
         
        RF_data=sum(RF_data,2);     %��ʱ����
        display(:,i)=RF_data(1:7952);       %ɨ���ߴ洢
           
        N=0:7932-1;
        p_display=N*(1e9/fs)*c/1000000;     %n*dt*c  ��λmm
        x_display(:,i)=p_display*sin(theta);     %ɨ���߸���������
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
       







