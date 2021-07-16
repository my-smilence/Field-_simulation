clc
clear
field_init(-1);
tic
%% %%%%%%%%%%%%%%%%%%%%
%  Example of use of the new Field II program running under Matlab
%
%  This example shows how a linear array B-mode system scans an image
%
%  This script assumes that the field_init procedure has been called
%
%  Example by Joergen Arendt Jensen, Version 2.0, March 22, 2011.
%  Generate the transducer apertures for send and receive
f0=7e6; %  Transducer center frequency [Hz]
fs=100e6; %  Sampling frequency [Hz]
c=1540; %  Speed of sound [m/s]
lambda=c/f0; %  Wave length [m]
width=lambda; %  Width of element
height=5/1000;  %  Height of element [m]
kerf=width/20; %  Kerf [m]
focus=[0 0 50]/1000; %  Fixed focal point [m]
N_elements=192; %  Number of elements in the transducer
N_active=64; %  Active elements in the transducer
%  Set the sampling frequency
set_sampling(fs);
Ts=1/fs;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%  此处表示发射变迹/加窗   激励和响应函数全部汉宁窗----改变发射信号波形！！！
LengthConv=(length(excitation) + 2*length(impulse_response)-2)*0.5;
%这是半个脉冲的长度   卷积：M+L-1=75+75-1

%% %%%%%%%%%%%%%%%%%%%%%    设定散射点
%  Load the computer phantom
[phantom_positions, phantom_amplitudes] = cyst_phantom(100000);  %计算散斑

% dz=5/1000;        %  散射点之间轴间距[m]
% dx=[repmat(-2/1000,1,5),repmat(2/1000,1,5)];        %  x轴方向投影为两个点
% z_start=35/1000;  %  开始深度 [m]
% Npoints=5;       %  共计10*2个点
% 
% positions = [dx; zeros(1,Npoints*2); repmat((1:Npoints)*dz+z_start,1,2)]';
% amp=ones(Npoints*2,1);
% %[-1 0 20] [1 0 20] ------- [-1 0 60][1 0 60]


%% %%%%%%%%%%%%%%%%%%%%%%%%%%
%  Do linear array imaging
no_lines=N_elements-N_active+1; %  Number of A-lines in image  线阵A扫扫描线：N-N—active+1
dx=width; %  Increment for image
z_focus=50/1000;
%  Pre-allocate some storage  成像数据（1*扫描线）
image_data=zeros(1,no_lines);

for i=1:no_lines%  Find position for imaging
x=(i-1-no_lines/2)*dx;
%  Set the focus for this direction
xdc_center_focus (Th, [x 0 0]);
xdc_focus (Th, 0, [x 0 z_focus]);
xdc_center_focus (Th2, [x 0 0]);
xdc_focus (Th2, 0, [x 0 z_focus]);

% %  Set the active elements using the apodization   接收变迹 xdc_apodization
apo=[zeros(1, i-1) hamming(N_active)' zeros(1, N_elements-N_active-i+1)];
 xdc_apodization (Th, 0, apo);
 xdc_apodization (Th2, 0, apo);
%发射和接收变迹，全部汉明窗  

%  Calculate the received response
[v, t1]=calc_scat(Th, Th2,phantom_positions,phantom_amplitudes);
%  Store the result
image_data(1:max(size(v)),i)=v;
times(i) = t1;

%  image_data=awgn( image_data, 60, 'measured'); %添加高斯白噪声；60=snr指定了每一个采样点信号与噪声的比率，单位为dB

end

%  Free space for apertures
xdc_free (Th)
xdc_free (Th2)
%  Adjust the data in time and display it as
%  a gray scale image
min_sample=min(times)*fs;
for i=1:no_lines        %希尔伯特变换取包络
rf_env=abs(hilbert([zeros(round(times(i)*fs-min_sample),1); image_data(:,i)])); %这里取最小时间部分作为起始
env(1:size(rf_env,1),i)=rf_env;
end
%  make logarithmic compression to a 60 dB dynamic range
%  with proper units on the axis
env_dB=20*log10(env);
env_dB=env_dB-max(max(env_dB));
env_gray=127*(env_dB+60)/60;  %60dB动态范围表示,且用灰度表示


%%  %%%%%%%%%%%%%%%%%  成像插值
depth=((0:size(env,1)-1)+min_sample)/fs*c/2; %差了0.8mm
x=((1:no_lines)-no_lines/2)*dx;

n=1:length(x);
ni=1:0.5:length(x);
z=1:length(depth);
zi=1:0.5:length(depth);
xi=interp1(n,x,ni,'linear');
depthi=interp1(z,depth,zi,'linear');    %拓展成像数据点数，四倍插值


env_gray_interpt=interp2(env_gray,2);  %zi=interp2(Z,n);做n次递归计算，在两元素间插入二维插值

%ZI = interp2(X,Y,Z,XI,YI,method)   用指定的算法method 计算二维插值：
%[X, Y] = meshgrid(x,y);
%% %%%%%%%%%%%%%%%%%%%%%%%
image(x*1000, depth*1000, env_gray)
xlabel('Lateral distance [mm]')
ylabel('Depth [mm]')
axis('image')
colormap(gray(128))
title('Image of cyst phantom (60 dB dynamic range)')

toc
