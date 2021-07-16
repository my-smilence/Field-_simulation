%  Example of use of the new Field II program running under Matlab
%
%  This example shows how a phased array B-mode system scans an image
%
%  This script assumes that the field_init procedure has been called
%
%  Example by Joergen Arendt Jensen, Nov. 28, 1995.
%  Generate the transducer apertures for send and receive
clc
clear
f0=3e6; %  Transducer center frequency [Hz]
fs=100e6; %  Sampling frequency [Hz]
c=1540; %  Speed of sound [m/s]
lambda=c/f0; %  Wavelength
height=5/1000;  %  Height of element [m]
kerf=0.1/1000; %  Kerf [m]
focus=[0 0 70]/1000; %  Fixed focal point [m]
N_elements=128;
width=lambda/2;
Ts=1/fs;
%% %%%%%%%%%%%%%%%%%%%%
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
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�趨��ɨ�Ƕ�-30-30��60���ܽǶȣ� sector=60��=60*pi/180
%�趨�������� 60��ɨ����  no_lines=60
%��ÿ�β����Ƕ�Ϊ  theta_d=sector/no_lines;
%theta=-sector/2:theta_d:sector/2; %�趨image_data=zeros(N,no_lines);�������ݴ洢

%  Do phased array imaging
point_position=[-1 0 50;1 0 50]/1000; %  Position of the point to be imaged

no_lines=60; %  Number of A-lines in image Aɨ����񷽷�  ����50��ɨ����
sector=60*pi/180; %  Size of image sector   ����Ƕ�/������С

d_theta=sector/no_lines; %  Increment in angle for 90 deg. image
%  Pre-allocate some storage
image_data=zeros(800,no_lines);

theta= -sector/2;
for i=1:no_lines
%  Set the focus for this direction
xdc_focus (Th, 0, [70*sin(theta) 0 70*cos(theta)]/1000);
xdc_center_focus (Th,[0 0 0]);
xdc_center_focus (Th2,[0 0 0]);
xdc_focus (Th2, 0, [70*sin(theta) 0 70*cos(theta)]/1000);
%  Calculate the received response
[v, t1]=calc_scat(Th, Th2, point_position, [1;1]);
%  Store the result
image_data(1:max(size(v)),i)=v';
times(i) = t1;
%  Steer in another angle
theta = theta + d_theta;
end
%  Here the display of the data is inserted
[n,m]=size(image_data);
%  Here the display of the data is inserted
image_data_plot=image_data/max(max(image_data));
for  i=1:no_lines
plot(image_data_plot(:,i)+i);hold on
end
hold off

%  Free space for apertures
xdc_free (Th)
xdc_free (Th2)
%% %%%%%%
min_sample=min(times)*fs;
for i=1:no_lines        %ϣ�����ر任ȡ����
rf_env=abs(hilbert([zeros(round(times(i)*fs-min_sample),1); image_data(:,i)])); %����ȡ��Сʱ�䲿����Ϊ��ʼ
env(1:size(rf_env,1),i)=rf_env;
end











