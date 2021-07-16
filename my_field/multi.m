clc
clear
%% %%%%%%%%%%%%%%%%%%%     calc_scat_multi��������
%���������calc_scat=sum(calc_scat_multi);
%  Set initial parameters
f0=3e6; %  Transducer center frequency [Hz]
fs=100e6; %  Sampling frequency [Hz]
c=1540; %  Speed of sound [m/s]
lambda=c/f0; %  Wavelength [m]
height=5/1000; %  Height of element [m]
width=1/1000; %  Width of element [m]
kerf=width/4; %  Distance between transducer elements [m]
N_elements=32; %  Number of elements
focus=[0 0 40]/1000;  %  Initial electronic focus  r������xdc_focus ��THʧЧ
Ts =1/fs;             % ����ʱ���� [s]
%% %%%%%%%%%%%%%%%%%%%%%%
field_init(-1);
%%  %%%%%%%%%%%%%%%%%

Th = xdc_linear_array (N_elements, width, height, kerf, 1, 1,focus);
Th2 = xdc_linear_array (N_elements, width, height, kerf, 1, 1,focus);
%����ͽ�������

BW=1*f0;
impulse_response=sin(2*pi*f0*(0:Ts:2/BW));
figure;
subplot 221
plot(impulse_response);
title('�Ӵ�ǰ������������');
impulse_response=impulse_response.*hanning(max(size(impulse_response)))';

%������Ӧ����

    xdc_impulse (Th, impulse_response);
    xdc_impulse (Th2, impulse_response);
    %�������������ͬ

excitation=sin(2*pi*f0*(0:Ts:2/BW));
excitation=excitation.*hanning(max(size(excitation)))';
xdc_excitation (Th, excitation);
%�����ź�  �����źź��˲�����ͬ��Ϊƥ���˲���ƥ���˲�Ƶ����ʱ����ƽ��t1

output_signal=conv(excitation,impulse_response);

subplot 222
plot(impulse_response);
title('���г����Ӧ����');
subplot 223
plot(excitation);
title('��������');
subplot 224
plot(output_signal);
title('ʵ�ʷ����źŲ���');


LengthConv=(length(excitation) + 2*length(impulse_response)-2)*0.5;
%���ǰ������ĳ���   �����M+L-1=75+75-1


%% %%%%%%%%%%%%%%%%%%%%%%%%%

%  Do the calculation
          xdc_center_focus (Th,[0 0 0]);
          xdc_focus (Th, 0, [0 0 50]/1000);      %����۽�
          xdc_center_focus (Th2,[0 0 0]);        %���վ۽�
          xdc_focus (Th2, 0, [0 0 50]/1000);
        
[v,t]=calc_scat_multi (Th, Th2, [0 0 20]/1000, 1); %ֻ��һ��ɢ��㣬ɢ��ǿ��Ϊ1
%  Plot the individual responses                   %ÿ��ͨ���ز�����
[v2,t1]=calc_scat(Th, Th2, [0 0 20]/1000, 1);
[N1,M1]=size(v2);
figure;
subplot(311)
[N,M]=size(v);
% v=v/max(max(v));
for i=1:N_elements
plot((0:N-1)/fs+t,v(:,i)+i), hold on
end
hold off
title('Individual traces')
xlabel('Time [s]')
ylabel('Normalized response')
subplot(312)
plot((0:N-1)/fs+t,sum(v'))
title('Summed response')
xlabel('Time [s]')
ylabel('Normalized response')

subplot(313)
plot((0:N1-1)/fs+t1,v2);
title('calc_scat response')
xlabel('Time [s]')
ylabel('Normalized response')











