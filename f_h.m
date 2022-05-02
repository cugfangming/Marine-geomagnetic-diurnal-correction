clear;clc;
%% 1D mesh model.sum
seawater=[10 15 20 40 50 60 90 100 200 200 300 300 400 400 500 500 600 600];%海水4385m
layer2 = [100.0 100.0 100.0 100.0 100.0 100.0 100.0 100.0 100.0 100.0];%沉积层1000m
layer3 = [100.0 100.0 100.0 100.0 100.0 100.0 100.0 100.0 100.0 100.0];%空隙充满海水的岩层1000m
layer4 = [100.0 100.0 100.0 100.0 100.0 100.0 100.0 100.0 100.0 100.0 100.0 100.0 100.0 100.0 100.0 100.0 100.0 100.0 100.0 100.0];
dz = [seawater layer2 layer3 layer4];
%电导率
sig(1:18) = 3.3;
sig(19:28) = 1.0;
sig(29:38) = 0.02;
sig(39:58) = 0.00002;
%把每层的深度转换成总深度
for k=1:length(dz)
    hm=0;
        for i=1:k
            hm=hm+dz(i);
        end
        h(k)=hm;%h(k)为深度
end
h=[0,h];
%%
%数据处理
dataH=readmatrix('磁总场强度.txt');%读入数据
x=dataH(:,4)';%台站坐标 北纬49.649879°；西经37.783464°；
L = length(x);%%数据采样点数
t=60:60:60*L;%单位s
figure(1)%原始数据时间域，频率域
subplot(2,1,1)
plot(t,x);title('原始数据时间域Amp');xlabel('时间t/s');ylabel('Amp');

Fs=1/60;%fft后有效的频段范围 2也行1/(60*L)-1/120
Y=fft(x);
f = Fs*(0:(L/2))/L;

P2 = abs(Y/L);
P1 = P2(1:L/2+1);%只取左半轴，...
P1(2:end-1) = 2*P1(2:end-1);
subplot(2,1,2)
plot(f,P1) 
semilogx(f,P1);title('原始数据频率域Amp');xlabel('频率f/Hz');ylabel('Amp');
%% 
%滤波部分与下面ifft不共存，滤波是为了画二维图
figure(2)
fss=find(f<4/(60*L));%这一段频率对应的幅值不准
subplot(2,1,1)
semilogx(f,P1);title('原始数据频率域Amp');xlabel('频率f/Hz');ylabel('Amp');

% Y(fss)=0;
% Y(L-fss+1)=0;

subplot(2,1,2)
f=f(fss(end)+1:end);
%P1(fss)=0;
P1=P1(fss(end)+1:end);
semilogx(f,P1);title('原始数据滤波频率域Amp');xlabel('频率f/Hz');ylabel('Amp');
Y=Y(fss(end)+1:L-fss(end));
%%
for i=1:length(f)
    hy(:,i) = mt1dtmfh(f(i),dz,sig,Y(i));%第i个频率f(i),不同深度对应的磁场Amp
end
%%
amp=abs(hy/L);%幅值
amp(:,2:end-1)=2*amp(:,2:end-1);
phs=atan2(imag(hy),real(hy))* 180.0 / pi;%相位
%%
%与figure(2)共存，因为前面低频也不准，所以省去
figure(3)%画不同深度下频率域的Amp图
i=4;
semilogx(f,amp(i,:));title('频率域Amp 深度: /m',h(i));xlabel('频率f/Hz');ylabel('Amp');
%%
%进行ifft
figure(4)%某一深度的ifft后时间域Amp
i=4;
Q=hy(i,:);
Q=[0*fss,Q];
Q(L/2+1:L)=Q(L/2:-1:1);
xfilt=ifft(Q);
plot(t,xfilt);title('时间域Amp 深度: /m',h(i));xlabel('时间t/s');ylabel('Amp');
%%
figure(5);%二维的图
pcolor(log10(f),h,log10(amp));shading flat
colorbar
title('不同频率不同深度下的磁场幅值');
xlabel('频率log10(f)/Hz');
set(gca,'XDir','reverse');
ylabel('深度h/m');
axis on;
grid off;