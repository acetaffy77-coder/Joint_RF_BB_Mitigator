clear; clc; close all;
%%  *%构建双音输入信号*
% 基本参数
fs=25e6;                            % 采样率25MHz
T=5e-4;                             % 时长1ms
N = round(fs * T);                  % 样本数
t = (0:N-1)' / fs;                  % 时间向量

% 输入信号构建
f1=2.4e6;                          % 中心2.6MHz
f2=2.8e6;
fc=(f1+f2)/2;


P_input_dBm=-25;              % 总输入功率
P_input=10^(P_input_dBm/10)/1000;

A=sqrt(P_input)/2;            % 每个音调的幅度
x_bb=A*(exp(1j*2*pi*f1*t)+exp(1j*2*pi*f2*t));                 %基带双音信号


%%  *%基带信号非线性建模*
% 非线性与 I/Q 参数
a1 = 5.62;
a2 = -(84351 + 1j*74391);
a3 = 3.16;
a4 = -1588.7;

g_m = 0.99;
phi_m = 0.0628;

k1=(1+g_m*exp(-1j*phi_m))/2;          % 公式(10)
k2=(1-g_m*exp(1j*phi_m))/2;


% RF非线性
y_rf=a1*x_bb+a2*x_bb.*abs(x_bb).^2;           % 公式(7)

% I/Q不平衡
y_tilde_I=real(y_rf);          % 公式(11)
y_tilde_Q=g_m*sin(phi_m)*real(y_rf) + g_m*cos(phi_m)*imag(y_rf);

% BB非线性
y_BB_I=a3*y_tilde_I+a4*y_tilde_I.^3;
y_BB_Q=a3*y_tilde_Q+a4*y_tilde_Q.^3;
y_BB=y_BB_I+1j*y_BB_Q;

% 添加 AWGN 
SNR_dB = 61;
y_BB = awgn(y_BB, SNR_dB, 'measured');

%freqz(y_BB(N-1023:N),1, 1024, 'whole');

%%  *%实/复带通/阻滤波器设计*

% 实低通滤波器设计
f_pass=0.5e6;        % 通带宽度1MHz
f_stop=1e6;
lpFilt =designfilt('lowpassfir','FilterOrder',500,'PassbandFrequency',0.3e6,'StopbandFrequency',0.5e6,'SampleRate',25000000);


%实带阻滤波器设计    %中心频率2.6MHz，通带宽度1MHz
bsFilt=designfilt('bandstopfir','FilterOrder',500,'PassbandFrequency1',2.1e6,'StopbandFrequency1',2.3e6,'StopbandFrequency2',2.9e6, ...
    'PassbandFrequency2',3.1e6,'SampleRate',25000000);



% 实低通频移获得复带通
N_filt1=length(lpFilt.Numerator);
n_filt1=0:N_filt1-1;
H_C_bpFilt=lpFilt.Numerator .* exp(1j * 2 * pi * fc/fs * n_filt1);

H_R_bsFilt=bsFilt.Numerator;
N_filt2=length(H_R_bsFilt);
gd_bp = (length(H_C_bpFilt)-1)/2; % 复带通滤波器延迟
gd_bs = (length(H_R_bsFilt)-1)/2; % 实带阻滤波器延迟
max_delay = gd_bp + gd_bs;
% freqz(H_R_bsFilt,1, 4096, 'whole');
% freqz(H_C_bpFilt,1, 4096, 'whole');

%%  *%缓解*

% 基函数
x_hat=filter(H_C_bpFilt,1,y_BB);
y_BB_delay1=[zeros(gd_bp, 1); y_BB(1:end-gd_bp)];          % HCbp延迟
d=y_BB_delay1-x_hat;

x_hat_conj=conj(x_hat);

z_hat=x_hat.*abs(x_hat).^2;
z_hat_filt=filter(H_R_bsFilt,1,z_hat);
z_hat_conj=conj(z_hat);

x_hat_pow3=x_hat.^3;
x_hat_pow3_I=real(x_hat_pow3);
x_hat_pow3_Q=imag(x_hat_pow3);
x_hat_pow3_I_filt=filter(H_R_bsFilt,1,x_hat_pow3_I);
x_hat_pow3_Q_filt=filter(H_R_bsFilt,1,x_hat_pow3_Q);


% S1 S2时间对齐

x_hat_conj_delay2 = [zeros(gd_bs, 1); x_hat_conj(1 : end-gd_bs)];          %补bs
z_hat_conj_delay2 = [zeros(gd_bs, 1); z_hat_conj(1 : end-gd_bs)];
d = [zeros(gd_bs, 1); d(1 : end-gd_bs)];



u = [0.1, 0.01, 0.01, 0.01, 0.01];


af = [1e-5, 1e-11, 1e-11, 1e-11, 1e-11];

% 自适应缓解
M=1;             % 无记忆
w=zeros(5*M,1);
W_trace = zeros(5, N);          % 记录权重变化
x_est=zeros(N,1);
e=zeros(N,1);


for idx=max_delay + 1 : N 

    % 构建s
    s=[ x_hat_conj_delay2(idx),...
        z_hat_conj_delay2(idx),...
        z_hat_filt(idx),...
        x_hat_pow3_I_filt(idx),...
        x_hat_pow3_Q_filt(idx)].';           %公式(21)

    % 计算归一化u 公式(27)
 for idx2 = 1:length(u)
    u_N(idx2) = u(idx2) / (af(idx2) + sum(abs(s(idx2)^2))); 
end



    e(idx)=w'*s;         % 公式(24)-(26)
    x_est(idx)=d(idx)-e(idx);
    w=w+diag(u_N)*conj(x_est(idx))*s;

    

W_trace(:, idx) = w;

end


x_hat_delay2 = [zeros(gd_bs, 1); x_hat(1:end-gd_bs)];
x_final = x_hat_delay2 + x_est;

%freqz(x_final(N-1023:N), 1, 1024, 'whole');
%freqz(y_BB,1, 40960, 'whole');
% freqz(x_last,1, 1024, 'whole');
% freqz(x_est(N-1023:N),1, 1024, 'whole');
% freqz(x_hat(N-1023:N),1, 1024, 'whole');
% freqz(y_BB(N-1023:N),1, 1024, 'whole');
%freqz(x_hat2,1, 40960, 'whole');freqz(x_est(N-1023:N),1, 1024, 'whole');

%%  *%性能评估*
% 
% figure('Color','w','Name','权重收敛曲线');
% plot((1:N)/fs*1e3, abs(W_trace'));
% grid on;
% xlabel('n');
% legend('w1 ', 'w2 ', 'w3 ', 'w4 ', 'w5 ')

% 截取最后一段
plot_len = 2048;
range = N - plot_len + 1 : N;


% 计算功率谱密度 
[P_y_BB, f] = periodogram(y_BB_delay1(range), blackmanharris(plot_len), plot_len, fs, 'centered');
[P_x_final, ~] = periodogram(x_final(range), blackmanharris(plot_len), plot_len, fs, 'centered');

% 归一化到主分量功率 (dBc)
P_y_BB_dB = 10*log10(P_y_BB / max(P_y_BB));
P_x_final_dB = 10*log10(P_x_final / max(P_x_final));


figure('Color', 'w', 'Name', '非线性缓解前后对比');
plot(f/1e6, P_y_BB_dB, 'r', 'LineWidth', 1, 'DisplayName', '缓解前');
hold on;
plot(f/1e6, P_x_final_dB, 'b', 'LineWidth', 1, 'DisplayName', '缓解后');
grid on;
xlabel('基带频率 (MHz)');
ylabel('相对功率 (dB)');
title('双音信号缓解对比');
legend('Location', 'south');
axis([-(fs/2e6) (fs/2e6) -100 5]);

% f_bandOut = (f < 2.3e6) | (f > 2.9e6); 
f_bandOut = ((f < -1.0e6 )| (f > 1.9e6 & f < 2.1e6) | (f > 3.1e6 & f < 3.3e6)| (f > 6.0e6 ));


avg_power_before = mean(P_y_BB(f_bandOut));
avg_power_after = mean(P_x_final(f_bandOut));




ANDR = 10 * log10(avg_power_before / avg_power_after);


text(-10, -20, sprintf('ANDR = %.2f dB', ANDR), 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');

