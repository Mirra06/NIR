%% Декодирование по алгоритму Витерби
% Канал: AWGN, Модуляция: BPSK
%% Исходные данные
clc;clear
snr0 = 0:10;
per0 = zeros(1,length(snr0));
bpsk = comm.BPSKModulator;
de_bpsk = comm.BPSKDemodulator;
M = 2; % Порядок модуляции
bl = 1e6; % Длина информационных символов
N = bl/2; % Длина блока 
ConstraintLength = 3; % Продолжительность ограничения
Poly_1 = 5;
Poly_2 = 7;
for cntr=1:length(snr0)
%% Создание сигнала
info_bits = randi([0 1],1,bl); %генерируем сигнал
% info_bits = [1 1 0 1 1 0 1 1]; %генерируем сигнал
%% Сверточное кодирование 
trellis = poly2trellis(ConstraintLength,[Poly_1 Poly_2]);
code = convenc(info_bits, trellis);
code_bit0 = zeros((length(code)/2),1);
code_bit1 = zeros((length(code)/2),1);
n=1;
for m=1:2:length(code)
 code_bit0(n,1) = code(m);
 n=n+1;
end
r=1;
for u=2:2:length(code)
 code_bit1(r,1)=code(u);
 r=r+1;
end
%% BPSK модулятор
mod_code_bit0=bpsk(code_bit0);
mod_code_bit1=bpsk(code_bit1);
%% Канал с шумом AWGN
 [y0, N0_0] = noise_awgn(mod_code_bit0, snr0(cntr), M, bl); % Добавление шума
 [y1, N0_1] = noise_awgn(mod_code_bit1, snr0(cntr), M, bl);
%% BPSK демодулятор
 y00 = demodulation(y0, N0_0);
 y11 = demodulation(y1, N0_0);
%% Сверточное декодирование функцией
q = 1;
d = 1;
y = zeros(1, 2*bl);
while (q < 2*bl)
    y(q) = y00(d);
    y(q+1) = y11(d);
    d = d + 1;
    q = q+2;
end
decodedout = vitdec(y,trellis, 10,'term','unquant');

[ERR1,ratio2] = biterr(info_bits,decodedout);
per0(cntr)=ratio2;
save('matlab0.mat', 'snr0', 'per0')
end
