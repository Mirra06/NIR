%% Декодирование по алгоритму Витерби
% Канал: AWGN, Модуляция: BPSK
%% Исходные данные
clc;clear
snr = 0:10;
per = zeros(1,length(snr));
bpsk = comm.BPSKModulator;
de_bpsk = comm.BPSKDemodulator;
M = 2; % Порядок модуляции
k = log2(M); %кол-во бит на символ
bl = 1e6; % Длина информационных символов
N = bl/2; % Длина блока 
ConstraintLength = 4; % Продолжительность ограничения
Poly_1 = 15;
Poly_2 = 17;
for cntr=1:length(snr)
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
%% Преобразование в целочисленные данные
dataSymbolsIn0 = bit2int(code_bit0, k);
dataSymbolsIn1 = bit2int(code_bit1, k);
%% BPSK модулятор
mod_code_bit0=bpsk(dataSymbolsIn0);
mod_code_bit1=bpsk(dataSymbolsIn1);
%% Канал с шумом AWGN
 [y0, N0_0] = noise_awgn(mod_code_bit0, snr(cntr), M, bl); % Добавление шума
 [y1, N0_1] = noise_awgn(mod_code_bit1, snr(cntr), M, bl);
%% BPSK демодулятор
 y00 = de_bpsk(y0);
 y11 = de_bpsk(y1);
 %% Преобразование в двоичные данные
 codedDataOut0 = int2bit(y00, k); 
 codedDataOut1 = int2bit(y11, k);
%% Сверточное декодирование функцией
q = 1;
d = 1;
y = zeros(1, 2*bl);
while (q < 2*bl)
    y(q) = codedDataOut0(d);
    y(q+1) = codedDataOut1(d);
    d = d + 1;
    q = q + 2;
end
decodedout = vitdec(y,trellis, 10,'trunc','hard');

[ERR1,ratio2] = biterr(info_bits,decodedout);
per(cntr)=ratio2;
save('matlab.mat', 'snr', 'per')
end


