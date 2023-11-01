%% Декодирование по алгоритму Витерби
% Канал: AWGN, Модуляция: BPSK
%% Исходные данные
clc;clear
snr1 = 0:10;
per1 = zeros(1,length(snr1));
bpsk = comm.BPSKModulator;
de_bpsk = comm.BPSKDemodulator;
M = 2; % Порядок модуляции
k = log2(M); %кол-во бит на символ
bl = 1e4; % Длина информационных символов
N = bl/2; % Длина блока 
ConstraintLength = 3; % Продолжительность ограничения
Poly_1 = 5;
Poly_2 = 7;
for cntr=1:length(snr1)
%% Создание сигнала
info_bits = randi([0 1],1,bl); %генерируем сигнал
% info_bits = [1 1 0 1 1 0]; %генерируем сигнал
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
mod_code_bit0=bpsk(code_bit0);
mod_code_bit1=bpsk(code_bit1);
fprintf('Implementing iteration_%d', cntr);fprintf('...\n');
%% Канал с шумом AWGN
 [y0, N0_0] = noise_awgn(mod_code_bit0, snr1(cntr), M, bl); % Добавление шума
 [y1, N0_1] = noise_awgn(mod_code_bit1, snr1(cntr), M, bl);
%% BPSK демодулятор
 y00 = de_bpsk(y0);
 y11 = de_bpsk(y1);
 %% Преобразование в двоичные данные
 codedDataOut0 = int2bit(y00, k); 
 codedDataOut1 = int2bit(y11, k);
%% Сверточное декодирование 
% Определение состояний декодера
R = trellis.numOutputSymbols/trellis.numInputSymbols;
len = length(0:2^(ConstraintLength-1)-1); % Количество состояний
binStr = de2bi(0:2^(ConstraintLength-1)-1); % Состояния
binStr = flip(binStr,2);
a = zeros(len, ConstraintLength);
b = zeros(len, ConstraintLength);
out_code_comb_1 = zeros(len, R);
out_code_comb_2 = zeros(len, R);
Poly_11 = oct2poly(Poly_1);
Poly_21 = oct2poly(Poly_2);
for p = 1:len
    a(p,:) = [0, binStr(p,:)];
    element1 = sum(a(p,:) .* Poly_11);
    element2 = sum(a(p,:) .* Poly_21);
    out_code_comb_1(p,:) = [module2(element1), module2(element2)];
end
a(:,ConstraintLength) = [];
% out_code_comb_1 = flip(out_code_comb_1, 2);
for p = 1:len
    b(p,:) = [1, binStr(p,:)];
    element1 = sum(b(p,:) .* Poly_11);
    element2 = sum(b(p,:) .* Poly_21);
    out_code_comb_2(p,:) = [module2(element1), module2(element2)];
end
b(:,ConstraintLength) = [];
% out_code_comb_2 = flip(out_code_comb_2, 2);
% Решетчатая диаграмма декодера
for p = 1:N:bl
    distance = zeros(len, N+1); % Состояния кодера, количество узлов N+1
    metric = zeros(2*len, N); % Количество метрик
    distance(1,1) = 0; % Первое состояние
    distance(2:end, 1) = Inf; % Остальные состояния

    for i = 1:N
    m = 1;
    s = 1;
        for j = 1:2:2*len
            metric(j,i)=sum(xor(out_code_comb_1(m,:),[y00(i+p-1),y11(i+p-1)])); % Метрика ветви
            m = m+1;
        end
        m = 1;
        for w = 1:4:2*len              
            distance(m,i+1) = min(distance(s,i)+metric(w,i),distance(s+1,i)+metric(w+2,i)); % Минимальная метрика
            s = s+2;
            m = m+1;
        end     
    
    m = 1;
    s = 1;
        for j = 2:2:2*len
            metric(j,i)=sum(xor(out_code_comb_2(m,:),[y00(i+p-1),y11(i+p-1)])); % Метрика ветви
            m = m+1;
        end
        m = 1;
        for w = 2:4:2*len
            distance(m+len/2,i+1) = min(distance(s,i)+metric(w,i),distance(s+1,i)+metric(w+2,i)); % Минимальная метрика
            s = s+2;
            m = m+1;
        end 
    end

[~,state]=min(distance(:,N+1)); % Выбираем минимальную метрику в последнем узле
for h = N:-1:1
    [state,decoded_bit]=prev_stage(state,distance(:,h),metric(:,h),len);
    decoded_bit_final(h+p-1)=decoded_bit; % Декодированная последовательность
end
end
%%
    [ERR1,ratio2] = biterr(info_bits,decoded_bit_final);
    per1(cntr)=ratio2;
end
save('matlab1.mat', 'snr1', 'per1')
