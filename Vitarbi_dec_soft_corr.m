%% Декодирование по алгоритму Витерби
% Канал: AWGN, Модуляция: BPSK
%% Исходные данные
clc;clear
snr2 = 0:0.5:10;
per2 = zeros(1,length(snr2));
bpsk = comm.BPSKModulator;
de_bpsk = comm.BPSKDemodulator('DecisionMethod', 'Log-likelihood ratio');
bl=1e6; % Длина информационных символов
N=bl/2; % Длина блока 
for cntr=1:length(snr2)
info_bits = randi([0 1],1,bl); %генерируем сигнал
% info_bits = [1 1 1 1 1 0 0 0 0 0];
%% Сверточное кодирование. 
trellis = poly2trellis(3,[5 7]);
code=convenc(info_bits, trellis);
code_bit0=zeros((length(code)/2),1);
code_bit1=zeros((length(code)/2),1);
n=1;
for m=1:2:length(code)
 code_bit1(n,1)=code(m);
 n=n+1;
end
r=1;
for u=2:2:length(code)
 code_bit0(r,1)=code(u);
 r=r+1;
end
%% BPSK модулятор
mod_code_bit0=bpsk(code_bit0);
mod_code_bit1=bpsk(code_bit1);
fprintf('Encoding completed...\n');
%% Канал с шумом AWGN
tb = 15;
% Добавление шума
 y0 = awgn(mod_code_bit0, snr2(cntr), 'measured');
 y1 = awgn(mod_code_bit1, snr2(cntr), 'measured');
% Демодулятор
 y00=de_bpsk(y0);
 y11=de_bpsk(y1);
% Сверточное декодирование
for p=1:N:bl
% Решетчатая диаграмма декодера
 distance=zeros(4,N+1); % Состояние кодера 4,количество узлов N+1
 metric=zeros(8,N); % Количество метрик
 distance(1,1)=0; % Первое состояние
 distance(2,1)=-Inf; % Второе состояние
 distance(3,1)=-Inf; % Третье состояние
 distance(4,1)=-Inf; % Четвертоесостояние
for i=1:N
 metric(1,i)=((y00(i+p-1)*(1)) +(y11(i+p-1)*(1)));  % Метрика ветви от 00 к 00
 metric(2,i)=((y00(i+p-1)*(-1))+(y11(i+p-1)*(-1))); % Метрика ветви от 00 к 10
 metric(3,i)=((y00(i+p-1)*(-1))+(y11(i+p-1)*(1)));  % Метрика ветви от 10 к 01
 metric(4,i)=((y00(i+p-1)*(1)) +(y11(i+p-1)*(-1))); % Метрика ветви от 10 к 11
 metric(5,i)=((y00(i+p-1)*(-1))+(y11(i+p-1)*(-1))); % Метрика ветви от 01 к 00
 metric(6,i)=((y00(i+p-1)*(1)) +(y11(i+p-1)*(1)));  % Метрика ветви от 01 к 10
 metric(7,i)=((y00(i+p-1)*(1)) +(y11(i+p-1)*(-1))); % Метрика ветви от 11 к 01
 metric(8,i)=((y00(i+p-1)*(-1))+(y11(i+p-1)*(1)));  % Метрика ветви от 11 к 11
distance(1,i+1)=max(distance(1,i)+metric(1,i),distance(3,i)+metric(5,i)); % Максимальная метрика в 00
distance(2,i+1)=max(distance(1,i)+metric(2,i),distance(3,i)+metric(6,i)); % Максимальная метрика в 10
distance(3,i+1)=max(distance(2,i)+metric(3,i),distance(4,i)+metric(7,i)); % Максимальная метрика в 01
distance(4,i+1)=max(distance(2,i)+metric(4,i),distance(4,i)+metric(8,i)); % Максимальная метрика в 11
end
[~,state]=max(distance(:,N+1)); % Выбораем Максимальную метрику в последнем узле
for j = N:-1:1
    [state,decoded_bit]=prev_stage(state,distance(:,j),metric(:,j));
    decoded_bit_final(j+p-1)=decoded_bit; % Декодированная последовательность
end
end
[ERR1,ratio1] = biterr(info_bits,decoded_bit_final);
per2(cntr)=ratio1;
end
save('matlab2.mat', 'snr2', 'per2')
%------------------------------------------------------------------------------------------------------
function [prev_state,decoded_bit]=prev_stage(curr_state,distance_prev,metric)
if(curr_state==1)
    if(distance_prev(1)+metric(1) >= distance_prev(3)+metric(5))
        prev_state=1;decoded_bit=0;
    else
        prev_state=3;decoded_bit=0;
    end
end
if(curr_state==2)
    if(distance_prev(1)+metric(2) >= distance_prev(3)+metric(6))
        prev_state=1;decoded_bit=1;
    else
        prev_state=3;decoded_bit=1;
    end
end
if(curr_state==3)
    if(distance_prev(2)+metric(3) >= distance_prev(4)+metric(7))
        prev_state=2;decoded_bit=0;
    else
        prev_state=4;decoded_bit=0;
    end
end
if(curr_state==4)
    if(distance_prev(2)+metric(4) >= distance_prev(4)+metric(8))
        prev_state=2;decoded_bit=1;
    else
        prev_state=4;decoded_bit=1;
    end
 end
end
