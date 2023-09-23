%% polar code
clc; clear
load('Q.mat')

snr = 0:0.5:20;
per = zeros(1,length(snr));

for cntr=1:length(snr)
    N = 16;
    K = N/2;
    msg = [0 0 0 1 0 0 0 1];
%     msg = randi([0 1], 1, K);
    n = log2(N);
    Q1 = Q(Q<=N);
    F = Q1(1:N-K);
    F_L = F(F<N/2);
    F_R = F(F>N/2);
    W = size(F_R);
    for i=1:1:W
        F_R(i) = i;
    end
    polar_code = zeros(1, N);
    polar_code(Q1(N-K+1:end)) = msg; 
    
    m = 1; %кол-во битов, которые нужно объединить
    for d = n-1:-1:0
        for i = 1:2*m:N
            a = polar_code(i:i+m-1);
            b = polar_code(i+m:i+2*m-1);
            polar_code(i:i+2*m-1) = [mod(a+b,2) b];
        end
        m = m*2;
    end
%% Channel 
    bpskModulator = comm.BPSKModulator;
    modSig = bpskModulator(polar_code'); %модуляция    
    r = awgn(modSig',snr(cntr)); % добавление шума
%% Polar decode
    L = zeros(n, N/2); %массив для всех сообщений
    L_soft = zeros(1, N);
    ucap_l = zeros(n+1, N); %массив для всех сообщений
    i = 0; %узел
    d = 0; %глубина
    state = zeros(1, 2*N-1); %вектор состояний узлов 
    L_soft(1,:) = r;
    mode = 0;%состояние работы декодера, работа завершена=1, не завершена=0
    while (mode == 0)
        if d == n %лист или нет
            if any(F_L == (i+1)) %лист
                ucap_l(n+1,i+1) = 0;
            else
                if L(n,i+1) >= 0
                    ucap_l(n+1,i+1) = 0;
                else
                    ucap_l(n+1,i+1) = 1;
                end
            end
            i = floor(i/2); d = d-1;
            else
                p_state = (2^d-1) + i + 1;%положение состояния узла в векторе состояний узла
            if (i == 0 && d == 0 && state(p_state) == 0)
                in_msg = 2^(n-d);
                Ln = L_soft(d+1, in_msg*i+1:in_msg*(i+1)); %входящие убеждения
                r1 = Ln(1:in_msg/2);
                r2 = Ln(in_msg/2+1:end);
                i = 2*i; d = d + 1;
                in_msg = in_msg / 2;
                L(d, in_msg*i+1:in_msg*(i+1)) = f(r1,r2);
                state(p_state) = 1;
            else
                if state(p_state) == 0 %не лист
                    in_msg = 2^(n-d);
                    Ln = L(d, in_msg*i+1:in_msg*(i+1)); %входящие убеждения
                    r1 = Ln(1:in_msg/2);
                    r2 = Ln(in_msg/2+1:end);
                    i = 2*i; d = d + 1;
                    in_msg = in_msg / 2;
                    L(d, in_msg*i+1:in_msg*(i+1)) = f(r1,r2);
                    state(p_state) = 1;
                else
                    if state(p_state) == 1
                        in_msg = 2^(n-d);
                        Ln = L(d, in_msg*i+1:in_msg*(i+1)); %входящие убеждения
                        r1 = Ln(1:in_msg/2);
                        r2 = Ln(in_msg/2+1:end);
                        l_i = 2*i; l_d = d + 1;
                        l_in_msg = in_msg / 2;
                        l_ucap = ucap_l(l_d+1, l_in_msg*l_i+1:l_in_msg*(l_i+1));
                        in_msg = in_msg / 2;
                        i = 2*i+1; d = d + 1;
                        L(d, in_msg*i+1:in_msg*(i+1)) = g(r1, r2, l_ucap);
                        state(p_state) = 2;
                    else
                        in_msg = 2^(n-d);
                        l_i = 2*i; r_i = 2*i+1; c_d = d + 1;
                        c_in_msg = in_msg / 2;
                        l_ucap = ucap_l(c_d+1, c_in_msg*l_i+1:c_in_msg*(l_i+1));
                        r_ucap = ucap_l(c_d+1, c_in_msg*r_i+1:c_in_msg*(r_i+1));
                        ucap_l(d+1, in_msg*i+1:in_msg*(i+1)) = [mod(l_ucap+r_ucap,2) r_ucap];
                        i = floor(i/2); d = d-1;
                        if (i == 0 && d == 0 && state(p_state) == 2)
                            mode = 1;  
                        end
                    end
                end
            end
        end
    end
    
    %-----------------------------------------------------------------------------------------------------
    state = [1 zeros(1, 2*N-2)]; %вектор состояний узлов
    mode = 0; %состояние работы декодера, работа завершена=1, не завершена=0
    ucap_r = zeros(n+1, N);
    while (mode == 0)
        if d == n %лист или нет
            if any(F_R == (i+1)) %лист
                ucap_r(n+1,i+1) = 0;
            else
                if L(n,i+1) >= 0
                    ucap_r(n+1,i+1) = 0;
                else
                    ucap_r(n+1,i+1) = 1;
                end
            end
            if i == floor((N-1)/2)
                mode = 1;
            else 
                i = floor(i/2); d = d-1;
            end
        else
            p_state = (2^d-1) + i + 1;%положение состояния узла в векторе состояний узла
            if (i == 0 && d == 0 && state(p_state) == 1)
                in_msg = 2^(n-d);
                Ln = L_soft(d+1, in_msg*i+1:in_msg*(i+1)); %входящие убеждения
                r1 = Ln(1:in_msg/2);
                r2 = Ln(in_msg/2+1:end);
                l_i = 2*i; l_d = d + 1;
                l_in_msg = in_msg / 2;
                l_ucap = ucap_l(l_d+1, l_in_msg*l_i+1:l_in_msg*(l_i+1));
                in_msg = in_msg / 2;
                d = d + 1;
                L(d, in_msg*i+1:in_msg*(i+1)) = g(r1, r2, l_ucap);
                state(p_state) = 2;
            else
                if state(p_state) == 0 %не лист
                    in_msg = 2^(n-d);
                    Ln = L(d, in_msg*i+1:in_msg*(i+1)); %входящие убеждения
                    r1 = Ln(1:in_msg/2);
                    r2 = Ln(in_msg/2+1:end);
                    i = 2*i; d = d + 1;
                    in_msg = in_msg / 2;
                    L(d, in_msg*i+1:in_msg*(i+1)) = f(r1,r2);
                    state(p_state) = 1;
                else
                    if state(p_state) == 1
                        in_msg = 2^(n-d);
                        Ln = L(d, in_msg*i+1:in_msg*(i+1)); %входящие убеждения
                        r1 = Ln(1:in_msg/2);
                        r2 = Ln(in_msg/2+1:end);
                        l_i = 2*i; l_d = d + 1;
                        l_in_msg = in_msg / 2;
                        l_ucap = ucap_r(l_d+1, l_in_msg*l_i+1:l_in_msg*(l_i+1));
                        in_msg = in_msg / 2;
                        i = 2*i+1; d = d + 1;
                        L(d, in_msg*i+1:in_msg*(i+1)) = g(r1, r2, l_ucap);
                        state(p_state) = 2;
                    else
                        in_msg = 2^(n-d);
                        l_i = 2*i; r_i = 2*i+1; c_d = d + 1;
                        c_in_msg = in_msg / 2;
                        l_ucap = ucap_r(c_d+1, c_in_msg*l_i+1:c_in_msg*(l_i+1));
                        r_ucap = ucap_r(c_d+1, c_in_msg*r_i+1:c_in_msg*(r_i+1));
                        ucap_r(d+1, in_msg*i+1:in_msg*(i+1)) = [mod(l_ucap+r_ucap,2) r_ucap];
                        i = floor(i/2); d = d-1;
                    end
                end
            end
        end
    end
    ucap = [ucap_l(n+1, 1:end/2) ucap_r(n+1, 1:end/2)];
    decoded_bit = ucap(1, Q1(N-K+1:end));
end

%-----------------------------------------------------------------------------------------------------
function f = f(r1,r2)
      f = sign(r1).*sign(r2).*min(abs(r1),abs(r2));
end

%-----------------------------------------------------------------------------------------------------
function g = g(r1, r2, b)
    g = r2 + (1-2*b).*r1;
end
