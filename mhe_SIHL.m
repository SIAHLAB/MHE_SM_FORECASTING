clc
clear all
close all

standardFormat = 'yyyy-MM-dd HH:mm:ss';

options = optimoptions('fmincon','Algorithm', 'sqp','MaxIterations',10000000,'MaxFunctionEvaluations',100000000,'ConstraintTolerance',1e-2,'OptimalityTolerance', 1e-5);

MAEAST0202= readtable('moisture_data.csv');
MAEAST0202.logged_at = datetime(MAEAST0202.logged_at, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ssXXX', 'TimeZone', 'UTC');
MAEAST0202.logged_at.Format = standardFormat;

WA0202 = readtable('weather_data.csv', 'HeaderLines', 5);
date_time_str = WA0202.Date_Time;
if ~iscell(date_time_str)
    date_time_str = cellstr(date_time_str);
end
date_time_str = cellfun(@(x) strrep(x, '/00', '/20'), date_time_str, 'UniformOutput', false);
WA0202.DateTime = datetime(date_time_str, 'InputFormat', 'MM/dd/yyyy hh:mm a', 'TimeZone', 'America/New_York');
WA0202.DateTime.Format = standardFormat;


% Specify the path to your Excel file
RIEGO = readtable('irrigation_data.xlsx');
RIEGO.Timestamp = datetime(RIEGO.Timestamp, 'InputFormat', 'yyyy-MM-dd HH:mm:ss Z', 'TimeZone', 'America/New_York');

ti = datetime('9-26-2023 15:00','InputFormat','MM-dd-yyyy HH:mm', 'TimeZone', 'America/New_York');
ti.Format = standardFormat;
tfi = datetime('11-9-2023 13:00','InputFormat','MM-dd-yyyy HH:mm', 'TimeZone', 'America/New_York');
tfi.Format = standardFormat;

tem_t=WA0202.DateTime; tem=WA0202.Temp__F; tem(isnan(tem))=0; tem=(tem-32)*(5/9); %temperature is in celcius
rain_t=WA0202.DateTime; rain=WA0202.Rain_In; 
%rain(isnan(rain))=(5/25.4);
moist_e_t=MAEAST0202.logged_at; moist_e_2=MAEAST0202.moisture_2;
riego_t = RIEGO.Timestamp; % Extracting the Timestamps from RIEGO table
riego_data = RIEGO.Sep_Fixed; % Sep_Fixed' is the data column

%% homogenious time
tau_s=15; t15mv=(ti:minutes(tau_s):tfi);%created time vector with 15m sampl
t15mv.TimeZone = 'America/New_York';

temp=interp1(tem_t,tem,t15mv);
rain=interp1(rain_t,rain,t15mv)*25.4; %milimeters
moist_e_2=interp1(moist_e_t,moist_e_2,t15mv);
riego_interp = interp1(riego_t, riego_data, t15mv)/4;

%% Create time vector acording irrigation period
tau_q=1;
td1samp=(ti+hours(tau_q):hours(tau_q):tfi); %time vector in steps of irrigation period
endpos=length(td1samp)+1; %this is the position of the last date in the time vector t15mv
delay=0;
%% initialconditions
beg=5; c1=0; c2=0; thetaant=[1 1 1 1]*0.00001; options = optimoptions('fmincon','Algorithm','sqp');
%begf=1; mr=8.77; thetwant=0.1;
begf=1; mr=9; thetwant=0.1;

%% Loop for aggregating data according to the irrigation period
ki = 2; % Starting position in td1samp for aggregation
k = ki; % Initialize k to start from the second position to ensure past data availability

% Initialize arrays to hold aggregated data for each irrigation period
%temps_d1 = zeros(length(td1samp)-1, 1); % Minus 1 because we start from position 2
%rains_d1 = zeros(length(td1samp)-1, 1);
%moist_e_2s_d1 = zeros(length(td1samp)-1, 1);

samp_rel=floor(hours(tau_q)/minutes(tau_s));

% Loop through each irrigation period defined in td1samp length(td1samp-2)
for k = 1:length(td1samp)
    % Find the start and end positions in t15mv for the current period
    startpos = samp_rel*(k-1) + 1;
    endpos = samp_rel*k + 1; % Use k+1 to get the end of the current period
    temps_d1(k,1) = mean(temp(startpos:endpos)); % Calculate mean temperature
    rains_d1(k,1) = sum(rain(startpos:endpos-1)); % Sum rainfall
    moist_e_2s_d1(k,1) = mean(moist_e_2(startpos:endpos)); % Calculate mean moisture
    riego_agg(k, 1) = sum(riego_interp(startpos:endpos-1));
 end



%% Preparation of data for MHE estimation
dco=10; %Maximum lenght of past components used
ts=tau_q; %used sampling time in hours
dmoist2=( moist_e_2s_d1(dco+1:end)-moist_e_2s_d1(dco:end-1) )/ts; %(m(k)-m(k-1))/ts
moistp1=moist_e_2s_d1(dco+1:end);   %m(k)
moistp2=moist_e_2s_d1(dco:end-1);   %m(k-1)
moistp3=moist_e_2s_d1(dco-1:end-2); %m(k-2)
tempp2= temps_d1(dco:end-1);        %c(k-1)
temppd1= temps_d1(dco-1:end-2);     %c(k-2)
temppd2= temps_d1(dco-2:end-3);     %c(k-3)
temppd3= temps_d1(dco-3:end-4);     %c(k-4)
temppd4= temps_d1(dco-4:end-5);     %c(k-5)
temppd5= temps_d1(dco-5:end-6);     %c(k-6)
temppd6= temps_d1(dco-6:end-7);     %c(k-7)
temppd7= temps_d1(dco-7:end-8);     %c(k-8)
temppd8= temps_d1(dco-8:end-9);     %c(k-9)
temppd9= temps_d1(dco-9:end-10);    %c(k-10)

time=   td1samp(dco+1:end);         %time(k)
rainp2= rains_d1(dco:end-1);        %r(k-1)
rainpd1=rains_d1(dco-1:end-2);      %r(k-2)
rainpd2=rains_d1(dco-2:end-3);      %r(k-3)
rainpd3=rains_d1(dco-3:end-4);      %r(k-4)
rainpd4=rains_d1(dco-4:end-5);      %r(k-5)
rainpd5=rains_d1(dco-5:end-6);      %r(k-6)
rainpd6=rains_d1(dco-6:end-7);      %r(k-7)
rainpd7=rains_d1(dco-7:end-8);      %r(k-8)
rainpd8=rains_d1(dco-8:end-9);      %r(k-9)
rainpd9=rains_d1(dco-9:end-10);     %r(k-10)

riegop2=riego_agg(dco:end-1);       %q(k-1)
riegopd1=riego_agg(dco-1:end-2);    %q(k-2)
riegopd2=riego_agg(dco-2:end-3);    %q(k-3)
riegopd3=riego_agg(dco-3:end-4);    %q(k-4)
riegopd4=riego_agg(dco-4:end-5);    %q(k-5)
riegopd5=riego_agg(dco-5:end-6);    %q(k-6)
riegopd6=riego_agg(dco-6:end-7);    %q(k-7)
riegopd7=riego_agg(dco-7:end-8);    %q(k-8)
riegopd8=riego_agg(dco-8:end-9);    %q(k-9)
riegopd9=riego_agg(dco-9:end-10);   %q(k-10)

%% contaminated prediction vectors
contper=0.2;

temps_d1_pred=temps_d1 + contper*std(temps_d1)*rand(size(temps_d1));

rain_chance = double(rains_d1 > 0); %this is a total certainty of the chance
rain_chance_cont = (rand(size(rains_d1)) > contper).*rain_chance + (rand(size(rains_d1)) <= contper).*(~rain_chance);
rain_noise=contper*std(rains_d1(rains_d1 > 0))*randn(size(rains_d1));
rain_noise_zero=rain_noise-mean(rain_noise);
rain_value_cont = rains_d1 + rain_noise_zero; %i contaminated with 10%
rain_cont= rain_value_cont.*rain_chance_cont;
rain_cont(rain_cont<0)=0;
rains_d1_pred=rain_cont;
%figure(1)
%stem(td1samp,rains_d1)
%hold on
%stem(td1samp,rains_d1_pred)


temppd1_pred= temps_d1_pred(dco-1:end-2);     %c(k-2)
temppd2_pred= temps_d1_pred(dco-2:end-3);     %c(k-3)
temppd3_pred= temps_d1_pred(dco-3:end-4);     %c(k-4)
temppd4_pred= temps_d1_pred(dco-4:end-5);     %c(k-5)
temppd5_pred= temps_d1_pred(dco-5:end-6);     %c(k-6)
temppd6_pred= temps_d1_pred(dco-6:end-7);     %c(k-7)
temppd7_pred= temps_d1_pred(dco-7:end-8);     %c(k-8)
temppd8_pred= temps_d1_pred(dco-8:end-9);     %c(k-9)
temppd9_pred= temps_d1_pred(dco-9:end-10);    %c(k-10)

rainp2_pred= rains_d1_pred(dco:end-1);        %r(k-1)
rainpd1_pred=rains_d1_pred(dco-1:end-2);      %r(k-2)
rainpd2_pred=rains_d1_pred(dco-2:end-3);      %r(k-3)
rainpd3_pred=rains_d1_pred(dco-3:end-4);      %r(k-4)
rainpd4_pred=rains_d1_pred(dco-4:end-5);      %r(k-5)
rainpd5_pred=rains_d1_pred(dco-5:end-6);      %r(k-6)
rainpd6_pred=rains_d1_pred(dco-6:end-7);      %r(k-7)
rainpd7_pred=rains_d1_pred(dco-7:end-8);      %r(k-8)
rainpd8_pred=rains_d1_pred(dco-8:end-9);      %r(k-9)
rainpd9_pred=rains_d1_pred(dco-9:end-10);     %r(k-10)

%%%%


kf=length(time);

esti_days=10; %number of days used for estimation
pred_days=3; %number of days used for predicition



kw= esti_days*24/tau_q; %size of estimation window
kp= pred_days*24/tau_q; %Size of predicition window

%thet2g2ant=0.0001*rand(24,1);
%thet2g2ant = 1.0e-03 * [0.4936, 0.0809, 0.5466, 0.2964, 0.2821, 0.4404, 0.3115, 0, 0, 0,...
                       % 0.0008, 0.0371,...  
                        %0.0142, 0.0015, 0.0016, 0.0244, 0.0168, 0.0030, 0.0029, 0, 0, 0,...
                        %0.8805,... 
                        %0.1723, 0.0455, 0.0185, 0.1063, 0.2763, 0.0082, 0.0085, 0, 0, 0,]';

thet2g2ant = rand * [0.0391, 0.0405, 0.0373, 0.0303, 0.0328, 0.0367, 0.0313, 0.0302, 0.0286, 0.0303,...
                        0.0012, 0.0110,...  
                        0.0007, 0.0006, 0.0004, 0.0003, 0.0003, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002,...
                        0.2790,... 
                        0.0073, 0.0070, 0.0069, 0.0067, 0.0064, 0.0062, 0.0061, 0.0062, 0.0063, 0.0066,]';

for kk=1:1:kf-kw-kp 
    tic
    %% Vectors for estimation
    riwant = rainp2(kk:kw+kk-1);    %r(k-1)
    rainwd1= rainpd1(kk:kw+kk-1);   %r(k-2)
    rainwd2= rainpd2(kk:kw+kk-1);   %r(k-3)
    rainwd3= rainpd3(kk:kw+kk-1);   %r(k-4)
    rainwd4= rainpd4(kk:kw+kk-1);   %r(k-5)
    rainwd5= rainpd5(kk:kw+kk-1);   %r(k-6)
    rainwd6= rainpd6(kk:kw+kk-1);   %r(k-7)
    rainwd7= rainpd7(kk:kw+kk-1);   %r(k-8)
    rainwd8= rainpd8(kk:kw+kk-1);   %r(k-9)
    rainwd9= rainpd9(kk:kw+kk-1);   %r(k-10)

    mpres=moistp1(kk:kw+kk-1);      %m(k) moisture vector
    mwant=moistp2(kk:kw+kk-1);      %m(k-1)moisture estimation vector
    mwantd=moistp3(kk:kw+kk-1);     %m(k-2)

    dmoiswant=dmoist2(kk:kw+kk-1);  %(m(k)-m(k-1))/ts
    
    tempwant=tempp2(kk:kw+kk-1);    %c(k-1)
    tempwd1 = temppd1(kk:kw+kk-1);  % c(k-2)
    tempwd2 = temppd2(kk:kw+kk-1);  % c(k-3)
    tempwd3 = temppd3(kk:kw+kk-1);  % c(k-4)
    tempwd4 = temppd4(kk:kw+kk-1);  % c(k-5)
    tempwd5 = temppd5(kk:kw+kk-1);  % c(k-6)
    tempwd6 = temppd6(kk:kw+kk-1);  % c(k-7)
    tempwd7 = temppd7(kk:kw+kk-1);  % c(k-8)
    tempwd8 = temppd8(kk:kw+kk-1);  % c(k-9)
    tempwd9 = temppd9(kk:kw+kk-1);  % c(k-10)

    riegowant=riegop2(kk:kw+kk-1);     %q(k-1)
    riegowd1 = riegopd1(kk:kw+kk-1);   % q(k-2)
    riegowd2 = riegopd2(kk:kw+kk-1);   % q(k-3)
    riegowd3 = riegopd3(kk:kw+kk-1);   % q(k-4)
    riegowd4 = riegopd4(kk:kw+kk-1);   % q(k-5)
    riegowd5 = riegopd5(kk:kw+kk-1);   % q(k-6)
    riegowd6 = riegopd6(kk:kw+kk-1);   % q(k-7)
    riegowd7 = riegopd7(kk:kw+kk-1);   % q(k-8)
    riegowd8 = riegopd8(kk:kw+kk-1);   % q(k-9)
    riegowd9 = riegopd9(kk:kw+kk-1);   % q(k-10)

    timeest=time(kk:kw+kk-1);       %time(k)
    %% Estimation
    
    rainvectors=[riwant rainwd1 rainwd2 rainwd3 rainwd4 rainwd5 rainwd6 rainwd7 rainwd8 rainwd9];
    moistvectors=[mwant mwantd];
    tempvectors=[tempwant tempwd1 tempwd2 tempwd3 tempwd4 tempwd5 tempwd6 tempwd7 tempwd8 tempwd9];
    riegovectors=[riegowant riegowd1 riegowd2 riegowd3 riegowd4 riegowd5 riegowd6 riegowd7 riegowd8 riegowd9];

    
    colrain=size(rainvectors,2);
    colmoist=size(moistvectors,2);
    coltemp=size(tempvectors,2);
    colriego=size(riegovectors,2);

    %alprain=2000*ones(1,colrain);
    %alpmoist=200000*ones(1,colmoist);
    %alptem=100*ones(1,coltemp);
    %alpinf=2000;
    %alprieg=100*ones(1,colriego);

    diagrain=(1:1:colrain);
    diagmoist=([1 2]);
    diatem=(1:1:coltemp);
    diaginf=1;
    diagrieg=(1:1:colriego);
        
    %Q=diag([alprain alpmoist alptem alpinf alprieg]);

    Q2=1*diag([diagrain diagmoist diatem diaginf diagrieg]);
    phi2 = [rainvectors -moistvectors -tempvectors ones(size(riwant)) riegovectors];
    Q=100*(phi2'*phi2);
    thet2g2cell(:,kk) = thet2g2ant;
    rootscell(:,kk) = roots([1 (-1+ts*thet2g2ant(11)) ts*thet2g2ant(12)]); 
    H = phi2'*phi2+Q+Q2;
    f = -dmoiswant'*phi2-thet2g2ant'*Q;
    %constraint stability
    Aest=[zeros(1,colrain) -ts -ts  zeros(1,coltemp) 0 zeros(1,colriego);
          zeros(1,colrain)  ts -ts  zeros(1,coltemp) 0 zeros(1,colriego);
          zeros(1,colrain)  0   ts  zeros(1,coltemp) 0 zeros(1,colriego);
          zeros(1,colrain)  0   -ts  zeros(1,coltemp) 0 zeros(1,colriego)];
    Best=[1 1 1 1]';
    Aeq = [];
    beq = [];
    lb = 1e-10 * ones(size(phi2, 2), 1); % Lower bound, slightly above zero to enforce theta > 0
    ub = []; %0.001 * ones(size(phi2, 2), 1); % No upper bound

    % Solve the quadratic programming problem
    thet2g2 = quadprog(H, f, Aest, Best, Aeq, beq, lb, ub);
    thet2g2ant = thet2g2;
    %% test estimation from kk to kk+kw
    moistestl(1,1)=mwantd(1); %moisture initial condition
    moistestl(2,1)=mwant(1); %moisture initial condition (you can estimat or use the same measurement)

    for k = 1:kw-1
            moistestl(k+2) = moistestl(k+1) + ts * (...
                  [riwant(k+1) rainwd1(k+1) rainwd2(k+1) rainwd3(k+1) rainwd4(k+1) rainwd5(k+1) rainwd6(k+1) rainwd7(k+1) rainwd8(k+1) rainwd9(k+1)] * thet2g2(1:10) ...
                - [moistestl(k+1) moistestl(k)] * thet2g2(11:12) ...
                - [tempwant(k+1) tempwd1(k+1) tempwd2(k+1) tempwd3(k+1) tempwd4(k+1) tempwd5(k+1) tempwd6(k+1) tempwd7(k+1) tempwd8(k+1) tempwd9(k+1)] * thet2g2(13:22) ...
                + thet2g2(23) ... % Constant term
                + [riegowant(k+1) riegowd1(k+1) riegowd2(k+1) riegowd3(k+1) riegowd4(k+1) riegowd5(k+1) riegowd6(k+1) riegowd7(k+1) riegowd8(k+1) riegowd9(k+1)] * thet2g2(24:33));
    end
        
    moistest=moistestl(2:end);
    moistestcell(:,kk)=moistest;
    moistRest(:,kk)=mwantd; %Measured moisture along the predicition stage
    eest(:,kk)=mwantd-moistest;
    meanerrest(:,kk)=mean(mwantd-moistest);
    varerrest(:,kk)=var(mwantd-moistest);
    stderrest(:,kk)=std(mwantd-moistest);
    timestcell(:,kk)=timeest;
    
    %% Define vectors used in prediction process from kw + kk -1 to kw + kk + kp -2    
    riwpost=rainp2_pred(kw + kk -1:kw + kk + kp -2); %rain prediction vector
    riwd1post=rainpd1_pred(kw + kk -1:kw + kk + kp -2);
    riwd2post=rainpd2_pred(kw + kk -1:kw + kk + kp -2);
    riwd3post=rainpd3_pred(kw + kk -1:kw + kk + kp -2);
    riwd4post=rainpd4_pred(kw + kk -1:kw + kk + kp -2);
    riwd5post=rainpd5_pred(kw + kk -1:kw + kk + kp -2);
    riwd6post=rainpd6_pred(kw + kk -1:kw + kk + kp -2);
    riwd7post=rainpd7_pred(kw + kk -1:kw + kk + kp -2);
    riwd8post=rainpd8_pred(kw + kk -1:kw + kk + kp -2);
    riwd9post=rainpd9_pred(kw + kk -1:kw + kk + kp -2);
    %tempost=(tempp2(kw + kk -1:kw + kk + kp -2) - mean(tewa))/std(tewa); %temperature estimation vector
    temppost = temppd2_pred(kw + kk - 1:kw + kk + kp - 2); % Direct temperature prediction vector (c(k-1))
    temppd1post = temppd1_pred(kw + kk - 1:kw + kk + kp - 2); % Temperature prediction vector (c(k-2))
    temppd2post = temppd2_pred(kw + kk - 1:kw + kk + kp - 2); % Temperature prediction vector (c(k-3))
    temppd3post = temppd3_pred(kw + kk - 1:kw + kk + kp - 2); % Temperature prediction vector (c(k-4))
    temppd4post = temppd4_pred(kw + kk - 1:kw + kk + kp - 2); % Temperature prediction vector (c(k-5))
    temppd5post = temppd5_pred(kw + kk - 1:kw + kk + kp - 2); % Temperature prediction vector (c(k-6))
    temppd6post = temppd6_pred(kw + kk - 1:kw + kk + kp - 2); % Temperature prediction vector (c(k-7))
    temppd7post = temppd7_pred(kw + kk - 1:kw + kk + kp - 2); % Temperature prediction vector (c(k-8))
    temppd8post = temppd8_pred(kw + kk - 1:kw + kk + kp - 2); % Temperature prediction vector (c(k-9))
    temppd9post = temppd9_pred(kw + kk - 1:kw + kk + kp - 2); % Temperature prediction vector (c(k-10))
    
    riegopost = riegop2(kw + kk - 1:kw + kk + kp - 2); % Direct irrigation prediction vector (q(k-1))
    riegopd1post = riegopd1(kw + kk - 1:kw + kk + kp - 2); % Irrigation prediction vector (q(k-2))
    riegopd2post = riegopd2(kw + kk - 1:kw + kk + kp - 2); % Irrigation prediction vector (q(k-3))
    riegopd3post = riegopd3(kw + kk - 1:kw + kk + kp - 2); % Irrigation prediction vector (q(k-4))
    riegopd4post = riegopd4(kw + kk - 1:kw + kk + kp - 2); % Irrigation prediction vector (q(k-5))
    riegopd5post = riegopd5(kw + kk - 1:kw + kk + kp - 2); % Irrigation prediction vector (q(k-6))
    riegopd6post = riegopd6(kw + kk - 1:kw + kk + kp - 2); % Irrigation prediction vector (q(k-7))
    riegopd7post = riegopd7(kw + kk - 1:kw + kk + kp - 2); % Irrigation prediction vector (q(k-8))
    riegopd8post = riegopd8(kw + kk - 1:kw + kk + kp - 2); % Irrigation prediction vector (q(k-9))
    riegopd9post = riegopd9(kw + kk - 1:kw + kk + kp - 2); % Irrigation prediction vector (q(k-10))

    %% test prediciton from kk+kw to kk+kw+kp
    moistpredl(1,1)=mpres(end-1); %moisture initial condition
    moistpredl(2,1)=mpres(end); %moisture initial condition
 
    for k=1:1:kp-1
        moistpredl(k+2,1)=moistpredl(k+1) + ts*( ...
                  [riwpost(k+1) riwd1post(k+1) riwd2post(k+1) riwd3post(k+1) riwd4post(k+1) riwd5post(k+1) riwd6post(k+1) riwd7post(k+1) riwd8post(k+1) riwd9post(k+1)] * thet2g2(1:10) ...
                - [moistpredl(k+1) moistpredl(k)] * thet2g2(11:12) ...
                - [temppost(k+1) temppd1post(k+1) temppd2post(k+1) temppd3post(k+1) temppd4post(k+1) temppd5post(k+1) temppd6post(k+1) temppd7post(k+1) temppd8post(k+1) temppd9post(k+1)] * thet2g2(13:22)...
                + thet2g2(23) ...
                + [riegopost(k+1) riegopd1post(k+1) riegopd2post(k+1) riegopd3post(k+1) riegopd4post(k+1) riegopd5post(k+1) riegopd6post(k+1) riegopd7post(k+1) riegopd8post(k+1) riegopd9post(k+1)] * thet2g2(24:33));
    end
    moistpred=moistpredl(2:end);
    moistpredcell(:,kk)=moistpred; %predicted moisture along the predicition stage
    moistRpred(:,kk)=moistp1(kk+kw-1:kk+kw+kp-2); %Measured moisture along the predicition stage
    epred(:,kk)=moistp1(kk+kw-1:kk+kw+kp-2)-moistpred; %predicition error along the predicition stage
    meanerrpred(kk)=mean(moistp1(kk+kw-1:kk+kw+kp-2)-moistpred); %mean of the predicition error
    varerrpred(kk)=var(moistp1(kk+kw-1:kk+kw+kp-2)-moistpred); %variance of the predicition error
    stderrpred(kk)=std(moistp1(kk+kw-1:kk+kw+kp-2)-moistpred); %standar deviation of the prediction error
    timepredcell(:,kk)=time(kk+kw-1:kk+kw+kp-2); %time vectos in the predicition
    TT(kk)=toc;  
end
%%% %

nitsekpi= ((meanerrpred(1:1:712).^2)*(1:1:712)')/712;
nisekpi= sum(abs(meanerrpred(1:1:712)))/712;

figure(1)
tl=tiledlayout(6,1);
tl.TileSpacing='none';
tl.Padding='compact';

nexttile
n=stem(td1samp,rains_d1_pred,'b--');
n.Marker = '.'; % Define marker shape
n.MarkerFaceColor = 'b'; % Color of the marker or the dot at each stem
n.LineWidth = 1;
hold on
h=stem(td1samp,rains_d1,'r');
h.Marker = '.'; % Define marker shape
h.MarkerFaceColor = 'r'; % Color of the marker or the dot at each stem
h.LineWidth = 2;
legend('Synthetic $\hat{r} (k)$','r (k)','Interpreter','latex', 'Orientation', 'horizontal')
ylabel({'r(k) (mm)'},'Interpreter','latex')
xticks([])
xlim([time(1) time(end)])

nexttile;
yyaxis right; % Activate the right y-axis
m = stem(time, riegop2, 'r', 'filled');
m.Marker = '.'; % Define marker shape
m.MarkerFaceColor = 'r'; % Color of the marker or the dot at each stem
m.LineWidth = 0.1;
legend('q (k)', 'Interpreter', 'latex', 'Orientation', 'horizontal');
ylabel({'q(k) (hours)'}, 'Interpreter', 'latex');
ax = gca; % Get current axis
ax.YAxis(2).Color = 'k'; % Set y-axis label and ticks color to black
ylim([0 1.1]);
xticks([]); % Remove x-axis ticks
yyaxis left; % Switch back to the left y-axis if needed
yticks([]); % Remove left y-axis ticks, if this is intended
xlim([time(1) time(end)]); % Set the x-axis limits



nexttile;
plot(td1samp, temps_d1_pred, 'b--', 'LineWidth', 3); % Plot the first dataset with specified properties
hold on; % Retain the current plot when adding new plots
plot(td1samp, temps_d1, 'r', 'LineWidth', 2); % Plot the second dataset with specified properties
hold off; % No longer retain the plot
legend('Synthetic $\hat{c} (k)$', 'c (k)', 'Interpreter', 'latex', 'Orientation', 'horizontal');
ylabel({'c(k) (C)'}, 'Interpreter', 'latex');
xticks([]);
xlim([time(1) time(end)]);



nexttile;
yyaxis right;
% Specify 'LineWidth' for the plot
plot(time, dmoist2, 'r', 'LineWidth', 2);
legend({'$\Delta{m}(k)/ \Delta{t}$'}, 'Interpreter', 'latex', 'Orientation', 'horizontal');
ylabel({'$\Delta{m}/ \Delta{t} (\% /t)$'}, 'Interpreter', 'latex');
ax = gca; % Get current axis
ax.YAxis(2).Color = 'k'; % Set y-axis label and ticks color to black
xticks([]);
yyaxis left;
yticks([]);
xlim([time(1) time(end)]);


nexttile;
hold on;
hActual = plot(time, moistp2, 'r', 'LineWidth', 1);
for kk = 1:1:kf-kw-kp + 1 - 6
    % Plot predictions in each loop iteration
    h = plot(timepredcell(:,kk), moistpredcell(:,kk), 'b', 'LineWidth', 0.1);
    % Optional: Update the plot handle for the actual data if it needs to be re-plotted in each iteration
    % hActual = plot(time, moistp2, 'r', 'LineWidth', 1); 
    pause(0.1); % Adjust as needed for visualization purposes
end
hActual = plot(time, moistp2, 'r', 'LineWidth', 2);
% Setting the legend once outside the loop if the legend entries do not change
legend([h, hActual], {'$\boldmath{\hat{m}(k_p \mid k)}$', '${m}(k)$'}, 'Interpreter', 'latex', 'Orientation', 'horizontal');

ylabel('m(k) (%)', 'Interpreter', 'latex');
hold off;
ylim([0, 40]);
xticks([]);
xlim([time(1), time(end)]);


nexttile
yyaxis right
timemedinp=find(time==timepredcell(1,1));
timemefinp=find(time==timepredcell(1,end));
plot(time,zeros(1,length(time)),'k')
hold on
plot(time(timemedinp:timemefinp)',(meanerrpred+varerrpred)','-','Color',[0.8 0.8 1],'LineWidth',0.001)
hold on
plot(time(timemedinp:timemefinp)',(meanerrpred-varerrpred)','-','Color',[0.8 0.8 1],'LineWidth',0.001)
hold on
plot(time(timemedinp:timemefinp)',meanerrpred,'b-','LineWidth',2)
legend({'','${\overline{e}(k_p \mid k)} + \sigma$','${\overline{e}(k_p \mid k)} - \sigma$','${\overline{e}(k_p \mid k)}$'},'Interpreter','latex', 'Orientation', 'horizontal')
ylabel({'${\overline{e}(k_p \mid k)} (\%)$'},'Interpreter','latex')
ax = gca; % Get current axis
ax.YAxis(2).Color = 'k'; % Set y-axis label and ticks color to black
ylim([-10 10 ])
yyaxis left
yticks([])
xlim([time(1) time(end)])





figure(2)
tl = tiledlayout(5, 1); % Setup for 5 subplots vertically aligned
tl.TileSpacing = 'none'; % No spacing between tiles
tl.Padding = 'compact'; % Minimal padding around the tiles

nexttile;
plot(timepredcell(1,:), thet2g2cell(1:10, :));
legend('$$\theta_{r_1}(k)$$', '$$\theta_{r_2}(k)$$', '$$\theta_{r_3}(k)$$', '$$\theta_{r_4}(k)$$', '$$\theta_{r_5}(k)$$', '$$\theta_{r_6}(k)$$', '$$\theta_{r_7}(k)$$', '$$\theta_{r_8}(k)$$', '$$\theta_{r_9}(k)$$', '$$\theta_{r_10}(k)$$', 'Interpreter', 'latex', 'Location', 'northwest', 'Orientation', 'horizontal');
ylabel('$$\theta_{r}(k)$$', 'Interpreter', 'latex');
ylim([0 max(max(thet2g2cell(1:10,:)))]);
xticks([]);
xlim([timepredcell(1,1) timepredcell(1,end)]);


% Plot for Tile 2: \theta_{m}(k), y-axis on the right
nexttile;
yyaxis right; % Activate the right y-axis
plot(timepredcell(1,:), thet2g2cell(11:12, :));
legend('$$\theta_{m_1}(k)$$', '$$\theta_{m_2}(k)$$', 'Interpreter', 'latex', 'Location', 'northwest', 'Orientation', 'horizontal');
ylabel('$$\theta_{m}(k)$$', 'Interpreter', 'latex');
ax = gca; % Get current axis
ax.YAxis(2).Color = 'k'; % Set y-axis label and ticks color to black
ylim([0 max(max(thet2g2cell(11:12,:)))]);
xticks([]); % Remove x-axis ticks
yyaxis left
yticks([])
xlim([timepredcell(1,1) timepredcell(1,end)]);


% Plot for Tile 3: \theta_{c}(k), y-axis on the left
nexttile;
plot(timepredcell(1,:), thet2g2cell(13:22,:));
ylabel('$$\theta_{c}(k)$$', 'Interpreter', 'latex');
legend('$$\theta_{c_1}(k)$$', '$$\theta_{c_2}(k)$$', '$$\theta_{c_3}(k)$$', '$$\theta_{c_4}(k)$$', '$$\theta_{c_5}(k)$$', '$$\theta_{c_6}(k)$$', '$$\theta_{c_7}(k)$$', '$$\theta_{c_8}(k)$$', '$$\theta_{c_9}(k)$$', '$$\theta_{c_10}(k)$$', 'Interpreter', 'latex', 'Location', 'northwest', 'Orientation', 'horizontal');
ylim([0 max(max(thet2g2cell(13:22,:)))]);
xticks([]);
xlim([timepredcell(1,1) timepredcell(1,end)]);


% Plot for Tile 4: \theta_{w}(k), y-axis on the right
nexttile;
yyaxis right;
plot(timepredcell(1,:), thet2g2cell(23,:));
legend('$$\theta_{w}(k)$$', 'Interpreter', 'latex', 'Location', 'northwest', 'Orientation', 'vertical');
ylabel('$$\theta_{w}(k)$$', 'Interpreter', 'latex');
ax = gca; % Get current axis
ax.YAxis(2).Color = 'k'; % Set y-axis label and ticks color to black
ylim([0 max(max(thet2g2cell(23,:)))]);
xticks([]); % Remove x-axis ticks
yyaxis left
yticks([])
xlim([timepredcell(1,1) timepredcell(1,end)]);

% Plot for Tile 5: \theta_{q}(k), y-axis on the left
nexttile;
plot(timepredcell(1,:), thet2g2cell(24:33,:));
ylabel('$$\theta_{q}(k)$$', 'Interpreter', 'latex');
legend('\(\theta_{q_1}(k)\)', '\(\theta_{q_2}(k)\)', '\(\theta_{q_3}(k)\)', '\(\theta_{q_4}(k)\)', '\(\theta_{q_5}(k)\)', '\(\theta_{q_6}(k)\)', '\(\theta_{q_7}(k)\)', '\(\theta_{q_8}(k)\)', '\(\theta_{q_9}(k)\)', '\(\theta_{q_10}(k)\)', 'Interpreter', 'latex', 'Location', 'northwest', 'Orientation', 'horizontal');
ylim([0 max(max(thet2g2cell(24:33,:)))]);
xlim([timepredcell(1,1) timepredcell(1,end)]);

figure(3)

% Define your poles
poles = rootscell;

% Get real and imaginary parts
real_poles = real(poles(:));
imag_poles = imag(poles(:));

% Create a 2D histogram (density of poles)
edges_real = linspace(-1.5, 1.5, 10);  % Increased resolution
edges_imag = linspace(-1.5, 1.5, 10);
[H, X, Y] = histcounts2(real_poles, imag_poles, edges_real, edges_imag);

% Define bin centers instead of edges
X_center = X(2:end) - diff(X)/2;
Y_center = Y(1:end-1) + diff(Y)/2;

% Plot the density of poles using contour plot
figure;
contourf(X_center, Y_center, H', 'LineColor', 'none'); % use X_center and Y_center
colorbar;
hold on;
plot(real_poles,imag_poles,'*','Markersize',8);  % plot poles over the contour plot using '.' marker

% Add a unit circle with its center at the origin
theta = 0:0.01:2*pi;
x = cos(theta);
y = sin(theta);
plot(x, y, 'Color', 'k', 'LineStyle', '--'); % Plot the unit circle in red

% Add a cross to represent the axis within the unit circle
line([-1, 1], [0, 0], 'Color', 'k', 'LineStyle', '--'); % x-axis
line([0, 0], [-1, 1], 'Color', 'k', 'LineStyle', '--'); % y-axis

hold off;

% Set the axes limits to match the edges of the histogram
%xlim([min(X_center) max(X_center)]);
%ylim([min(Y_center) max(Y_center)]);

% Add labels for real and complex axis
xlabel('Real');
ylabel('Imaginary');
