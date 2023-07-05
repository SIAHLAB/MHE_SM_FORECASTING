clc
close all
clear all

strain=5; %sampling time used for the weather station
stmoist=30; %sampling time used for the SM sensor

%% Read data from moisture and rain temperature CSV files
Tabmoist = readtable('moisttime.csv');
Tabrain= readtable('raintemp.csv');
%% Define the date and time range for the analysis
g2i = datetime('3-11-2022 15:00','InputFormat','MM-dd-yyyy HH:mm');
g2f = datetime('3-30-2022 13:00','InputFormat','MM-dd-yyyy HH:mm');
%% Data extracting and labeling
dhrain=Tabrain.dhrain; rain=Tabrain.Rain_In; tempf=(Tabrain.Temp__F); %data from time rain ad temperature station
dhmoist=Tabmoist.dhmoist; moisture_2=Tabmoist.moisture_2; %data from time and moisture sensor
prg2i=find(dhrain==g2i); prg2f=find(dhrain==g2f);
dhrg2=dhrain(prg2i:prg2f); raing2=rain(prg2i:prg2f); tempfg2=tempf(prg2i:prg2f); %original time, rain vector into selected time
t5sg2=(g2i:minutes(strain):g2f);%created time vector with 5m sampl (important)
%% repair data from weather station
raing25s=interp1(dhrg2,raing2,t5sg2)'; %interpolated rain, fix mised data
tempfg25s=(interp1(dhrg2,tempfg2,t5sg2)'); %interpolated temperature, fix mised data
%% repair data from SM sensor
pmg2i=find(dhmoist==g2i); pmg2f=find(dhmoist==g2f); 
dhmg2=dhmoist(pmg2i:pmg2f); %original moisture time vector
m2g2=moisture_2(pmg2i:pmg2f);
m2g25s=interp1(dhmg2,m2g2,t5sg2)'; %interpolated moisture data

%% Preparation of data for MHE estimation
ts=60*5; %used sampling time in seconds
dmoist2=(m2g25s(8:end)-m2g25s(7:end-1))/ts; %moisture variation
moistp1=m2g25s(8:end-1); %Previous moisture
moistp2=m2g25s(7:end-1); %Previous moisture
moistp3=m2g25s(6:end-2); %Previous moisture
tempp2=tempfg25s(7:end-1); %previous temperature
time=t5sg2(7:end-1);
rainp2=raing25s(7:end-1); %previous rain
rainpd1=raing25s(6:end-2); %previous rain
rainpd2=raing25s(5:end-3); %previous rain
rainpd3=raing25s(4:end-4); %previous rain
rainpd4=raing25s(3:end-5); %previous rain
rainpd5=raing25s(2:end-6); %previous rain
rainpd6=raing25s(1:end-7); %previous rain
kf=length(t5sg2);
kkw=7;
apri=34;
apos=3.5*10^11;
kw=kkw*24*60/5; %size of estimation window change days
kp= 3*24*60/5; %Size of predicition window
thet2g2ant=(1/300)*rand(10,1);
thetantf=(1/300)*rand;

for kk=1:1:kf-kw-kp + 1-7
    tic
    riwant=rainp2(kk:kw+kk-1); %rain estimation vector
    rainwd1=rainpd1(kk:kw+kk-1);
    rainwd2=rainpd2(kk:kw+kk-1);
    rainwd3=rainpd3(kk:kw+kk-1);
    rainwd4=rainpd4(kk:kw+kk-1);
    rainwd5=rainpd5(kk:kw+kk-1);
    rainwd6=rainpd6(kk:kw+kk-1);
    mpres=moistp1(kk:kw+kk-1); %m(k+1) moisture vector
    mwant=moistp2(kk:kw+kk-1); %m(k)moisture estimation vector
    mwantd=moistp3(kk:kw+kk-1); %m(k-1)
    dmoiswant=dmoist2(kk:kw+kk-1); %moisture variation estimation vector
    alphaprior=apri;
    alphapost=apos;
    Q=alphaprior*eye(10);
    tewa=tempp2(kk:kw+kk-1); %temperature estimation vector
    tewant=(tewa-mean(tewa))/std(tewa);
    timeest=time(kk:kw+kk-1);

    %% Priory estimation stage    
    phi2=[riwant rainwd1 rainwd2 rainwd3 rainwd4 rainwd5 rainwd6 -mwant -mwantd tewant];
    Aest=[0 0 0 0 0 0 0  -ts -ts 0;
          0 0 0 0 0 0 0 ts -ts 0;
          0 0 0 0 0 0 0  0  ts 0];
    Best=[0 0 1]';
    ppp=-0.000;
    thet2g2cell(:,kk)=thet2g2ant;
    rootscell(:,kk)=roots([1 (-1+ts*thet2g2ant(8)) ts*thet2g2ant(9)]); %ojo you need to review that the discrete polynimial is in function of the variation (1)
    thet2g2 = quadprog(phi2'*phi2+Q,-dmoiswant'*phi2-thet2g2ant'*Q,Aest,Best,[],[],[0 [0 0 0 0 0 0] 0 0 -100],[10 10 10 10 10 10 10 10 10 0],[],[]);
    thet2g2ant=thet2g2;
    %% Subsequent estimation stage
    const(kk)=thetantf; %theta_q
    options = optimoptions('fmincon','Algorithm','sqp');
    [thetf,fval,exitflag,output]=fmincon(@(thetf)  costfunctestff(mpres,thetf,mwant,mwantd,kw,thet2g2,riwant,rainwd1,rainwd2,rainwd3,rainwd4,rainwd5,rainwd6,tewant,ts,thetantf,alphapost),thetantf,[],[],[],[],0,0.1,[],options);
    thetantf=thetf;
    %% test estimation from kk to kk+kw
    moistest(1)=mwantd(1); %moisture initial condition
    moistest(2)=mwant(1); %moisture initial condition

    for k=2:1:kw-1
        moistest(k+1)=moistest(k) + ts*(thet2g2(1)*riwant(k) + thet2g2(2) *rainwd1(k) ...
        + thet2g2(3) *rainwd2(k) + thet2g2(4) *rainwd3(k) + thet2g2(5) *rainwd4(k) + thet2g2(6) *rainwd5(k)...
        + thet2g2(7) *rainwd6(k) - thet2g2(8)*moistest(k) - thet2g2(9)*moistest(k-1) + thet2g2(10)*tewant(k)...
        +thetf);
    end
    moistestcell(:,kk)=moistest';
    moistRest(:,kk)=mwantd; %Measured moisture along the predicition stage
    eest(:,kk)=mwantd-moistest';
    meanerrest(:,kk)=mean(mwantd-moistest');
    varerrest(:,kk)=var(mwantd-moistest');
    stderrest(:,kk)=std(mwantd-moistest');
    timestcell(:,kk)=timeest;
    %% Define vectors used in prediction process from kw + kk -1 to kw + kk + kp -2    
    riwpost=rainp2(kw + kk -1:kw + kk + kp -2); %rain prediction vector
    riwd1post=rainpd1(kw + kk -1:kw + kk + kp -2);
    riwd2post=rainpd2(kw + kk -1:kw + kk + kp -2);
    riwd3post=rainpd3(kw + kk -1:kw + kk + kp -2);
    riwd4post=rainpd4(kw + kk -1:kw + kk + kp -2);
    riwd5post=rainpd5(kw + kk -1:kw + kk + kp -2);
    riwd6post=rainpd6(kw + kk -1:kw + kk + kp -2);
    tempost=(tempp2(kw + kk -1:kw + kk + kp -2) - mean(tewa))/std(tewa); %temperature estimation vector

    %% test prediciton from kk+kw to kk+kw+kp
    moistpred(1)=moistp3(kk+kw-1); %moisture initial condition
    moistpred(2)=moistp2(kk+kw-1); %moisture initial condition
    for k=2:1:kp-1
        moistpred(k+1)=moistpred(k) + ts*(thet2g2(1)*riwpost(k) + thet2g2(2)*riwd1post(k) + thet2g2(3)*riwd2post(k)... 
        + thet2g2(4)*riwd3post(k)+ thet2g2(5)*riwd4post(k)+ thet2g2(6)*riwd5post(k)+ thet2g2(7)*riwd6post(k)...
        - thet2g2(8)*moistpred(k) - thet2g2(9)*moistpred(k-1) + thet2g2(10)*tempost(k) +thetf);
    end
    moistpredcell(:,kk)=moistpred'; %predicted moisture along the predicition stage
    moistRpred(:,kk)=moistp1(kk+kw-1:kk+kw+kp-2); %Measured moisture along the predicition stage
    epred(:,kk)=moistp1(kk+kw-1:kk+kw+kp-2)-moistpred'; %predicition error along the predicition stage
    meanerrpred(kk)=mean(moistp1(kk+kw-1:kk+kw+kp-2)-moistpred'); %mean of the predicition error
    varerrpred(kk)=var(moistp1(kk+kw-1:kk+kw+kp-2)-moistpred'); %variance of the predicition error
    stderrpred(kk)=std(moistp1(kk+kw-1:kk+kw+kp-2)-moistpred'); %standar deviation of the prediction error
    timepredcell(:,kk)=time(kk+kw-1:kk+kw+kp-2); %time vectos in the predicition
    TT(kk)=toc;  
end

tl=tiledlayout(4,1);
tl.TileSpacing='none';
tl.Padding='compact';

nexttile
plot(time,rainp2,'r')
legend('r (k)')
ylabel({'r(k) (mm)'},'Interpreter','latex')
xticks([])

nexttile
plot(time,dmoist2,'r')
legend({'$\Delta{m}(k)/ \Delta{t}$'},'Interpreter','latex')
ylabel({'$\Delta{m}/ \Delta{t} (\% /t)$'},'Interpreter','latex')
xticks([])

nexttile
for kk=1:30:kf-kw-kp + 1 - 6
h=plot(timepredcell(:,kk),moistpredcell(:,kk),time,moistp2);
set(h(1),'color','b','LineWidth',0.1)
set(h(2),'color','r','LineWidth',1)
hold on
legend({'$\boldmath{\hat{m}(k_p \mid k)}$','${m}(k)$'},'Interpreter','latex')
pause(0.1)
end
ylabel({'m(k) (\%)'},'Interpreter','latex')
hold off 
ylim([0 50])
xticks([])

nexttile
timemedinp=find(time==timepredcell(1,1));
timemefinp=find(time==timepredcell(1,end));

plot(time,zeros(1,length(time)),'k')
hold on
plot(time(timemedinp:timemefinp),meanerrpred+varerrpred,'Color',[0.8 0.8 1],'LineWidth',0.001)
hold on
plot(time(timemedinp:timemefinp),meanerrpred-varerrpred,'Color',[0.8 0.8 1],'LineWidth',0.001)
hold on
plot(time(timemedinp:timemefinp),meanerrpred,'b','LineWidth',2)
legend({'','${\overline{e}(k_p \mid k)} + \sigma$','${\overline{e}(k_p \mid k)} - \sigma$','${\overline{e}(k_p \mid k)}$'},'Interpreter','latex')
ylabel({'${\overline{e}(k_p \mid k)} (\%)$'},'Interpreter','latex')
ylim([-20 20])
