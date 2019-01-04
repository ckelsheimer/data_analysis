ADD_TO_STARTUP;
cd('~/Box/CJ Kelsheimer/PEGDA_HGS')
load('kin.mat') %remove if different source of kinetics data
options.t2_array = 10; %units of ps! 
options.damping = '1exp';
options.pnames = {'Delta1 (cm-1)','tau1 (ps)','anh (cm-1)'};

options.dt = 0.10; %ps
options.n_t = 64;
options.order = 3;

Delta1_cm  = 4; %cm-1
tau1 = 20; %ps
anh_cm =25; %doesn't matter for linear
mean_w_0 = 2340;
w1 = [2200:2:2500];
w3 = w1; %not used
T2 = 3;

p0 = [Delta1_cm tau1 anh_cm]; %[delta1_cm delta2_cm delta3_cm delta4_cm tau1_ps tau2_ps tau3_ps anh_cm]

% not a fitting parameter right now
options.w_0_cm = mean_w_0;
options.p0 = p0; %save where we started
options.w_nu2_cm = -12; %put this somewhere better later - shift to 

[S,extra] = analyticalResponseFunctionsFun_cjk(p0,w1,w3,options);
t = extra.t;
% S1 = extra.S1;
g = extra.g;
c2 = extra.c2;

figure(1),
%plot(w1,S)
figure(1),clf,my2dPlot(w1,w3,S)

options.order = 1;
S = analyticalResponseFunctionsFun(p0,w1,w3,options);

figure(2),
subplot(3,1,1)
plot(t,c2)
ylabel('c_2(t)') 
subplot(3,1,2)
plot(t,g)
ylabel('g(t)')
subplot(3,1,3)
plot(w1,S)
ylabel('S(\omega)')


Delta1 = Delta1_cm*wavenumbersToInvPs*2*pi;
Lambda1 = 1/tau1;
ratio = Delta1/Lambda1;
disp('ratio Delta1/Lambda1')
disp(ratio)
%%

%%
cd('~/Box/CJ_Kelsheimer/PEGDA_HGS')
%cd('/home/sgrlab/Documents/CJK/PEGDA_HGS/')
load('PEGDA_2D_Expt.mat')
dataMatrix = prepareGlobalFitData(PEGDA_2d_crop);

%% 1exp1fast
ADD_TO_STARTUP;
%cd('~/Box/CJ_Kelsheimer/PEGDA_HGS') %MAC
%cd /home/sgrlab/Documents/CJK/PEGDA_HGS/ %LINUX
load('kin.mat') %remove if different (or no) source of  kinetics data 
options.t2_array = [PEGDA_2d_crop.t2]/1000; %units of ps! %put entire t2_array CHECK UNITS
options.damping = '1exp1fast'; %see caRFF to select damping
options.pnames = {'Delta1 (cm-1)','tau1 (ps)','T2 (ps)','anh (cm-1)','phi (deg)','mu12_2','k1','k2','k3'};

options.dt = 0.20; % fft stuff
options.n_t = 64;  % fft stuff
options.order = 3; %3 - 2D, 1 - FTIR
options.kin = kin; %comment if no kinetics data 

%data updated 1.3.19
Delta1_cm  = 1.75; %cm-1 
tau1 = 110; %ps
anh_cm =23.5; %doesn't matter for linear
mean_w_0 = 2337.5;
% fix calibration of experimental spectrum

w1 = PEGDA_2d_crop.w1;
w3 = PEGDA_2d_crop.w3; %not used
T2 = 7.6;
mu12_2 = 1.51;
phi_deg = 0;

dE = 667;%cm^1
%assumption: kT in wavenumbers at 298K = 207cm^-1
k_B_cm = 0.69503476;
RT = 25;
T = 273+RT;
Keq = 2*exp(-dE/(T*k_B_cm));

k1 = .0028;%0.05
k2 =  0.00156; %if this goes above 1 intensities switch
k3 = k2/Keq;


p0 = [Delta1_cm tau1 T2 anh_cm phi_deg mu12_2 k1 k2 k3]; %[delta1_cm delta2_cm delta3_cm delta4_cm tau1_ps tau2_ps tau3_ps anh_cm]

% not a fitting parameter right now
options.w_0_cm = mean_w_0;

options.p0 = p0; %save where we started
options.w_nu2_cm = -12; %put tfprintf(fid,'\n');his somewhere better later - shift to bend/strech diag



[S,extra] = analyticalResponseFunctionsFun_cjk(p0,w1,w3,options);
t = extra.t;
% S1 = extra.S1;
g = extra.g;
c2 = extra.c2;


for ii = 1:size(S,3)
    for jj = ii
    if(ii>1),pause,end
    n1 = abs(min(min(dataMatrix(:,:,jj))));
    n2 = abs(min(min(S(:,:,ii))));
    
figure(1),
    clf,my2dPlot(w1,w3,S(:,:,ii)./n2,'n_contours',20)

figure(2)
    my2dPlot(w1,w3,dataMatrix(:,:,jj)./n1,'n_contours',20)

figure(3)
    my2dPlot(w1,w3,dataMatrix(:,:,jj)./n1-S(:,:,ii)./n2,'n_contours',20);
    end
    
end

%options.order = 1;
%S = analyticalResponseFunctionsFun(p0,w1,w3,options);

%Delta1 = Delta1_cm*wavenumbersToInvPs*2*pi;
%Lambda1 = 1/tau1;
%ratio = Delta1/Lambda1;
%disp('ratio Delta1/Lambda1')
%disp(ratio) 1

%% 2exp1fast

%%
ADD_TO_STARTUP;
cd('~/Box/CJ_Kelsheimer/PEGDA_HGS')
load('kin.mat') %remove if different source of kinetics data
options.t2_array = [PEGDA_2d_crop(1).t2/1000]; %units of ps! 
options.damping = '2exp1fast';
% options.pnames = {'Delta1 (cm-1)','tau1 (ps)','anh (cm-1)'};
options.pnames = {'Delta1 (cm-1)','tau1 (ps)','Delta2 (cm-1)','tau2(ps)','T2 (ps)','anh (cm-1)','k1','k2','k3'};

options.dt = 0.20; %ps
options.n_t = 64;
options.order = 3;
%options.kin = kin;


anh_cm =25.5; %doesn't matter for linear
mean_w_0 = 2338;
% fix calibration of experimental spectrum

Delta1_cm  = 0.7; %cm-1
tau1 = 5; %ps
Delta2_cm = 0.5;
tau2 = 150;

w1 = PEGDA_2d_crop.w1;
w3 = PEGDA_2d_crop.w3; %not used
T2 = 20;

dE = 667;%cm^1
%assumption: kT in wavenumbers at 298K = 207cm^-1
Keq = 2*exp(-dE/207);

k1 = .005;%0.05
k2 =  0.008; %if this goes above 1 intensities switch
k3 = k2/Keq;

p0 = [Delta1_cm tau1 Delta2_cm tau2 T2 anh_cm k1 k2 k3]; %[delta1_cm delta2_cm delta3_cm delta4_cm tau1_ps tau2_ps tau3_ps anh_cm]

% not a fitting parameter right now
options.w_0_cm = mean_w_0;
options.p0 = p0; %save where we started
options.w_nu2_cm = -12; %put this somewhere better later - shift to bend/strech diag



[S,extra] = analyticalResponseFunctionsFun_cjk(p0,w1,w3,options);
t = extra.t;
% S1 = extra.S1;
g = extra.g;
c2 = extra.c2;

for ii = 1%:size(S,3)
    n1 = abs(min(min(dataMatrix(:,:,ii))));
    n2 = abs(min(min(S(:,:,ii))));
    
    figure(1),
    clf,my2dPlot(w1,w3,S(:,:,ii)./n2,'zlimit',1,'n_contours',20)


figure(2)
    my2dPlot(w1,w3,dataMatrix(:,:,ii)./n1,'zlimit',1,'n_contours',20)

    figure(3)
    my2dPlot(w1,w3,dataMatrix(:,:,ii)./n1-S(:,:,ii)./n2,'n_contours',20);
end

% options.order = 1;
% S = analyticalResponseFunctionsFun(p0,w1,w3,options);

%Delta1 = Delta1_cm*wavenumbersToInvPs*2*pi;
%Lambda1 = 1/tau1;
%ratio = Delta1/Lambda1;
%disp('ratio Delta1/Lambda1')
%disp(ratio)

%%
%p0 = [Delta1_cm tau1 Delta2_cm tau2 T2 anh_cm k1 k2 k3] 2exp1fast
% ub = [1  10 1 250  5 30 1 1 5];
% lb = [0 0.5 0 10  0.1 20 0.0001 0.0001 0.0001];
% added mu12_2
% ub = [10 500  15 30 1.66 1 1 5];
% lb = [0 10  0.1 20 1.4 0.0001 0.0001 0.0001];
% added phi_deg
ub = [10 500  15 30 40 1.6 1 1 5];
lb = [0 10  0.1 20 -40 1.4 0.0001 0.0001 0.0001];

[pfit,err,gfstruct] = globalFit_cjk(dataMatrix,w1,w3,p0,options,lb, ub);
log_gfit;
%%

[S2,extra2] = analyticalResponseFunctionsFun_cjk(pfit,w1,w3,options);
t2 = extra.t;

% testing phase real fast
% p_uh =pfit;
% p_uh(5) = -20;
% options.t2 = 0.2;
% [S2,extra2] = analyticalResponseFunctionsFun_cjk(p_uh,w1,w3,options);
% t2 = extra.t;

%%
for ii = 1:size(S2,3)
    for jj = ii
    if(ii>1),pause,end
    n1 = abs(min(min(dataMatrix(:,:,jj))));
    n2 = abs(min(min(S(:,:,ii))));
    
figure(1),
    clf,my2dPlot(w1,w3,S2(:,:,ii)./n2,'n_contours',20)
drawnow
figure(2)
    my2dPlot(w1,w3,dataMatrix(:,:,jj)./n1,'n_contours',20)

figure(3)
    my2dPlot(w1,w3,dataMatrix(:,:,jj)./n1-S(:,:,ii)./n2,'n_contours',20);
    end
    
end
% 
% for ii = 1:size(S2,3)
%     if(ii>1),pause,end
%     figure(1)
%     my2dPlot(w1,w3,S(:,:,ii),'zlimit',0.1,'n_contours',20)
% 
%     figure(2)
%     my2dPlot(w1,w3,S3(:,:,ii),'zlimit',0.1,'n_contours',20)
% 
%     figure(3)
%     my2dPlot(w1,w3,dataMatrix(:,:,ii),'zlimit',0.1,'n_contours',20)
%     
% end

%%
ub = [10 500  15 30 1.66 1 1 5];
lb = [0 10  0.1 20 1.45 0.0001 0.0001 0.0001];


pfit2 = globalFit_cjk(dataMatrix,w1,w3,pfit,options,lb, ub)

%%
ub = [10 500  15 30 1.66 1 1 5];
lb = [0 10  0.1 20 1.0 0.0001 0.0001 0.0001];

pfit3 = globalFit_cjk(dataMatrix,w1,w3,pfit2,options,lb, ub)
%%
[S3,extra3] = analyticalResponseFunctionsFun_cjk(pfit2,w1,w3,options);
t3 = extra.t;

%%
ADD_TO_STARTUP;
cd('~/Box/CJ_Kelsheimer/PEGDA_HGS')
load('kin.mat') %remove if different source of kinetics data
options.t2_array = PEGDA_2d_crop.t2; %units of ps! 
options.damping = '2exp1fast';
% options.pnames = {'Delta1 (cm-1)','tau1 (ps)','anh (cm-1)'};
options.pnames = {'Delta1 (cm-1)','tau1 (ps)','Delta2 (cm-1)','tau2(ps)','T2 (ps)','anh (cm-1)','k1','k2','k3'};

options.dt = 0.20; %ps
options.n_t = 64;
options.order = 3;
options.kin = kin;

Delta1_cm  = 0.1; %cm-1
tau1 = 5; %ps
Delta2_cm = 0.1;
tau2 = 85;
anh_cm =25; %doesn't matter for linear
mean_w_0 = 2337.5;

w1 = PEGDA_2d_crop.w1;
w3 = PEGDA_2d_crop.w3; %not used
T2 = 3;

dE = 667;%cm^1
%assumption: kT in wavenumbers at 298K = 207cm^-1
Keq = 2*exp(-dE/207);

k1 = .004;%0.05
k2 =  0.008; %if this goes above 1 intensities switch
k3 = k2/Keq;

p0 = [Delta1_cm tau1 Delta2_cm tau2 T2 anh_cm k1 k2 k3]; %[delta1_cm delta2_cm delta3_cm delta4_cm tau1_ps tau2_ps tau3_ps anh_cm]

% not a fitting parameter right now
options.w_0_cm = mean_w_0;
options.p0 = p0; %save where we started
options.w_nu2_cm = -12; %put this somewhere better later - shift to bend/strech diag



[S,extra] = analyticalResponseFunctionsFun_cjk(p0,w1,w3,options);
t = extra.t;
% S1 = extra.S1;
g = extra.g;
c2 = extra.c2;

% figure(1),
% plot(w1,S)
% for ii = 1:size(S,3)
%     figure(ii),clf,my2dPlot(w1,w3,S(:,:,ii),'zlimit',1,'n_contours',20)
% end

%options.order = 1;
%S = analyticalResponseFunctionsFun(p0,w1,w3,options);

%Delta1 = Delta1_cm*wavenumbersToInvPs*2*pi;
%Lambda1 = 1/tau1;
%ratio = Delta1/Lambda1;
%disp('ratio Delta1/Lambda1')
%disp(ratio)
%% 2exp1fast

ADD_TO_STARTUP;
cd('~/Box/CJ_Kelsheimer/PEGDA_HGS')
load('kin.mat') %remove if different source of kinetics data
options.t2_array = PEGDA_2d_crop.t2; %units of ps! 
%options.t2_array = 10; %units of ps! 
options.damping = '2exp1fast';
options.pnames = {'Delta1 (cm-1)','tau1 (ps)','Delta2 (cm-1)','tau2','T2 (ps)','anh (cm-1)'};

options.dt = 0.200; %ps
options.n_t = 128;
options.order = 3;
options.kin = kin;

Delta1_cm  = 3; %cm-1
tau1 = 5; %ps
Delta2_cm  = 2; %cm-1
tau2 = 100; %ps
anh_cm =25; %doesn't matter for linear
mean_w_0 = 2340;

w1 = PEGDA_2d_crop.w1;
T2 = 0.5;
w3 = PEGDA_2d_crop.w3;

dE = 667;%cm^1
%assumption: kT in wavenumbers at 298K = 207cm^-1
Keq = 2*exp(-dE/207);

k1 = .004;%0.05
k2 =  0.008; %if this goes above 1 intensities switch
k3 = k2/Keq;

p0 = [Delta1_cm tau1 Delta2_cm tau2 T2 anh_cm k1 k2 k3]

options.w_0_cm = mean_w_0;
options.p0 = p0; %save where we started
options.w_nu2_cm = -12; %put this somewhere better later - shift to bend/strech diag

[S,extra] = analyticalResponseFunctionsFun_cjk(p0,w1,w3,options);
t = extra.t;
% S1 = extra.S1;
g = extra.g;


figure(1),
%plot(w1,S)
for ii = 1:size(S,3)
    figure(ii),clf,my2dPlot(w1,w3,S(:,:,ii),'zlimit',1,'n_contours',20)
end


% options.order = 1;
% S = analyticalResponseFunctionsFun(p0,w1,w3,options);



% Delta1 = Delta1_cm*wavenumbersToInvPs*2*pi;
% Lambda1 = 1/tau1;
% ratio = Delta1/Lambda1;
% disp('ratio Delta1/Lambda1')
% disp(ratio)

%%
% options.t2_array = 0/1000; %units of ps! 
% options.damping = 'voigt';
% options.pnames = {'Delta1 (cm-1)','tau1 (ps)','anh (cm-1)'};
% options.dt = 0.005; %ps
% options.n_t = 256;
% options.order = 1; %linear only!
% 
% Delta1_cm  = 30; %cm-1
% tau1 = .5; %ps
% anh_cm = 25; %doesn't matter for linear
% mean_w_0 = 2050;
% w1 = [1850:2:2250];
% w3 = []; %not used
% 
% p0 = [Delta1_cm tau1 anh_cm]; %[delta1_cm delta2_cm delta3_cm delta4_cm tau1_ps tau2_ps tau3_ps anh_cm]
% 
% % not a fitting parameter right now
% options.w_0_cm = mean_w_0;
% options.p0 = p0; %save where we started
% 
% [S,extra] = analyticalResponseFunctionsFun(p0,w1,w3,options);
% t = extra.t;
% S1 = extra.S1;
% g = extra.g;
% c2 = extra.c2;
% 
% figure(1),
% plot(w1,S)
% 
% figure(2),
% subplot(3,1,1)
% plot(t,c2)
% ylabel('c2(t)')
% subplot(3,1,2)
% plot(t,g)
% ylabel('g(t)')
% subplot(3,1,3)
% plot(t,S1)
% ylabel('S(t)')
% 
% Delta1 = Delta1_cm*wavenumbersToInvPs*2*pi;
% Lambda1 = 1/tau1;
% ratio = Delta1/Lambda1;
% disp('ratio Delta1/Lambda1')
% disp(ratio)
% 
