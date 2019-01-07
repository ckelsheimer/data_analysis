%% Load Experimental Data for Reference/Fit
    % data needs to be in format that comes out of
    % load2DIRdata/sort2DIRdata/calibrate2DIRdata/cropData - CALIBRATE DATA FIRST
    
    switch computer
        case 'MACI64'
            cd('~/Box/CJ_Kelsheimer/PEGDA_HGS') %MAC (for CJ's data)
        case 'GLNXA64'
            cd('/home/sgrlab/Documents/CJK/data/') %LINUX MACHINE (for CJ's data)
    end
    
load('pegda_022018.mat') 
%this will get your data into the format alt-globalFit needs.
dataMatrix = prepareGlobalFitData(crop_data); 
%% CHI BY EYE - use caRFF to ensure that your data 
ADD_TO_STARTUP;
%load('kin.mat') %remove if different source of  kinetics data
crosspeaks_kinetics_matlabfun;
%%
options.t2_array = [crop_data.t2]/1000; %units of ps! 
options.damping = '1exp1fast'; %select in aRFF (or whatever version you're using)
% The damping you choose WILL IMPACT EVERYTHING BELOW. Ensure you have all
% required variables

options.pnames = {'Delta1 (cm-1)','tau1 (ps)','T2 (ps)','anh (cm-1)','phi (deg)','mu12_2','k1','k2','k3'};
%These are your fitting parameters

Delta1_cm  = 1.75; %cm-1 - from correlation function
tau1 = 110; %ps - from correlation function
anh_cm = 23.5; %anharmonicity
mean_w_0 = 2337.5; %CENTER FREQUENCY

w1 = crop_data.w1;  %set as open range or crop to match data being fit
w3 = crop_data.w3; %same as above
T2 = 10; % dephasing time (1/value scales with antidiagonal width)
mu12_2 = 2; %intensity of dipole for 1-2 transition (can be removed and left at default value)
phi_deg = -9; %phase contributions

%KINETICS INFORMATION
dE = 667;%cm^1 %energy gap 
k_B_cm = 0.69503476; %boltzmann's constant in wavenumbers
RT = 23; % temperature
T = 273+RT; % conversion to kelvin
Keq = 2*exp(-dE/(T*k_B_cm)); %calculate equlibirum constant

k1 = .0028;%0.05 %growth of HGS peak
k2 =  0.0015; %K_up
k3 = k2/Keq; %K_down

%STARTING VALUES/PARAMETERS

p0 = [Delta1_cm tau1 T2 anh_cm phi_deg mu12_2 k1 k2 k3]; %[delta1_cm delta2_cm delta3_cm delta4_cm tau1_ps tau2_ps tau3_ps anh_cm]

%options 

options.dt = 0.20; %ps fft stuff
options.n_t = 64; %fft stuff
options.order = 3; %1 - ftir, 3 - 2DIR
options.kin = kin; %comment if no kinetics data 
options.w_0_cm = mean_w_0; %centerfrequency
options.p0 = p0; %save where we started
options.w_nu2_cm = -12; %shift for diag peaks

%%
[S,extra] = analyticalResponseFunctionsFun_cjk(p0,w1,w3,options); %generates calculated spectra from input parameters
t = extra.t; %honestly no idea what this does 
%%
% this plots generated spectra(S) vs dataMatrix (real data) - minimize
% residuals as much as possible.
% COMMENT THESE OUT ONCE THE FIT (by eye) IS A GOOD STARTING POINT
% % 
% for ii = 1:size(S,3)
%     for jj = ii
%     if(ii>1),pause,end
%     n1 = abs(min(min(dataMatrix(:,:,jj)))); %normalization
%     n2 = abs(min(min(S(:,:,ii))));
%     
% %calculated spectrum
% figure(1),
%     clf,my2dPlot(w1,w3,S(:,:,ii)./n2,'n_contours',20)
% %experimental spectrum
% figure(2)
%     my2dPlot(w1,w3,dataMatrix(:,:,jj)./n1,'n_contours',20)
% %residual
% figure(3)
%     my2dPlot(w1,w3,dataMatrix(:,:,jj)./n1-S(:,:,ii)./n2,'n_contours',20);
%     end
% end
%% GLOBAL FITTING
%p0 = [Delta1_cm tau1 T2 anh_cm phi_deg mu12_2 k1 k2 k3]
ub = [5 500  15  30  40  2.5      1      1      1];
lb = [0  10  0.1  20 -40 1.4 0.0001 0.0001 0.0001];

[pfit,err,gfstruct] = globalFit_cjk(dataMatrix,w1,w3,p0,options,lb, ub);%actually fits the spectrum
log_gfit; %records some information about the fit
%%
[S2,extra] = analyticalResponseFunctionsFun_cjk(pfit,w1,w3,options); %generates calculated spectra from input parameters
t = extra.t; %honestly no idea what this does 
%%
% this plots generated spectra(S) vs dataMatrix (real data) - minimize
% residuals as much as possible.
% COMMENT THESE OUT ONCE THE FIT (by eye) IS A GOOD STARTING POINT
% 
for ii = 1:size(S2,3)
    for jj = ii
    if(ii>1),pause,end
    n1 = abs(min(min(dataMatrix(:,:,jj)))) %normalization
    n2 = abs(min(min(S2(:,:,ii))))
    

%calculated spectrum
figure(1),
    clf,my2dPlot(w1,w3,S2(:,:,ii)./n2,'n_contours',20)
    
%experimental spectrum
figure(2)
    my2dPlot(w1,w3,dataMatrix(:,:,jj)./n1,'n_contours',20)
%residual
figure(3)
    my2dPlot(w1,w3,dataMatrix(:,:,jj)./n1-S2(:,:,ii)./n2,'n_contours',20);
    end
end