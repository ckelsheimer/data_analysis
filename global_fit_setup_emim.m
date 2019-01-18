%% Load Experimental Data for Reference/Fit
    % data needs to be in format that comes out of
    % load2DIRdata/sort2DIRdata/calibrate2DIRdata/cropData - CALIBRATE DATA FIRST
    
    switch computer
        case 'MACI64'
            cd('~/Box/SGRLAB RESEARCH/Globalfitting_data') %MAC (for CJ's data)
        case 'GLNXA64'
            cd('/home/sgrlab/Documents/CJK/data/') %LINUX MACHINE (for CJ's data)
    end
    
load('emimtf2n_pol_gfit.mat') 
%this will get your data into the format alt-globalFit needs.
dataMatrix = prepareGlobalFitData(crop_par); 
%% CHI BY EYE - use caRFF to ensure that your data 
ADD_TO_STARTUP;
%load('kin.mat') %remove if different source of  kinetics data
crosspeaks_kinetics_matlabfun;
%%
options.t2_array = [crop_par.t2]/1000; %units of ps! 
options.damping = '1exp1fast'; %select in aRFF (or whatever version you're using)
% The damping you choose WILL IMPACT EVERYTHING BELOW. Ensure you have all
% required variables

options.pnames = {'Delta1 (cm-1)','tau1 (ps)','T2 (ps)','anh (cm-1)','phi (deg)','mu12_2','k1','k2','k3'};
%These are your fitting parameters

Delta1_cm  = 2.6; %cm-1 - from correlation function
tau1 = 90; %ps - from correlation function
anh_cm = 24.75; %anharmonicity
mean_w_0 = 2340.5; %CENTER FREQUENCY

w1 = crop_par.w1;  %set as open range or crop to match data being fit
w3 = crop_par.w3; %same as above
T2 = 6.5; % dephasing time (1/value scales with antidiagonal width)
mu12_2 = 1.45; %intensity of dipole for 1-2 transition (can be removed and left at default value)
phi_deg = -12; %phase contributions

%KINETICS INFORMATION
dE = 667;%cm^1 %energy gap 
k_B_cm = 0.69503476; %boltzmann's constant in wavenumbers 
RT = 25.4; % temperature
T = 273+RT; % conversion to kelvin
Keq = 2*exp(-dE/(T*k_B_cm)); %calculate equlibirum constant

k1 = 0.0015;%0.05 %growth of HGS peak
k2 = 0.1;%;0.025; %K_up
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
%options.poptrans = 1;
%options.hgs = 1

%%
[S,extra] = analyticalResponseFunctionsFun_cjk(p0,w1,w3,options); %generates calculated spectra from input parameters
t = extra.t; %honestly no idea what this does 

%%
% this plots generated spectra(S) vs dataMatrix (real data) - minimize
% residuals as much as possible.
% COMMENT THESE OUT ONCE THE FIT (by eye) IS A GOOD STARTING POINT
% % 
for ii = 1:size(S,3)
    for jj = ii
    if(ii>1),pause,end
    n1 = abs(min(min(dataMatrix(:,:,jj)))); %normalization
    n2 = abs(min(min(S(:,:,ii))));
    
%calculated spectrum
figure(1),
    clf,my2dPlot(w1,w3,S(:,:,ii)./n2,'n_contours',20,'zlimit',.1)
%experimental spectrum
figure(2)
    my2dPlot(w1,w3,dataMatrix(:,:,jj)./n1,'n_contours',20,'zlimit',.1)
%residual
figure(3)
    my2dPlot(w1,w3,dataMatrix(:,:,jj)./n1-S(:,:,ii)./n2,'n_contours',20,'zlimit',0.1);
    end
end
 %% GLOBAL FITTING
%p0 = [Delta1_cm tau1 T2 anh_cm phi_deg mu12_2 k1 k2 k3]
ub = [10 250  17  30  40  2.1  1      0.49      4];
lb = [0  0  0.01  20 -40 1.4 0.000000001 0.00001 0.5];

[pfit,err,gfstruct] = globalFit_cjk(dataMatrix,w1,w3,p0,options,lb, ub);%actually fits the spectrum
log_gfit; %record s some information about the fit