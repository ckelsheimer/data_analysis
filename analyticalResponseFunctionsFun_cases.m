function [out,extra] = analyticalResponseFunctionsFun_cjk(p,w1_in,w3_in,options)
%Create simulated spectral data using an analytical form for the nth-order
%response function (see Hamm and Zanni Ch. 7 in particular).
%  [OUT,EXTRA] = analyticalResponseFunctionsFun(P,W1_IN,W3_IN,OPTIONS) will
%  generate a matrix with the simulated response in W1 and W3, based on the
%  starting point (P) and options (OPTIONS) you provide.
%
%   OPTIONS should be a structure with the fields:
%
%       't2_array'              t2 values in ps, in array format, e.g.:
%                               [1 2 3 4 5]
%       'dt'                    Size of the timestep (0.4 is the default)
%       'n_t'                   Number of time points (64 is the default)
%       'noise' (optional)
%       'order' (optional)      1 -- First order (linear)
%                               3 -- Third order (2D)
%                               5 -- Fifth order (3D)
%       'w_0_cm'                Now a fitting parameter -- allows variation 
%                               of w_0 for the fit -- Leave this [].
%       'bootstrap'             No bootstrapping for now (set as []).
%       ----------------------
%       'damping'               Gives the form of the correlation function
%       'pnames'                Gives names for the inputs from P
%       'p0'                    Starting point.
%
%       Damping determines the elements in pnames and p0. We'll
%       have to make this a little more user friendly over time. The
%       general form for OPTIONS.pnames is ~~:
%           options.pnames = {'Delta (cm-1)','tau (ps)','anh (cm-1)',...
%                               'mu12_2','w0 (cm-1)','phi (rad)'};
%       ----------------------
%
%
%   ALL TIME PARAMETERS ARE IN UNITS OF PICOSECONDS (ps). By default, this
%   function currently simulates 3rd order (2D-IR) spectra.


damping = options.damping;
t2_array = options.t2_array; %can be a vector of times

%for 5th order also define t4 times (not used for 2DIR)
t4_array = t2_array;
n_t2_array = length(t2_array);

flag_print = 0; %1 => figures or 0 => no figures
flag_plot = 0;
order = 3; %order of spectroscopy to calculate. 3 = 2DIR, 5 = 3DIR

if isfield(options,'dt')
    dt = options.dt;
else
    dt = 0.200;%added_cjk
end

if isfield(options,'n_t')
    n_t = options.n_t;
else
    n_t = 128;%added_cjk
end

w_0_cm = options.w_0_cm;% %center frequency
w_nu2_cm = options.w_nu2_cm; %shift to bend/stretch diagonal peak
phi = 0; %phase shift (radians) exp(1i*phi)
mu01_2 = 1; %default
mu12_2 = 2; %default

%look to see if parameters from inputs
nparams = length(options.pnames);
for ii = 1:nparams
    switch options.pnames{ii}
        case 'w0 (cm-1)'
            w_0_cm = p(ii);
        case 'mu01_2'
            mu01_2 = p(ii);
        case 'mu12_2'
            mu12_2 = p(ii);
        case 'phi (rad)'
            phi = p(ii);
        case 'phi (deg)'
            phi = p(ii)*pi/180;
    end   
end

%check for what order (third 
if isfield(options,'order')
    order = options.order;
end

flag_rotating_frame = true;
if isfield(options,'flag_rotating_frame')
    flag_rotating_frame = options.flag_rotating_frame;
end

two_level_system = false;
%two_level_system = true;

flag_bootstrap = false;
if isfield(options,'bootstrap')
    if ~isempty(options.bootstrap)
        flag_bootstrap = true;
        bootstrap_index = options.bootstrap;
    end
end

%details of fft

fft_type = 'sgrsfft';
n_zp = 2*n_t; %total length after zeropadding 
            %(make n_zp=n_t for no zero padding)
n_under = 0;%2;

%type of projection to calculate
%projection_type = 'window';
projection_type = 'all';

% Body of the calculation

extra = [];

%-------------------------------------------------------------
%
% start calculation
%
%-------------------------------------------------------------
c = 2.9979e10;
wavenumbersToInvPs=c*1e-12;
invPsToWavenumbers=1/wavenumbersToInvPs;
c_cmfs = c*1e-15;
t=0:dt:(n_t-1)*dt;

w = fftFreqAxis(t,...
  'time_units','ps',...
  'freq','wavenumbers',...
  'shift','on',...
  'zeropad',n_zp,...
  'undersampling',n_under);
if flag_rotating_frame
    w = w + w_0_cm;
end

dw = w(2)-w(1);
w_0 = w_0_cm*2*pi*wavenumbersToInvPs; %convert to radians
w_nu2 = w_nu2_cm*2*pi*wavenumbersToInvPs;
  
c2 = [];
g = [];
switch damping,

    case 'voigt'
        Delta1_cm = p(1);%linewidth (sigma) in wavenumbers of one motion
        T2 = p(2); %first timescale (ps)
        anh_cm = p(3);

        Delta1 = Delta1_cm*wavenumbersToInvPs*2*pi;
        anh = anh_cm*wavenumbersToInvPs*2*pi;

        c2 = @(t) (t==0)/T2 + Delta1^2;
        g = @(t) t./T2 + Delta1^2/2.*t.^2;
        
    case {'overdamped', '1exp'}
        %overdamped exp(-t/tau)
        Delta1_cm = p(1);%linewidth (sigma) in wavenumbers of one motion
        tau1 = p(2); %first timescale (ps)
        anh_cm = p(3);

        Delta1 = Delta1_cm*wavenumbersToInvPs*2*pi;
        Lambda1 = 1/tau1;
        anh = anh_cm*wavenumbersToInvPs*2*pi;

        c2 = @(t) Delta1^2.*exp(-Lambda1.*t);
        g = @(t) Delta1^2/Lambda1^2.*(exp(-Lambda1.*t)-1+Lambda1.*t);

    case {'1exp1fast'}
        %overdamped exp(-t/tau)
        Delta1_cm = p(1);%linewidth (sigma) in wavenumbers of one motion
        tau1 = p(2); %first timescale (ps)
        T2 = p(3); %T2 time (fast process)
        anh_cm = p(4);

        Delta1 = Delta1_cm*wavenumbersToInvPs*2*pi;
        Lambda1 = 1/tau1;
        anh = anh_cm*wavenumbersToInvPs*2*pi;

        c2 = @(t) (t==0)/T2 + Delta1^2.*exp(-Lambda1.*t);
        g = @(t) t./T2 + Delta1^2/Lambda1^2.*(exp(-Lambda1.*t)-1+Lambda1.*t);

    case {'1exp1slow'}
        Delta1_cm = p(1);%linewidth (sigma) in wavenumbers of one motion
        tau1 = p(2);%first motion's timescale (ps)
        Delta2_cm = p(3);%linewidth (sigma) in wavenumbers of the inhomogeneous component
        anh_cm = p(4); %Anharmonicity in wavenumbers

        Delta1 = Delta1_cm*wavenumbersToInvPs*2*pi;
        Lambda1 = 1/tau1;
        Delta2 = Delta2_cm*wavenumbersToInvPs*2*pi;
        anh = anh_cm*wavenumbersToInvPs*2*pi;

        c2 = @(t) Delta1^2.*exp(-Lambda1.*t) ...
            + Delta2^2;
        g = @(t) Delta1^2/Lambda1^2.*(exp(-Lambda1.*t)-1+Lambda1.*t) ...
        + (Delta2^2).*t.^2/2;
        
    case {'1exp1fast1slow'}
        Delta1_cm = p(1);%linewidth (sigma) in wavenumbers of one motion
        tau1 = p(2);%first motion's timescale (ps)
        Delta2_cm = p(3);%linewidth (sigma) in wavenumbers of the inhomogeneous component
        T2 = p(4); %Homogeneous dephasing time
        anh_cm = p(5); %Anharmonicity in wavenumbers

        Delta1 = Delta1_cm*wavenumbersToInvPs*2*pi;
        Lambda1 = 1/tau1;
        Delta2 = Delta2_cm*wavenumbersToInvPs*2*pi;
        anh = anh_cm*wavenumbersToInvPs*2*pi;

        c2 = @(t) (t==0)/T2 + Delta1^2.*exp(-Lambda1.*t) ...
            + Delta2^2;
        g = @(t) t./T2 + Delta1^2/Lambda1^2.*(exp(-Lambda1.*t)-1+Lambda1.*t) ...
        + (Delta2^2).*t.^2/2;
  
    case {'2exp1fast'}
        % Developed for use with Zhe's SCN- data in ILs, where he has two
        % inhomogeneous timescales, one longer, and one shorter.
        % Expanded form of '1exp1fast' above --Tom
        Delta1_cm = p(1);%linewidth (sigma) in wavenumbers of one motion
        tau1 = p(2);%first timescale (ps)
        Delta2_cm = p(3); %linewidth of second motion
        tau2 = p(4); %second timescale (ps)
        T2 = p(5); %T2 time (fast / homogeneous processes)
        anh_cm = p(6);
        
        Delta1 = Delta1_cm*wavenumbersToInvPs*2*pi;
        Lambda1 = 1/tau1;
        Delta2 = Delta2_cm*wavenumbersToInvPs*2*pi;
        Lambda2 = 1/tau2;
        anh = anh_cm*wavenumbersToInvPs*2*pi;

        c2 = @(t) (t==0)/T2 + Delta1^2.*exp(-Lambda1.*t) ...
            + Delta2^2.*exp(-Lambda2.*t);
        g = @(t) t./T2 + Delta1^2/Lambda1^2.*(exp(-Lambda1.*t)-1+Lambda1.*t) ...
            + Delta2^2/Lambda2^2.*(exp(-Lambda2.*t)-1+Lambda2.*t);
    
  case {'2exp1fast1slow'}
        Delta1_cm = p(1);%linewidth (sigma) in wavenumbers of one motion
        tau1 = p(2);%first motion's timescale (ps)
        Delta2_cm = p(3);%linewidth (sigma) in wavenumbers of the other motion
        tau2 = p(4);%second motion's timescale (ps)
        Delta3_cm = p(5);%linewidth (sigma) in wavenumbers of the inhomogeneous component
        T2 = p(6); %Homogeneous dephasing time
        anh_cm = p(7); %Anharmonicity in wavenumbers
      
        Delta1 = Delta1_cm*wavenumbersToInvPs*2*pi;
        Lambda1 = 1/tau1;
        Delta2 = Delta2_cm*wavenumbersToInvPs*2*pi;
        Lambda2 = 1/tau2;
        Delta3 = Delta3_cm*wavenumbersToInvPs*2*pi;
        anh = anh_cm*wavenumbersToInvPs*2*pi;

        c2 = @(t) (t==0)/T2 + Delta1^2.*exp(-Lambda1.*t) ...
            + Delta2^2.*exp(-Lambda2.*t) + Delta3^2;
        g = @(t) t./T2 + Delta1^2/Lambda1^2.*(exp(-Lambda1.*t)-1+Lambda1.*t) ...
        + Delta2^2/Lambda2^2.*(exp(-Lambda2.*t)-1+Lambda2.*t) + (Delta3^2).*t.^2/2;
    
    case {'2exp1slow'}
        Delta1_cm = p(1);%linewidth (sigma) in wavenumbers of one motion
        tau1 = p(2);%first motion's timescale (ps)
        Delta2_cm = p(3);%linewidth (sigma) in wavenumbers of the other motion
        tau2 = p(4);%second motion's timescale (ps)
        Delta3_cm = p(5);%linewidth (sigma) in wavenumbers of the inhomogeneous component
        anh_cm = p(6); %Anharmonicity in wavenumbers
      
        Delta1 = Delta1_cm*wavenumbersToInvPs*2*pi;
        Lambda1 = 1/tau1;
        Delta2 = Delta2_cm*wavenumbersToInvPs*2*pi;
        Lambda2 = 1/tau2;
        Delta3 = Delta3_cm*wavenumbersToInvPs*2*pi;
        anh = anh_cm*wavenumbersToInvPs*2*pi;

        c2 = @(t) Delta1^2.*exp(-Lambda1.*t) ...
            + Delta2^2.*exp(-Lambda2.*t) + Delta3^2;
        g = @(t) Delta1^2/Lambda1^2.*(exp(-Lambda1.*t)-1+Lambda1.*t) ...
        + Delta2^2/Lambda2^2.*(exp(-Lambda2.*t)-1+Lambda2.*t) + (Delta3^2).*t.^2/2;
  
    case {'3exp1fast'}
        Delta1_cm = p(1);%linewidth (sigma) in wavenumbers of one motion
        tau1 = p(2);%first timescale (ps)
        Delta2_cm = p(3); %linewidth of second motion
        tau2 = p(4); %second timescale (ps)
        Delta3_cm = p(5);
        tau3 = p(6);
        T2 = p(7); %T2 time (fast / homogeneous processes)
        anh_cm = p(8);
        
        Delta1 = Delta1_cm*wavenumbersToInvPs*2*pi;
        Lambda1 = 1/tau1;
        Delta2 = Delta2_cm*wavenumbersToInvPs*2*pi;
        Lambda2 = 1/tau2;
        Delta3 = Delta3_cm*wavenumbersToInvPs*2*pi;
        Lambda3 = 1/tau3;
        anh = anh_cm*wavenumbersToInvPs*2*pi;

        c2 = @(t) (t==0)/T2 + Delta1^2.*exp(-Lambda1.*t) ...
            + Delta2^2.*exp(-Lambda2.*t) + Delta3^2.*exp(-Lambda3.*t);
        g = @(t) t./T2 + Delta1^2/Lambda1^2.*(exp(-Lambda1.*t)-1+Lambda1.*t) ...
            + Delta2^2/Lambda2^2.*(exp(-Lambda2.*t)-1+Lambda2.*t) ...
            + Delta3^2/Lambda3^2.*(exp(-Lambda3.*t)-1+Lambda3.*t);
    
    
    case 'critical'
    %critically damped (1+2t/tau)exp(-2t/tau)
    Delta1_cm = p(1);%linewidth (sigma) in wavenumbers of one motion
    tau1 = p(2); %first timescale (ps)
    anh_cm = p(3);

    Delta1 = Delta1_cm*wavenumbersToInvPs*2*pi;
    Lambda1 = 1/tau1;
    anh = anh_cm*wavenumbersToInvPs*2*pi;

    c2 = @(t) Delta1^2.*(1+2*Lambda1.*t).*exp(-2*Lambda1.*t);
    g = @(t) Delta1^2/4/Lambda1^2.*exp(-2.*Lambda1.*t) ...
      .*(3 + 2*Lambda1*t + exp(2.*Lambda1.*t).*(4*Lambda1.*t - 3));
  
  case '2expcrit'
    Delta1_cm = p(1);%linewidth (sigma) in wavenumbers of one motion
    Delta2_cm = p(2);%linewidth (sigma) in wavenumbers of the other motion
    tau1 = p(3); %first timescale (ps)
    tau2 = p(4); %second timescale (ps)
    anh_cm = p(5);

    Delta1 = Delta1_cm*wavenumbersToInvPs*2*pi;
    Delta2 = Delta2_cm*wavenumbersToInvPs*2*pi;
    Lambda1 = 1/tau1;
    Lambda2 = 1/tau2;
    anh = anh_cm*wavenumbersToInvPs*2*pi;
    g = @(t) Delta1^2/4/Lambda1^2 ...
      .*(3.*exp(-2.*Lambda1.*t) + 2*Lambda1*t.*exp(-2.*Lambda1.*t) + (4*Lambda1.*t - 3)) ...
      + Delta2^2/4/Lambda2^2 ...
      .*(3.*exp(-2.*Lambda2.*t) + 2*Lambda2*t.*exp(-2.*Lambda2.*t) + (4*Lambda2.*t - 3));
  
  case '3expcrit'
    Delta1_cm = p(1);%linewidth (sigma) in wavenumbers of one motion
    Delta2_cm = p(2);%linewidth (sigma) in wavenumbers of the other motion
    Delta3_cm = p(3);%linewidth (sigma) in wavenumbers of the other motion
    tau1 = p(4); %first timescale (ps)
    tau2 = p(5); %second timescale (ps)
    tau3 = p(6); %second timescale (ps)
    anh_cm = p(7);

    Delta1 = Delta1_cm*wavenumbersToInvPs*2*pi;
    Delta2 = Delta2_cm*wavenumbersToInvPs*2*pi;
    Delta3 = Delta3_cm*wavenumbersToInvPs*2*pi;
    Lambda1 = 1/tau1;
    Lambda2 = 1/tau2;
    Lambda3 = 1/tau3;
    anh = anh_cm*wavenumbersToInvPs*2*pi;
    g = @(t) Delta1^2/4/Lambda1^2 ...
      .*(3.*exp(-2.*Lambda1.*t) + 2*Lambda1*t.*exp(-2.*Lambda1.*t) + (4*Lambda1.*t - 3)) ...
      + Delta2^2/4/Lambda2^2 ...
      .*(3.*exp(-2.*Lambda2.*t) + 2*Lambda2*t.*exp(-2.*Lambda2.*t) + (4*Lambda2.*t - 3)) ...      
      + Delta3^2/4/Lambda3^2 ...
      .*(3.*exp(-2.*Lambda3.*t) + 2*Lambda3*t.*exp(-2.*Lambda3.*t) + (4*Lambda3.*t - 3));

  case '3exp1offcrit'
    Delta1_cm = p(1);%linewidth (sigma) in wavenumbers of one motion
    Delta2_cm = p(2);%linewidth (sigma) in wavenumbers of the other motion
    Delta3_cm = p(3);%linewidth (sigma) in wavenumbers of the other motion
    Delta4_cm = p(4);%linewidth (sigma) in wavenumbers of static component
    tau1 = p(5); %first timescale (ps)
    tau2 = p(6); %second timescale (ps)
    tau3 = p(7); %second timescale (ps)
    anh_cm = p(8);

    Delta1 = Delta1_cm*wavenumbersToInvPs*2*pi;
    Delta2 = Delta2_cm*wavenumbersToInvPs*2*pi;
    Delta3 = Delta3_cm*wavenumbersToInvPs*2*pi;
    Delta4 = Delta4_cm*wavenumbersToInvPs*2*pi;
    Lambda1 = 1/tau1;
    Lambda2 = 1/tau2;
    Lambda3 = 1/tau3;
    anh = anh_cm*wavenumbersToInvPs*2*pi;
    g = @(t) Delta1^2/4/Lambda1^2 ...
      .*(3.*exp(-2.*Lambda1.*t) + 2*Lambda1*t.*exp(-2.*Lambda1.*t) + (4*Lambda1.*t - 3)) ...
      + Delta2^2/4/Lambda2^2 ...
      .*(3.*exp(-2.*Lambda2.*t) + 2*Lambda2*t.*exp(-2.*Lambda2.*t) + (4*Lambda2.*t - 3)) ...      
      + Delta3^2/4/Lambda3^2 ...
      .*(3.*exp(-2.*Lambda3.*t) + 2*Lambda3*t.*exp(-2.*Lambda3.*t) + (4*Lambda3.*t - 3)) ...
      + Delta4^2.*t.^2/2;

  case '2exp'
    Delta1_cm = p(1);%linewidth (sigma) in wavenumbers of one motion
    Delta2_cm = p(2);%linewidth (sigma) in wavenumbers of the other motion
    tau1 = p(3); %first timescale (ps)
    tau2 = p(4); %second timescale (ps)
    anh_cm = p(5);

    Delta1 = Delta1_cm*wavenumbersToInvPs*2*pi;
    Delta2 = Delta2_cm*wavenumbersToInvPs*2*pi;
    Lambda1 = 1/tau1;
    Lambda2 = 1/tau2;
    anh = anh_cm*wavenumbersToInvPs*2*pi;
    g = @(t) Delta1^2/Lambda1^2 ...
      .*(exp(-Lambda1.*t) - 1 + Lambda1*t) ...
      + Delta2^2/Lambda2^2 ...
      .*(exp(-Lambda2.*t) - 1 + Lambda2*t);

  case '3exp'
    Delta1_cm = p(1);%linewidth (sigma) in wavenumbers of one motion
    Delta2_cm = p(2);%linewidth (sigma) in wavenumbers of the other motion
    Delta3_cm = p(3);%linewidth (sigma) in wavenumbers of the other motion
    tau1 = p(4); %first timescale (ps)
    tau2 = p(5); %second timescale (ps)
    tau3 = p(6); %second timescale (ps)
    anh_cm = p(7);

    Delta1 = Delta1_cm*wavenumbersToInvPs*2*pi;
    Delta2 = Delta2_cm*wavenumbersToInvPs*2*pi;
    Delta3 = Delta3_cm*wavenumbersToInvPs*2*pi;
    Lambda1 = 1/tau1;
    Lambda2 = 1/tau2;
    Lambda3 = 1/tau3;
    anh = anh_cm*wavenumbersToInvPs*2*pi;
    g = @(t)Delta1^2/Lambda1^2 ...
      .*(exp(-Lambda1.*t) - 1 + Lambda1*t) ...
      + Delta2^2/Lambda2^2 ...
      .*(exp(-Lambda2.*t) - 1 + Lambda2*t) ...
      + Delta3^2/Lambda3^2 ...
      .*(exp(-Lambda3.*t) - 1 + Lambda3*t);

    case '5exp'
    Delta1_cm = p(1);%linewidth (sigma) in wavenumbers of one motion
    Delta2_cm = p(2);%linewidth (sigma) in wavenumbers of the other motion
    Delta3_cm = p(3);%linewidth (sigma) in wavenumbers of the other motion
    Delta4_cm = p(4);
    Delta5_cm = p(5);
    tau1 = p(6); %first timescale (ps)
    tau2 = p(7); %second timescale (ps)
    tau3 = p(8); %second timescale (ps)
    tau4 = p(9);
    tau5 = p(10);
    anh_cm = p(11);

    Delta1 = Delta1_cm*wavenumbersToInvPs*2*pi;
    Delta2 = Delta2_cm*wavenumbersToInvPs*2*pi;
    Delta3 = Delta3_cm*wavenumbersToInvPs*2*pi;
    Delta4 = Delta4_cm*wavenumbersToInvPs*2*pi;
    Delta5 = Delta5_cm*wavenumbersToInvPs*2*pi;
    Lambda1 = 1/tau1;
    Lambda2 = 1/tau2;
    Lambda3 = 1/tau3;
    Lambda4 = 1/tau4;
    Lambda5 = 1/tau5;
    anh = anh_cm*wavenumbersToInvPs*2*pi;
    g = @(t)Delta1^2/Lambda1^2.*(exp(-Lambda1.*t) - 1 + Lambda1*t) ...
      + Delta2^2/Lambda2^2.*(exp(-Lambda2.*t) - 1 + Lambda2*t) ...
      + Delta3^2/Lambda3^2.*(exp(-Lambda3.*t) - 1 + Lambda3*t) ...
      + Delta4^2/Lambda4^2.*(exp(-Lambda4.*t) - 1 + Lambda4*t) ...
      + Delta5^2/Lambda5^2.*(exp(-Lambda5.*t) - 1 + Lambda5*t);
  
  
    case 'toms_5exp_T1_or_rlx'
    Delta1_cm = p(1);%linewidth (sigma) in wavenumbers of one motion
    Delta2_cm = p(2);%linewidth (sigma) in wavenumbers of the other motion
    Delta3_cm = p(3);%linewidth (sigma) in wavenumbers of the other motion
    Delta4_cm = p(4);
    Delta5_cm = p(5);
    tau1 = p(6); %first timescale (ps)
    tau2 = p(7); %second timescale (ps)
    tau3 = p(8); %second timescale (ps)
    tau4 = p(9);
    tau5 = p(10);
    T_or = p(11);
    T1 = p(12);
    anh_cm = p(13);

    Delta1 = Delta1_cm*wavenumbersToInvPs*2*pi;
    Delta2 = Delta2_cm*wavenumbersToInvPs*2*pi;
    Delta3 = Delta3_cm*wavenumbersToInvPs*2*pi;
    Delta4 = Delta4_cm*wavenumbersToInvPs*2*pi;
    Delta5 = Delta5_cm*wavenumbersToInvPs*2*pi;
    Lambda1 = 1/tau1;
    Lambda2 = 1/tau2;
    Lambda3 = 1/tau3;
    Lambda4 = 1/tau4;
    Lambda5 = 1/tau5;
    anh = anh_cm*wavenumbersToInvPs*2*pi;
    g = @(t)Delta1^2/Lambda1^2.*(exp(-Lambda1.*t) - 1 + Lambda1*t) ...
      + Delta2^2/Lambda2^2.*(exp(-Lambda2.*t) - 1 + Lambda2*t) ...
      + Delta3^2/Lambda3^2.*(exp(-Lambda3.*t) - 1 + Lambda3*t) ...
      + Delta4^2/Lambda4^2.*(exp(-Lambda4.*t) - 1 + Lambda4*t) ...
      + Delta5^2/Lambda5^2.*(exp(-Lambda5.*t) - 1 + Lambda5*t) ...
      + t./(3*T_or) + t./(2*T1);
  
  
    case '3exp1off'
    Delta1_cm = p(1);%linewidth (sigma) in wavenumbers of one motion
    Delta2_cm = p(2);%linewidth (sigma) in wavenumbers of the other motion
    Delta3_cm = p(3);%linewidth (sigma) in wavenumbers of the other motion
    Delta4_cm = p(4);%linewidth (sigma) in wavenumbers of static component
    tau1 = p(5); %first timescale (ps)
    tau2 = p(6); %second timescale (ps)
    tau3 = p(7); %second timescale (ps)
    anh_cm = p(8);

    Delta1 = Delta1_cm*wavenumbersToInvPs*2*pi;
    Delta2 = Delta2_cm*wavenumbersToInvPs*2*pi;
    Delta3 = Delta3_cm*wavenumbersToInvPs*2*pi;
    Delta4 = Delta4_cm*wavenumbersToInvPs*2*pi;
    Lambda1 = 1/tau1;
    Lambda2 = 1/tau2;
    Lambda3 = 1/tau3;
    anh = anh_cm*wavenumbersToInvPs*2*pi;
    g = @(t)Delta1^2/Lambda1^2 ...
      .*(exp(-Lambda1.*t) - 1 + Lambda1*t) ...
      + Delta2^2/Lambda2^2 ...
      .*(exp(-Lambda2.*t) - 1 + Lambda2*t) ...
      + Delta3^2/Lambda3^2 ...
      .*(exp(-Lambda3.*t) - 1 + Lambda3*t) ...
      + Delta4^2.*t.^2/2;

  case 'hynesform'
  %fit c2 to hynes type function, double integrate with Mathematica
  %to find g(t)
  a1 = 0.3232;
  k11 = 30.01; %ps-1
  k12 = 17.41; %ps-1
  a2 = 0.3378;
  k2 = 8.270; %ps-1
  a3 = 0.3455;
  k3 = 1.897; %ps-1
  
  %yuck
  g = @(t) Delta^2*exp(-t.*(k12+k2+k3))/(k2^2*k3^2*(k11^2+k12^2)^2) ...
      .*(a1*k2^2*k3^2.*exp((k2+k3).*t).*(cos(k11.*t).*(k12^2-k11^2)-2*k11*k12*sin(k11*t) ...
					 + exp(k12.*t).*(k11^2*(k12.*t+1)+k12^2*(k12.*t-1))) ... 
	 + (k11^2+k12^2)^2.*exp(k12.*t).*(a3*k2^2*exp(k2.*t).*(exp(k3.*t).*(k3.*t-1)+1) ...
					  +a2*k3^2*exp(k3.*t).*(exp(k2.*t).*(k2*t-1)+1)));
  otherwise
    error('damping value is unknown');
end

if order==1
  %S = exp(-g(t)).*cos(w_0.*t);
  S1 = exp(-g(t));
  %S2 = exp(-g(t)).*cos(w_0.*t).*(2-2*exp(-sqrt(-1)*anh.*t));
  
  S  = fftshift(real(sgrsfft(S1,n_zp)));

  
 s.c2 = c2(t);
 s.g = g(t);
 s.S1 = S1;
 s.t = t;
 extra = s;
 out = interp1(w,S,w1_in,'pchip');
    return
end

if order==3
    P = zeros(length(w3_in),length(w1_in)); %signle time step
    out = zeros(length(w3_in),length(w1_in),n_t2_array); %array for output
    
  [T1,T3] = meshgrid(t,t);
  for i=1:n_t2_array
    t2 = t2_array(i);
      
   %blue peaks (main band)
    Br=exp(-g(T1)+g(t2)-g(T3)-g(T1+t2)-g(t2+T3)+g(T1+t2+T3)).*(2*mu01_2); %need these too
         Bn=exp(-g(T1)-g(t2)-g(T3)+g(T1+t2)+g(t2+T3)-g(T1+t2+T3)).*(2*mu01_2);
    %red peaks (main band)
    Rr=exp(-g(T1)+g(t2)-g(T3)-g(T1+t2)-g(t2+T3)+g(T1+t2+T3)).*(-mu12_2.*exp(-1i*anh.*T3)); %need these too
        Rn=exp(-g(T1)-g(t2)-g(T3)+g(T1+t2)+g(t2+T3)-g(T1+t2+T3)).*(-mu12_2.*exp(-1i*anh.*T3));
        
%  P1=exp(-g(T1)+g(t2)-g(T3)-g(T1+t2)-g(t2+T3)+g(T1+t2+T3)).*(2*mu01_2-mu12_2.*exp(-1i*anh.*T3)); %need these too
%         P2=exp(-g(T1)-g(t2)-g(T3)+g(T1+t2)+g(t2+T3)-g(T1+t2+T3)).*(2*mu01_2-mu12_2.*exp(-1i*anh.*T3));
%   %HGS - Red Peak Only      
%  P9r=exp(-g(T1)+g(t2)-g(T3)-g(T1+t2)-g(t2+T3)+g(T1+t2+T3)).*(-mu12_2.*exp(-1i*anh.*T3)); %need these too
%         P9n=exp(-g(T1)-g(t2)-g(T3)+g(T1+t2)+g(t2+T3)-g(T1+t2+T3)).*(-mu12_2.*exp(-1i*anh.*T3));

%     %bend/stretch diag peak
%     P3 = P1.*exp(1i*w_nu2.*(-T1+T3));
%     P4 = P2.*exp(1i*w_nu2.*(T1+T3));

        %coupled bend/stretch diagonal peak blue
     Brd = Br.*exp(1i*w_nu2.*(-T1+T3));
     Bnd = Bn.*exp(1i*w_nu2.*(T1+T3));
        %coupled bend/stretch diagonal peak red
     Rrd = Rr.*exp(1i*w_nu2.*(-T1+T3));
     Rnd = Rn.*exp(1i*w_nu2.*(T1+T3));
     
%     %bend/stretch cp1
%     P5 = P1.*exp(1i*w_nu2.*(T3));
%     P6 = P2.*exp(1i*w_nu2.*(T3));

        %bend/stretch w1 -> w3 peak blue (Blue rephasing/nonrephasing righthand)
     Brr = Br.*exp(1i*w_nu2.*(T3));
     Bnr = Bn.*exp(1i*w_nu2.*(T3));
        %bend/stretch w1 --> w3 peak red (Red rephasing/nonrephasing righthand)
     Rrr = Rr.*exp(1i*w_nu2.*(T3));
     Rnr = Rn.*exp(1i*w_nu2.*(T3));
     
%     %bend/stretch cp2
%     P7 = P1.*exp(1i*w_nu2.*(-T1));
%     P8 = P2.*exp(1i*w_nu2.*(T1));

        %bend/stretch w3-->w1 peak blue (Blue rephasing/nonrephasing lefthand)
     Brl = Br.*exp(1i*w_nu2.*(-T1));
     Bnl = Bn.*exp(1i*w_nu2.*(T1));
        %bend/stretch w3 --> w1 peak red (Red rephasing/nonrephasing lefthand)
     Rrl = Rr.*exp(1i*w_nu2.*(-T1));
     Rnl = Rn.*exp(1i*w_nu2.*(T1));
     
%     %HGS - needs to be shifted above the original ESA peak so -T3 (idk it
%     %worked)
%     P9r = P9r.*exp(1i*w_nu2.*(-T3));
%     P9n = P9n.*exp(1i*w_nu2.*(-T3));
     
        %HGS
     Hr = Rr.*exp(1i*w_nu2.*(-T3));
     Hn = Rn.*exp(1i*w_nu2.*(-T3));
%Adding back in phase information
        Br = exp(1i*phi).*Br;
        Bn = exp(-1i*phi).*Bn;
        Rr = exp(1i*phi).*Rr;
        Rn = exp(-1i*phi).*Rn;
        Brd = exp(1i*phi).*Brd;
        Bnd = exp(-1i*phi).*Bnd;
        Rrd = exp(1i*phi).*Rrd;
        Rnd = exp(-1i*phi).*Rnd;
        Brr = exp(1i*phi).*Brr;
        Bnr = exp(-1i*phi).*Bnr;
        Rrr = exp(1i*phi).*Rrr;
        Rnr = exp(-1i*phi).*Rnr;
        Brl = exp(1i*phi).*Brl;
        Bnl = exp(-1i*phi).*Bnl;
        Rrl = exp(1i*phi).*Rrl;
        Rnl = exp(-1i*phi).*Rnl;
        Hr = exp(1i*phi).*Hr;
        Hn = exp(-1i*phi).*Hn;
        
    if flag_rotating_frame == false
        P1 = exp(1i*w_0.*(-T1+T3)).*P1;
        P2 = exp(1i*w_0.*(T1+T3)).*P2;
    end
      
    %do fft
%         P1=sgrsfft2(P1,n_zp);
%         P2=sgrsfft2(P2,n_zp);
% 
%         P3=sgrsfft2(P3,n_zp);
%         P4=sgrsfft2(P4,n_zp);
% %         
%         P5=sgrsfft2(P5,n_zp);
%         P6=sgrsfft2(P6,n_zp);
%         
%         P7=sgrsfft2(P7,n_zp);
%         P8=sgrsfft2(P8,n_zp);
%        
%         P9r=sgrsfft2(P9r,n_zp);
%         P9n=sgrsfft2(P9n,n_zp);
       
% new broken up peaks
        Br=sgrsfft2(Br,n_zp);
        Bn=sgrsfft2(Bn,n_zp);

        Rr=sgrsfft2(Rr,n_zp);
        Rn=sgrsfft2(Rn,n_zp);
        
        Brd=sgrsfft2(Brd,n_zp);
        Bnd=sgrsfft2(Bnd,n_zp);
       
        Rrd=sgrsfft2(Rrd,n_zp);
        Rnd=sgrsfft2(Rnd,n_zp);
        
        Brl=sgrsfft2(Brl,n_zp);
        Bnl=sgrsfft2(Bnl,n_zp);
        
        Rrl=sgrsfft2(Rrl,n_zp);
        Rnl=sgrsfft2(Rnl,n_zp);
        
        Brr=sgrsfft2(Brr,n_zp);
        Bnr=sgrsfft2(Bnr,n_zp);
       
        Rrr=sgrsfft2(Rrr,n_zp);
        Rnr=sgrsfft2(Rnr,n_zp);
        
        Hr=sgrsfft2(Hr,n_zp);
        Hn=sgrsfft2(Hn,n_zp);


    %change how data packaged
%     P=-fftshift(real(fliplr(circshift(P1,[0 -1]))+P2));
%     P = P - fftshift(real(fliplr(circshift(P3,[0 -1]))+P4));
%     P = P - fftshift(real(fliplr(circshift(P5,[0 -1]))+P6));
%     P = P - fftshift(real(fliplr(circshift(P7,[0 -1]))+P8));
%     P = P - fftshift(real(fliplr(circshift(P9r,[0 -1]))+P9n));

    k1 = p(end-2);
    k2 = p(end-1);
    k3 = p(end);
    
    
    %added 1/4/19 - prevents code from breaking if kinetics are not present
    
if isfield(options,'kin')
    
    a1 = options.kin{1};
    a2 = options.kin{2};
    a3 = options.kin{3};
    a4 = options.kin{4};
    a5 = options.kin{5};
    a6 = options.kin{6};
    a7 = options.kin{7};
    a8 = options.kin{8};
end

if isfield(options,'kin')
    
    P=-fftshift(real(fliplr(circshift(Br,[0 -1]))+Bn))*double(a1(k2,k3)*a3(k1,k2,k3,t2)); %coeff between - & fft
    P = P - fftshift(real(fliplr(circshift(Rr,[0 -1]))+Rn))*double(a1(k2,k3)*a3(k1,k2,k3,t2)*a8(k1,k2,k3,t2)); %or here ;

else
   
     P=-fftshift(real(fliplr(circshift(Br,[0 -1]))+Bn)); %coeff between - & fft
    P = P - fftshift(real(fliplr(circshift(Rr,[0 -1]))+Rn)); %or here ;
end    

    [W1,W3] = meshgrid(w,w);
    P = P./abs(min(P(:))); %normalize to the 01 band (negative)
    out(:,:,i) = interp2(W1,W3,P,w1_in,w3_in','*linear');
    
  end %end t2_array loop
  
  if flag_bootstrap
      out = out(bootstrap_index);
  end
  
  try
      s.c2 = c2(t);
      s.g = g(t);
      s.t = t;
      extra = s;
  end
end %end if order 3


%WORKS as of 12/6/18
