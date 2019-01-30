<<<<<<< HEAD
function s=construct2dPP(varargin)

switch nargin
    case 0
        s = struct('freq',[],...
            'time',[],...
            't2',[],...
            't3',0,...
            'w1',[],...
            'w3',[],...
            'R',[],...
            'PP',[],...
            ...% 'R2',[],...
            ...%'hene_x',[],...
            ...%'hene_y',[],...
            'igram',[],...
            'phase',0,...
            't0_bin',[],...
            't_bin_shift',0,...
            'PP_noise',[],...
            ...%'R2_noise',[],...
            'basename',[],...
            'undersampling',0,...
            'centerfreq',[],...
            'resolution',[],...
            'zeropad',[],...
            'fft_type','sgrsfft',...
            'time_units','fs',...
            'freq_units','wavenumbers',...
            'spec_calib',[],...
            'pump_probe',[],...
            'pump_probe_freq',[],...
            'comment',[],...
            'time_stamp',[],...
            'PARAMS',[],...
            'noise',[],...
            'background',[],...
            'apod',struct('name',[],'params',[]));
    case 1
        if isa(varargin{1},'struct')
            s = varargin{1};
        else
            error('Construct2dPP wrong argument type');
        end
end
=======
function s=construct2dPP(varargin)

switch nargin
  case 0
    s = struct('freq',[],...
      'time',[],...
      't2',[],...
      't3',0,...
      'w1',[],...
      'w3',[],...
      'R',[],...
      'PP',[],...
      ...% 'R2',[],...
      ...%'hene_x',[],...
      ...%'hene_y',[],...
      'igram',[],...
      'phase',0,...
      't0_bin',[],...
      't_bin_shift',0,...
      'PP_noise',[],...
      ...%'R2_noise',[],...
      'basename',[],...
      'undersampling',0,...
      'centerfreq',[],...
      'resolution',[],...
      'zeropad',[],...
      'fft_type','sgrsfft',...
      'time_units','fs',...
      'freq_units','wavenumbers',...
      'spec_calib',[],...
      'pump_probe',[],...
      'pump_probe_freq',[],...
      'comment',[],...
      'time_stamp',[],...
      'PARAMS',[],...
      'noise',[],...
      'background',[],...
      'apod',struct('name',[],'params',[]));
  case 1
    if isa(varargin{1},'struct')
      s = varargin{1};
    else
      error('Construct2dPP wrong argument type');
    end
end
>>>>>>> 94a2492ef59acd231e6ae002a4046a19b8ffdcbd
