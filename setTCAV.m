function S = setTCAV(Efocus, Vrf)
% For nominal 2 bunch:  E_witness = 9.925, E_drive = 10.1
%                       E_both ~ 10.0

global BEAMLINE  PS

if nargin == 0
    Efocus = 10;
    Vrf = 20;
    disp('Setting defaults - Vrf = 20 MeV and Efocus = 10 GeV');
elseif nargin == 1
    Vrf = 20;
    disp('Setting Vrf = 20 MeV');    
end

xtc=findcells(BEAMLINE,'Name','XTCAVF');   % tcav, split into 2

% Opt settings 2/6/19
quadPS = [ 362   363   364   367   368   369   370   371   372];
PSvals = [ -7.4420  -12.2894   13.7426   17.4449  -27.6751    7.6582  -11.0168   17.3447  -11.0178];


% Set and scale quads to energy
for ii=1:length(quadPS)
    % Ensure B = 0.5 in in each quad
    for jj = PS(quadPS(ii)).Element
        BEAMLINE{jj}.B = 0.5;
    end
    PS(quadPS(ii)).Ampl = PSvals(ii)*Efocus/10;
end

% Turn tcav into an LCAV
for i = 1:length(xtc)
  iele = xtc(i);
  BEAMLINE{iele}.Class='TCAV';
  BEAMLINE{iele}.Volt=Vrf/length(xtc);
  BEAMLINE{iele}.Phase=90;
  BEAMLINE{iele}.Tilt=0;
end


% Streak stength
S = 13.4004;

1


end
