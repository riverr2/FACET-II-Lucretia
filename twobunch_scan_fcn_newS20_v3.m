% Parallel scan function Lucretia
function [data]=twobunch_scan_fcn_newS20_v3(scanparams)

global BEAMLINE PS KLYSTRON WF

load 2bunch
load beam_10_fac.mat
Beam = end_inj;

%istart = 45;

%de = 1.2;
%de = 1.7; % for use with gpt2
%BEAMLINE{68}.Volt = 7.6273 - de/cos(-2.5*pi/180);%5.9256;
%BEAMLINE{68}.Egain = 7.62 - de;
%BEAMLINE{68}.Phase = - acos((7.62 - de) / 7.6273)*180/pi;

S = setTCAV(10,20);
% Physics effects on/off control
ISRon=1; % 1= Incoherent Sync. Rad. on, 0= ISR off
CSRon=0; % 1= Coherent Sync. Rad. on, 0= CSR off
LSCon=1; % 1= Longitudinal space charge on, 0= LSC off
% Decimate beam (input beam from IMPACT-T has 2e5 macro-particles)
decBeam=20; % use integer>1 to sample 1:decBeam:end of beam

% Set initial beam parameters

Initial.Q=Initial.Q*(1+scanparams.qi); % Initial charge offset
Beam.Bunch.x(1,:)=Beam.Bunch.x(1,:)+scanparams.dx;% Initial offset in x
Beam.Bunch.x(3,:)=Beam.Bunch.x(3,:)+scanparams.dy;% Initial offset in y

% --- Set Linac Phases / degrees off-crest and amplitudes

idfbL2_1=findcells(BEAMLINE,'Name','K14_4*'); % don't touch E feedback stations
idfbL2_2=findcells(BEAMLINE,'Name','K14_5*'); 
idfbL3_1=[findcells(BEAMLINE,'Name','K19_1*') findcells(BEAMLINE,'Name','K19_2*') findcells(BEAMLINE,'Name','K19_3*')];
idfbL3_2=[findcells(BEAMLINE,'Name','K19_4*') findcells(BEAMLINE,'Name','K19_5*') findcells(BEAMLINE,'Name','K19_6*')];
for ip=[findcells(BEAMLINE,'Name','BEGL1F'):findcells(BEAMLINE,'Name','L1XFBEG') ...
    findcells(BEAMLINE,'Name','L1XFEND'):findcells(BEAMLINE,'Name','BEGBC11')]
  if isfield(BEAMLINE{ip},'Phase')      
    BEAMLINE{ip}.Phase=scanparams.P1;
   
  end
    if isfield(BEAMLINE{ip},'Volt')      
    BEAMLINE{ip}.Volt=BEAMLINE{ip}.Volt*cos(-19.2*pi/180)/cos(scanparams.P1*pi/180);%*(1+scanparams.V1);% Relative offset
  end
end
for ip=findcells(BEAMLINE,'Name','BEGL2F'):findcells(BEAMLINE,'Name','BEGBC14')
  if isfield(BEAMLINE{ip},'Phase') && ~ismember(ip,[idfbL2_1 idfbL2_2])% Change the linac phase for all but the E feedback stations
    BEAMLINE{ip}.Phase=scanparams.P2;  
  end
  if isfield(BEAMLINE{ip},'Volt') && ~ismember(ip,[idfbL2_1 idfbL2_2])% Change the linac phase for all but the E feedback stations
    BEAMLINE{ip}.Volt=BEAMLINE{ip}.Volt*cos(-38.35*pi/180)/cos(scanparams.P2*pi/180);%(1+scanparams.V2);% Relative offset  
  end
end
for ip=findcells(BEAMLINE,'Name','BEGL3F'):findcells(BEAMLINE,'Name','BEGBC20')%findcells(BEAMLINE,'Name','FBEG')
  if isfield(BEAMLINE{ip},'Phase') && ~ismember(ip,[idfbL3_1 idfbL3_2])% Change the linac phase for all but the E feedback stations
      BEAMLINE{ip}.Phase=scanparams.P3;% Don't scan L3 for now
  end
  if isfield(BEAMLINE{ip},'Volt') && ~ismember(ip,[idfbL3_1 idfbL3_2])
    BEAMLINE{ip}.Volt=BEAMLINE{ip}.Volt*cos(45*pi/180)/cos(scanparams.P3*pi/180);
  end
end


% Collective effects
%-- Sync. Radiation
eleSR=findcells(BEAMLINE,'Class','SBEN');
for iele=eleSR
  BEAMLINE{iele}.TrackFlag.SynRad=2*ISRon;
  BEAMLINE{iele}.TrackFlag.CSR=-1*CSRon;
  BEAMLINE{iele}.TrackFlag.CSR_SmoothFactor=0;
  BEAMLINE{iele}.TrackFlag.CSR_DriftSplit=25;
  BEAMLINE{iele}.TrackFlag.Split=25;
  BEAMLINE{iele}.TrackFlag.Diagnostics=0;
  % Turn off 2d CSR
        BEAMLINE{iele}.TrackFlag.CSR_2D = 0;
        BEAMLINE{iele}.TrackFlag.CSR_USEGPU = 0;
end

% -- Space charge
for iele=findcells(BEAMLINE,'TrackFlag')
  BEAMLINE{iele}.TrackFlag.LSC=LSCon;
  BEAMLINE{iele}.TrackFlag.LSC_storeData=0;
  % Set NBPM on LCAV elements to ensure 0.1m drift sections for
  % application of LSC
  if strcmp(BEAMLINE{iele}.Class,'LCAV')
    BEAMLINE{iele}.NBPM=LSCon*BEAMLINE{iele}.L/0.1;
    BEAMLINE{iele}.GetSBPMData=LSCon;
    BEAMLINE{iele}.GetInstData=LSCon;
  end
end
%
% -- Wakefields 
SetTrackFlags('SRWF_Z',1,1,findcells(BEAMLINE,'Name','BEGBC20'));
SetTrackFlags('SRWF_T',0,1,findcells(BEAMLINE,'Name','BEGBC20'));
%% -- Tracking
% Define indices
iLH=findcells(BEAMLINE,'Name','HTRUNDF');
iL1A=findcells(BEAMLINE,'Name','BEGL1F');
SYAG1=findcells(BEAMLINE,'Name','BPM11333');% Location of fake syag in the middle of BC11
iBC1B=findcells(BEAMLINE,'Name','ENDBC11'); % start of L2
SYAG2=findcells(BEAMLINE,'Name','BPM14801');% Location of fake syag in the middle of BC14
iBC2B=findcells(BEAMLINE,'Name','BEGL3F'); % start of L3
BC20BEG =findcells(BEAMLINE,'Name','BEGBC20');
BC20END =findcells(BEAMLINE,'Name','ENDBC20');
SYAG3=findcells(BEAMLINE,'Name','SYAG');% Location of SYAG in BC20
iFF=findcells(BEAMLINE,'Name','MFFF'); % start of FFS
ip=findcells(BEAMLINE,'Name','PENT'); % IP
pmdmp=findcells(BEAMLINE,'Name','PDUMP'); % Find spectrometer screen
xtc = findcells(BEAMLINE,'Name','XTCAVF');

%TwissPlot(1,pmdmp,Initial,[1 1 0]);
%return
% Dispersion at SYAG is 8.4 cm
% Decimate beam?
if decBeam>1
  decBeam=floor(decBeam);
  Beam.Bunch.x=Beam.Bunch.x(:,1:decBeam:end);
  Beam.Bunch.stop=Beam.Bunch.stop(1:decBeam:end);
  Beam.Bunch.Q=Beam.Bunch.Q(1:decBeam:end).*decBeam;
end

% Track to end of each bending section and re-center beam in each case to
% take care of phase-slip with respect to RF and orbit excursions due to SR
% energy losses. In reality this is done with RF phasing and orbit feedbacks.
%LH=0; % laser heater energy spread / keV
%T=Track(Beam);
% The code below tracks without changing the beta function at the TCAV
%beam=T.beamOut;
%beam.Bunch.x(6,~beam.Bunch.stop)=beam.Bunch.x(6,~beam.Bunch.stop)+randn(size(beam.Bunch.x(6,~beam.Bunch.stop))).*LH.*1e-6;
% Get the phase advance between tcav and screen
% In=TwissToInitial(Twiss,BC20END,Initial); 
% [~,Twiss]=GetTwiss(BC20END,length(BEAMLINE),In.x.Twiss,In.y.Twiss);
% pa = (Twiss.nux(pmdmp-BC20END+1) - Twiss.nux(xtc(1)-BC20END+1) )*2*pi;
% disp(['Phase adv = ',num2str(pa/pi),'pi'])
%% Tracking from LH to IP in steps dumping data after each BC
beam = Beam;
%indici = [istart iLH 189 190 iBC1B SYAG2 iBC2B SYAG3 BC20END,ip,pmdmp];
indici = [istart iLH iBC1B SYAG2 iBC2B SYAG3 BC20END,ip,pmdmp];

%zloc = ["iLH", "prel1", "postl1","BC11end","SYAG14","BC14end","SYAG20","BC20END","IP","Phosphor"];
zloc = ["iLH", "BC11end","SYAG14","BC14end","SYAG20","BC20END","IP","Phosphor"];
data = struct();
data.scanparams = scanparams;
for i=1:length(indici)-1
disp(zloc(i));
disp(indici(i));
T = Track(beam);
T.startInd=indici(i);
% These two commented lines were the difference btween my sims and Doug's -
% this makes all the difference to the output LPS
T.centerZInd=[iL1A iBC1B iBC2B iFF];
T.centerTInd=[iL1A iBC1B iBC2B iFF];
T.centerZInd = [BC20END];
T.finishInd=indici(i+1);
T.trackThru();
beam=T.beamOut;
id = find(beam.Bunch.stop==0);

data(i).beamparams = beamImage_noplot(beam,0);
[data(i).beamparams.nx,data(i).beamparams.ny] = GetNEmitFromBeam(beam,1) ;

data(i).beamparams.centroidx = mean(beam.Bunch.x(1,id));% X centroid in the chicanes is a proxy for the energy
data(i).zpos = zloc(i);
% Dump current after BC11 and BC14
switch zloc(i)
    case "iLH"        
        % Add laser heater energy spread        
        LH=0; % laser heater energy spread / keV
    beam.Bunch.x(6,~beam.Bunch.stop)=beam.Bunch.x(6,~beam.Bunch.stop)+randn(size(beam.Bunch.x(6,~beam.Bunch.stop))).*LH.*1e-6;
    
    case "soon dead"        
        % Add laser heater energy spread        
        LH=0; % laser heater energy spread / keV
        data(i).beam = beam;
                   
    case "dead"        
        % Add laser heater energy spread        
        LH=0; % laser heater energy spread / keV
        data(i).beam = beam;
     
    case "BC11end"        
        % Calculate current
        ebins = 0.335*[(1-5e-2):1e-3:(1+4e-2)];% GeV
        zbins = ([-750:5:2000]*1e-6);% m
        %[~,~,~,data(i).I,data(i).Eprof]=MakeBeamLPS(beam,zbins,ebins,0);     
        nbins = 500;
        data(i).beam = beam;
        %beamImage(beam,0,mean(beam.Bunch.x(6,id)),'false',nbins);
        
    case "BC14end"                
        ebins = 4.48*[(1-5e-2):1e-3:(1+7e-2)];% GeV
        zbins = ([-250:1:400]*1e-6);% m
        %[~,~,~,data(i).I,data(i).Eprof]=MakeBeamLPS(beam,zbins,ebins,1);  
        data(i).beam = beam;
        
    case "SYAG20"      
        xbins = -.25e-3:0.01e-3:3.25e-3;ybins = -1e-3:0.01e-3:1e-3;
        %[~,~,~,data(i).Xvals,~] =MakeSYAGImage(beam,xbins,ybins,1);
        data(i).beam = beam;        
        % Calculate current         
        ebins = [9.8:0.00046:10.4];% GeV
        zbins = ([-250:0.25:200]*1e-6);% m
        %[~,~,~,data(i).I,data(i).Eprof]=MakeBeamLPS(beam,zbins,ebins,1);    
        %beamImage(beam,0,mean(beam.Bunch.x(6,id)),'false');
                
    case "BC20END"
        % Re-center manually in case the T.centerZInd doesn't work  
        beam.Bunch.x(5,id) = beam.Bunch.x(5,id)-mean(beam.Bunch.x(5,id));     
                nbins = 500;
          % Calculate current         
        ebins = [9.8:0.00046:10.4];% GeV
        %zbins = ([-250:0.25:150]*1e-6);% m
        zbins = ([-350:1:250]*1e-6);% m
        %[~,~,~,data(i).I,data(i).Eprof]=MakeBeamLPS(beam,zbins,ebins,1);                       
        data(i).beam = beam;    
        
    case "IP"
        % Calculate current         
        ebins = [9.8:0.00046:10.4];% GeV
        zbins = ([-250:0.25:150]*1e-6);% m
        data(i).dz = zbins(2)-zbins(1);
        %[~,~,~,data(i).I,data(i).Eprof]=MakeBeamLPS(beam,zbins,ebins,1);                     
        %beamImage(beam,0,mean(beam.Bunch.x(6,id)),'false');
end

%{
    if contains(zloc(i),"BC11end") 
        nbins = 500;
        % Calculate current
        ebins = 0.335*[(1-5e-2):1e-3:(1+4e-2)];% GeV
        zbins = ([-750:5:2000]*1e-6);% m
        [data(i).Zbins,data(i).Ebins,N,data(i).I,data(i).Eprof]=MakeBeamLPS(beam,zbins,ebins,1);     
        %beamImage(beam,0,mean(beam.Bunch.x(6,id)),'false',nbins);
    end
    
    if contains(zloc(i),"BC14end")
        nbins = 500;
        ebins = 4.48*[(1-5e-2):1e-3:(1+7e-2)];% GeV
        zbins = ([-250:1:400]*1e-6);% m
        [data(i).Zbins,data(i).Ebins,N,data(i).I,data(i).Eprof]=MakeBeamLPS(beam,zbins,ebins,1);               
    end
    
    if contains(zloc(i),"SYAG20")      
        %beamImage(beam,0,mean(beam.Bunch.x(6,id)),'false',500);
        data(i).beam = beam;
    end
    
    if contains(zloc(i),"BC20END");   beam.Bunch.x(5,id) = beam.Bunch.x(5,id)-mean(beam.Bunch.x(5,id)); end    % Re-center manually in case the T.centerZInd doesn't work  
    
    if contains(zloc(i),"IP")        
        % Calculate current         
        ebins = [9.8:0.00046:10.4];% GeV
        zbins = ([-250:0.25:150]*1e-6);% m
        [data(i).Zbins,data(i).Ebins,N,data(i).I,data(i).Eprof]=MakeBeamLPS(beam,zbins,ebins,1);                        
    end
%}
end

ndumps = length(indici)-1;
[data(ndumps).Xedges,data(ndumps).Yedges,data(ndumps).N_cam]=MakeImage(beam,'OTR');

% Dump the beam structure at phosphor to compare with what you get from MakeImage
data(ndumps).beam = beam;

% Decimate the XTCAV image because it's too large
scalefactor = 4;
data(ndumps).N_cam = data(ndumps).N_cam(1:scalefactor:end,1:scalefactor:end);
data(end).Xedges = data(end).Xedges(1:scalefactor:end);
data(end).Yedges = data(end).Yedges(1:scalefactor:end);

% Calculate current from XTCAV image 
%S = 13.4004;%um/um - get this from Dougs_beamline_2bunch.m script OR setTCAV
N_sum = sum(data(ndumps).N_cam,2);
id = find(beam.Bunch.stop==0);
charge = sum(beam.Bunch.Q(id));
zaxis = data(ndumps).Xedges*1e3/S;% in um

N_int = trapz(zaxis*1e-6/3e8,N_sum);
conv_factor = charge/N_int;% From N electron to I [kA]
data(ndumps).xtcav_current = N_sum*conv_factor*1e-3;
data(ndumps).xtcav_zaxis = zaxis;

% Crop the LPS image if you want before storing the output
% xmin = 2000; xmax = 4000;
% ymin = 400; ymax = 1750;
% 
% data(ndumps).Xedges = data(ndumps).Xedges(1,xmin:xmax);
% data(ndumps).Yedges = data(ndumps).Yedges(1,ymin:ymax);
% data(ndumps).N_cam = data(ndumps).N_cam(xmin:xmax,ymin:ymax);



% Scan beta at IP and screen
%{ 
     betaIP = 1;
     betascreen = 1;
     tcavdata = scan_tcav_betas(betaIP,betascreen,beambc20end);  
%     data.x(1,:) = tcavdata.Beam.Bunch.x(1,:);
%     data.x(2,:) = tcavdata.Beam.Bunch.x(2,:);
%     data.x(5,:) = tcavdata.Beam.Bunch.x(5,:);
%     data.x(6,:) = tcavdata.Beam.Bunch.x(6,:);
%     MakeImage(beam,'OTR');
%}
