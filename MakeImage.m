function [Xedges,Yedges,N_cam] = MakeImage(beam, imtype);
set(0,'defaulttextinterpreter','latex')
% Makes an OTR or YAG image from beam structure
% Storey 8/8/18
% imtype = 'YAG' or 'OTR', or left blank defaults to YAG

if nargin > 1
    if  strcmpi(imtype, 'YAG')
        imtype = 1;
    elseif strcmpi(imtype, 'OTR')
        imtype = 0;
    else
        warning('Unknown image type')
        return
    end
else
    imtype = 0;  % Default is OTR image
end


E = 10000; %MeV
%E = mean(beam.Bunch.x(6,:))*0.511;%MeV

% Set pixel size and image dimensions
px = 3.25/1000;

Xedges = -3.5:px:3;

Yedges = -2.5:px:2.5;

% Condition beam structure
id=find(beam.Bunch.stop==0);
rays=beam.Bunch.x(:,id)';
Q=beam.Bunch.Q(id);
x=1000*rays(:,1);
y=1000*rays(:,3);
  idx=true(1,length(x))';
  idy=true(1,length(x))';

[~,u,BIN]=histcounts(y(idy), Yedges);
u=u(1:end-1); npixsqrY=length(u);
idy(BIN==0)=false; BIN(BIN==0)=[];

[~,u,BIN]=histcounts(x(idx), Xedges);
u=u(1:end-1); npixsqrX=length(u);
idx(BIN==0)=false; BIN(BIN==0)=[];

id = idx&idy;

% Make 2D hist
xb=linspace(min(Xedges),max(Xedges),npixsqrX);
yb=linspace(min(Yedges),max(Yedges),npixsqrY);

xr=interp1(xb,1:npixsqrX,x(id),'nearest');
yr=interp1(yb,1:npixsqrY,y(id),'nearest');
N = accumarray([xr(:) yr(:)], Q(id), [npixsqrX npixsqrY]);

Xedges = xb; Yedges = yb;

%% Rescale bunch charge 
%  N = N/sum(sum(N))*0.5e-9

% Estimate the real LPS in similar way
%%% [N,Xedges,Yedges] = histcounts2(beam.Bunch.x(5,~beam.Bunch.stop)*1000*15,(beam.Bunch.x(6,~beam.Bunch.stop)-10)/.6*4, -3.5:px:3.5,-4.15:px:4.15);


%% Test case
% %             Lmat = 1000;
% %             sigx = 4.69e-6*1000;  sigy = 255e-6*1000;
% %             N = zeros(Lmat*2-1);
% %             for ii = 1:Lmat*2-1
% %                 for jj = 1:Lmat*2-1
% %                     N(ii,jj) = exp( -(ii-Lmat)^2*px^2/2/sigx^2 - (jj-Lmat)^2*px^2/2/sigy^2 );
% %                 end
% %             end
% %             N = N/sum(sum(N))*2e-9;


%% Estimate light output
if imtype
    %% Do YAG screen
    sigYAG = 30e-6;
    LO = 8000/4/pi;  % Light output of YAG in gamma/MeV/st
    Edep = 2.1*4.56*.002; % 2.1 MeV.cm/g, 4.56 g/cm3, 20µm,  PER ELECTRON
    Omega = 0.046264; % st
    imname = 'YAG '; fprintf('Making YAG image\n')
    
    Q_dens = max(max(N))*1e9/px^2;
    if Q_dens > 10
        fprintf('Caution: Peak charge density %.0f nC/mm^2 could saturate YAG (>10 nC/mm^2)\n', Q_dens)
    end
    
    N_photons = N/1.602e-19*Edep*LO*Omega;
    
else
    %% Do OTR foil
    sigYAG = 2e-6;
    alpha = 1/137;
    R = 0.65;
    imname = 'OTR '; fprintf('Making OTR image\n')
    
    % Assum all photons collected
    N_photons = N/1.602e-19 * alpha/pi*2*log(E/0.511-1);
   
end

%% Convolve with guassians to estimate bluring from YAG and aberations/diffraction

sigY = 1*sigYAG*1000;
Lmat = ceil(4*sigY/px);
Myag = zeros(Lmat*2-1);

for ii = 1:Lmat*2-1
    for jj = 1:Lmat*2-1
        Myag(ii,jj) = exp( -(ii-Lmat)^2*px^2/2/sigY^2 - (jj-Lmat)^2*px^2/2/sigY^2 );
    end
end
% Normalize
Myag = Myag/sum(sum(Myag));

% Diffraction resolution
sigDiff = 4/1000;
Lmat = ceil(4*sigDiff/px);
MDiff = zeros(Lmat*2-1);

for ii = 1:Lmat*2-1
    for jj = 1:Lmat*2-1
        MDiff(ii,jj) = exp( -(ii-Lmat)^2*px^2/2/sigDiff^2 - (jj-Lmat)^2*px^2/2/sigDiff^2 );
    end
end
MDiff = MDiff/sum(sum(MDiff));

% Convolve with YAG res
N_YAG = conv2(N_photons,Myag, 'same');

% Convolve with Diffraction res
N_phot_cam = conv2(N_YAG,MDiff, 'same');

% Photons to electrons in CCD
QE = 0.4; % Quantum efficienc
FullWell = 30000; % Saturation intensity

N_cam = N_phot_cam*QE/FullWell;

% Change the x and y values on the screen to z E values
dispersion = 6.96e-2;% From Twiss.etay(1480)
dgamma = Yedges*1e-3/dispersion;
dgamma = (dgamma-median(dgamma))*100;
%% Make plots if you want

figure(3)
subplot(2,2,1)
imagesc(Xedges,Yedges,N_cam')
colorbar
set(gca,'YDir','normal','LineWidth',2,'FontSize',16)
axis tight
xlabel('x position [mm]','interpreter','latex') ; 
ylabel('y position [mm]','interpreter','latex') ;
title( [imname 'Image On Camera $\sigma_{screen}$ = ' num2str(sigY*1000) ' $\mu m$, $\sigma_{optics}$ = ' num2str(sigDiff*1000) ' $\mu m$'])
%xlim([-4,4])
% Plot projections
hold on
plot(Xedges, sum(N_cam,2)*max(ylim)*0.4/max(sum(N_cam,2))+min(ylim), 'w', 'LineWidth', 2)
plot(sum(N_cam,1)*max(xlim)*0.4/max(sum(N_cam,1))+min(xlim), Yedges, 'w', 'LineWidth', 2)
hold off
subplot(2,2,3)
plot(dgamma,sum(N_cam,1))
xlabel('dE/E [\%]')
subplot(2,2,2)
plot(Xedges,sum(N_cam,2))
xlabel('x [mm]')

%}


