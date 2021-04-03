function [Xedges,Yedges,N,I,Eprof] = MakeBeamLPS(beam,zbins,ebins,makeplot)

% Makes a LPS image from beam structure
% Returns energy projection and current
% Cemma 02-25-2019
id=find(beam.Bunch.stop==0);
E0 = mean(beam.Bunch.x(6,:))*1000; %MeV

% Set pixel size and image dimensions
px = 3.25/1000;
Xedges = zbins;
Yedges = ebins;

% Condition beam structure
rays=beam.Bunch.x(:,id)';
Q=beam.Bunch.Q(id);
z=rays(:,5);
E=rays(:,6);
idz=true(1,length(z))';
idE=true(1,length(z))';

[~,u,BIN]=histcounts(E(idE), Yedges);
u=u(1:end-1); npixsqrY=length(u);
idE(BIN==0)=false; BIN(BIN==0)=[];

[~,u,BIN]=histcounts(z(idz), Xedges);
u=u(1:end-1); npixsqrX=length(u);
idz(BIN==0)=false; BIN(BIN==0)=[];

id = idz&idE;

% Make 2D hist
xb=linspace(min(Xedges),max(Xedges),npixsqrX);
yb=linspace(min(Yedges),max(Yedges),npixsqrY);

xr=interp1(xb,1:npixsqrX,z(id),'nearest');
yr=interp1(yb,1:npixsqrY,E(id),'nearest');
N = accumarray([xr(:) yr(:)], Q(id), [npixsqrX npixsqrY]);

Xedges = xb; Yedges = yb;

% Calculate current and energy profile
sliceq = sum(N,2);
I=1e-3.*sliceq.*(1/(((zbins(2)-zbins(1)))/299792458)); % y-axis Q->I (kA)
Eprof = sum(N,1)*1e12;% in pC
%% Plot stuff
if makeplot
figure
subplot(2,2,1)
imagesc(xb*1e6,yb,N')
%colormap jetvar
set(gca,'YDir','normal')
xlabel('z [$\mu$m]','interpreter','latex') ; 
ylabel('Energy [GeV]','interpreter','latex') ;
enhance_plot
%xlim([-200,200])

subplot(2,2,2)
ax = gca;
plot(yb,Eprof,'LineWidth',2);
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
xlabel('Energy [GeV]','interpreter','latex','FontSize',20,'FontName','Times') ;
ylabel('Q [pC]','interpreter','latex','FontSize',20,'FontName','Times') ;
view([-90,90])

subplot(2,2,3)
plot(xb*1e6,I)
enhance_plot
xlabel('z [$\mu$m]','interpreter','latex') ; 
ylabel('I [kA]','interpreter','latex') ;
%xlim([-200,200])
end
