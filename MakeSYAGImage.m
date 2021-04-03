function [Xedges,Yedges,N,Xvals,Yvals,dgamma] = MakeSYAGImage(beam,xbins,ybins,makeplot)


% Makes a LPS image from beam structure
% Returns energy projection and current
% Dispersion at SYAG is 5.5 cm

% Cemma 02-25-2019
id=find(beam.Bunch.stop==0);
E0 = mean(beam.Bunch.x(6,:))*1000; %MeV

% Set pixel size and image dimensions
px = 3.25/1000;
Xedges = xbins;
Yedges = ybins;

% Condition beam structure
id=find(beam.Bunch.stop==0);
rays=beam.Bunch.x(:,id)';
Q=beam.Bunch.Q(id);
x=rays(:,1);
y=rays(:,3);
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

% Calculate current and energy profile
Xvals = sum(N,2)*1e12;% Q in C
Yvals = sum(N,1)*1e12;% in pC
%% Plot stuff
if makeplot
figure
subplot(2,2,1)
imagesc(xb*1e6,yb*1e6,N')
set(gca,'YDir','normal')
xlabel('x [$\mu$m]','interpreter','latex') ; 
ylabel('y [$\mu$m]','interpreter','latex') ;
enhance_plot
subplot(2,2,2)
ax = gca;
plot(yb*1e6,Yvals,'LineWidth',2);
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
xlabel('y [$\mu$ m]','interpreter','latex','FontSize',20,'FontName','Times') ;
ylabel('Q [pC]','interpreter','latex','FontSize',20,'FontName','Times') ;
view([-90,90])
xlim([yb(1),yb(end)]*1e6)
subplot(2,2,3)
plot(xb*1e6,Xvals)
enhance_plot
xlabel('x [$\mu$m]','interpreter','latex') ; 
ylabel('Q [pC]','interpreter','latex') ;
xlim([xb(1),xb(end)]*1e6)

dispersion = 8.4e-2;% from Twiss.etax(SYAG3)
dgamma = xb/dispersion;
dgamma = (dgamma-median(dgamma))*100;

figure
subplot(2,1,1)
imagesc(Xedges*1e3,yb*1e3,N')
set(gca,'YDir','normal')
%xlabel('\Delta E [%]','interpreter','tex') ; 
xlabel('x [mm]','interpreter','tex') ; 
ylabel('y [mm]','interpreter','latex') ;
enhance_plot
xlim([Xedges(1),Xedges(end)]*1e3)

subplot(2,1,2)
plot(dgamma,Xvals)
enhance_plot
xlabel('\Delta E [%]','interpreter','tex') ; 
ylabel('Q [pC]','interpreter','latex') ;

xlim([dgamma(1),dgamma(end)])

% Calculate energy separation and charge ratio
% Define the 'drive beam' as being everything above the spike + 1 FWHM 
% And the witness beam below 1 fwhm
deltaE = round(fwhm(dgamma,Xvals)/(dgamma(end)-dgamma(end-1)));
deltaE = 0;% Otherwise just use the max as the cutoff between dr/wit bunches

[maxi,ind] = max(Xvals);
Qw = trapz(Xvals(1:ind));

Qd = trapz(Xvals(ind+1:end));

[maxi,ind]=max(Xvals);
wb = Xvals(1:ind);
db = Xvals(ind:end);
wb_com = round(sum(wb'.*(1:length(wb)))/sum(wb));
db_com = ind+round(sum(db'.*(1:length(db)))/sum(db));

energysep = (db_com-wb_com)*(dgamma(end)-dgamma(end-1));

disp(wb_com)
disp(db_com)
disp(energysep)
chargeratio = Qd/Qw;
%disp(Qd*1e-12)
%disp(Qw*1e-12)
%disp(chargeratio)
end
