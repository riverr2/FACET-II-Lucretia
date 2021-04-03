function [] = LPSCurrent(new_data, ming, zbound, zsteps, name)

    new_data.Bunch.x(5,:) = new_data.Bunch.x(5,:) - min(new_data.Bunch.x(5,:));
    idx = new_data.Bunch.stop==0;%(new_data.Bunch.x(5,:)<zbound)&(new_data.Bunch.x(6,:)>ming);
    new_data.Bunch.x = new_data.Bunch.x(:,idx);
    new_data.Bunch.Q = new_data.Bunch.Q(idx);
    new_data.Bunch.x(5,:) = new_data.Bunch.x(5,:) - min(new_data.Bunch.x(5,:));
    zmin = min(new_data.Bunch.x(5,:));
    zmax = max(new_data.Bunch.x(5,:));
    dz = (zmax - zmin) / zsteps;
    emits = [];
    curs = [];
    for i = 1:zsteps
       idx = (new_data.Bunch.x(5,:)<zmin+(i+1)*dz)&(new_data.Bunch.x(5,:)>=zmin+i*dz);
       xrms = std(new_data.Bunch.x(1,idx));
       xprms = std(new_data.Bunch.x(2,idx));
       xmean = mean(new_data.Bunch.x(1,idx));
       xpmean = mean(new_data.Bunch.x(2,idx));
       gavg = mean(new_data.Bunch.x(6,idx)/0.000511);
       emits(i) = gavg * sqrt(xrms^2*xprms^2 - mean((new_data.Bunch.x(1,idx)-xmean).*(new_data.Bunch.x(2,idx)-xpmean))^2);
       curs(i) = sum(new_data.Bunch.Q(idx)) / (dz/3e8);
    end
    zs = (1:zsteps)*dz;

    figure();
    
    subplot(2,1,1);
    
    plot(1e6*zs, curs)
    ylabel('Current (A)')
    xlabel('Longitudinal Position (um)')
    %xlim([390, 425]);
    %ylim([0, 3e3]);
    subplot(2,1,2);
    
    plot(1e6*new_data.Bunch.x(5,:), new_data.Bunch.x(6,:), '.')
    ylabel('Energy (GeV)')
    xlabel('Longitudinal Position (um)')
    %xlim([390, 425]);
    sgtitle(name);
end

