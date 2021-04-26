clear all
close all 

load 2bunch 
load end_inj_100pC2

orig = Beam;
new = end_inj_100pC2;

% bunches have changes 0.05 and 0.15 e-13 and e-14

Q = 0%0.15e-13;

idx_orig = (orig.Bunch.Q ~= Q)&(orig.Bunch.stop==0);
idx_new = (new.Bunch.Q ~= Q)&(new.Bunch.stop==0);

arr1 = 5;
arr2 = 6;

figure();
fig = plot(orig.Bunch.x(arr1,idx_orig),orig.Bunch.x(arr2,idx_orig),'.',new.Bunch.x(arr1,idx_new),new.Bunch.x(arr2,idx_new),'.');
legend("Original", "GPT");

xlabel("Longitudinal Position, z (m)")
ylabel("Energy (GeV)")

disp(mean(orig.Bunch.x(1, idx_orig).*orig.Bunch.x(2,idx_orig)))
disp(mean(new.Bunch.x(1, idx_new).*new.Bunch.x(2,idx_new)))

disp(std(orig.Bunch.x(2, idx_orig)))
disp(std(new.Bunch.x(2, idx_new)))

figure()
histogram(orig.Bunch.x(arr1,idx_orig), 100, 'Normalization', 'probability');
hold on
histogram(new.Bunch.x(arr1, idx_new), 100, 'Normalization', 'probability', 'FaceColor', 'red');
legend("Original", "GPT");
xlabel("Longitudinal Position, z (m)")
sgtitle('Witness');

figure()
histogram(orig.Bunch.x(arr2,idx_orig), 100, 'Normalization', 'probability');
hold on
histogram(new.Bunch.x(arr2, idx_new), 100, 'Normalization', 'probability', 'FaceColor', 'red');
legend("Original", "GPT");
xlabel("Transverse Angle, x' (rad)")
