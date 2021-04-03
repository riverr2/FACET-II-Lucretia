clear all % Including global variables
close all
% Launch single Lucretia simulation

% -5 and -2.
% -6 and -2.4
% -7 and -2.8

P1 = 0;
P2 = 0;

linacVals = struct('P1',-19.2+P1,'P2',-38.35+P2,'P3',-45,'V1',0,'V2',0,'qi',0,'dx',0,'dy',0);% The defaults for full compression case 
tic
dati = twobunch_scan_fcn_newS20_v3(linacVals);
toc

%LPSCurrent(dati(2).beam, 0, 10e-3, 500, dati(2).zpos)
%LPSCurrent(dati(4).beam, 0, 10e-3, 500, dati(4).zpos)
%LPSCurrent(dati(6).beam, 0, 10e-3, 500, dati(6).zpos)
%LPSCurrent(dati(8).beam, 0, 10e-3, 500, dati(8).zpos)

bc11 = dati(2).beam.Bunch;
save('bc11.mat', 'bc11');
bc14 = dati(4).beam.Bunch;
save('bc14.mat', 'bc14');
bc20 = dati(6).beam.Bunch;
save('bc20.mat', 'bc20');
