%% Main Script to define the Geometrical simulator, Receivers, and Interference
clear; clc;close all hidden;
rng(42); % For reproducibility
%% Define Parameters and Ground stations
P01_Parameters
P02_GStations
%% Create Scenario
fprintf('Creating satellite scenario...\n');
sc = satelliteScenario(startTime, stopTime, sampleTime);
%% Satellite and GS creation
P03_GeometricSimulatoion
% play(sc,PlaybackSpeedMultiplier=100);
% save('Geometric',"gsAntenna","fc","c","leoAntenna","ElLEO","RhoLEO","eff", ...
%     "geoAntenna","ElGEO", "RhoGEO", "NumLeoUser", "NumGeoUser", "ts", "GS", ...
%     "leoNum", "geoNum", "GSLEOFilter", "GSGEOFilter", "leoPower", "geoPower");
% save('Geometric');
%% Power Simulation
% clear; clc;close all hidden
% load("Geometric")
% P01_Parameters
P04_RxSimulation
% save('Power',"PrxLEO","PrxGEO");
%% Channel allocation and Serving satellites
% clear; clc;close all hidden
% load("Geometric")
% load('Power');
P05_ChannelAllocation
% save('Data');
%% Evaluation => Interference, throughput, spectral efficiency
P06_Intf_Eval
