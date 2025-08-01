%% P04_RxSimulation
%% Gain Calculation
%% Rx Gain GS
LEOGSrx = 10* log10((pi * leogsAntenna *fc /c)^2 * eff);
GEOGSrx = 10* log10((pi * geogsAntenna *fc /c)^2 * eff);
ThermalNoisedBm = 10 * log10(kb * TempK * ChannelBW) +30; % Noise in dBm
ThermalNoise = 10^(ThermalNoisedBm / 10);
%% LEO Tx Gain
% Define a realistic 1D sinc-squared approximation pattern shape antenna case (2)
LeoGainMax = 10* log10((pi * leo.Antenna *fc /c)^2 * eff);
GtxLEO = LeoGainMax + 10 * log10( (sinc(0.2 * leotheta / leo.psi)).^2 );
%% GEO Tx Gain
GtxGEO = 10* log10((pi * geoAntenna *fc /c)^2 * eff);
%% LEO Power calculations
RhoLEO(ElLEO<0) = Inf;
PathLoss = 20*log10(fc) + 20*log10(RhoLEO) -147.55;
AtmoLLEO = F01_ComputeAtmosphericLoss(fc, ElLEO, Att);
FadingLEO = F02_MultipathFadingLoss(FadingModel, ElLEO);
% PrxLEO = leoPower + GtxLEO + LEOGSrx - PathLoss - AtmoLLEO - FadingLEO;
PrxLEO = leoPower + GtxLEO + LEOGSrx - PathLoss;
SNRLEO = PrxLEO - ThermalNoisedBm;
%% GEO Power calculations
RhoGEO(ElGEO<0) = Inf;
PathLoss = 20*log10(fc) + 20*log10(RhoGEO) -147.55;
AtmoLGEO = F01_ComputeAtmosphericLoss(fc, ElGEO , Att);
FadingGEO = F02_MultipathFadingLoss(FadingModel, ElGEO);
PrxGEO = geoPower + GtxGEO + GEOGSrx - PathLoss - AtmoLGEO - FadingGEO;
% PrxGEO = geoPower + GtxGEO + Grx - PathLoss;
SNRGEO = PrxGEO - ThermalNoisedBm;
