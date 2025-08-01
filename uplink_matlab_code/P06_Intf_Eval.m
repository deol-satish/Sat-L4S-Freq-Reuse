fprintf('Uplink interference calculation step...\n');
T = length(ts);
SINR = NaN(NumGS, T);
Intf = NaN(NumGS, T);
Thrpt = NaN(NumGS, T);
SINR_mW_dict = NaN(NumGS, T);
Intf_mW_dict = NaN(NumGS, T);
Noise_mW = 10^(ThermalNoisedBm / 10);
berQPSK = NaN(NumGS, T);
berMQAM = NaN(NumGS, T);

for t = 1:T
    PrxLEOt = PrxLEO(:, :, t);              % [NumGS x LEO] — uplink received power at LEO
    PrxGEOt = PrxGEO(:, :, t);              % [NumGS x GEO]
    ChannelListLeot = ChannelListLeo(:, :, t);  % GS's uplink channel per LEO
    ChannelListGeot = ChannelListGeo(:, :, t);  % GS's uplink channel per GEO
    PservLEOt = PservLEO(:, t);
    Serv_idxLEOt = Serv_idxLEO(:, t);
    PservGEOt = PservGEO(:, t);
    Serv_idxGEOt = Serv_idxGEO(:, t);
    
    for userIdx = 1:NumGS
        isLEOUser = GSLEOFilter(userIdx);
        isGEOUser = GSGEOFilter(userIdx);

        %% Identify serving satellite
        if isLEOUser
            s_serv = Serv_idxLEOt(userIdx);
            if s_serv == 0 || isnan(s_serv), continue; end
            ch_user = ChannelListLeot(userIdx, s_serv);
            Psig_dBm = PservLEOt(userIdx);  % uplink received power at LEO from GS
            Psig_mW = 10^(Psig_dBm / 10);
        elseif isGEOUser
            s_serv = Serv_idxGEOt(userIdx);
            if s_serv == 0 || isnan(s_serv), continue; end
            ch_user = ChannelListGeot(userIdx, s_serv);
            Psig_dBm = PservGEOt(userIdx);
            Psig_mW = 10^(Psig_dBm / 10);
        else
            continue;
        end

        %% Uplink Interference Calculation
        Pint_mW = 0;
        interferers = [];

        for otherIdx = 1:NumGS
            if otherIdx == userIdx, continue; end
            % Must be same type (LEO or GEO) and use same sat
            if isLEOUser && GSLEOFilter(otherIdx)
                ch_other = ChannelListLeot(otherIdx, s_serv);
                if ch_other == ch_user
                    Pint_dBm = PrxLEOt(otherIdx, s_serv);  % interference from this GS
                    if ~isnan(Pint_dBm) && ~isinf(Pint_dBm)
                        testInterf_mW = Pint_mW + 10^(Pint_dBm/10);
                        SINR_test = Psig_mW / (testInterf_mW + Noise_mW);
                        if 10*log10(SINR_test) < SINRThreshold
                            Pint_mW = testInterf_mW;
                            interferers(end+1) = otherIdx;
                        end
                    end
                end
            elseif isGEOUser && GSGEOFilter(otherIdx)
                ch_other = ChannelListGeot(otherIdx, s_serv);
                if ch_other == ch_user
                    Pint_dBm = PrxGEOt(otherIdx, s_serv);
                    if ~isnan(Pint_dBm) && ~isinf(Pint_dBm)
                        testInterf_mW = Pint_mW + 10^(Pint_dBm/10);
                        SINR_test = Psig_mW / (testInterf_mW + Noise_mW);
                        if 10*log10(SINR_test) < SINRThreshold
                            Pint_mW = testInterf_mW;
                            interferers(end+1) = otherIdx;
                        end
                    end
                end
            end
        end

        %% Final SINR and Throughput
        Pint_totaldB = 10 * log10(Pint_mW + eps);
        SINR_mW = Psig_mW / (Pint_mW + Noise_mW);
        SINR(userIdx, t) = 10 * log10(SINR_mW);
        SINR_mW_dict(userIdx, t) = SINR_mW;
        Intf_mW_dict(userIdx, t) = Pint_mW;
        Intf(userIdx, t) = Pint_totaldB;

        Thrpt(userIdx, t) = log2(1 + SINR_mW);

        % BERs
        berQPSK(userIdx, t) = qfunc(sqrt(2 * SINR_mW));
        M = 16;
        berMQAM(userIdx, t) = (4 / log2(M)) * (1 - 1 / sqrt(M)) * qfunc(sqrt(3 * SINR_mW / (M - 1)));

        %% Debug Print
        
        fprintf('[Uplink t=%d] GS %d → Sat %d (Ch %d): Psig=%.2f dBm, Intf=%.2f dBm, SINR=%.2f dB, Thrpt=%.2f, BER(QPSK)=%.4f, BER(MQAM)=%.4f\n', ...
            t, userIdx, s_serv, ch_user, Psig_dBm, Pint_totaldB, SINR(userIdx, t), ...
            Thrpt(userIdx, t), berQPSK(userIdx, t), berMQAM(userIdx, t));

        
        if ~isempty(interferers)
            fprintf('    ↳ Interferers: %s\n', mat2str(interferers));
        end
    end
end
