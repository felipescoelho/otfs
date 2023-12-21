% otfs_vs_ofdm.m
%
% Script to compare OTFS and OFDM
%
% luizfelipe.coelho@smt.ufrj.br
% Dec. 19, 2023
%


clc
clear
close all

addpath('functions/ofdm')
addpath('functions/utils')

% Simple CP-OFDM:
% Definitions:
modSize = 4;  % Number of symbols in QAM constellation
channelMode = 'EVA';
numSubcarriers = 256;
numBlocks = 16;  % Number of information blocks in a frame
carrierFreq = 4*1e9;
samplingPeriod = 200*1e-9;
blockDuration = numSubcarriers*samplingPeriod;
subcarrierSpacing = 1/blockDuration;
snrdB = linspace(-15, 50, 25);
maxSpeed = linspace(0, 300, 10);

% Simulation settings + memory allocation:
numChannels = 1;
ensembleMonteCarlo = 100;
ofdmSER = zeros(length(snrdB), length(maxSpeed));
otfsSER = zeros(length(snrdB), length(maxSpeed));

for speedIdx = 1:length(maxSpeed)
    for channelIdx = 1:numChannels
        [gi, ki, li] = channel_parameters(channelMode, ...
            maxSpeed(speedIdx+1), carrierFreq, samplingPeriod, ...
            subcarrierSpacing/numBlocks);
        cpLength = max(li);  % Enough CP
        [H, h] = channel_matrix_ofdm(gi, ki, li, numSubcarriers, ...
            numBlocks, cpLength);
        serOFDM_chan = zeros(length(snrdB), 1);
        for snrIdx = 1:length(snrdB)
            serOFDM_mc = 0;
            serOTFS_mc = 0;
            for mcIdx = 1:ensembleMonteCarlo
                % CP-OFDM
                txSymbols = qam_random_symbols(modSize, numSubcarriers, ...
                    numBlocks);
                txOFDM = add_redundancy_matrix(numSubcarriers, cpLength) ...
                    * unit_dftmtx(numSubcarriers)*txSymbols;
                u = txOFDM(:);
                y = H*u;
                [y, ~] = add_wgn(y, snrdB(end));
                rxSymbols = conj(unit_dftmtx(numSubcarriers)) ...
                    * rm_redundancy_matrix(numSubcarriers, cpLength) ...
                    * reshape(y, numSubcarriers+cpLength, numBlocks);
                % Equalization:
                hEstimated = rxSymbols(:, 1)./txSymbols(:, 1);
                rxSymbolsEqualized = rxSymbols(:, 2:end)./repmat(hEstimated, 1, numBlocks-1);
                serOFDM_mc = serOFDM_mc + estimate_ser(rxSymbolsEqualized, ...
                    txSymbols(:, 2:end), modSize);
            end
            serOFDM_chan(snrIdx) = serOFDM_mc/ensembleMonteCarlo;
        end
    end
end

