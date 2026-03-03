function getFigureAllFreqSvHystogram(azfpData)
%GETFIGUREALLFREQSVHYSTOGRAM Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    azfpData    
end

figure

for i = 1:4
    subplot(2,2,i)
    sv_data = azfpData(i).Sv(:);
    histogram(sv_data, 100, 'Normalization', 'probability')
    xlabel('Sv (dB)')
    ylabel('Probability')
    title(sprintf('%d kHz', azfpData(i).Freq))
    grid on
end

end