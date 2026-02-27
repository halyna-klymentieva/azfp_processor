function getTestFigure(src)
%GETTESTFIGURE Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    src struct    
end

figure%(1)
imagesc(src(3).Sv')
colormap('jet');
clim([-110 -40]);
xlabel('Ping Number')
ylabel('Range')
h = colorbar;
set(get(h,'label'),'string','Sv (dB scattering per unit volume)');

end