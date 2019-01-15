data = importdata('cubic_raman_raw.txt');
freq = data(:,1);
intensity = data(:,2);
xgrid = 5:0.01:160;
sigma = 0.8;
sizef = length(freq);
sizex = length(xgrid);
gauss = zeros(sizef,sizex);
broaden = zeros(1,sizex);
for i=1:sizef
    for j=1:sizex
        gauss(i,j) = intensity(i)*1/sqrt(2*pi)/sigma*exp(-(xgrid(j) - freq(i))^2/2/sigma^2);
    end
end

for i=1:sizef
    broaden = broaden + gauss(i,:);
end


s = plot(xgrid,broaden);
title('Dynamical MAPbI_3 Raman Shift (Orthorhombic)')
xlabel('Frequency/cm^{-1}')
s.LineWidth = 1.5;
set(gca,'FontSize',24)
