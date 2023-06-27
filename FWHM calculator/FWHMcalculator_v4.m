clear
excelFile = 'BandwidthCalcs.xlsx'; % This Excel file contains example 
% bandwidth measurement data for single cell protein blots in gel (sheets g)
% and nitrocellulose post transfer (sheets n) for 1, 2, 3, 4 minutes
sheetName = 'n1'; % nitrocellulose data, 1 minute. Change to other sheets
% to see what gel and nitrocellulose data looks like at 1, 2, 3, 4 minutes
data = xlsread(excelFile, sheetName);

% Create output table
T = table();

% Extract x and y data
% x data is in the first column, microns
% y data is all subsequent columns, represents EGFP fluorescence measured at
% each point x
x = data(:, 1);
xx = linspace(x(1), x(end), 100);
y = data(:, 2:end);
numColumns = size(y, 2);
disp(['Number of columns: ' num2str(numColumns)]);

n = numColumns; % number of columns of y data to smooth, plot, and analyze
for i = 1:n
    z = smooth(y(:, i));
    yy = interp1(x, z, xx, 'spline');
    plot(xx, yy); % Changed 'plot' to 'plot(xx, yy)'
    hold on
    
    % Find the half height - midway between min and max y values.
    halfHeight = (min(yy) + max(yy)) / 2;
    yline(halfHeight, 'Color', 'g', 'LineWidth', 2);
    
    % Find left edge
    index1 = find(yy >= halfHeight, 1, 'first');
    x1 = xx(index1);
    line([x1, x1], [0, yy(index1)], 'Color', 'r', 'LineWidth', 2);
    text(x1 + 1, 5, sprintf('x1 = %.2f', x1), 'FontSize', 12, 'Color', 'r');
    
    % Find right edge
    index2 = find(yy >= halfHeight, 1, 'last');
    x2 = xx(index2);
    line([x2, x2], [0, yy(index2)], 'Color', 'r', 'LineWidth', 2);
    text(x2 + 1, 5, sprintf('x2 = %.2f', x2), 'FontSize', 12, 'Color', 'r');
    
    % Compute the full width, half max.
    fwhm = x2 - x1;
    fwhm_t = table(fwhm(1));
    T = [T; fwhm_t];
    
    caption = sprintf('Full Width, Half Max = %.2f', x2 - x1);
    title(caption, 'FontSize', 12);
end

% Display the output table
disp(T);
