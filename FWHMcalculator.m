
x = x;
xx= linspace(x(1,:), x(end,:), 100);
T = table; %data to ebe analyzed imported from Excel file as a table
n = 13; %number of columns of y data to smooth, plot and analyze
for i=1:n

   % Smooth noisy y data if needed as it can affect FWHM estiamtion
    y = smooth(N(:,i));
    yy = interp1(x, y, xx, 'spline');
    plot = line(xx,yy);
    hold on
    % Find the half height - midway between min and max y values
    halfHeight = (min(yy) + max(yy)) / 2;
    hold on;
    yline(halfHeight, 'Color', 'g', 'LineWidth', 2);
    % Find left edge of the plot
    index1 = find(yy >= halfHeight, 1, 'first');
    x1 = xx(index1);
    line([x1, x1], [0, yy(index1)], 'Color', 'r', 'LineWidth', 2);
    text(x1+1, 5, sprintf('x1 = %.2f', x1), 'FontSize', 12, 'Color', 'r');
    % Find right edge of the plot
    index2 = find(yy >= halfHeight, 1, 'last');
    x2 = xx(index2);
    line([x2, x2], [0, yy(index2)], 'Color', 'r', 'LineWidth', 2);
    text(x2+1, 5, sprintf('x2 = %.2f', x2), 'FontSize', 12, 'Color', 'r');
    % Compute the full width, half max of the plot
    fwhm = x2 - x1
    fwhm_t=table(fwhm(1));
    T = [T; fwhm_t];
   
    
end
