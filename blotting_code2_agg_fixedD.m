clear all;
close all;

% Initialize parameters
nx = 500; % number of spatial points
L = 1000; % length of domain in µm
dx = L/nx; % spatial step size in µm
D = 89.8; % diffusion coefficient for GFP in free solution (Stokes-Einstein), in µm^2 s^-1, see Damkohler_calcs_fixedD
%dt = dx^2/(2*D); % time step size
%dt= 0.1; %time in s
dt= 0.02; %time in s
T = 300; % total simulation time in s
%nt = floor(T/dt); % number of time steps to plot
nt = 50;
%immobilization rate (s^-1)
%0.1-0.00001
%k=0.1
FWHM_holder_holder=[];
percent_bound_holder_holder=[];

%looking at band spreading over time for different immobilization rates k

for k=[0.000897604,0.008976043,0.017952086,0.02692813,0.044880216,0.062832302,0.089760432,0.179520864,0.44880216];
%for k=[0.02692813,0.44880216],

% Initialize position and time vectors
x = linspace(0, L, nx);
t = 0:dt:T;

% Initialize Gaussian initial condition
A = 1; % amplitude

%say starting band width is 100 µm FWHM, then sigma = 100/2.355 = 42 µm
band_width_start=100; %starting band width is 100 µm FWHM
sigma_start = 42; % standard deviation
%initial distribution of protein
u0 = A*exp(-(x-L/2).^2/(2*sigma_start^2));
%initial bound protein distribution
u2_0=zeros(1,length(x));

% Define Gaussian function for FWHM fit
gaussiann = @(params,x1)(params(1)*exp(-(x-params(2)).^2/(2*params(3)^2)));

%params(1).*exp(-((x1-params(2))/params(3)).^2)

% Time-stepping loop that simulates diffusion and immobilization pde's
FWHM_holder=[];
percent_bound_holder=[];
u = u0;
u2=u2_0;
% figure
% xlabel('Position')
% ylabel('Concentration')
% axis([0 L 0 1])
% hold on;
% plot(x, u2)
for i = 1:length(t)
%for i = 1,
    
    % Forward difference in time for unbound protein (FTCS scheme)
    u = u + D*dt/dx^2*([u(2:end), u(end)] + [u(1), u(1:end-1)] - 2*u) - k.*u.*dt;
    
%     %try second type of scheme
%     % Compute spatial derivative
%     spatial_derivative = (circshift(u, [-1 0]) - circshift(u, [1 0]))/(2*dx);
%     % Compute time derivative
%     time_derivative = D*spatial_derivative - k*u;
%     % Update concentration
%     u = u + time_derivative*dt;

    %try third type of scheme
    %u = u + dt*([0 D*(diff(u,2)) 0]);

    %conc of bound protein 
    u2=u2+k.*u.*dt;

%     % Plot concentration at some time points
%     if any([0:T/nt:T]==t(i)),
%     plot(x, u2)
%     else,
%     end;
%     %makes a nice dynamic plot
%     %pause(0.001)

    %now fit a Gaussian to each iteration and determine FWHM

    % Fit conc of bound protein to Gaussian function
    params=[0.5,L/2,sigma_start];
    param = nlinfit(x,u2,gaussiann,params);

    % Compute FWHM
    sigma = param(3);
    FWHM = 2*sqrt(2*log(2))*sigma;

    [FWHM_holder]=[FWHM_holder;FWHM];

    %compute % protein bound
    percent_bound=100*(trapz(u2)/(trapz(u)+trapz(u2)));
    percent_bound_holder=[percent_bound_holder;percent_bound];
end

% figure
% plot([0:dt:T],FWHM_holder);
% %axis equal
% axis([0 T 0 L/3])
% 
% %color points by % immobilized
% figure
% scatter([0:dt:T],FWHM_holder,20,[zeros(size(FWHM_holder,1),1) (percent_bound_holder-.01)/100 zeros(size(FWHM_holder,1),1)],'filled');
% %axis equal
% axis([0 T 0 L/3])
% 
% figure
% hold on
% plot(x, u2)
% plot(x,gaussiann(param,x))

FWHM_holder_holder=[FWHM_holder_holder FWHM_holder];
percent_bound_holder_holder=[percent_bound_holder_holder percent_bound_holder];

end;

%plot FWHM of immobilized protein fraction for different k values
figure
hold on
for j=1:size(FWHM_holder_holder,2),
    plot([0:dt:T]/60,FWHM_holder_holder(:,j)-band_width_start);
end;
%plot until 4.5 min
axis([0 4.5 0 80]);

%color points by % immobilized
fig=figure;
hold on
for j=1:size(FWHM_holder_holder,2),
    scatter([0:dt:T]/60,FWHM_holder_holder(:,j)-band_width_start,120,[zeros(size(FWHM_holder,1),1) (percent_bound_holder_holder(:,j)-.01)/100 zeros(size(FWHM_holder,1),1)],'filled');
end;
%plot until 4.5 min
axis([0 4.5 0 80]);
set(gca, 'TickDir', 'out');
fig.Position=[10 10 1000 800];

% save to transparent image
set(gcf, 'color', 'none');    
set(gca, 'color', 'none');

% saving to transparent requires some convoluted scheme starting with export with two background colors
exportgraphics(fig, 'test1.png', 'ContentType', 'image', 'BackgroundColor', 'k'); % black background
exportgraphics(fig, 'test2.png', 'ContentType', 'image', 'BackgroundColor', 'w'); % white background

% load exported images back in and scale to [0,1]
u = imread('test1.png');
u = double(u) / 255;
v = imread('test2.png');
% have to add this sometimes since for some reason v can have one more row than u
%v = v(1:(size(v,1)-1),:,:);
v = double(v) / 255;

% recover transparency as a
a = 1 - v + u;
a = mean(a, 3);
a = max(0, min(1, a));
m = a > eps;

% recover rgb
c = zeros(size(u));
for i = 1 : 3
    ui = u(:, :, i);
    ci = c(:, :, i);
    ci(m) = ui(m) ./ a(m);
    c(:, :, i) = ci;
end
c = max(0, min(1, c));

% store again
imwrite(uint8(c*255), 'Damkohler_calcs_overlayplot_matlab.png', 'Alpha', a);

%output csv of change in band FWHM vs. time
FWHM_holder_holder_diff=[([0:dt:T]/60)', FWHM_holder_holder-band_width_start];
FWHM_holder_holder_diff_sub=FWHM_holder_holder_diff(1:5:size(FWHM_holder_holder_diff,1),:);
csvwrite('spread_dist.csv',FWHM_holder_holder_diff_sub);