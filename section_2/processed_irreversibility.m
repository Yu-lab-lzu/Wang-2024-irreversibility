%% 
clear;close all;clc;

%% data
figure(1)
subplot(2,2,1)
load('data_500_entropy_260.mat')
plot(a,b,'k.-','LineWidth',3,'MarkerSize',30,'MarkerFaceColor','k')
set(gca,'LineWidth',2,'FontSize',18,'FontName','Arial');
hold on
%
set(gca,'LineWidth',2,'FontSize',18,'FontName','Arial');
hold on
P = polyfit(a,b,3);
xi = 0.32:0.01:0.63;
yi = polyval(P,xi);
hold on
plot(xi,yi,'-r','LineWidth',4)
%
x = a;
y = b;
clear b

% Generate a range of x values for a regular interval,
% from the minimum to the maximum value of the original x
x_new = linspace(min(x), max(x), length(x));

% Use linear interpolation to obtain the corresponding y value
y_new = interp1(x, y, x_new, 'linear');  

% Estimate sampling frequency based on interpolated x_new 
delta_x_new = diff(x_new);
Fs = 1 / mean(delta_x_new); 

% Set the cutoff frequency of the low-pass filter
f_c = 0.5; 
Wn = f_c / (Fs / 2);  

% Design a FIR low-pass filter (using Hanning window)
N = 2;  
b = fir1(N, Wn, 'low', hanning(N+1));  

% Filter the interpolated data
filtered_y_new = filtfilt(b, 1, y_new);
hold on
plot(a,filtered_y_new,'-b','LineWidth',4)

legend('Data','Fit','Filter','location','northwest')
legend('box','off')

%%
subplot(2,2,2)
clear a
load('data_500_entropy_260');
 for i=1:length(b)-1
     d(i)=(b(i+1)-b(i))/(a(i+1)-a(i));
     i=i+1;
 end
for j=1:length(a)-1
    xd(j)=a(j);
    j=j+1;
end

plot(xd,d,'g.-','LineWidth',3,'MarkerSize',30,'MarkerFaceColor','k')
set(gca,'LineWidth',2,'FontSize',18,'FontName','Arial');
hold on
ylabel('$$ \rm {{d \dot{I}}}$$','Interpreter','latex');%B
xlabel('mean synchronization');%B
ylim([-0.01 0.22])
xlim([0.3 0.65])

%%
subplot(2,2,3)
clear a xd d
 for i=1:length(xi)-1
     d(i)=(yi(i+1)-yi(i))/(xi(i+1)-xi(i));
     i=i+1;
 end
for j=1:length(xi)-1
    xd(j)=xi(j);
    j=j+1;
end
plot(xd,d,'r-','LineWidth',4,'MarkerSize',30,'MarkerFaceColor','k')
set(gca,'LineWidth',2,'FontSize',18,'FontName','Arial');
hold on
ylabel('$$ \rm {{d \dot{I}}}$$','Interpreter','latex');%B
xlabel('mean synchronization');%B
xlim([0.3 0.65])
ylim([-0.05 0.17])

%%
subplot(2,2,4)
clear a d xd
load('data_500_entropy_260');

for i=1:length(filtered_y_new)-1
     d(i)=(filtered_y_new(i+1)-filtered_y_new(i))/(x_new(i+1)-x_new(i));
     i=i+1;
 end
for j=1:length(x_new)-1
    xd(j)=x_new(j);
    j=j+1;
end
plot(xd,d,'k.-','LineWidth',3,'MarkerSize',30,'MarkerFaceColor','k')
xlim([0.3 0.65]);ylim([0.04 0.158])
ylabel('$$ \rm {{d \dot{I}}}$$','Interpreter','latex');%B
xlabel('mean synchronization');%B
set(gca,'LineWidth',2,'FontSize',18,'FontName','Arial');


