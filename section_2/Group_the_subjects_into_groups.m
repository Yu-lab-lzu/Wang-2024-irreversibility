%% Group the subjects into groups

clear;close all;clc;
subjects = 260;
load('ms_se_rest.mat')
[sort_x index_x] = sort(x);

subplot(3,4,1)
sub_260_1 = index_x(1:subjects);
% plot(x,y,'*');
% hold on
% plot(x(sub_260_1),y(sub_260_1),'*r')
% set (gcf,'unit','centimeters','position', [10 10 20 15]);
% ylim([3.8 4.9])
% xlim([0.2 0.9])
% xlabel('rest all brain ms');ylabel('rest all brain se')
% set(gca,'LineWidth',2,'FontSize',18,'FontName','Arial');

subplot(3,4,2)
sub_260_2 = index_x((subjects/2)+1:3/2*subjects);
% plot(x,y,'*');
% hold on
% plot(x(sub_260_2),y(sub_260_2),'*r')
% set (gcf,'unit','centimeters','position', [10 10 20 15]);
% ylim([3.8 4.9])
% xlim([0.2 0.9])
% xlabel('rest all brain ms');ylabel('rest all brain se')
% set(gca,'LineWidth',2,'FontSize',18,'FontName','Arial');

subplot(3,4,3)
sub_260_3 = index_x((subjects)+1:2*subjects);
% plot(x,y,'*');
% hold on
% plot(x(sub_260_3),y(sub_260_3),'*r')
% set (gcf,'unit','centimeters','position', [10 10 20 15]);
% ylim([3.8 4.9])
% xlim([0.2 0.9])
% xlabel('rest all brain ms');ylabel('rest all brain se')
% set(gca,'LineWidth',2,'FontSize',18,'FontName','Arial');

subplot(3,4,4)
sub_260_4 = index_x((3/2*subjects)+1:5/2*subjects);
% plot(x,y,'*');
% hold on
% plot(x(sub_260_4),y(sub_260_4),'*r')
% set (gcf,'unit','centimeters','position', [10 10 20 15]);
% ylim([3.8 4.9])
% xlim([0.2 0.9])
% xlabel('rest all brain ms');ylabel('rest all brain se')
% set(gca,'LineWidth',2,'FontSize',18,'FontName','Arial');

subplot(3,4,5)
sub_260_5 = index_x((2*subjects)+1:3*subjects);
% plot(x,y,'*');
% hold on
% plot(x(sub_260_5),y(sub_260_5),'*r')
% set (gcf,'unit','centimeters','position', [10 10 20 15]);
% ylim([3.8 4.9])
% xlim([0.2 0.9])
% xlabel('rest all brain ms');ylabel('rest all brain se')
% set(gca,'LineWidth',2,'FontSize',18,'FontName','Arial');


subplot(3,4,6)
sub_260_6 = index_x((5/2*subjects)+1:7/2*subjects);
% plot(x,y,'*');
% hold on
% plot(x(sub_260_6),y(sub_260_6),'*r')
% set (gcf,'unit','centimeters','position', [10 10 20 15]);
% ylim([3.8 4.9])
% xlim([0.2 0.9])
% xlabel('rest all brain ms');ylabel('rest all brain se')
% set(gca,'LineWidth',2,'FontSize',18,'FontName','Arial');

subplot(3,4,7)
sub_260_7 = index_x((3*subjects)+1:4*subjects);
% plot(x,y,'*');
% hold on
% plot(x(sub_260_7),y(sub_260_7),'*r')
% set (gcf,'unit','centimeters','position', [10 10 20 15]);
% ylim([3.8 4.9])
% xlim([0.2 0.9])
% xlabel('rest all brain ms');ylabel('rest all brain se')
% set(gca,'LineWidth',2,'FontSize',18,'FontName','Arial');

subplot(3,4,8)
sub_260_8 = index_x((7/2*subjects)+1:9/2*subjects);
% plot(x,y,'*');
% hold on
% plot(x(sub_260_8),y(sub_260_8),'*r')
% set (gcf,'unit','centimeters','position', [10 10 20 15]);
% ylim([3.8 4.9])
% xlim([0.2 0.9])
% xlabel('rest all brain ms');ylabel('rest all brain se')
% set(gca,'LineWidth',2,'FontSize',18,'FontName','Arial');

subplot(3,4,9)
sub_260_9 = index_x((4*subjects)+1:5*subjects);
% plot(x,y,'*');
% hold on
% plot(x(sub_260_9),y(sub_260_9),'*r')
% set (gcf,'unit','centimeters','position', [10 10 20 15]);
% ylim([3.8 4.9])
% xlim([0.2 0.9])
% xlabel('rest all brain ms');ylabel('rest all brain se')
% set(gca,'LineWidth',2,'FontSize',18,'FontName','Arial');


subplot(3,4,10)
sub_260_10 = index_x((9/2*subjects)+1:11/2*subjects);
% plot(x,y,'*');
% hold on
% plot(x(sub_260_10),y(sub_260_10),'*r')
% set (gcf,'unit','centimeters','position', [10 10 20 15]);
% ylim([3.8 4.9])
% xlim([0.2 0.9])
% xlabel('rest all brain ms');ylabel('rest all brain se')
% set(gca,'LineWidth',2,'FontSize',18,'FontName','Arial');


subplot(3,4,11)
sub_260_11 = index_x((5*subjects)+1:6*subjects);
% plot(x,y,'*');
% hold on
% plot(x(sub_260_11),y(sub_260_11),'*r')
% set (gcf,'unit','centimeters','position', [10 10 20 15]);
% ylim([3.8 4.9])
% xlim([0.2 0.9])
% xlabel('rest all brain ms');ylabel('rest all brain se')
% set(gca,'LineWidth',2,'FontSize',18,'FontName','Arial');

subplot(3,4,12)
sub_260_12 = index_x(end-subjects+1:end);
% plot(x,y,'*');
% hold on
% plot(x(sub_260_12),y(sub_260_12),'*r')
% set (gcf,'unit','centimeters','position', [10 10 20 15]);
% ylim([3.8 4.9])
% xlim([0.2 0.9])
% xlabel('rest all brain ms');ylabel('rest all brain se')
% set(gca,'LineWidth',2,'FontSize',18,'FontName','Arial');


save('sub_middle_sup_260.mat','sub_260_1','sub_260_2','sub_260_3',...
    'sub_260_4','sub_260_5','sub_260_6','sub_260_7','sub_260_8',...
    'sub_260_9','sub_260_10','sub_260_11','sub_260_12')


