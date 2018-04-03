Wnd        = 2048;
Freq       = 250*(10:150)/Wnd; %(1:1024)/Wnd; 
load('Part_1_4');
%COH_A  = COH;
COH_A = COH_F;

load('Part_2_4');
%COH_B  = COH;
COH_B = COH_F;

load('Part_3_4');
%COH_C  = COH;
COH_C = COH_F;


subplot(3,1,1)
imagesc((1:size(COH_A,1))*Wnd*8/(250*3600),Freq,COH_A');
colorbar
colormap('jet')
xlabel('hour')
ylabel('Freq')
set(gca,'fontsize',24)
title('PNAS')

subplot(3,1,2)
imagesc((1:size(COH_C,1))*Wnd/(250*3600),Freq,COH_C');
colorbar
colormap('jet')
xlabel('hour')
ylabel('Freq')
set(gca,'fontsize',24)
title('Method 1')
subplot(3,1,3)
imagesc((1:size(COH_B,1))*Wnd/(250*3600),Freq,COH_B');
colorbar
colormap('jet')
xlabel('hour')
ylabel('Freq')
title('Method 2')
set(gca,'fontsize',24)
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 48 36])
print('-dpng',['data_4.png']);



subplot(3,1,1)
imagesc((1:size(COH_A,1))*Wnd*8/(250*3600),Freq,COH_A');
colorbar
colormap('jet')
xlabel('hour')
ylabel('Freq')
set(gca,'fontsize',24)
title('PNAS')
ylim([8 12])

subplot(3,1,2)
imagesc((1:size(COH_C,1))*Wnd/(250*3600),Freq,COH_C');
colorbar
colormap('jet')
xlabel('hour')
ylabel('Freq')
set(gca,'fontsize',24)
title('Method 1')
ylim([8 12])
subplot(3,1,3)
imagesc((1:size(COH_B,1))*Wnd/(250*3600),Freq,COH_B');
colorbar
colormap('jet')
xlabel('hour')
ylabel('Freq')
title('Method 2')
set(gca,'fontsize',24)
ylim([8 12])
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 48 36])
print('-dpng',['data_5.png']);
