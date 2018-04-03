Wnd = 2048;
sample_rate = 250; 
Freq      = [10 150]*sample_rate/Wnd;
data_det  =  data_det(2.5e4:size(data_det,1)-2.5e4,:);
ch_no       = size(data_det,2);
if ch_no==32
    for i=1:4
        for j=1:8
            subplot(4,8,(i-1)*8+j)
            spectrogram(data_det(:,(i-1)*8+j),2048,0,2048,250,'yaxis');
            colorbar('off')
            ylim([Freq(1)  Freq(2)]);
            title(['ch'  num2str((i-1)*8+j)])
        end
    end
else
    % if ch_no==64
    for i=1:8
        for j=1:8
            subplot(8,8,(i-1)*8+j)
            spectrogram(data_det(:,(i-1)*8+j),2048,0,2048,250,'yaxis');
            colorbar('off')
            ylim([Freq(1)  Freq(2)]);
            title(['ch'  num2str((i-1)*8+j)])
            colormap('jet');
        end
    end
end
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 64 64])
print('-dpng',['data_spec.png']);

% spectrogram(data_det(:,2),2048,0,2048,250,'yaxis');
%             colorbar('off')
%             ylim([Freq(1)  Freq(2)]);
%             title(['ch'  num2str(2)])
%             colormap('jet');