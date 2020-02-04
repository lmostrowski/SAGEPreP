function time_series_image(data,times,gcapos)

try 
    icadefs
catch
	GUIBUTTONCOLOR   = [0.66 0.76 1]; 
end

data=squeeze(data); axis off;
times=times(:);
plot_height=0.2;
scroll_height=0.1;

% time trace axes:
ax1=axes('Position', ...
    [gcapos(1) gcapos(2)+scroll_height*gcapos(4) gcapos(3) plot_height*gcapos(4)]);
plot(times,data(:,1),'k','LineWidth',2)
xlabel('Time (ms)')
ylabel('Voltage (\mum)')

% spectrogram axes:
ax2=axes('Position',...
    [gcapos(1) gcapos(2)+(plot_height+scroll_height)*gcapos(4) ...
    gcapos(3) (1-plot_height-scroll_height)*gcapos(4)]);
srate = 1000*(length(times)-1)/(times(end)-times(1));
dt=1/srate; fNQ=srate/2; T=round(max(times)/1000); df=1/T; n=T*srate;
H=hann(n);
for j=1:size(data,2)
    x=data(:,j); x=x(:);
    x = x-mean(x); % Re-reference EEG data to common average
    x = H.*x; % Apply Hanning taper (aka. "raised cosine" taper)
    xf = fft(x); % Compute power using Fast Fourier Transform
    Sh = 2 * dt^2 * 1/n * abs(xf).^2; % store periodogram for the j window
    Sh2 = Sh(1:n/2+1); % One-sided spectrum (i.e. only positive frequencies)
    powerSpectra(:,j)=Sh2;
end
faxis = 0:df:fNQ;
tt=1:1:size(data,2);
imagesc(tt,faxis,10*log10(powerSpectra), [-60 0])
axis xy; set(gca,'XTick',[]); ylim([0 50]); ylabel('Frequency (Hz)'); hold on;
scatter(0,50,'d','filled','MarkerEdgeColor','r','MarkerFaceColor','r','LineWidth',4)
xlim([0 size(data,2)])
ylim([0 50])
hold on
plot([1 1], [0 50], 'r', 'LineWidth',2)

% scrollbar
hScrollBar = uicontrol(gcf,'Style','Slider','Min',1,'Max',size(data,2),'Value',1,...
    'Units','Normalized','BackgroundColor',GUIBUTTONCOLOR,'Position',...
    [gcapos(1) gcapos(2)-gcapos(4)*scroll_height/1.2 gcapos(3) gcapos(4)*scroll_height*.3],...
    'Callback',@slidecommand,'SliderStep',[1/size(data,2) 0.2]);

function slidecommand(hObject,eventdata)
scrollPlace=round(hObject.Value);
plot1trace(scrollPlace);
end

function plot1trace(i)
% plot time series data in time window
axes(ax1)
plot(times,data(:,i),'k','LineWidth',2)
xlabel('Time (ms)')
ylabel('Voltage (\mum)')

% plot spectrogram
axes(ax2)
children = get(gca, 'children'); delete(children(2)); delete(children(1));
imagesc(tt,faxis,10*log10(powerSpectra), [-60 0])
axis xy;set(gca,'XTick',[]);ylim([0 50]);ylabel('Frequency (Hz)'); hold on;
scatter(i,50,'d','filled','MarkerEdgeColor','r','MarkerFaceColor','r','LineWidth',4)
xlim([0 size(data,2)])
ylim([0 50])
hold on
plot([i i], [0 50], 'r', 'LineWidth',2)
end

end