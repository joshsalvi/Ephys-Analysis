j=1;
h = gcf; %current figure handle
axesObjs = get(h, 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
objTypes = get(dataObjs, 'Type');  %type of low-level graphics object
xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
ydata = get(dataObjs, 'YData');

ydataall = [];xdataall=[];

for j = 1:length(ydata)
    ydataall = [ydataall ydata{j}];
    xdataall = [xdataall  xdata{j}];
end
xdu = unique(xdataall);

for j = 1: length(xdu)
    qr = find(xdataall==xdu(j));
    VSmeanall(j) = mean(ydataall(qr));
    VSsemall(j) = std(ydataall(qr))./sqrt(length(ydata));
end

%%
figure;
ciplot(VSmeanall(VSsemall~=0)-VSsemall(VSsemall~=0)./sqrt(1),VSmeanall(VSsemall~=0)+VSsemall(VSsemall~=0)./sqrt(1),xdu(VSsemall~=0),'g')
hold on;
h=plot(xdu(VSsemall~=0),VSmeanall(VSsemall~=0),'k');
ha = get(h,'Parent');
set(ha,'XScale','log')
set(ha,'XLim',[min(xdu)-0.2*min(xdu) max(xdu)+0.2*max(xdu)])
xlabel('Stimulus frequency (Hz)');ylabel('Vector strength');

%%
figure;
ciplot(VSmeanall()-VSsemall()./sqrt(1),VSmeanall()+VSsemall()./sqrt(1),xdu(),'g')
hold on;h=plot(xdu(),VSmeanall(),'k');
ha = get(h,'Parent');
set(ha,'XScale','log')
set(ha,'XLim',[min(xdu)-0.2*min(xdu) max(xdu)+0.2*max(xdu)])
xlabel('Stimulus frequency (Hz)');ylabel('Vector strength');