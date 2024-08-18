% Make figures for Fig 5C
% saving images showing colony contour & cooperative regions (Fig 5C)
figure(2)
set(gcf,'position',[680   158   450   420]);
cooperReg = double(N <= N_upper & N >= N_lower); % region of cooperative swarming
pcolor(xx, yy, cooperReg);
shading interp; axis equal; axis([-L/2 L/2 -L/2 L/2]); 
colormap(flipud(gray)); caxis([0 5])
hold on
contour(xx, yy, totalt,[0.01 0.01],'k')
set(gca,'YTick',[], 'XTick',[], 'Visible', 'off')
fig = gcf; fig.PaperPositionMode = 'auto';
print(['results\Fig5C_' prefix '_t=' num2str(i*dt) '.jpg'],'-dtiff','-r800')