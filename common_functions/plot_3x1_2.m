function plot_3x1_2(x1, y1, x2, y2, ttl)

figure;

for i = 1:3
    subplot(3, 1, i)
    plot(x1, y1(i,:));
    hold on;
    plot(x2, y2(i,:));
    hold off;
end

subplot(3, 1, 1)
title(ttl, 'interpreter', 'latex');

end