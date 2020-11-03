function plot_3x1_1(x1, y1, ttl, line)

for i = 1:3
    subplot(3, 1, i)
    plot(x1, y1(i,:), line, 'LineWidth', 1.0);
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
    hold on;
end
xlabel('t (s)');

subplot(3, 1, 2)
ylabel(ttl, 'interpreter', 'latex');
end