function plot_3x1_2e(x, y1, y2, ttl)

figure;
subplot(3, 1, 1)
plot(x, y1(1,:), 'k--', 'LineWidth', 1.0);
hold on;
plot(x, y2(1,:), 'r', 'LineWidth', 0.5);
% legend('True', 'Measured', 'Estimated');
% xlim([0, 10]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
hold off;

subplot(3, 1, 2)
plot(x, y1(2,:), 'k--', 'LineWidth', 1.0);
hold on;
plot(x, y2(2,:), 'r', 'LineWidth', 0.5);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
ylabel(ttl, 'interpreter', 'latex');
% xlim([0, 10]);
hold off;

subplot(3, 1, 3)
plot(x, y1(3,:), 'k--', 'LineWidth', 1.0);
hold on;
plot(x, y2(3,:), 'r', 'LineWidth', 0.5);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
% ylim([-1, 1]);
% xlim([0, 10]);
xlabel('t (s)')
hold off;

end