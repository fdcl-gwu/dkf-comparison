function plot_3x3_3(t, true, measured, estimated, ttl)

figure;
% suptitle(ttl);

for i = 1:3
    for j = 1:3
        k = 3 * (i - 1) + j;
        subplot(3, 3, k)
        plot(t, squeeze(true(i,j,:)), 'b--', 'LineWidth', 1.5)
        hold on;
        plot(t, squeeze(measured(i,j,:)), 'k', 'LineWidth', 0.2);
        plot(t, squeeze(estimated(i,j,:)), 'r', 'LineWidth', 0.5);
        hold off;
        ylim([-1 1]);
        set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
    end
end
ylabel('t (s)', 'interpreter', 'latex')

subplot(3, 3, 4);
ylabel('$R$', 'interpreter', 'latex')

subplot(3, 3, 8);
xlabel('t (s)', 'interpreter', 'latex')

% legend('True', 'Estimated')
end