function plot_3x1_sigma(t, true, estimated, sigma, ttl, lim)

figure;

for i = 1:3
    error = true(i,:) - estimated(i,:);
    
    subplot(3, 1, i);
    pl1 = [t, fliplr(t)];
    pl2 = [error+sigma(i,:), fliplr(error-sigma(i,:))];
    pl3 = fill(pl1, pl2, 'b', 'LineWidth', 0.01, 'edgealpha', 0.0);
    hold on;
    alpha(pl3, 0.1);
    plot(t, error, 'r');
    hold off;
    
    ylim([-lim, lim]);
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
end
xlabel('t (s)');

% subplot(3, 1, 1)
% title(ttl, 'interpreter', 'latex');

subplot(3, 1, 2)
ylabel(ttl, 'interpreter', 'latex');
end