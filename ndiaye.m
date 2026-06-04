% =========================================================
% TP - Interpolation de Lagrange et intégration numérique
% Objectif : approximation polynomiale de f(z) sur [-2, 6]
%            puis calcul numérique de l'intégrale
% ===========================================================
clear; clc; close all;

% ---------------------------------------------------------
% DONNÉES DE BASE
% Points de collocation (xj, fj) 
% ---------------------------------------------------------
x = [-2, 0, 2, 6];
f = [3, 5, 8, 5];
n = length(x) - 1; % degré du polynôme de Lagrange

I_exacte = 448/9; % valeur exacte de l'intégrale

% =========================================================
% QUESTION 0 — Créer donnees.dat et tracer les points
% =========================================================

writematrix([x', f'], 'donnees.dat', 'Delimiter', 'space');

figure;
plot(x, f, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
xlabel('x_j'); ylabel('f_j');
title('Question 0 — Points de mesure');
grid on;
saveas(gcf, 'q0_points.png');

% =========================================================
% QUESTION 1 — Calcul de P(ztest) par interpolation de Lagrange
% =========================================================

ztest = 1.3; % valeur quelconque dans [-2, 6]
Ptest = 0;

for i = 1:n+1
    Li = 1;
    for j = 1:n+1
        if j ~= i
            Li = Li * (ztest - x(j)) / (x(i) - x(j));
        end
    end
    Ptest = Ptest + f(i) * Li;
end

fprintf('--- Question 1 ---\n');
fprintf('ztest = %.4f => P(ztest) = %.6f\n\n', ztest, Ptest);

% =========================================================
% QUESTION 2 — Construction de P point par point (m = 10)
% =========================================================

m  = 10;
x0 = x(1);   % borne gauche = -2
xn = x(end); % borne droite =  6
h  = (xn - x0) / m;

z = x0 + (0:m) * h; % abscisses z_k = x0 + k*h, k = 0,...,m
P = zeros(1, m+1);  % vecteur P à m+1 composantes

for k = 1:m+1
    for i = 1:n+1
        Li = 1;
        for j = 1:n+1
            if j ~= i
                Li = Li * (z(k) - x(j)) / (x(i) - x(j));
            end
        end
        P(k) = P(k) + f(i) * Li;
    end
end

fprintf('--- Question 2 --- (m = %d, h = %.4f)\n', m, h);
fprintf('  k\t  z_k\t\t  P_k\n');
for k = 1:m+1
    fprintf('  %d\t  %.4f\t\t  %.6f\n', k-1, z(k), P(k));
end
fprintf('\n');

% =========================================================
% QUESTION 3 — Fichier lagrange.dat 
% =========================================================

writematrix([z', P'], 'lagrange.dat', 'Delimiter', 'space');

% Tracé MATLAB ()
figure;
plot(z, P, 'b-o', 'MarkerSize', 5, 'LineWidth', 1.5); hold on;
plot(x, f, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
xlabel('z'); ylabel('P(z)');
title('Question 3 — Polynôme de Lagrange + points de collocation');
legend('P(z_k)', 'Points (x_j, f_j)');
grid on;
saveas(gcf, 'figures/q3_lagrange.png');

% =========================================================
% QUESTIONS 4, 5, 6 — Intégration numérique + étude de l'erreur
% =========================================================

% Valeurs de h à tester (question 6)
h_vals = [0.1, 0.01, 0.001, 0.0001, 0.00001];
E_trap = zeros(size(h_vals));
E_simp = zeros(size(h_vals));

fprintf('--- Questions 4, 5, 6 ---\n');
fprintf('%-10s %-15s %-15s %-15s %-15s\n', 'h', 'Ih_trap', 'E_trap', 'Ih_simp', 'E_simp');

for idx = 1:length(h_vals)
    h_k = h_vals(idx);
    m_k = round((xn - x0) / h_k);
    z_k = linspace(x0, xn, m_k + 1);

    % Évaluation du polynôme de Lagrange sur z_k
    Pk = zeros(1, m_k + 1);
    for kk = 1:m_k+1
        for i = 1:n+1
            Li = 1;
            for j = 1:n+1
                if j ~= i
                    Li = Li * (z_k(kk) - x(j)) / (x(i) - x(j));
                end
            end
            Pk(kk) = Pk(kk) + f(i) * Li;
        end
    end

    % --- Question 4 : Trapèzes composites ---
    Ih_trap = h_k * (Pk(1)/2 + sum(Pk(2:end-1)) + Pk(end)/2);
    E_trap(idx) = abs(Ih_trap - I_exacte);

    % --- Question 5 : Simpson composite (m doit être pair) ---
    if mod(m_k, 2) ~= 0
        m_k = m_k + 1;
        z_k = linspace(x0, xn, m_k + 1);
        Pk  = zeros(1, m_k + 1);
        for kk = 1:m_k+1
            for i = 1:n+1
                Li = 1;
                for j = 1:n+1
                    if j ~= i
                        Li = Li * (z_k(kk) - x(j)) / (x(i) - x(j));
                    end
                end
                Pk(kk) = Pk(kk) + f(i) * Li;
            end
        end
    end

    Ih_simp = (h_k/3) * (Pk(1) + 4*sum(Pk(2:2:end-1)) + 2*sum(Pk(3:2:end-2)) + Pk(end));
    E_simp(idx) = abs(Ih_simp - I_exacte);

    fprintf('%-10.5f %-15.6f %-15.2e %-15.6f %-15.2e\n', ...
        h_k, Ih_trap, E_trap(idx), Ih_simp, E_simp(idx));
end

% --- Fichier erreur.dat ---
writematrix([h_vals', E_trap', E_simp'], 'erreur.dat', 'Delimiter', 'space');

% --- Tracé log-log ---
figure;
loglog(h_vals, E_trap, 'b-o', 'LineWidth', 1.5); hold on;
loglog(h_vals, E_simp, 'r-s', 'LineWidth', 1.5);
loglog(h_vals, h_vals.^2 * (E_trap(1)/h_vals(1)^2), 'b--'); % pente 2 théorique
loglog(h_vals, h_vals.^4 * (E_simp(1)/h_vals(1)^4), 'r--'); % pente 4 théorique
xlabel('h'); ylabel('E(h)');
title('Question 6 — Erreur en fonction du pas (log-log)');
legend('Trapèzes', 'Simpson', 'Ordre 2 théorique', 'Ordre 4 théorique');
grid on;
saveas(gcf, 'figures/q6_erreur_loglog.png');

fprintf('\nFichiers générés : donnees.dat, lagrange.dat, erreur.dat\n');
fprintf('Figures sauvegardées dans le dossier figures/\n');