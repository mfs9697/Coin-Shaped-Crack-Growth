function fig_TI_relaxation_moduli(cfg)
% fig_TI_relaxation_moduli
% ------------------------------------------------------------
% Plot relaxation moduli used for TI/transversal viscoelastic material:
%   M(t) = M_inf + (M_inst - M_inf)*exp(-t/tau)
% where M is {Et, Ep, Gt}.
%
% Expected fields:
%   cfg.mat.Et_inst, cfg.mat.Ep_inst, cfg.mat.Gt_inst
%   cfg.mat.tau_k
%   cfg.mat.ratio_E2_relax, cfg.mat.ratio_E3_relax, cfg.mat.ratio_G_relax

    m = cfg.mat;

    % --- time grid (log is convenient for viscoelasticity)
    tau = m.tau_k;
    tmin = max(1e-2*tau, 1e-8);
    tmax = 1e2*tau;
    t = logspace(log10(tmin), log10(tmax), 600);

    % --- instantaneous and relaxed values
    Et0 = m.Et_inst;
    Ep0 = m.Ep_inst;
    Gt0 = m.Gt_inst;

    EtInf = Et0 / m.ratio_E2_relax;
    EpInf = Ep0 / m.ratio_E3_relax;
    GtInf = Gt0 / m.ratio_G_relax;

    % --- relaxation laws
    Et = EtInf + (Et0 - EtInf).*exp(-t./tau);
    Ep = EpInf + (Ep0 - EpInf).*exp(-t./tau);
    Gt = GtInf + (Gt0 - GtInf).*exp(-t./tau);

    % --- plot
    figure(1); clf(1);
    plot(t, Et, 'LineWidth', 1.7); hold on;
    plot(t, Ep, 'LineWidth', 1.7);
    plot(t, Gt, 'LineWidth', 1.7);
    grid on;
    set(gca,'XScale','log');

    xlabel('time, t');
    ylabel('modulus (same units as cfg.mat)');
    title('Relaxation moduli (TI / transversal): single-exponential');

    legend({'E_t(t)','E_p(t)','G_t(t)'}, 'Location','best');

    % --- annotate key parameters/levels
    txt = sprintf(['\\tau = %.4g\n' ...
                   'E_t: %.4g \\rightarrow %.4g (ratio %.4g)\n' ...
                   'E_p: %.4g \\rightarrow %.4g (ratio %.4g)\n' ...
                   'G_t: %.4g \\rightarrow %.4g (ratio %.4g)'], ...
                   tau, Et0, EtInf, m.ratio_E2_relax, ...
                   Ep0, EpInf, m.ratio_E3_relax, ...
                   Gt0, GtInf, m.ratio_G_relax);
    text(0.02, 0.02, txt, 'Units','normalized', ...
         'VerticalAlignment','bottom', 'Interpreter','tex');
end
