function plot_Fig2_live(cfg, ctx, state, tnow, sigz, aFront, kdt, outDir)
% plot_Fig2_live
% Live plotting for Fig.2:
%   Subplot(1): profiles 2u_z/Delta_max (left axis) and sigma_z/sigma_max (right axis)
%   Subplot(2): time-radius points a(t)
%
% Mesh coordinate is in cm: x = state.yb(:,2)
% x-limits use cfg.geom.cr_a2 (cm)

persistent fig axProf axTR
persistent hU hSig hTip hLabel hRc
persistent last_t last_a
persistent colors styleIdx

% -------------------- init (first call / figure closed) --------------------
needInit = isempty(fig) || ~isvalid(fig) || ...
           isempty(axProf) || ~isgraphics(axProf) || ...
           isempty(axTR)   || ~isgraphics(axTR);

if needInit
    fig = figure(2); clf(fig);
    set(fig,'Color','w');

    axProf = subplot(1,2,1,'Parent',fig); hold(axProf,'on'); grid(axProf,'on'); box(axProf,'on');
    axTR   = subplot(1,2,2,'Parent',fig); hold(axTR,'on');   grid(axTR,'on'); box(axTR,'on');

    % LaTeX everywhere
    set(fig, 'DefaultTextInterpreter','latex');
    set(fig, 'DefaultAxesTickLabelInterpreter','latex');
    set(fig, 'DefaultLegendInterpreter','latex');

    % Subplot 1 labels/titles
    title(axProf, '$\mathbf{Profiles:}\; 2u_z/\Delta_{\max}\;\&\;\sigma_z/\sigma_{\max}$');
    xlabel(axProf, '$x\;[\mathrm{cm}]$');

    yyaxis(axProf,'left');
    ylabel(axProf, '$2u_z/\Delta_{\max}$');

    yyaxis(axProf,'right');
    ylabel(axProf, '$\sigma_z/\sigma_{\max}$');

    % Subplot 2 (report-ready title)
    title(axTR, '$\mathbf{Time\!-\!radius\ curve:}\; a(t)$');
    xlabel(axTR,'$t\;[\mathrm{s}]$');
    ylabel(axTR,'$a\;[\mathrm{cm}]$');

    % enforce x-limits in cm
    xlim(axProf, [0, 1.2*cfg.geom.cr_a2]);

    % handles
    hU     = gobjects(0);
    hSig   = gobjects(0);
    hTip   = gobjects(0);
    hLabel = gobjects(1);
    hRc    = gobjects(1);

    colors   = lines(10);
    styleIdx = 0;

    last_t = NaN;
    last_a = NaN;
end

% -------------------- profile data (mesh is cm) --------------------
x_all = state.yb(:,2);   % cm

% normalized opening: 2u_z / Delta_max
U_norm = 2*state.U(state.yb(:,1)*ctx.ncoord) / ctx.Dym;

% sigma_z profile: normalize by sigma_max
ci = state.ci;
x_sig     = state.yb(ci:end,2);                  % cm
% sigma_z at those nodes
sig_norm = sigz(state.yb(ci:end,1));
sig_norm = max(0, min(1, sig_norm));

yyaxis(axProf,'left');
ylim(axProf,[-.02, 1.02*max(U_norm)])

yyaxis(axProf,'right');
ylim(axProf, [0 1.02]);

% -------------------- draw Rc marker in Subplot 1 (LaTeX label) --------------------
if ~isempty(state.ci0) && state.ci0 >= 1 && state.ci0 <= size(state.yb,1)
    xRc = state.yb(state.ci0,2);  % cm
    if isempty(hRc) || ~isgraphics(hRc)
        % IMPORTANT: explicitly parent to axProf so it never lands on subplot 2
        hRc = xline(axProf, xRc, ':', 'LineWidth', 1.0);
        hRc.Label = '$R_c$';
        hRc.LabelVerticalAlignment   = 'bottom';
        hRc.LabelHorizontalAlignment = 'left';
        hRc.Interpreter = 'latex';
    end
end

% -------------------- choose style for this time step --------------------
styleIdx = styleIdx + 1;
icol = 1 + mod(styleIdx-1, size(colors,1));
col  = colors(icol,:);

% -------------------- de-emphasize older curves --------------------
for i = 1:numel(hU)
    if isgraphics(hU(i))
        hU(i).LineWidth = 0.75;
        hU(i).Color     = 0.65*hU(i).Color + 0.35*[1 1 1];
    end
end
for i = 1:numel(hSig)
    if isgraphics(hSig(i))
        hSig(i).LineWidth = 0.75;
        hSig(i).Color     = 0.65*hSig(i).Color + 0.35*[1 1 1];
    end
end
for i = 1:numel(hTip)
    if isgraphics(hTip(i))
        hTip(i).MarkerSize = 4;   % <-- enlarged (was ~2)
        hTip(i).Color      = 0.6*hTip(i).Color + 0.4*[1 1 1];
        if isprop(hTip(i),'MarkerEdgeColor'); hTip(i).MarkerEdgeColor = hTip(i).Color; end
        if isprop(hTip(i),'MarkerFaceColor'); hTip(i).MarkerFaceColor = hTip(i).Color; end
    end
end


% -------------------- add newest curves (highlight) --------------------
% Left axis: 2u_z/Delta_max (solid)
yyaxis(axProf,'left');
hU(end+1) = plot(axProf, x_all, U_norm, '-', 'Color', col, 'LineWidth', 1.5);

% marker at crack-front node (ci) on U_norm curve
if ~isempty(ci) && ci>=1 && ci<=numel(x_all)
    hTip(end+1) = plot(axProf, x_all(ci), U_norm(ci), 'o', ...
        'Color', col, 'MarkerEdgeColor', col, 'MarkerFaceColor', col, 'MarkerSize', 4);
end

% Right axis: sigma_z/sigma_max (SOLID, as requested)
yyaxis(axProf,'right');
hSig(end+1) = plot(axProf, x_sig, sig_norm, '-', 'Color', col, 'LineWidth', 1.5);

% -------------------- LaTeX title + newest label (no legend) --------------------
yyaxis(axProf,'left');  % place label in left-axis coordinates

if ~isempty(hLabel) && isgraphics(hLabel)
    delete(hLabel);
end

xlab = min(1.2*cfg.geom.cr_a2, x_all(end));
ylab = U_norm(end);
hLabel = text(axProf, xlab, ylab, sprintf('$t=%.4g$', tnow), ...
    'HorizontalAlignment','right', ...
    'VerticalAlignment','bottom', ...
    'Color', col);

title(axProf, sprintf('$\\mathbf{Profiles}\\; (k=%d):\\; t=%.4g$', kdt, tnow));
xlim(axProf, [0, 1.2*cfg.geom.cr_a2]);  % lock again for safety

% -------------------- Subplot 2: time-radius points --------------------
plot(axTR, tnow, aFront, 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
if ~isnan(last_t) && ~isnan(last_a)
    plot(axTR, [last_t, tnow], [last_a, aFront], '-', 'Color', 0.2*[1 1 1], 'LineWidth', 0.9);
end
last_t = tnow;
last_a = aFront;

drawnow limitrate


% ---- export Fig.2 as PNG (600 dpi) ----
fname = sprintf('Fig2_k%03d_t%.4g.png', kdt, tnow);
fpath = fullfile(outDir, fname);

exportgraphics(fig, fpath, 'Resolution', 600);
end
