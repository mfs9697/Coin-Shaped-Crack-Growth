function Crack3Dintern05_Step2
% Crack3Dintern05_Step2
% Driver with detailed dt-finding diagnostics (scan + fzero iterations).
% Adds time-radius subplot (from Step1) and enforces fsolve precision.

outDir = fullfile(pwd,'Fig2_frames');
if ~exist(outDir,'dir')
    mkdir(outDir);
end

cfg = cfg_CoinCrack_TI_legacy();

fig_TI_relaxation_moduli(cfg);

[ctx, state] = init_problem(cfg);

% --- user caps / diagnostics knobs ---
Tmax   = cfg.solve.Tmax;
dt_eps = cfg.solve.dt_eps;        % you set dt_eps = 2 (OK)
dt_max = cfg.solve.dt_max;        % e.g. Tmax
Nscan  = cfg.solve.Nscan;

% --- figures ---
figure(1); clf(1)
figure(3); clf(3)
figure(4); clf(4)   % dt-scan curve

nt  = cfg.solve.nt;
t1  = zeros(nt+1,1);
t1(1) = 0;
kdt = 0;
p   = cfg.solve.p;

% arrays for stress averaging (kept local)
nnodes  = size(ctx.coords,2);
sigz    = zeros(nnodes,1);
sigzip  = zeros(ctx.nelem*ctx.nip3,1);
sigzt   = zeros(nnodes,30);
nodelmn = zeros(nnodes,30);
indU    = zeros(ctx.eldf,1);

% ---- one-time diagnostics ----
fprintf('--- Mesh/DOF ---\n');
fprintf('nelem=%d, nnodes=%d, ndof=%d, nelnodes=%d, nip3=%d\n', ...
    ctx.nelem, nnodes, ctx.ndof, ctx.nelnodes, ctx.nip3);
fprintf('connect size: %dx%d (expect 10 x nelem)\n', size(ctx.connect,1), size(ctx.connect,2));
fprintf('mc=%d, me=%d, |jc|=%d, |je|=%d, crf=%dx%d\n', ...
    ctx.mc, ctx.me, numel(ctx.jc), numel(ctx.je), size(ctx.crf,1), size(ctx.crf,2));

fprintf('--- Material/CZM/Time ---\n');
tau = cfg.mat.tau_k(1);
fprintf('tau=%g s\n', tau);
fprintf('sigma_applied=%.6g, sc1=%.6g, Dym=%.6g\n', state.sig, ctx.sc1, ctx.Dym);
fprintf('TSL: a1=%.6g, a2=%.6g, c=%.6g\n', ctx.a1, ctx.a2, (3-2*ctx.a1+3*ctx.a2)/6);
fprintf('Tmax=%g, dt_eps=%g, dt_max=%g, Nscan=%d\n\n', Tmax, dt_eps, dt_max, Nscan);

% header
fprintf('%-4s %-10s %-10s %-7s %-6s %-8s %-12s %-12s %-10s\n', ...
    'kdt','t','dt','kD','ci','aFront','sig_eps','sig_out','bDmax');

% continuation guess (context)
dt_prev = cfg.solve.Tmax;

while kdt <= nt && t1(kdt+1) <= Tmax

    dt = 0;
    sig_out = NaN;
    sig_eps = NaN;

    if kdt == 0
        % ---- initial elastic step ----
        ctxK = ctx;
        ctxK.K = StifFnd(ctx.Eo, ctx);

        stateA = struct('U', state.U, 'sig', state.sig, ...
                        'Fpsig', state.Fpsig, 'Ftsig', state.Ftsig);

        fun = @(du) F0c_core(du, ctxK, stateA);
        state.dU = full(fsolve(fun, state.dU0, cfg.solve.fsolve_F0c));

        state.K  = ctxK.K;
        state.kD = 2*state.dU(state.vc(state.ci)) / ctx.Dym;

        t1(kdt+1) = 0;

    else
        % ---- visco step: find dt so that sig_out(dt) == sigma_applied ----
        fprintf('\n---- dt solve: step kdt=%d (kD=%.4f, ci=%d, t=%.4g) ----\n', ...
            kdt, state.kD, state.ci, t1(kdt));

        dt_cap = min(dt_max, Tmax - t1(kdt));
        dt_cap = min(dt_cap, 5*tau);  % optional physics cap

        if dt_cap <= dt_eps
            warning('dt_cap <= dt_eps; stopping. dt_cap=%g', dt_cap);
            break;
        end

        g = @(dtn) sig_only(dtn, cfg, ctx, state) - state.sig;

        sig_eps = sig_only(dt_eps, cfg, ctx, state);
        g_eps   = sig_eps - state.sig;

        if ~isfinite(g_eps)
            error('g(dt_eps) is NaN/Inf at dt_eps=%g. Inner solve unstable.', dt_eps);
        end
        fprintf('g(dt_eps=%g): sig_eps=%.6g, target=%.6g => g=%.6g\n', ...
            dt_eps, sig_eps, state.sig, g_eps);

        % --- log scan for FIRST sign change ---
        dts = logspace(log10(dt_eps), log10(dt_prev), Nscan);
        gvals = zeros(size(dts));
        gvals(1) = g(dts(1));
        idx = [];

        fprintf('Scan (logspace dt):\n');
        fprintf('  %3s  %-12s  %-12s\n', 'i', 'dt', 'g(dt)');
        fprintf('  %3d  %-12.4g  %-12.4g\n', 1, dts(1), gvals(1));

        for ii = 2:Nscan
            gvals(ii) = g(dts(ii));
            if ~isfinite(gvals(ii))
                error('g(dt) became NaN/Inf during scan at dt=%g.', dts(ii));
            end
            fprintf('  %3d  %-12.4g  %-12.4g\n', ii, dts(ii), gvals(ii));

            % IMPORTANT: break at first sign change
            if gvals(ii-1) * gvals(ii) <= 0
                idx = ii-1;
                fprintf('  -> FIRST sign change between i=%d and i=%d\n', ii-1, ii);
                break;
            end
        end

        % Plot scan curve for this step
        figure(4); clf(4);
        semilogx(dts, gvals, '-o'); grid on;
        xlabel('dt'); ylabel('g(dt)=sig(dt)-sig_{applied}');
        title(sprintf('dt scan at kdt=%d (kD=%.3f, ci=%d)', kdt, state.kD, state.ci));

        if isempty(idx)
            fprintf('No sign change found up to dt_cap=%g.\n', dt_cap);
            break;
        end

        dt_lo = dts(idx);   g_lo = gvals(idx);
        dt_hi = dts(idx+1); g_hi = gvals(idx+1);
        fprintf('Bracket: [%.6g, %.6g], g=[%.6g, %.6g]\n', dt_lo, dt_hi, g_lo, g_hi);

        optsZ = cfg.solve.fzero_dt;
        optsZ = optimset(optsZ, ...
            'OutputFcn', @(x,ov,st) fzero_outfun(x,ov,st,state.sig));

        dt = fzero(@(x) g(x), [dt_lo, dt_hi], optsZ);

        fprintf('fzero result: dt=%.6g\n', dt);

        [sig_out, state] = Ucurr_core(dt, cfg, ctx, state);
        dt_prev = dt;
        t1(kdt+1) = t1(kdt) + dt;
    end

    % ---- advance displacement ----
    state.U = state.U + state.dU;

    % ---- stress update & averaging ----
    sigzt(:) = 0;

    for lmn = 1:ctx.nelem
        nod = ctx.connect(:,lmn);
        for i = 1:ctx.nelnodes
            indU((3*i-2):(3*i)) = (3*nod(i)-2):(3*nod(i));
        end

        for j = 1:ctx.nip3
            B1    = B(ctx.xip3(:,j), lmn, ctx);
            deps  = B1*state.dU(indU);
            lmnip = (lmn-1)*ctx.nip3 + j;

            if kdt
                tsig = 0;
                for m = 1:ctx.M
                    tsig = tsig + sum((1-state.zetamn(:,:,m)).*state.S(:,:,lmnip,m),2);
                end
                sigzip(lmnip) = sigzip(lmnip) - tsig(3) + state.tE(3,:)*deps;
            else
                sigzip(lmnip) = ctx.Eo(3,:)*deps;
            end

            deps1 = repmat(deps',6,1);
            for m = 1:ctx.M
                if kdt
                    state.S(:,:,lmnip,m) = state.zetamn(:,:,m).*state.S(:,:,lmnip,m) + ...
                        state.etamn(:,:,m).*(1-state.zetamn(:,:,m)).*deps1;
                else
                    state.S(:,:,lmnip,m) = ctx.Em(:,:,m).*deps1;
                end
            end
        end

        sigzv = ctx.Nextr * sigzip((lmnip-ctx.nip3+1):lmnip);
        for j = 1:ctx.nelnodes
            k1 = find(sigzt(nod(j),:)==0,1,'first');
            sigzt(nod(j),k1) = sigzv(j);
            if kdt==0, nodelmn(nod(j),k1) = lmn; end
        end
    end

    for j = 1:nnodes
        k1 = length(find(sigzt(j,:)));
        sigz(j) = sigzt(j,1:k1) * ctx.velem(nodelmn(j,1:k1)) / ...
                  sum(ctx.velem(nodelmn(j,1:k1))) / ctx.sc1;
        if sigz(j) > 1, sigz(j) = 1; 
        elseif sigz(j) < 0, sigz(j) = 0; 
        end
    end

    % ---- Fpsig update ----
    if kdt
        state.Fpsig = state.K*state.dU + state.Fpsig - state.Ftsig;
    else
        state.Fpsig = state.K*state.dU;
    end

    % ---- cohesive opening diagnostics ----
    bD   = 2*state.U(3*ctx.crf) / ctx.Dym;
    bDmax = max(bD(:));

    aFront = state.yb(state.ci,2);

    fprintf('%-4d %-10.3e %-10.3e %-7.4f %-6d %-8.4f %-12.6g %-12.6g %-10.3g\n', ...
        kdt, t1(kdt+1), dt, state.kD, state.ci, aFront, sig_eps, sig_out, bDmax);

    % ---- visualization ----
    coords1 = ctx.coords;
    for i = 1:nnodes
        coords1(:,i) = ctx.coords(:,i) + state.U((i-1)*3+(1:3))*cfg.solve.fact;
    end
    Mesh.Points   = coords1;
    Mesh.Elements = ctx.connect;

    figure(1); clf(1)
    colormap(jet(10));
    plot3D(Mesh,'ColorMapData',sigz); hold on
    plot3D(Mesh,'FaceAlpha',0,'EdgeColor',.6*[1,1,1]); hold on
    view(26,-26); 
    zoom(4);
    axis equal

    plot_Fig2_live(cfg, ctx, state, t1(kdt+1), sigz, aFront, kdt, outDir);


    % Figure(3): crack-face map (unchanged)
    figure(3); clf(3)
    xp = ctx.coords(1,:);
    yp = ctx.coords(2,:);
    patch(xp(ctx.crf), yp(ctx.crf), bD, 'EdgeColor', .6*[1 1 1]);
    axis equal; xlim([0 ctx.d(1)]); ylim([0 ctx.d(2)]); colorbar()

    % ---- crack-front bookkeeping ----
    if kdt==0
        state.kD = ceil(p*state.kD)/p + 1/p;
    elseif state.kD < 1-1e-8
        state.kD = state.kD + 1/p;
    else
        state.ci = state.ci + 1;
    end

    kdt = kdt + 1;
    pause(.05)
end

end

% ============================================================
function sig = sig_only(dtn, cfg, ctx, state)
% Return only sigma(dt) (first output of Ucurr_core), without committing.
[sig, ~] = Ucurr_core(dtn, cfg, ctx, state);
end

function stop = fzero_outfun(x, optimValues, stateFlag, sigTarget)
stop = false;
switch stateFlag
    case 'init'
        fprintf('fzero iterations:\n');
        fprintf('  %-12s  %-12s\n','dt','g(dt)');
    case 'iter'
        fprintf('  %-12.6g  %-12.6g\n', x, optimValues.fval);
    case 'done'
end
end
