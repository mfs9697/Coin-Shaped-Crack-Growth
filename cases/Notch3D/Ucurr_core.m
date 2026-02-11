function [sig_out, state] = Ucurr_core(dtn, cfg, ctx, state)
%UCURR_CORE  Time-increment operator (no globals).
% Returns:
%   sig_out : scalar (the solved sigma)
%   state   : updated (K,dU,Ftsig,tE,zetamn,etamn, etc.)



% --- shorthand ---
M   = ctx.M;
rhom= ctx.rhom;
Em  = ctx.Em;

% 1) internal variables
for m1 = 1:M
    state.zetamn(:,:,m1) = exp(-dtn ./ rhom(:,:,m1));
    state.etamn(:,:,m1)  = Em(:,:,m1) .* rhom(:,:,m1) / dtn;
end

% 2) assemble history force Ftsig
ndof = ctx.ndof;
Ftsig = zeros(ndof,1);

indU = zeros(ctx.eldf,1);
for lmn1 = 1:ctx.nelem
    nod = ctx.connect(:,lmn1);
    for j1 = 1:ctx.nelnodes
        indU((3*j1-2):(3*j1)) = (3*nod(j1)-2):(3*nod(j1));
    end

    Ftsige = zeros(ctx.eldf,1);
    for j1 = 1:ctx.nip3
        B3 = B(ctx.xip3(:,j1), lmn1, ctx);
        tsig1 = 0;
        lmnip1 = (lmn1-1)*ctx.nip3 + j1;

        for m1 = 1:M
            tsig1 = tsig1 + sum((1-state.zetamn(:,:,m1)) .* state.S(:,:,lmnip1,m1), 2);
        end

        Ftsige = Ftsige + ctx.w3(j1) * (B3' * tsig1) * ctx.velem(lmn1);
    end

    Ftsig(indU) = Ftsig(indU) + Ftsige;
end
Ftsig(ctx.fixvars) = 0;
state.Ftsig = Ftsig;

% 3) tangent modulus and stiffness
tE = ctx.tEi;
for m1 = 1:M
    tE = tE + state.etamn(:,:,m1) .* (1 - state.zetamn(:,:,m1));
end
state.tE = tE;
state.K  = StifFnd(tE, ctx);


% 4) solve coupled system (F0a_core is already globals-free)
ctxA = ctx;            % needs current K
ctxA.K = state.K;

stateA = struct();
stateA.U     = state.U;
stateA.vc = state.vc;
stateA.ci    = state.ci;
stateA.kD    = state.kD;
stateA.Fpsig = state.Fpsig;
stateA.Ftsig = state.Ftsig;

fun = @(us) F0a_core(us, ctxA, stateA);
Us  = fsolve(fun, [state.dU0; 0.01*state.sig], cfg.solve.fsolve_F0c);

state.dU  = full(Us(1:ndof));
sig_out   = Us(end);
end
