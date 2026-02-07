function Notch3D02

global rwi cli st dU K tE Ftsig Fpsig S

clear var
opt=optimset('Display','off','Jacobian','on','TolFun',1e-18); %'TolX',1e-50
opt2=optimoptions('fsolve','Display','iter','TolFun',1e-18,'MaxIterations',30);

phi1=150;

a0=0.001; a1=0.001; a2=a1; c=(3-2*a1+3*a2)/6;
%sc1=5e-3;
sc1=4/1000;
Dym=phi1/c/sc1*1e-7;

%{
figure(1); clf(1)
nx=300; x=linspace(-.001,1,nx)';
bT=G(x);
subplot(1,2,1); plot(x,bT(:,1),'b'); ylim([-.1,1.1]); hold on
%subplot(1,2,2); plot(x,bT(:,2),'b'); ylim([-.1,1.1]); hold on
%}

% Quadrature Integration Points and Weights for Volume Integrals
nip3=4; al=0.58541020; be=0.13819660;
xip3=[al be be be
    be al be be
    be be al be];
w3=ones(1,4)/4;

% {
al1=(1-be)/(al-be); be1=-be/(al-be); al2=(.5-be)/(al-be);
xip=[al1 be1 be1 be1 al2 be1 al2 al2 be1 be1
    be1 al1 be1 be1 al2 al2 be1 be1 al2 be1
    be1 be1 al1 be1 be1 al2 al2 be1 be1 al2];
Nextr=[xip' 1-sum(xip)'];
%}

R=[6 -1 -1  0 -4  0
    -1  6 -1  0  0 -4
    -1 -1  6 -4  0  0
    0  0 -4 32 16 16
    -4  0  0 16 32 16
    0 -4  0 16 16 32]/180;

coords=readmatrix('Geometry\crd05.txt','Delimiter',' ');
connect4=readmatrix('Geometry\con05.txt','Delimiter',' '); 
d=max(coords,[],2); 

model = createpde();
geometryFromMesh(model,coords,connect4);

figure(2); clf(2); view([55,-25])
subplot(1,2,1); 
pdegplot(model,"FaceAlpha",0.9); hold on
view([120,30])
subplot(1,2,2); 
pdemesh(model,"FaceAlpha",0.9)
view([120,30])
%zoom(6)

% numeration of the tetrahedron base is counterclockwise
connect4=connect4([1,3,2,4],:); me=max(connect4,[],[1,2]);

nelem=size(connect4,2); 
nelnodes=10; %nodes per element

[coords,connect,~]=T4to10(coords(:,1:me),connect4);

nnodes=size(coords,2);
ncoord=3; %problem dimension
nfnodes=6; %number of nodes per a face
eldf=ncoord*nelnodes; %degrees of freedom for element
ndof=ncoord*nnodes; %total number of degrees of freedom
faces=[1 2 3 5 6 7; 1 4 2 8 9 5; 2 4 3 9 10 6; 3 4 1 10 8 7]; %numbering of element faces
indq0=reshape(reshape(1:36,6,6)',36,1);


% Volumes of Elements
velem=zeros(nelem,1);
for i=1:nelem
    coordc=coords(:,connect(1:4,i));
    velem(i)=det([coordc;1 1 1 1]')/6;
end

% on the continuation of the crack plain
e=1e-8; b1=1.5; lcoh=3;
cind=find(coords(3,:)<e & ...
          coords(2,:)<b1+lcoh+e);

fixed1=find(coords(3,:)<e & ...
            coords(2,:)>b1+lcoh-e);
fixed2=find(coords(2,:)<e);
fixed3=find(coords(1,:)<e);

rwa=[ncoord*fixed1'-0;ncoord*fixed2'-1;ncoord*fixed3'-2]; nfix=length(rwa);

cohsym=find(coords(3,:)<e & coords(1,:)<e & coords(2,:)<b1+lcoh+e)';
cohsym=sortrows([cohsym,coords(2,cohsym)'],2); cohsym=cohsym(:,1);
ci0=find(cohsym>b1+lcoh-e,1,'first');
vc=cohsym*ncoord;

eind=find(abs(coords(3,:)-d(3))<e);

ne=0; er=zeros(1e3*(nfnodes-3),1); ej=er;
nc=0; % number of cohesive faces
cj=zeros(1e3*nfnodes,1);
cfar=zeros(1e3,1);  %areas of cohesive faces
%figure(1); clf(1)
for lmn=1:nelem
    [eC,ea]= intersect(connect(:,lmn),eind);
    [cC,ca]= intersect(connect(:,lmn),cind);
    if length(eC)==nfnodes
        ne=ne+1;
        for face=1:4
            if length(intersect(ea,faces(face,:)))==nfnodes
                nodelist=faces(face,:);
                lmncoords=coords(1:2,connect(nodelist(1:3),lmn));
                ind=(ne-1)*3+(1:3);
                ej(ind)=connect(nodelist(4:6),lmn)*ncoord;
                er(ind)=polyarea(lmncoords(1,:),lmncoords(2,:))/3;
                break;
            end
        end
    end
    if length(cC)==nfnodes
        nc=nc+1;
        for face=1:4
            if length(intersect(ca,faces(face,:)))==nfnodes
                nodelist=faces(face,:); lmncon=connect(nodelist,lmn);
                lmncoords=coords(1:2,lmncon(1:3));
                ind=(nc-1)*6+(1:6);
                cj(ind)=lmncon*ncoord;
                cfar(nc)=polyarea(lmncoords(1,:),lmncoords(2,:));
                %{
                cloads(:,nc)=[lmn;face];
                patch('Faces',1:3,'Vertices',lmncoords',...
                    'EdgeColor','blue','FaceColor',.85*[1,1,1]);
                for i=1:3
                    text(coords(1,lmncon(i)),coords(2,lmncon(i)),...
                        sprintf('%d',lmncon(i))); hold on
                end
                %}
                break;
            end
        end
    end
end
er=er(1:ne*(nfnodes-3)); ej=ej(1:ne*(nfnodes-3));
cj=cj(1:nc*nfnodes); cfar=cfar(1:nc);

ci=1; % cohesive index

cohlsym=find(abs(coords(1,:))<e & abs(coords(3,:))<e)';
tmp=sortrows([cohlsym,coords(2,cohlsym)'],2); cohlsym=tmp(:,1);
vcls=cohlsym*ncoord;

%cohlext=find(abs(coords(1,:)-d(1))<e & abs(coords(3,:))<e)';
%tmp=sortrows([cohlext,coords(2,cohlext)'],2); cohlext=tmp(:,1);
%vcle=cohlext*ncoord;

[M,Eo,Em,rhom,tEi]=ViscMod();  
zetamn=rhom; etamn=rhom;

U=zeros(ndof,1); dU=U; Fpsig=U; Ftsig=U;
S=zeros(6,6,nip3*nelem,M);

% Formation of matrix S for encremental viscoelastic scheme
sigz=zeros(nnodes,1); % Stress in Nodes
sigzip=zeros(nelem*nip3,1); % Stress in Integration Points
sigzt=zeros(nnodes,30); % Stress in Nodes Calculated for All Adjacent Elements
nodelmn=zeros(nnodes,30); % Adjacent Elements for Node

indU=zeros(eldf,1);

figure(3); clf(3)
sig=sc1*.45; 
nt=100; t1=zeros(nt+1,1); kdt=0; p=5; dU0=zeros(ndof,1);
while kdt<=nt
    
    if kdt==0
        
        K=Stif(Eo);
        
        dU=fsolve(@F0c,dU0,opt); dU=full(dU);
        kD=2*dU(vc(1))/Dym;
    else
        dt0=1e-4;
        sig0=Ucurr(dt0,kD,ci);
        fprintf('%.4f %.4f\n',sig0*1e3,sig*1e3)
        if sig0>sig
            dt=fsolve(@(dt) Ucurr(dt,kD,ci)-sig,.1,opt2);
        else, break;
        end
    end
    
    U=U+dU;
    sigzt(:)=0;
    for lmn=1:nelem
        nod=connect(:,lmn);
        for i=1:nelnodes, indU((3*i-2):(3*i))=(3*nod(i)-2):(3*nod(i)); end
        
        % Incremental viscoelastic scheme
        for j=1:nip3
            B=BN(xip3(:,j),lmn);
            deps=B*dU(indU); lmnip=(lmn-1)*nip3+j;
            
            if kdt
                tsig=0; % tilda sigma
                for m=1:M
                    tsig=tsig+sum((1-zetamn(:,:,m)).*S(:,:,lmnip,m),2);
                end
                sigzip(lmnip)=sigzip(lmnip)-tsig(3)+tE(3,:)*deps;
            else
                sigzip(lmnip)=Eo(3,:)*deps;
            end
            
            deps1=repmat(deps',6,1);
            for m=1:M
                if kdt
                    S(:,:,lmnip,m)=zetamn(:,:,m).*S(:,:,lmnip,m)+...
                        etamn(:,:,m).*(1-zetamn(:,:,m)).*deps1;
                else
                    S(:,:,lmnip,m)=Em(:,:,m).*deps1;
                end
            end
        end
        
        %Stress in the Node Calculated for all the Elements Containing This Node
        sigzv=Nextr*sigzip((lmnip-nip3+1):lmnip);
        for j=1:nelnodes
            k1=find(sigzt(nod(j),:)==0,1,'first');
            sigzt(nod(j),k1)=sigzv(j);
            if kdt==0,  nodelmn(nod(j),k1)=lmn; end
        end
    end
    
    % Averaging over the Elements Containing Node
    for j=1:nnodes
        k1=length(find(sigzt(j,:)));
        sigz(j)=sigzt(j,1:k1)*velem(nodelmn(j,1:k1))/...
            sum(velem(nodelmn(j,1:k1)))/sc1;
        if sigz(j)>1, sigz(j)=1; end
    end
    
    %integral of previous sigma
    if kdt, Fpsig=K*dU+Fpsig-Ftsig;
    else,   Fpsig=K*dU;
    end
    
    if kdt==1, t1(2)=dt;
    elseif kdt>1, t1(kdt+1)=t1(kdt)+dt;
    end
    
    fact=100; 
    Mesh.Points=coords+fact*reshape(U,3,nnodes);
    Mesh.Elements=connect;
    figure(1); clf(1)
    
    plot3D(Mesh,'ColorMapData',sigz); hold on
    plot3D(Mesh,'FaceAlpha',0,'EdgeColor',.6*[1,1,1]); hold on
    
    col=(linspace(-.1,1,12))';
    cmap=jet(length(col)-1);
    clim([col(1) col(end)]); colormap(cmap);
    colorbar('Ticks',col); hold on  
    
    view(60,-30); axis equal
    zoom(6)
    str=sprintf('stress%d.png',kdt);
    %exportgraphics(gcf,str,'Resolution',600);

    figure(3)
    subplot(1,3,1:2); 
    
    xlim([0,b1+lcoh]); hold on; ylim([-.1,3]); hold on
    Ds=2*U(vcls)/Dym; 
    %De=2*U(vcle)/Dym; plot(coords(2,cohlext),De,'b--','LineWidth',.75); hold on
    plot(coords(2,cohlsym),Ds,'k','LineWidth',.75); hold on
    plot(coords(2,cohsym(2*ci-1)),2*U(vc(2*ci-1))/Dym,'ko','MarkerSize',2); hold on
    
    subplot(1,3,3);
    plot(t1(kdt+1),coords(2,cohsym(2*ci-1)),'ko','MarkerSize',3); hold on
    if kdt>0, plot(t1(kdt:kdt+1),coords(2,cohsym(2*[ci0,ci]-1)),'k','LineWidth',.75); hold on;  end
        
    ci0=ci;
    if kdt==0, kD=ceil(p*kD)/p+1/p;
    elseif kD<1-e, kD=kD+1/p;
    else, ci=ci+1;
    end
    
    kdt=kdt+1;
    pause(.1)
end
figure(3); exportgraphics(gcf,'Disp2.eps');
%figure(5); clf(5); colormap(cmap); colorbar()
%close(v)

    function f=Ucurr(dtn,kD,ci)
        
        for m1=1:M
            zetamn(:,:,m1)=exp(-dtn./rhom(:,:,m1));
            etamn(:,:,m1)=Em(:,:,m1).*rhom(:,:,m1)/dtn;
        end
        
        Ftsig=zeros(ndof,1); % Global viscoelastic residual vector
        for lmn1=1:nelem
            nod=connect(:,lmn1);
            for j1=1:nelnodes 
                indU((3*j1-2):(3*j1))=(3*nod(j1)-2):(3*nod(j1)); 
            end
            
            % Viscoelastic residual vector of the element
            Ftsige=zeros(eldf,1);
            for j1=1:nip3
                B1=BN(xip3(:,j1),lmn1);
                tsig1=0; lmnip1=(lmn1-1)*nip3+j1;
                for m1=1:M
                    tsig1=tsig1+sum((1-zetamn(:,:,m1)).*S(:,:,lmnip1,m1),2);
                end
                Ftsige=Ftsige+w3(j1)*B1'*tsig1*velem(lmn1);
            end
            % Global viscoelastic residual vector
            Ftsig(indU)=Ftsig(indU)+Ftsige;
        end
        % Displacement boundary conditions
        Ftsig(rwa)=0;
        
        tE=tEi;
        for m1=1:M, tE=tE+etamn(:,:,m1).*(1-zetamn(:,:,m1)); end
        K=Stif(tE);
        
        f=F1(kD,ci);
        
    end

    function sign=F1(kD,ci)
        Us=fsolve(@(us) F0a(us,kD,ci),[dU0;sig/4],opt);
        sign=Us(ndof+1); dU=full(Us(1:ndof));
    end

    function B=BN(xi,lmn)
        
        xi4=1-sum(xi); 
        
        dNdxi=zeros(nelnodes,ncoord);
        dNdxi(1,1) = 4*xi(1)-1;
        dNdxi(2,2) = 4*xi(2)-1;
        dNdxi(3,3) = 4*xi(3)-1;
        dNdxi(4,1) =-(4*xi4-1);
        dNdxi(4,2) =-(4*xi4-1);
        dNdxi(4,3) =-(4*xi4-1);
        dNdxi(5,1) = 4*xi(2);
        dNdxi(5,2) = 4*xi(1);
        dNdxi(6,2) = 4*xi(3);
        dNdxi(6,3) = 4*xi(2);
        dNdxi(7,1) = 4*xi(3);
        dNdxi(7,3) = 4*xi(1);
        dNdxi(8,1) = 4*(xi4-xi(1));
        dNdxi(8,2) =-4*xi(1);
        dNdxi(8,3) =-4*xi(1);
        dNdxi(9,1) =-4*xi(2);
        dNdxi(9,2) = 4*(xi4-xi(2));
        dNdxi(9,3) =-4*xi(2);
        dNdxi(10,1)=-4*xi(3);
        dNdxi(10,2)=-4*xi(3);
        dNdxi(10,3)= 4*(xi4-xi(3));
        
        dxdxi=coords(:,connect(:,lmn))*dNdxi;
        dNdx=dNdxi*inv(dxdxi);
        
        B=zeros(6,eldf);
        for i1=1:nelnodes
            b=dNdx(i1,:);
            B(:,ncoord*(i1-1)+1:ncoord*i1)=...
                [b(1) 0    0
                0    b(2) 0
                0    0    b(3)
                b(2) b(1) 0
                b(3) 0    b(1)
                0    b(3) b(2)];
        end
    end

    function f=Stif(D)
        %disp('Assemble the global stiffness matrix'); tic
        
        m1=round(.1*ndof^2);
        rwi=zeros(m1,1); cli=zeros(m1,1); st=zeros(m1,1); ki=0;
        for lmn1=1:nelem    % Loop over all the elements
            
            % Set up the stiffness for the current element
            K1=zeros(eldf);
            for j1=1:nip3
                B1=BN(xip3(:,j1),lmn1);
                K1=K1+w3(j1)*B1'*D*B1*velem(lmn1);
            end
            
            % Add the current element stiffness to the global stiffness
            for i1=1:nelnodes
                for ii=1:ncoord
                    rw=ncoord*(connect(i1,lmn1)-1)+ii;
                    if ~any(rwa==rw)
                        for j1=1:nelnodes
                            for jj=1:ncoord
                                cl=ncoord*(connect(j1,lmn1)-1)+jj;
                                ki=ki+1; rwi(ki)=rw; cli(ki)=cl;
                                st(ki)=K1(ncoord*(i1-1)+ii,ncoord*(j1-1)+jj);
                            end
                        end
                    end
                end
            end
        end
        rwi=rwi(1:ki); cli=cli(1:ki); st=st(1:ki);
        
        % Modify the global stiffness and residual to include constraints
        %for i1=1:nfix, st(rwi==rwa(i1))=0; end
        rwi=[rwi;rwa]; cli=[cli;rwa]; st=[st;ones(nfix,1)];
        
        f=sparse(rwi,cli,st,ndof,ndof);
        %toc
        %disp('* * * * * * * * * * * * ');
    end

    function [f,jac]=F0a(us,kD,ci)
        
        du=us(1:ndof); sigi=us(end);
        
        D1=2*(U(cj)+du(cj))/Dym;
        tmp1=-sc1*G(D1); t=tmp1(:,1);
        
        cr=zeros(6*nc,1);
        p1=zeros(36*nc,1); q1=p1; jac=p1; ta=tmp1(:,2)*2/Dym;
        
        for i4=1:nc
            ind3=6*(i4-1)+(1:6);
            R1=cfar(i4)*R;
            
            cr(ind3)=R1*t(ind3);
            
            % Jacobi matrix
            indp=36*(i4-1)+(1:36); indq=36*(i4-1)+indq0;
            p1(indp)=repmat(cj(ind3),6,1); q1(indq)=p1(indp);
            jac(indp)=-reshape(R1.*repmat(ta(ind3)',6,1),36,1);
        end
        
        r=sparse([ej;cj],ones(3*ne+6*nc,1),[sigi*er;cr],ndof,1);
        
        f=K*du-(r-Fpsig+Ftsig);
        f=[f;U(vc(2*ci-1))+du(vc(2*ci-1))-kD*Dym/2];
        
        rw1=[rwi;p1;ej;ndof+1];
        cl1=[cli;q1; (ndof+1)*ones(ne*(nfnodes-3),1);vc(2*ci-1)];
        st1=[st; jac;-er;1];
        jac=sparse(rw1,cl1,st1,ndof+1,ndof+1);
        
    end

    function [f,jac]=F0c(du) %
        
        D1=2*du(cj)/Dym;
                
        tmp1=-sc1*G(D1); t=tmp1(:,1);
        
        cr=zeros(6*nc,1);
        p1=zeros(36*nc,1); q1=p1; jac=p1; ta=tmp1(:,2)*2/Dym;
        
        for i4=1:nc
            ind3=6*(i4-1)+(1:6);
            R1=cfar(i4)*R;
            
            cr(ind3)=R1*t(ind3);
            
            % Jacobi matrix
            indp=36*(i4-1)+(1:36); indq=36*(i4-1)+indq0;
            p1(indp)=repmat(cj(ind3),6,1); q1(indq)=p1(indp);
            jac(indp)=-reshape(R1.*repmat(ta(ind3)',6,1),36,1);
        end
        
        r=sparse([ej;cj],ones(3*ne+6*nc,1),[sig*er;cr],ndof,1);
        jac=sparse([rwi;p1],[cli;q1],[st;jac],ndof,ndof);
        
        f=K*du-(r-Fpsig+Ftsig);
        
    end

    function bT=G(bD)
        nbD=length(bD); bT=zeros(nbD,2);
        for n1=1:nbD
            y=bD(n1);
            if y<0
                bT(n1,1)=y/a0*(2-y/a0);
                bT(n1,2)=2/a0*(1-y/a0);
            elseif y<a1
                bT(n1,1)=y/a1*(2-y/a1);
                bT(n1,2)=2/a1*(1-y/a1);
            elseif y<a2, bT(n1)=1;
            elseif y<=1
                bT(n1,1)=((1+2*y-3*a2).*(1-y).^2)/(1-a2)^3;
                bT(n1,2)=-6*(1-y)*(y-a2)/(1-a2)^3;
            end
        end
    end

    function [M,Eo,Em,rhom,tEi]=ViscMod()
        % Viscoelastic parameters
        ro=20;
        %{
        E20=8;  E2i=1;  E21=E20-E2i;
        E30=24; E3i=6;  E31=E30-E3i;
        G0 =3;  Gi =.4; G1 =G0 -Gi;
        nu21=.3; nu32=.3;
        %}
        % {
        E20=4;  E2i=1;  E21=E20-E2i;
        E30=4;  E3i=1;  E31=E30-E3i;
        %G0 =3;  Gi =.4; G1 =G0 -Gi;
        nu21=.3; nu32=.3;
        G0=E30/2/(1+nu32); Gi=E3i/2/(1+nu32); G1=G0-Gi;
        %}
        
        g=nu32^2; hp=1+nu21; hm=1-nu21;
        
        aa0=E30*hm-2*E20*g;
        aa1=E31*hm-2*E21*g; bb1=E21*(E31-E21*g)/hp; d1=E21*(E31*nu21+E21*g)/hp;
        ai =E3i*hm-2*E2i*g; bi=E2i*(E3i-E2i*g)/hp; di=E2i*(E3i*nu21+E2i*g)/hp;
        
        cc=hm*g*(E2i*E30-E20*E3i)^2/aa0/aa1/ai;
        E=[ bi/ai             bb1/aa1             -cc
            di/ai             d1/aa1             -cc
            nu32*E2i*E3i/ai   nu32*E21*E31/aa1   -cc*2*nu32
            hm*E3i^2/ai       hm*E31^2/aa1       -cc*4*g
            Gi                G1               0
            E2i/hp/2          E21/hp/2         0];
        
        ro2=aa0/ai*ro;
        
        M=2;
        
        %{
        lta=0; ltb=2.2; lt1=20;
        lt=linspace(lta,ltb,lt1); t1=10.^lt;
        %dtn=1; t0=1; tn=50; t1=t0:dtn:tn; lt1=length(t1);

        e=exp(-t1/ro);
        E3=E3i+E31*e;
        E2=E2i+E21*e;
        G =Gi +G1 *e;
        nu23=E2./E3*nu32;
        EG=[E3/10;E2/10;G/10;nu23];
        %}
        
        %{
        rho=[ro,ro2
            ro,ro2
            ro,ro2
            ro,ro2
            ro,0
            ro,0];
        
        Del=hp*(1-nu21-2*nu23*nu32);
        
        D=[E2.*(1-nu23*nu32)./Del
            E2.*(nu21+nu23*nu32)./Del
            nu32*E2*hp./Del
            E3.*(1-nu21^2)./Del
            G
            E2/2/hp];
        
        
        function f=Cm(k,t)
        %{
        D11=Cm(1,t1);
        D12=Cm(2,t1);
        D13=Cm(3,t1);
        D33=Cm(4,t1);
        D44=Cm(5,t1);
        D66=Cm(6,t1);
        %}
            f=E(k,1)*ones(size(t));
            for m1=1:M
                f=f+E(k,m1+1)*exp(-t/rho(k,m1));
            end
        end
        
        %{
        figure(2); clf(2)
        for k1=1:6
            D1=Cm(k1,t1);
            subplot(1,2,1)
            plot(lt,D1);  hold on
            plot(lt,D(k1,:),'k--');
            subplot(1,2,2)
            plot(lt,D1./D(k1,:)-1,'k--');
        end
        %}
        %}
        
        rhom=zeros(6,6,M); Em=rhom;
        
        rhom(:,:,1)=...
            [ro ro ro 0 0 0
            ro ro ro 0 0 0
            ro ro ro 0 0 0
            0 0 0 ro 0 0
            0 0 0 0 ro 0
            0 0 0 0 0 ro];
        rhom(:,:,2)=...
            [ro2 ro2 ro2 0 0 0
            ro2 ro2 ro2 0 0 0
            ro2 ro2 ro2 0 0 0
            0 0 0 0 0 0
            0 0 0 0 0 0
            0 0 0 0 0 0];
        
        for m1=2:M+1
            Em(:,:,m1-1)=...
                [E(1,m1) E(2,m1) E(3,m1) 0 0 0
                E(2,m1) E(1,m1) E(3,m1) 0 0 0
                E(3,m1) E(3,m1) E(4,m1) 0 0 0
                0 0 0 E(6,m1) 0 0
                0 0 0 0 E(5,m1) 0
                0 0 0 0 0 E(5,m1)];
        end
        
        % Instantaneous Moduli
        Eonz=sum(E,2);
        Eo=[Eonz(1) Eonz(2) Eonz(3) 0 0 0
            Eonz(2) Eonz(1) Eonz(3) 0 0 0
            Eonz(3) Eonz(3) Eonz(4) 0 0 0
            0 0 0 Eonz(6) 0 0
            0 0 0 0 Eonz(5) 0
            0 0 0 0 0 Eonz(5)];
        
        % Long-time Moduli
        tEi=[E(1,1) E(2,1) E(3,1) 0 0 0
            E(2,1) E(1,1) E(3,1) 0 0 0
            E(3,1) E(3,1) E(4,1) 0 0 0
            0 0 0 E(6,1) 0 0
            0 0 0 0 E(5,1) 0
            0 0 0 0 0 E(5,1)];
        %}
        
        %{
        lt=linspace(0,2,10); ex=exp(-10.^lt/ro);
        E1=E1i+(E10-E1i)*ex;
        E2=E2i+(E20-E2i)*ex;
        G12=Gi+(G0-Gi)*ex;
        nu12=nu21*E1./E2;
        ro_o=sqrt(E1.*E2)/2./G12-sqrt(nu21*nu12);
        figure(2); clf(2)
        plot(lt,ro_o)
        %}
        
    end

    function [coords,connect,edget]=T4to10(coords4,connect4)
        edge=[1,2;2,3;3,1;1,4;2,4;3,4];
        midedge=5:10;

        connect=zeros(nelnodes,nelem); connect(1:4,:)=connect4;
        edges=zeros(6*nelem,2); edgessort=edges;
        for lmn1=1:nelem
            for i1=1:6
                ind1=(lmn1-1)*6+i1;
                edges(ind1,:)=connect4(edge(i1,:),lmn1);
                edgessort(ind1,:)=sort(edges(ind1,:));
            end
        end
        [~,~,iedge]=unique(edgessort,'rows'); ne=length(unique(iedge));

        n=max(connect4,[],[1,2]);
        coords=zeros(3,n+ne);
        coords(:,1:n)=coords4;
        edget=zeros(ne,3); % all the edges
        for i1=1:ne
            edgei=find(iedge==i1);
            edgeilmn=floor((edgei-1)/6)+1;
            edgeinum=mod(edgei,6); edgeinum(edgeinum==0)=6;
            for m1=1:length(edgei)
                nn=n+1;
                connect(midedge(edgeinum(m1)),edgeilmn(m1))=nn;
                if m1==1, edget(i1,:)=[edges(edgei(1),1) nn edges(edgei(1),2)]; end
            end
            ac=coords(:,edges(edgei(1),1));
            bc=coords(:,edges(edgei(1),2));
            coords(:,nn)=(ac+bc)/2;
            n=n+1;
        end
    end
    
end