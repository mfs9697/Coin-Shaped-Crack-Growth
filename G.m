function bT = G(bD,a1,a2)
% G  Dimensionless trapezoidal cohesive shape function
%     g(bD) in [0,1]

nbD=length(bD); bT=zeros(nbD,2);
for n1=1:nbD
    y=bD(n1);
    if y<a1
        bT(n1,1)=y/a1*(2-y/a1);
        bT(n1,2)=2/a1*(1-y/a1);
    elseif y<a2, bT(n1)=1;
    elseif y<=1
        bT(n1,1)=((1+2*y-3*a2).*(1-y).^2)/(1-a2)^3;
        bT(n1,2)=-6*(1-y)*(y-a2)/(1-a2)^3;
    end
end

end
