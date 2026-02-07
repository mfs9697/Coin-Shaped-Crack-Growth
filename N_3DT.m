 function N=N_3DT(xi,nelnodes)
        x1=xi(1); x2=xi(2); x3=xi(3); x4=xi(4);
        if nelnodes==10
            N=[2*x1*(x1-.5)
                2*x2*(x2-.5)
                2*x3*(x3-.5)
                2*x4*(x4-.5)
                4*x1*x2
                4*x2*x3
                4*x3*x1
                4*x1*x4
                4*x2*x4
                4*x3*x4];
        else
            xp1=x1*(x1-1/3);
            xp2=x2*(x2-1/3);
            xp3=x3*(x3-1/3);
            xp4=x4*(x4-1/3);
            
            N=zeros(20,1);
            
            N(1:4)=4.5*[xp1*(x1-2/3)
                xp2*(x2-2/3)
                xp3*(x3-2/3)
                xp4*(x4-2/3)];
            
            N(5:16)=13.5*[xp1*x2
                xp2*x1
                xp2*x3
                xp3*x2
                xp3*x1
                xp1*x3
                xp1*x4
                xp4*x1
                xp2*x4
                xp4*x2
                xp3*x4
                xp4*x3];
            
            N(17:20)=27*[x1*x2*x3
                x2*x1*x4
                x3*x4*x2
                x4*x3*x1];
        end
    end