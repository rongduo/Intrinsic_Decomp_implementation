function [R,D,O,C] =solution1_Ninner_pcg(n,s,I,t1,t2,t3,t4,t5,t6,W, Q1,Q2,Q3,Q4,Q5,R,D,O,C) 



K=[1 1 1];

mu1 = 0.2;
gamma1 = 1.2;
mu2 = 0.1;
gamma2 = 1.2;

factor_1=0.6;
factor_2=0.7;

[m1, n1 ] = size(R);
[u1, v1 ] = size(Q1);
[m2, n2 ] = size(D);
[u2, v2 ] = size(Q2);
Y1 = zeros(u1, n1);
Q1t = Q1';
Q1tQ1 = Q1t*Q1;
Y2 = zeros(u2, n2);

Q2t = Q2';
Q2tQ2 = Q2t*Q2;
Q3tQ3 = Q3'*Q3;
Q4tQ4 = Q4'*Q4;
Q5tQ5 = Q5'*Q5;
WtW = W'*W;

maxitr=3; % number of iterations
for ii = 1:maxitr
    Fr=I-C-D*K-O*K;
    if(ii>1)
       t1=t1*factor_1;
       t2=t2*factor_2;
    end  
    % ==================
    %   Shrinkage Step
    % ==================
    Q1R = Q1*R;
    Z1 = Q1R - Y1/mu1;
    B1 = max(abs(Z1) - (t1/mu1), 0).*sign(Z1);  
    Cr= (mu1*Q1tQ1 + 2*t6*WtW);
    dr=(2*t6*WtW*Fr + Q1t*Y1+ mu1*Q1t*B1);

    R(:,1)=pcg(Cr,dr(:,1));
    R(:,2)=pcg(Cr,dr(:,2));
    R(:,3)=pcg(Cr,dr(:,3));

    Y1 = Y1 + mu1*(B1 - Q1*R);
    mu1 = gamma1*mu1;

    %compute D,l1 norm
    Fd=I-R-C-O*K;
    Q2D = Q2*D;
    Z2 = Q2D - Y2/mu2;
    B2 = max(abs(Z2) - (t2/mu2), 0).*sign(Z2);

    Cd=6*t6*WtW+mu2*Q2tQ2;
    dd=2*t6*WtW*Fd*K'+Q2t*Y2+mu2*Q2t*B2;
    D=pcg(Cd,dd);
    Y2 = Y2 + mu2*(B2 - Q2*D);
    mu2 = gamma2*mu2;

    % compute O ,l2 norm    
    Fo=I-R-C-D*K;
    Co=3*t6*WtW+t3*Q3tQ3+t4*Q4tQ4;
    do=t6*WtW*Fo*K';
    N=pcg(Co,do);

    %compute C,l2 norm   
    Fc=I-R-D*K-O*K;
    Cc=t6*WtW+t5*Q5tQ5;
    dc=t6*WtW*Fc;
    C(:,1)=pcg(Cc,dc(:,1));
    C(:,2)=pcg(Cc,dc(:,2));
    C(:,3)=pcg(Cc,dc(:,3));

end  
end


