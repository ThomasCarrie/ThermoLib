function Z = EOS(P,T,z,q,Component)
    % Constant;
    R = 8.314;
    c = length(z);
    z = z./sum(z);
    Pc = ones(c,1);
    Tc = ones(c,1);
    omega = ones(c,1);
    Mw = ones(c,1);
    for i = 1:c
        Pc(i) = Component(i).Pc;
        Tc(i) = Component(i).Tc;
        omega(i) = Component(i).omega;
        Mw(i) = Component(i).Mw;
    end
    
    % PR 方程系数
    Omega_a = 0.45724;
    Omega_b = 0.0778;
    m = [0.37464, 1.54226, -0.2699];
    
    alpha = (1+ (m(1) + m(2).*omega + m(3).*omega.^2).*(1-(T./Tc).^0.5) ).^2;
    a_i = Omega_a.*R.^2.*Tc.^2./Pc .* alpha; 
    b_i = Omega_b*R.*Tc./Pc;
    a = MixRule(z,a_i,0,2);
    b = MixRule(z,b_i,0,1);
    
    % PR 方程关系式
    B = b*P/R/T;
    delta = 2*B;
    varepsilon = -B^2;
    eta = B;
    A = a*P/(R*T)^2;
    
    EOS_coeff = [1, delta-B-1, A+varepsilon-delta*(B+1), -(varepsilon*(B+1)+A*eta)];
    Z = roots(EOS_coeff);
    Z = Z(imag(Z) == 0);
    
    % q, 液相分率, [0,1]
    if q == 1
        Z = min(Z);
    elseif q == 0
        Z = max(Z);
    end
    
    function Q = MixRule(z,Q_i,kij,rule)
        % z, 为摩尔组成, 1*c
        % Q_i, 为物性参数量, c*1
        % rule, 混合规则
        % kij, 二元交互参数, c*c
        % 一般情况下, 默认i=j时, kij = 0;
        if nargin == 3
            kij = 0;
        end
        
        switch rule
            case 1
                Q = z*Q_i;
            case 2
                Q = z*( (Q_i*Q_i').^0.5 ).*(1-kij)*z';
        end
        
    end
    
end