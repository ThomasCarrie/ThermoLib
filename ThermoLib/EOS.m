function [V,phi,H] = EOS(P,T,z,q,Component,flag)
    % Input:
        % P: 压力, Pa, 1*1;
        % T: 温度, K, 1*1;
        % z: 摩尔组成, -, 1*c;
        % q: 相态, 0,气相;1,液相, 1*1;
        % Component: 组分物性参数, -, 1*c;
        % flag: EOS方程, 'RK,SRK,PR', 1*1;
    % Output:
        % V: 摩尔体积, m^3/kmol, 1*1;
        % phi, 组分的逸度系数, -, c*1;
        % H: 焓, J/kmol, 1*1;
        % S: 熵, J/koml/K, 1*1;
        
    % Constant;
    % 理想气体常数: R = 8314, J/koml/K
    R = Component.R;
    % 基准温度: Tref = 298.15, K
    Tref = Component.Tref;
    % 基准压力: Pref = 101325, Pa;
    Pref = Component.Pref;
    
    % 组分数c
    c = length(z);
    % 摩尔组成归一化
    z = z./sum(z);
    
    % 物性参数
    Pc_i = ones(c,1);
    Tc_i = ones(c,1);
    omega_i = ones(c,1);
    Mw_i = ones(c,1);
    Hideal_i = ones(c,1);
    S_temp_i = ones(c,1);
    for i = 1:c
        Pc_i(i) = Component(i).Pc;
        Tc_i(i) = Component(i).Tc;
        omega_i(i) = Component(i).omega;
        Mw_i(i) = Component(i).Mw;
        Hideal_i(i) = integral(Component(i).IdealGasHeatCapacityCp_func, Tref, T);
        S_temp_i(i) = integral(@(T) Component(i).IdealGasHeatCapacityCp_func(T)/T, Tref, T, 'ArrayValued',true);
    end
    
    % 对比温度
    Tr_i = T./Tc_i;
    % 对比压力
    Pr_i = P./Pc_i;
    switch flag
        case 'RK'
            sigma = 1;
            epsilon = 0;
            Omega = 0.08664;
            Psi = 0.42748;
            m_i = (Tr_i.^(-0.25) - 1)./(1-Tr.^0.5);
        case 'SRK'
            sigma = 1;
            epsilon = 0;
            Omega = 0.08664;
            Psi = 0.42748;
            m = [0.480, 1.574, -0.176];
            m_i = [ones(c,1) omega_i omega_i.^2]*m';
        case 'PR'
            sigma = 1+sqrt(2);
            epsilon = 1-sqrt(2);
            Omega = 0.07780;
            Psi = 0.45724;
            m = [0.37464, 1.54226, -0.26992];
            m_i = [ones(c,1) omega_i omega_i.^2]*m';
    end
    alpha_i = (1+m_i.*(1-Tr_i.^0.5) ).^2;
    
    a_i = Psi.*R.^2.*Tc_i.^2./Pc_i.* alpha_i; 
    b_i = Omega*R.*Tc_i./Pc_i;
    
    % k,二元交互参数, c*c矩阵
    k_ij = 0;
    
    % 混合规则的选择
    a_ij = (a_i*a_i').^0.5.*(1-k_ij);
    b_ij = (b_i+b_i)/2;
    
    a = z*a_ij*z';
    b = z*b_ij;
    
    switch q
        case 0
            % 气相体积根
            f = @(V) R*T/P + b - a/P*(V-b)/( (V+epsilon*b) * (V+sigma*b) ) - V;
            V = fsolve(f, R*T/P);
        case 1
            % 液相体积根
            f = @(V) b + (V+epsilon*b) * (V+sigma*b) * (R*T + b*P - V*P)/a - V;
            V = fsolve(f, b);
    end
    
    % 混合物的压缩因子Z
    Z = P*V/R/T;
    
    % 定义q,beta,I
    q = a/b/R/T;
    beta = b*P/R/T;
    I = 1/(sigma-epsilon)*log( (Z+sigma*beta)/(Z+epsilon*beta) );
    
    
    % 计算混合物的焓
    
    % 混合理想气体在研究态工况下的焓
    Hideal = z*Hideal_i;
    % 偏离焓的计算
    dalpha_idT = -m_i.*sqrt(Tr_i./alpha_i);
    HR_i = R.*T.*(Z-1+q.*I.*(dalpha_idT-1));
    HR = z*HR_i;
    H = Hideal + HR;

    % 熵的计算目前还存在问题
    % 计算混合物的熵
    SR_i = log(Z-beta) + dalpha_idT.*q.*I;
    SR = z*SR_i;

    % 计算混合物中组分i的系数系数
    db_idn = b_i;
    dq_idn = q.*(2*a_ij*z'./a-b_i./b);
    phi = exp(db_idn./b.*(Z-1) - log(Z-beta) - dq_idn.*I); 
end
