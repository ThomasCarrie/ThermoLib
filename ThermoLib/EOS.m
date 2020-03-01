function [V,phi,H] = EOS(P,T,z,q,Component,flag)
    % Input:
        % P: ѹ��, Pa, 1*1;
        % T: �¶�, K, 1*1;
        % z: Ħ�����, -, 1*c;
        % q: ��̬, 0,����;1,Һ��, 1*1;
        % Component: ������Բ���, -, 1*c;
        % flag: EOS����, 'RK,SRK,PR', 1*1;
    % Output:
        % V: Ħ�����, m^3/kmol, 1*1;
        % phi, ��ֵ��ݶ�ϵ��, -, c*1;
        % H: ��, J/kmol, 1*1;
        % S: ��, J/koml/K, 1*1;
        
    % Constant;
    % �������峣��: R = 8314, J/koml/K
    R = Component.R;
    % ��׼�¶�: Tref = 298.15, K
    Tref = Component.Tref;
    % ��׼ѹ��: Pref = 101325, Pa;
    Pref = Component.Pref;
    
    % �����c
    c = length(z);
    % Ħ����ɹ�һ��
    z = z./sum(z);
    
    % ���Բ���
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
    
    % �Ա��¶�
    Tr_i = T./Tc_i;
    % �Ա�ѹ��
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
    
    % k,��Ԫ��������, c*c����
    k_ij = 0;
    
    % ��Ϲ����ѡ��
    a_ij = (a_i*a_i').^0.5.*(1-k_ij);
    b_ij = (b_i+b_i)/2;
    
    a = z*a_ij*z';
    b = z*b_ij;
    
    switch q
        case 0
            % ���������
            f = @(V) R*T/P + b - a/P*(V-b)/( (V+epsilon*b) * (V+sigma*b) ) - V;
            V = fsolve(f, R*T/P);
        case 1
            % Һ�������
            f = @(V) b + (V+epsilon*b) * (V+sigma*b) * (R*T + b*P - V*P)/a - V;
            V = fsolve(f, b);
    end
    
    % ������ѹ������Z
    Z = P*V/R/T;
    
    % ����q,beta,I
    q = a/b/R/T;
    beta = b*P/R/T;
    I = 1/(sigma-epsilon)*log( (Z+sigma*beta)/(Z+epsilon*beta) );
    
    
    % �����������
    
    % ��������������о�̬�����µ���
    Hideal = z*Hideal_i;
    % ƫ���ʵļ���
    dalpha_idT = -m_i.*sqrt(Tr_i./alpha_i);
    HR_i = R.*T.*(Z-1+q.*I.*(dalpha_idT-1));
    HR = z*HR_i;
    H = Hideal + HR;

    % �صļ���Ŀǰ����������
    % �����������
    SR_i = log(Z-beta) + dalpha_idT.*q.*I;
    SR = z*SR_i;

    % �������������i��ϵ��ϵ��
    db_idn = b_i;
    dq_idn = q.*(2*a_ij*z'./a-b_i./b);
    phi = exp(db_idn./b.*(Z-1) - log(Z-beta) - dq_idn.*I); 
end
