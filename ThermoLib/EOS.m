
function [V,phi,H] = EOS(P,T,z,q,Component)
    % Constant;
    % �������峣��,J/koml/K
    R = 8.314 * 1000;
    % ��׼�¶�Tref = 298.15,K
    % ��׼ѹ��Pref = 101325,Pa
    Tref = 298.15;
    Pref = 101325;
    
    c = length(z);
    z = z./sum(z);
    Pc = ones(c,1);
    Tc = ones(c,1);
    omega = ones(c,1);
    Mw = ones(c,1);
    Hideal_i = ones(c,1);
    S_temp_i = ones(c,1);
    Sideal_i = ones(c,1);
    for i = 1:c
        Pc(i) = Component(i).Pc;
        Tc(i) = Component(i).Tc;
        omega(i) = Component(i).omega;
        Mw(i) = Component(i).Mw;
        Hideal_i(i) = integral(Component(i).IdealGasHeatCapacityCp_func, Tref, T);
        S_temp_i(i) = integral(@(T) Component(i).IdealGasHeatCapacityCp_func(T)/T, Tref, T, 'ArrayValued',true);
        Sideal_i(i) = Component(i).Sform;
    end
    
    % PR ����ϵ��
    Omega_a = 0.45724;
    Omega_b = 0.0778;
    m = [0.37464, 1.54226, -0.2699];
    
    fw_i = m(1) + m(2).*omega + m(3).*omega.^2;
    alpha = (1+ fw_i.*(1-(T./Tc).^0.5) ).^2;
    ac_i = Omega_a.*R.^2.*Tc.^2./Pc;
    a_i = ac_i.* alpha; 
    b_i = Omega_b*R.*Tc./Pc;
    
    % k,��Ԫ��������, c*c����
    kij = 0;
    
    % ��Ϲ����ѡ��
    a_ij = (a_i*a_i').^0.5.*(1-kij);
    b_ij = (b_i+b_i)/2;
    
    a = z*a_ij*z';
    b = z*b_ij;
    
    % PR ����delta,varepsilon,eta��a,b�Ĺ�ϵʽ
    delta = 2*b;
    DELTA = delta*P/R/T;
    
    varepsilon = -b^2;
    VAREPSILON = varepsilon*(P/R/T)^2;
    
    eta = b;
    ETA = eta*P/R/T;
    
    A = a*P/(R*T)^2;
    B = b*P/R/T;
    
    EOS_coeff = [1, DELTA-B-1, A+VAREPSILON-DELTA*(B+1), -(VAREPSILON*(B+1)+A*ETA)];
    Z = roots(EOS_coeff);
    Z = Z(imag(Z) == 0);
    
    % q == 1, Һ��
    % q == 0, ����
    if q == 1
        Z = min(Z);
    elseif q == 0
        Z = max(Z);
    end
    
    % ���������
    V = Z*R*T/P;
    
    % PR ���̵�dadT
    dadT = -z *( fw_i.*sqrt(a_i.*ac_i./T./Tc) );
    
    % �����������
    
    % ��������������о�̬�����µ���
    Hideal = z*Hideal_i;
    % ƫ���ʵļ���
    Temp1 = T*dadT - a;
    Temp2 = R*T*(delta^2-4*varepsilon)^0.5;
    Temp3 = 2*V+delta+(delta^2-4*varepsilon)^0.5;
    Temp4 = 2*V+delta-(delta^2-4*varepsilon)^0.5;
    Temp5 = log(Temp3/Temp4);
    HR = R*T*( Temp1 / Temp2*Temp5 +Z-1 );
    % �������� = �����������+ƫ����
    H = Hideal + HR;

    % �صļ���Ŀǰ����������
    % �����������
    
    % ��������������о�̬�����µ���
    Sideal = z*(Sideal_i - R*log(P/Pref)-R*log(z') + S_temp_i);
    % ƫ���صļ���
    Temp1 =  R*(delta^2-4*varepsilon)^0.5;
    Temp2 = log(Temp3/Temp4);
    Temp5 = log(P*(V-b)/R/T);
    SR = R*(Temp5 + dadT/Temp1 * Temp2);
    S = Sideal + SR;

    % �������������i��ϵ��ϵ��
    Temp1 = b_i./b.*(Z-1)-log(Z-B);
    Temp2 = A./2^1.5/B.*(2*a_ij*z'./a - b_i./b);
    Temp3 = log( (Z+(1+sqrt(2))*B)./(Z+(1-sqrt(2))*B) );
    phi = exp(Temp1-Temp2.*Temp3);
    
end