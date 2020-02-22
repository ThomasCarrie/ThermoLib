function Property = Component(CAS)
    CAS = [CAS, '.json'];
    PropertyParameter = jsondecode(fileread(CAS));

    Property.CAS = PropertyParameter.CAS;
    Property.Ename = PropertyParameter.CompoundID;
    Property.StructureFormula = PropertyParameter.StructureFormula;

    % Pc, �ٽ�ѹ��, Pa
    Property.Pc = str2double(PropertyParameter.CriticalPressure(1));
    % Tc, �ٽ��¶�, K    
    Property.Tc = str2double(PropertyParameter.CriticalTemperature(1));
    % Vc, �ٽ����, m3/kmol
    Property.Vc = str2double(PropertyParameter.CriticalVolume(1));
    % Zc, �ٽ�ѹ������, -
    Property.Zc = str2double(PropertyParameter.CriticalCompressibility(1));
    % Tb, ��ѹ�е��¶�, K
    Property.Tb = str2double(PropertyParameter.NormalBoilingPointTemperature(1));
    % Ttriple, ������¶�, K
    Property.Ttriple = str2double(PropertyParameter.TriplePointTemperature(1));
    % Ptriple, �����ѹ��, Pa
    Property.Ptriple = str2double(PropertyParameter.TriplePointPressure(1));
    % Mw, Ħ������, kg/kmol
    Property.Mw = str2double(PropertyParameter.MolecularWeight(1));
    % Omega, ƫ������, -
    Property.Omega = str2double(PropertyParameter.AcentricityFactor(1));
    % Hform, ��׼Ħ��������, J/kmol
    Property.Hform = str2double(PropertyParameter.HeatOfFormation(1));
    % Gform, ��׼Ħ������������, J/kmol
    Property.Gform = str2double(PropertyParameter.GibbsEnergyOfFormation(1));
    
    % LiquidDensity, Һ����ܶ�, kmol/m3
    % LiquidDensity_func, DIPPR����
    PropertyName = 'LiquidDensity';
    group = PropertyParameter.LiquidDensity(1);
    equationNo = group{1}.eqno;
    [Property.LiquidDensity, Property.LiquidDensity_func] = DIPPR(equationNo);
    
    % VaporPressure, ��������ѹ, Pa
    % VaporPressure_func, DIPPR����
    PropertyName = 'VaporPressure';
    group = PropertyParameter.VaporPressure(1);
    equationNo = group{1}.eqno;
    [Property.VaporPressure, Property.VaporPressure_func] = DIPPR(equationNo);
    
    % HeatOfVaporization, ��׼Ħ��������, J/kmol
    % HeatOfVaporization_func, DIPPR����
    PropertyName = 'HeatOfVaporization';
    group = PropertyParameter.HeatOfVaporization(1);
    equationNo = group{1}.eqno;
    [Property.HeatOfVaporization, Property.HeatOfVaporization_func] = DIPPR(equationNo);
    
    % LiquidHeatCapacityCp, Һ���ѹ������, J/kmol/K
    % LiquidHeatCapacityCp_func, DIPPR����
    PropertyName = 'LiquidHeatCapacityCp';
    group = PropertyParameter.LiquidHeatCapacityCp(1);
    equationNo = group{1}.eqno;
    [Property.LiquidHeatCapacityCp, Property.LiquidHeatCapacityCp_func] = DIPPR(equationNo);
    
    % IdealGasHeatCapacityCp, ���������ѹ������, J/kmol/K
    % IdealGasHeatCapacityCp_func, DIPPR����
    PropertyName = 'IdealGasHeatCapacityCp';
    group = PropertyParameter.IdealGasHeatCapacityCp(1);
    equationNo = group{1}.eqno;
    [Property.IdealGasHeatCapacityCp, Property.IdealGasHeatCapacityCp_func] = DIPPR(equationNo);
    
    % LiquidViscosity, Һ��ճ��, Pa.s
    % LiquidViscosity, DIPPR����
    PropertyName = 'LiquidViscosity';
    group = PropertyParameter.LiquidViscosity(1);
    equationNo = group{1}.eqno;
    [Property.LiquidViscosity, Property.LiquidViscosity_func] = DIPPR(equationNo);
    
    % VaporViscosity, ����ճ��, Pa.s
    % VaporViscosity_func, DIPPR����
    PropertyName = 'VaporViscosity';
    group = PropertyParameter.VaporViscosity(1);
    equationNo = group{1}.eqno;
    [Property.VaporViscosity, Property.VaporViscosity_func] = DIPPR(equationNo);
    
    % LiquidThermalConductivity, Һ�嵼��ϵ��, W/m/K
    % LiquidThermalConductivity_func, DIPPR����
    PropertyName = 'LiquidThermalConductivity';
    group = PropertyParameter.LiquidThermalConductivity(1);
    equationNo = group{1}.eqno;
    [Property.LiquidThermalConductivity, Property.LiquidThermalConductivity_func] = DIPPR(equationNo);
    
    % VaporThermalConductivity, ���嵼��ϵ��, W/m/K
    % VaporThermalConductivity_func, DIPPR����
    PropertyName = 'VaporThermalConductivity';
    group = PropertyParameter.VaporThermalConductivity(1);
    equationNo = group{1}.eqno;
    [Property.VaporThermalConductivity, Property.VaporThermalConductivity_func] = DIPPR(equationNo);
    
    % AntoineVaporPressure, ����������ѹ, Pa
    % AntoineVaporPressure_func, DIPPR����
    PropertyName = 'AntoineVaporPressure';
    group = PropertyParameter.AntoineVaporPressure(1);
    equationNo = group{1}.eqno;
    [Property.AntoineVaporPressure, Property.AntoineVaporPressure_func] = DIPPR(equationNo);
    
    
    % SurfaceTension, ��������, N/m
    % SurfaceTension_func, DIPPR����
    PropertyName = 'SurfaceTension';
    group = PropertyParameter.SurfaceTension(1);
    equationNo = group{1}.eqno;
    [Property.SurfaceTension, Property.SurfaceTension_func] = DIPPR(equationNo);
    
    function [para, func] = DIPPR(equationNo)
        equation = PropertyParameter.(PropertyName)(1);
        Tmin = str2double(equation{1}.Tmin);
        Tmax = str2double(equation{1}.Tmax);
        switch equationNo
            case '1'
                A = str2double(equation{1}.A);
                para = [A,Tmin,Tmax];
                func = @(T) A;
            case '2'
                A = str2double(equation{1}.A);
                B = str2double(equation{1}.B);
                para = [A,B,Tmin,Tmax];
                func = @(T) A + B.*T;
            case '3'
                A = str2double(equation{1}.A);
                B = str2double(equation{1}.B);
                C = str2double(equation{1}.C);
                para = [A,B,C,Tmin,Tmax];
                func = @(T) A + B.*T + C.*T.^2;
            case '4'
                A = str2double(equation{1}.A);
                B = str2double(equation{1}.B);
                C = str2double(equation{1}.C);
                D = str2double(equation{1}.D);
                para = [A,B,C,D,Tmin,Tmax];
                func = @(T) A + B.*T + C.*T.^2 + D.*T.^3;
            case '5'
                A = str2double(equation{1}.A);
                B = str2double(equation{1}.B);
                C = str2double(equation{1}.C);
                D = str2double(equation{1}.D);
                E = str2double(equation{1}.E);
                para = [A,B,C,D,E,Tmin,Tmax];
                func = @(T) A + B.*T + C.*T.^2 + D.*T.^3 + E.*T.^4;
            case '6'
                A = str2double(equation{1}.A);
                B = str2double(equation{1}.B);
                C = str2double(equation{1}.C);
                D = str2double(equation{1}.D);
                E = str2double(equation{1}.E);
                para = [A,B,C,D,E,Tmin,Tmax];
                func = @(T) A + B.*T + C.*T.^2 + D.*T.^3 + E./T.^2;
            case '10'
                A = str2double(equation{1}.A);
                B = str2double(equation{1}.B);
                C = str2double(equation{1}.C);
                para = [A,B,C,Tmin,Tmax];
                func = @(T) exp(A - B./(C+T)); 
            case '16'
                A = str2double(equation{1}.A);
                B = str2double(equation{1}.B);
                C = str2double(equation{1}.C);
                D = str2double(equation{1}.D);
                E = str2double(equation{1}.E);
                para = [A,B,C,D,E,Tmin,Tmax];
                func = @(T) A + exp(B./T + C + D*T + E.*T.^2);
            case '100'
                A = str2double(equation{1}.A);
                B = str2double(equation{1}.B);
                C = str2double(equation{1}.C);
                D = str2double(equation{1}.D);
                E = str2double(equation{1}.E);
                para = [A,B,C,D,E,Tmin,Tmax];
                func = @(T) A + B.*T + C.*T.^2 + D.*T.^3 + E.*T.^4;
            case '101'
                A = str2double(equation{1}.A);
                B = str2double(equation{1}.B);
                C = str2double(equation{1}.C);
                D = str2double(equation{1}.D);
                E = str2double(equation{1}.E);
                para = [A,B,C,D,E,Tmin,Tmax];
                func = @(T) exp(A + B./T + C.*log(T) + D.*T.^E);
            case '102'
                A = str2double(equation{1}.A);
                B = str2double(equation{1}.B);
                C = str2double(equation{1}.C);
                D = str2double(equation{1}.D);
                para = [A,B,C,D,Tmin,Tmax];
                func = @(T) A.*T.^B./(1 + C./T + D./T.^2);
            case '103'
                A = str2double(equation{1}.A);
                B = str2double(equation{1}.B);
                C = str2double(equation{1}.C);
                D = str2double(equation{1}.D);
                para = [A,B,C,D,Tmin,Tmax];
                func = @(T) A + B.*exp(-C./T.^D);
            case '104'
                A = str2double(equation{1}.A);
                B = str2double(equation{1}.B);
                C = str2double(equation{1}.C);
                D = str2double(equation{1}.D);
                E = str2double(equation{1}.E);
                para = [A,B,C,D,E,Tmin,Tmax];
                func = @(T) A + B./T + C*10^6./T.^3 + D*10^16./T.^8 + E*10^18./T.^9;
            case '105'
                A = str2double(equation{1}.A);
                B = str2double(equation{1}.B);
                C = str2double(equation{1}.C);
                D = str2double(equation{1}.D);
                para = [A,B,C,D,Tmin,Tmax];
                func = @(T) A./B.^(1+(1-T/C).^D);
            case '106'
                A = str2double(equation{1}.A);
                B = str2double(equation{1}.B);
                C = str2double(equation{1}.C);
                D = str2double(equation{1}.D);
                E = str2double(equation{1}.E);
                para = [A,B,C,D,E,Tmin,Tmax];
                func = @(T) A.*(1-T./Property.Tc).^(B + C.*(T./Property.Tc) + D.*(T./Property.Tc).^2 + E.*(T./Property.Tc).^3);
            case '107'
                A = str2double(equation{1}.A);
                B = str2double(equation{1}.B);
                C = str2double(equation{1}.C);
                D = str2double(equation{1}.D);
                E = str2double(equation{1}.E);
                para = [A,B,C,D,E,Tmin,Tmax];
                func = @(T) A + B.*(C./T./sinh(C./T)).^2 + E.*(D./T./cosh(D./T)).^2;
            case '114'
                A = str2double(equation{1}.A);
                B = str2double(equation{1}.B);
                C = str2double(equation{1}.C);
                D = str2double(equation{1}.D);
                para = [A,B,C,D,Tmin,Tmax];
                func = @(T) A.*T + B.*T.^2/2 + C.*T.^3/3 + D.*T.^4/4;
            case '117'
                A = str2double(equation{1}.A);
                B = str2double(equation{1}.B);
                C = str2double(equation{1}.C);
                D = str2double(equation{1}.D);
                E = str2double(equation{1}.E);
                para = [A,B,C,D,E,Tmin,Tmax];
                func = @(T) A.*T + B.*(C./T)./tanh(C./T) - D.*(E./T)./tanh(E./T);
        end   
    end
    
end