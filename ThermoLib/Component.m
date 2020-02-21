function Property = Component(CAS)
    CAS = [CAS,'.json'];
    PropertyParameter = jsondecode(fileread(CAS));

    Property.CAS = PropertyParameter.CAS;
    Property.Ename = PropertyParameter.CompoundID;
    Property.StructureFormula = PropertyParameter.StructureFormula;

    % Pc,�ٽ�ѹ��, Pa
    Property.Pc = str2double(PropertyParameter.CriticalPressure(1));
    % Tc,�ٽ��¶�, K    
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
    equation = PropertyParameter.LiquidDensity(1);
    A = str2double(equation{1}.A);
    B = str2double(equation{1}.B);
    C = str2double(equation{1}.C);
    D = str2double(equation{1}.D);
    Tmin = str2double(equation{1}.Tmin);
    Tmax = str2double(equation{1}.Tmax);
    Property.LiquidDensity = [A, B, C, D, Tmin, Tmax];
    Property.LiquidDensity_func = @(T) A ./ B.^(1+ (1-T./C).^D);
end