function Property = Component(CAS)
    CAS = [CAS,'.json'];
    PropertyParameter = jsondecode(fileread(CAS));

    Property.CAS = PropertyParameter.CAS;
    Property.Ename = PropertyParameter.CompoundID;
    Property.StructureFormula = PropertyParameter.StructureFormula;

    % Pc,临界压力, Pa
    Property.Pc = str2double(PropertyParameter.CriticalPressure(1));
    % Tc,临界温度, K    
    Property.Tc = str2double(PropertyParameter.CriticalTemperature(1));
    % Vc, 临界体积, m3/kmol
    Property.Vc = str2double(PropertyParameter.CriticalVolume(1));
    % Zc, 临界压缩因子, -
    Property.Zc = str2double(PropertyParameter.CriticalCompressibility(1));
    % Tb, 常压沸点温度, K
    Property.Tb = str2double(PropertyParameter.NormalBoilingPointTemperature(1));
    % Ttriple, 三相点温度, K
    Property.Ttriple = str2double(PropertyParameter.TriplePointTemperature(1));
    % Ptriple, 三相点压力, Pa
    Property.Ptriple = str2double(PropertyParameter.TriplePointPressure(1));
    % Mw, 摩尔质量, kg/kmol
    Property.Mw = str2double(PropertyParameter.MolecularWeight(1));
    % Omega, 偏心因子, -
    Property.Omega = str2double(PropertyParameter.AcentricityFactor(1));
    % Hform, 标准摩尔生成焓, J/kmol
    Property.Hform = str2double(PropertyParameter.HeatOfFormation(1));
    % Gform, 标准摩尔生成自由能, J/kmol
    Property.Gform = str2double(PropertyParameter.GibbsEnergyOfFormation(1));
    
    % LiquidDensity, 液体的密度, kmol/m3
    % LiquidDensity_func, DIPPR函数
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