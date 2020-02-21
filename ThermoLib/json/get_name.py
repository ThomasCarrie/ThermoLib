import os 
import json
import numpy as np 

def file_name(file_dir):   
	for root, dirs, files in os.walk(file_dir):  
		# print(root) #当前目录路径  
		# print(dirs) #当前路径下所有子目录  
		name_list = files[0:-1] #当前路径下所有非目录子文件
		name = []
		for f in name_list:
			name_tmp = os.path.splitext(f)
			name.append(name_tmp[0])
	return name

path = os.getcwd()
name = file_name(path)
n = 0 
for i in name:
	CAS = i
	with open(CAS + '.json') as f:
		data = json.load(f)
		f.close()

	print(data['CAS'])
	print(data['CompoundID'])
	print(data['StructureFormula'])
	print(data['CriticalTemperature'][0])
	print(data['CriticalPressure'][0])
	print(data['CriticalVolume'][0])
	print(data['CriticalCompressibility'][0])
	print(data['NormalBoilingPointTemperature'][0])
	print(data['NormalMeltingPointTemperature'][0])
	print(data['TriplePointTemperature'][0])
	print(data['TriplePointPressure'][0])
	print(data['MolecularWeight'][0])
	print(data['LiquidVolumeAtNormalBoilingPoint'][0])
	print(data['AcentricityFactor'][0])
	print(data['SolubilityParameter'][0])
	print(data['DipoleMoment'][0])
	print(data['HeatOfFormation'][0])
	print(data['GibbsEnergyOfFormation'][0])
	print(data['AbsEntropy'][0])
	print(data['HeatOfFusionAtMeltingPoint'][0])
	print(data['HeatOfCombustion'][0])

	Tmin = float(data['LiquidDensity'][0]['Tmin'])
	Tmax = float(data['LiquidDensity'][0]['Tmax'])

	property_name = 'LiquidDensity'
	# property_name = 'VaporPressure'
	# property_name = 'HeatOfVaporization'
	# property_name = 'LiquidHeatCapacityCp'
	# property_name = 'IdealGasHeatCapacityCp'
	# property_name = 'SecondVirialCoefficient'
	# property_name = 'LiquidViscosity'
	# property_name = 'VaporViscosity'
	# property_name = 'LiquidThermalConductivity'
	# property_name = 'VaporThermalConductivity'
	# property_name = 'RPPHeatCapacityCp'
	# property_name = 'RelativeStaticPermittivity'
	# property_name = 'LiquidDensity'
	# property_name = 'LiquidViscosityRPS'

	T = np.arange(Tmin,Tmax,1)
	Value = Thermo_Prop(property_name,T)