import re
from rdkit import Chem

def get_monoisotopic_mass(formula):
	parts = re.findall("[A-Z][a-z]?|[0-9]+", formula)
	print (parts)
	mass = 0

	for index in range(len(parts)):
		if parts[index].isnumeric():
			continue

		multiplier = int(parts[index + 1]) if len(parts) > index + 1 and parts[index + 1].isnumeric() else 1
		isotopeMass = Chem.PeriodicTable.GetMostCommonIsotopeMass(Chem.GetPeriodicTable(), parts[index])
		mass += isotopeMass * multiplier

	return mass


