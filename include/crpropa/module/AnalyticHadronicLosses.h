#pragma once
#include "crpropa/Module.h"
#include "crpropa/Vector3.h"
#include "crpropa/Common.h"
#include "crpropa/massDistribution/Density.h"

#include <functional>
#include <string>
#include <vector>
#include <unordered_map>

namespace crpropa {

/**
 * Continuous hadronic cooling for pions and kaons in gas (pion/kaon + p).
 * dE/dx = n_target * sigma_inel(E) * kappa(E) * E
 *
 * Tables are species-specific and may be provided in millibarn.
 * File format per line (comments with '#'):
 *   log10(E/GeV)   sigma_inel   kappa
 * Units:
 *   - Energies: J, tables use log10(E/GeV)
 *   - n_target: m^-3
 *   - sigma_inel: mb
 *   - kappa: dimensionless inelasticity per interaction
 */
class AnalyticHadronicLosses : public Module {
public:
	AnalyticHadronicLosses(ref_ptr<Density> dens, double limit = 0.1);

	// Number density (of target protons)
	double getDensityAtPosition(Vector3d &pos) const;

	// Load sigma_inel tables for particles.
	void setSigmaTableForPDG(int pdgAbs, const std::string& path);

	void setLimit(double limit); 
	double getLimit() const;

	// Module API
	void process(Candidate* candidate) const;
	std::string getDescription() const;

private:
	struct Table {
		std::vector<double> logEGeV, sigma_mb;
		bool valid() const { return logEGeV.size() >= 2 && sigma_mb.size() == logEGeV.size(); }
	};


    ref_ptr<Density> density; //< target field (proton number density)
	double limit;

	// Per-PDG (abs) tables; family fallback used if exact not found.
	std::unordered_map<int, Table> tableByPDG;

	// helpers
	static bool isPionOrKaon(int pdgAbs);
	static Table loadTableFile(const std::string& absPath);

	// read the table
	const Table* resolveTable(int pdgAbs) const;

};

} // namespace crpropa
