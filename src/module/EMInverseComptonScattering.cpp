#include "crpropa/module/EMInverseComptonScattering.h"
#include "crpropa/Units.h"
#include "crpropa/Random.h"
#include "crpropa/Common.h"

#include <algorithm>
#include <vector>
#include <cmath>
#include <stdexcept>

namespace crpropa {

static const double mec2 = mass_electron * c_squared;

EMInverseComptonScattering::EMInverseComptonScattering(ref_ptr<PhotonField> photonField, bool havePhotons, double thinning, double limit, ref_ptr<Surface> surface) {
	setSurface(surface);
	setPhotonField(photonField);
	setHavePhotons(havePhotons);
	setLimit(limit);
	setThinning(thinning);
}

EMInverseComptonScattering::EMInverseComptonScattering(ref_ptr<PhotonField> photonField, ref_ptr<InteractionRates> interactionRates, bool havePhotons, double thinning, double limit) {
	if (!photonField)
		throw std::invalid_argument("EMInverseComptonScattering: the photon field must not be null");
	if (!interactionRates)
		throw std::invalid_argument("EMInverseComptonScattering: the interaction rates must not be null");
	// as setPhotonField(), but with the given rates instead of tables loaded by the field name
	this->photonField = photonField;
	setDescription("EMInverseComptonScattering: " + photonField->getFieldName());
	setInteractionRates(interactionRates);
	setHavePhotons(havePhotons);
	setLimit(limit);
	setThinning(thinning);
}

void EMInverseComptonScattering::setPhotonField(ref_ptr<PhotonField> photonField) {
	
	this->photonField = photonField;
	std::string fname = photonField->getFieldName();
	setDescription("EMInverseComptonScattering: " + fname);
	
	if (!this->photonField->hasPositionDependence()) {
		this->interactionRates = new InteractionRatesHomogeneous(
			getDataPath("EMInverseComptonScattering/rate_" + fname + ".txt"),
			getDataPath("EMInverseComptonScattering/cdf_" + fname + ".txt")
		);
	} else {
		this->interactionRates = new InteractionRatesPositionDependent(
			getDataPath("EMInverseComptonScattering/" + fname + "/Rate/"),
			getDataPath("EMInverseComptonScattering/" + fname + "/CumulativeRate/"),
			this->surface
		);
	}
}

void EMInverseComptonScattering::setHavePhotons(bool havePhotons) {
	this->havePhotons = havePhotons;
}

void EMInverseComptonScattering::setLimit(double limit) {
	this->limit = limit;
}

void EMInverseComptonScattering::setThinning(double thinning) {
	this->thinning = thinning;
}

void EMInverseComptonScattering::setSurface(ref_ptr<Surface> surface) {
	this->surface = surface;
}

ref_ptr<Surface> EMInverseComptonScattering::getSurface() const {
	return this->surface;
}

void EMInverseComptonScattering::setInteractionRates(ref_ptr<InteractionRates> intRates) {
	this->interactionRates = intRates;
}

ref_ptr<InteractionRates> EMInverseComptonScattering::getInteractionRates() const {
	return this->interactionRates;
}

void EMInverseComptonScattering::initRate(std::string path) {
	this->interactionRates->initRate(path);
}

void EMInverseComptonScattering::initCumulativeRate(std::string path) {
	this->interactionRates->initCumulativeRate(path);
}

// Energy distribution of the electron after the scattering, and sampling from it
ICSSecondariesEnergyDistribution::ICSSecondariesEnergyDistribution() {
	Nrer = 1000;
	s_min = mec2 * mec2;
	// bins of the original table (1000 bins from s_min to 2e23 eV^2), continued up to at least 1e30 eV^2, so that
	// the bins below 2e23 eV^2 are unchanged
	const double dls = (log(2e23 * eV * eV) - log(s_min)) / 1000;
	Ns = static_cast<size_t>(std::ceil(log(1e30 * eV * eV / s_min) / dls));
	data = std::vector< std::vector<double> >(Ns, std::vector<double>(Nrer));
	std::vector<double> data_i(Nrer);

	// tabulate s bin borders
	s_values = std::vector<double>(Ns + 1);
	for (size_t i = 0; i < Ns + 1; ++i)
		s_values[i] = s_min * exp(i*dls);

	// for each s tabulate cumulative differential cross section
	for (size_t i = 0; i < Ns; i++) {
		// r = s_min / s at the bin centre; x0 = (1 - beta) / (1 + beta) = r exactly
		double r = exp(-(i + 0.5) * dls);
		double x0 = r;
		double dlx = -log(x0) / Nrer;

		// cumulative midpoint integration
		data_i[0] = differentialCrossSection(x0, r) * expm1(dlx);
		for (size_t j = 1; j < Nrer; j++) {
			double x = x0 * exp((j+0.5) * dlx);
			double dx = exp((j+1) * dlx) - exp(j * dlx);
			data_i[j] = differentialCrossSection(x, r) * dx;
			data_i[j] += data_i[j-1];
		}
		data[i] = data_i;
	}
}

// differential cross-section, see Lee '96 (arXiv:9604098), eq. 23 for x = Ee'/Ee
double ICSSecondariesEnergyDistribution::dSigmadE(double x, double beta) {
	return differentialCrossSection(x, (1 - beta) / (1 + beta));
}

double ICSSecondariesEnergyDistribution::differentialCrossSection(double x, double r) const {
	// with beta = (1 - r) / (1 + r): (1 - beta) / beta = 2 r / (1 - r) and (1 + beta) / beta = 2 / (1 - r)
	double q = (2 * r / (1 - r)) * (1 - 1./x);
	return (2 / (1 - r)) * (x + 1./x + 2 * q + q * q);
}

// draw random energy of the electron after the scattering, Ee'(Ee, s)
double ICSSecondariesEnergyDistribution::sample(double Ee, double s) {
	// bin idx holds s_values[idx] <= s < s_values[idx + 1]; below s_min the first bin, above s_max the last one
	size_t idx = std::upper_bound(s_values.begin(), s_values.end(), s) - s_values.begin();
	idx = (idx == 0) ? 0 : std::min(idx - 1, data.size() - 1);
	std::vector<double> s0 = data[idx];
	Random &random = Random::instance();
	size_t j = random.randBin(s0) + 1; // draw random bin (upper bin boundary returned)
	double x0 = std::min(1., s_min / s); // (1 - beta) / (1 + beta)
	double dlx = -log(x0) / Nrer;
	double binWidth = x0 * (exp(j * dlx) - exp((j-1) * dlx));
	double Ep = (x0 * exp((j-1) * dlx) + binWidth) * Ee;
	return std::min(Ee, Ep); // prevent Ep > Ee from numerical inaccuracies
}

double EMInverseComptonScattering::getRate(double E, const Vector3d &position, double z) const {
	return this->interactionRates->getProcessRate(E, position) * pow_integer<2>(1 + z) * photonField->getRedshiftScaling(z);
}

void EMInverseComptonScattering::performInteraction(Candidate *candidate) const {

	// scale the particle energy instead of background photons
	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy() * (1 + z);
	
	Vector3d position = candidate->current.getPosition();
	
	std::vector<double> tabE;
	std::vector<double> tabs;
	std::vector<std::vector<double>> tabCDF;
	
	this->interactionRates->loadPerformInteractionTabs(position, tabE, tabs, tabCDF);
	
	if (E < tabE.front() or E > tabE.back())
		return;
	
	// sample the value of s
	Random &random = Random::instance();
	size_t i = closestIndex(E, tabE);
	size_t j = random.randBin(tabCDF[i]);
	double s_kin = pow(10, log10(tabs[j]) + (random.rand() - 0.5) * 0.1);
	double s = s_kin + mec2 * mec2;
	
	// sample electron energy after scattering
	static ICSSecondariesEnergyDistribution distribution;
	double Enew = distribution.sample(E, s);
	
	// add up-scattered photon
	double Esecondary = E - Enew;
	double f = Enew / E;
	if (havePhotons) {
		if (random.rand() < pow(1 - f, thinning)) {
			double w = 1. / pow(1 - f, thinning);
			Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
			candidate->addSecondary(22, Esecondary / (1 + z), pos, w, interactionTag);
		}
	}
	
	// update the primary particle energy; do this after adding the secondary to correctly set the secondary's parent
	candidate->current.setEnergy(Enew / (1 + z));
}

void EMInverseComptonScattering::process(Candidate *candidate) const {
	// check if electron / positron
	int id = candidate->current.getId();
	if (abs(id) != 11)
		return;
	
	// scale the particle energy instead of background photons
	double z = candidate->getRedshift();
	double time = candidate->getTime();
	double E = candidate->current.getEnergy() * (1 + z);
	Vector3d position = candidate->current.getPosition();
	
	// interaction rate
	double rate = getRate(E, position, z);
	if (rate < 0)
		return;
	
	rate *= pow_integer<2>(1 + z) * photonField->getRedshiftScaling(z);

	// time scaling
	rate *= photonField->getTimeScaling(time);

	// run this loop at least once to limit the step size
	double step = candidate->getCurrentStep();
	Random &random = Random::instance();
	do {
		double randDistance = -log(random.rand()) / rate;
		
		// check for interaction; if it doesn't ocurr, limit next step
		if (step < randDistance) {
			candidate->limitNextStep(limit / rate);
			return;
		}
		performInteraction(candidate);
		
		// repeat with remaining step
		step -= randDistance;
	} while (step > 0);
}

void EMInverseComptonScattering::setInteractionTag(std::string tag) {
	interactionTag = tag;
}

std::string EMInverseComptonScattering::getInteractionTag() const {
	return interactionTag;
}

} // namespace crpropa
