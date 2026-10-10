#include "crpropa/module/EMPairProduction.h"
#include "crpropa/Units.h"
#include "crpropa/Random.h"
#include "crpropa/Common.h"
#include "kiss/logger.h"

#include <algorithm>
#include <vector>
#include <cmath>
#include <stdexcept>

namespace crpropa {

static const double mec2 = mass_electron * c_squared;

EMPairProduction::EMPairProduction(ref_ptr<PhotonField> photonField, bool haveElectrons, double thinning, double limit, ref_ptr<Surface> surface) {
	setSurface(surface);
	setPhotonField(photonField);
	setThinning(thinning);
	setLimit(limit);
	setHaveElectrons(haveElectrons);
}

EMPairProduction::EMPairProduction(ref_ptr<PhotonField> photonField, ref_ptr<InteractionRates> interactionRates, bool haveElectrons, double thinning, double limit) {
	if (!photonField)
		throw std::invalid_argument("EMPairProduction: the photon field must not be null");
	if (!interactionRates)
		throw std::invalid_argument("EMPairProduction: the interaction rates must not be null");
	// as setPhotonField(), but with the given rates instead of tables loaded by the field name
	this->photonField = photonField;
	setDescription("EMPairProduction: " + photonField->getFieldName());
	setInteractionRates(interactionRates);
	setThinning(thinning);
	setLimit(limit);
	setHaveElectrons(haveElectrons);
}

void EMPairProduction::setPhotonField(ref_ptr<PhotonField> photonField) {
		
	this->photonField = photonField;
	std::string fname = photonField->getFieldName();
	setDescription("EMPairProduction: " + fname);
		
	if (!this->photonField->hasPositionDependence()){	
		this->interactionRates = new InteractionRatesHomogeneous(
			getDataPath("EMPairProduction/rate_" + fname + ".txt"),
			getDataPath("EMPairProduction/cdf_" + fname + ".txt")
		);
	} else {
		this->interactionRates = new InteractionRatesPositionDependent(
			getDataPath("EMPairProduction/"+fname+"/Rate/"),
			getDataPath("EMPairProduction/"+fname+"/CumulativeRate/"),
			this->surface
		);
	}
}

void EMPairProduction::setHaveElectrons(bool haveElectrons) {
	this->haveElectrons = haveElectrons;
}

void EMPairProduction::setLimit(double limit) {
	this->limit = limit;
}

void EMPairProduction::setThinning(double thinning) {
	this->thinning = thinning;
}

void EMPairProduction::setSurface(ref_ptr<Surface> surface) {
	this->surface = surface;
}

ref_ptr<Surface> EMPairProduction::getSurface() const {
	return this->surface;
}

void EMPairProduction::setInteractionRates(ref_ptr<InteractionRates> intRates) {
	this->interactionRates = intRates;
}

ref_ptr<InteractionRates> EMPairProduction::getInteractionRates() const {
	return this->interactionRates;
}

void EMPairProduction::initRate(std::string path) {
	this->interactionRates->initRate(path);
}

void EMPairProduction::initCumulativeRate(std::string path) {
	this->interactionRates->initCumulativeRate(path);
}

// Hold an data array to interpolate the energy distribution on
PPSecondariesEnergyDistribution::PPSecondariesEnergyDistribution() {
	N = 1000;
	const double s_min = 4 * mec2 * mec2;
	// bins of the original table (1000 bins from s_min to 1e23 eV^2), continued up to at least 1e30 eV^2, so that the
	// bins below 1e23 eV^2 are unchanged
	const double dls = log(1e23 * eV * eV / s_min) / 1000;
	const size_t Ns = static_cast<size_t>(std::ceil(log(1e30 * eV * eV / s_min) / dls));
	data = std::vector< std::vector<double> >(Ns, std::vector<double>(N));
	tab_s = std::vector<double>(Ns + 1);
	
	for (size_t i = 0; i < Ns + 1; ++i)
		tab_s[i] = s_min * exp(i*dls); // tabulate s bin borders
	
	for (size_t i = 0; i < Ns; i++) {
		// with r = s_min / s exactly: beta^2 = 1 - r, 1 - beta = r / (1 + beta) (no cancellation for s >> s_min)
		double r = exp(-(i + 0.5) * dls);
		double beta = sqrt(1 - r);
		double x0 = r / (1 + beta) / 2; // (1 - beta) / 2
		double dx = (2 * log1p(beta) - log(r)) / N; // ln((1 + beta) / (1 - beta)) / N
		
		// cumulative midpoint integration
		std::vector<double> data_i(N);
		data_i[0] = differentialCrossSection(x0, r) * expm1(dx);
		for (size_t j = 1; j < N; j++) {
			double x = x0 * exp(j*dx + 0.5*dx);
			double binWidth = exp((j+1)*dx)-exp(j*dx);
			data_i[j] = differentialCrossSection(x, r) * binWidth + data_i[j-1];
		}
		data[i] = data_i;
	}
}

// differential cross section for pair production for x = Epositron/Egamma, compare Lee 96 arXiv:9604098
double PPSecondariesEnergyDistribution::dSigmadE_PPx(double x, double beta) {
	return differentialCrossSection(x, 1 - beta * beta);
}

double PPSecondariesEnergyDistribution::differentialCrossSection(double x, double y) const {
	// y = 1 - beta^2 = s_min / s
	double A = (x / (1. - x) + (1. - x) / x );
	double B =  (1. / x + 1. / (1. - x) );
	return A + y * B - y * y / 4 * B * B;
}

// sample the energy fraction of the softer lepton from cdf(s)
double PPSecondariesEnergyDistribution::sampleSoftFraction(double s) {
	// get distribution for given s: bin idx holds tab_s[idx] <= s < tab_s[idx + 1], and s at the upper end the last
	// bin (data has one bin less than tab_s has borders)
	if (!(s >= tab_s.front() && s <= tab_s.back()))
		return NAN;
	size_t idx = std::upper_bound(tab_s.begin(), tab_s.end(), s) - tab_s.begin() - 1;
	idx = std::min(idx, data.size() - 1);
	
	std::vector<double> s0 = data[idx];
	
	// draw random bin; x = x0 exp(j dx) is its upper edge (j = N: x = 1 - x0)
	Random &random = Random::instance();
	size_t j = random.randBin(s0) + 1;
	
	// with r = s_min / s exactly, as in the constructor
	double r = 4. * mec2 * mec2 / s;
	double beta = sqrt(1. - r);
	double x0 = r / (1. + beta) / 2.;
	double dx = (2. * log1p(beta) - log(r)) / N;
	double x = x0 * exp(j * dx);
	// 1 - x without cancellation: x0 exp(N dx) = 1 - x0, so 1 - x = x0 + (1 - x0) (1 - exp(-(N - j) dx))
	double oneMinusX = x0 - (1. - x0) * expm1(-double(N - j) * dx);
	return std::min(x, oneMinusX);
}

// sample the energy of one lepton (positron or electron at random)
double PPSecondariesEnergyDistribution::sample(double E0, double s) {
	double soft = sampleSoftFraction(s);
	if (!std::isfinite(soft))
		return NAN;
	// the distribution is symmetric in x and 1 - x: the softer or the harder lepton with equal probability
	return (Random::instance().rand() < 0.5) ? E0 * soft : E0 * (1. - soft);
}

double EMPairProduction::getRate(double E, const Vector3d &position, double z) const {
	return this->interactionRates->getProcessRate(E, position) * pow_integer<2>(1 + z) * photonField->getRedshiftScaling(z);
}

void EMPairProduction::performInteraction(Candidate *candidate) const {
	
	// scale particle energy instead of background photon energy
	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy() * (1 + z);
	Vector3d position = candidate->current.getPosition();
	
	// cosmic ray photon is lost if interaction
	candidate->setActive(false);
	
	// check if secondary electron pair needs to be produced
	if (not haveElectrons)
		return;
	
	std::vector<double> tabE;
	std::vector<double> tabs;
	std::vector<std::vector<double>> tabCDF;
	
	this->interactionRates->loadPerformInteractionTabs(position, tabE, tabs, tabCDF);
	
	// check if in tabulated energy range
	if (E < tabE.front() or (E > tabE.back())) {		
		KISS_LOG_WARNING 
			<< "EMPairProduction: Energy " 
			<< E / eV << " eV is not in tabulated range";
		return;
	}

	// sample the value of s
	Random &random = Random::instance();
	size_t i = closestIndex(E, tabE);  // find closest tabulation point
	size_t j = random.randBin(tabCDF[i]);
	// bin j of the cumulative rate holds tabs[j-1] < s <= tabs[j]; bin 0 holds the rate from the threshold
	// (2 me c^2)^2 to the first tabulated value tabs[0], which lies above the threshold in the tables
	const double sMin = 4 * mec2 * mec2;
	double lo = (j == 0) ? sMin : std::max(sMin, tabs[j-1]);
	double hi = std::max(lo, tabs[j]);
	double s = lo + random.rand() * (hi - lo);
	
	// sample electron / positron energy: the softer lepton (fraction computed without cancellation) is the electron
	// or the positron with equal probability, the harder one gets the rest (no cancellation, energy conserved)
	static PPSecondariesEnergyDistribution interpolation;
	double soft = interpolation.sampleSoftFraction(s);
	double Ee, Ep;
	if (random.rand() < 0.5) {
		Ee = E * soft;
		Ep = E - Ee;
	} else {
		Ep = E * soft;
		Ee = E - Ep;
	}
	double f = Ep / E;
	
	// s outside the tabulated range of the energy distribution gives NaN
	if (not std::isfinite(Ee) || not std::isfinite(Ep)) {
		KISS_LOG_WARNING
			<< "EMPairProduction: Sampled energies are not finite for primary energy "
			<< E / eV << " eV and s = " << s / (eV * eV) << " eV^2 (maximum tabulated s = "
			<< tabs.back() / (eV * eV) << " eV^2).";
		return;
	}

	// sample random position along current step
	Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
	// apply sampling
	if (random.rand() < pow(f, thinning)) {
		double w = 1. / pow(f, thinning);
		candidate->addSecondary(11, Ep / (1 + z), pos, w, interactionTag);
	}
	if (random.rand() < pow(1 - f, thinning)) {
		double w = 1. / pow(1 - f, thinning);
		candidate->addSecondary(-11, Ee / (1 + z), pos, w, interactionTag);
	}
}

void EMPairProduction::process(Candidate *candidate) const {
	
	// check if photon
	if (candidate->current.getId() != 22)
		return;
	
	// scale particle energy instead of background photon energy
	double z = candidate->getRedshift();
	double time = candidate->getTime();
	double E = candidate->current.getEnergy() * (1 + z);
	Vector3d position = candidate->current.getPosition();
	
	double rate = getRate(E, position, z);
	// geometric rate
	rate *= photonField->getTimeScaling(time);

	if (rate < 0)
		return;
	
	// run this loop at least once to limit the step size
	double step = candidate->getCurrentStep();
	Random &random = Random::instance();
	do {
		double randDistance = -log(random.rand()) / rate;
		// check for interaction; if it doesn't occur, limit next step
		if (step < randDistance) {
			candidate->limitNextStep(limit / rate);
		} else {
			performInteraction(candidate);
			return;
		}
		step -= randDistance;
	} while (step > 0.);
}

void EMPairProduction::setInteractionTag(std::string tag) {
	interactionTag = tag;
}

std::string EMPairProduction::getInteractionTag() const {
	return interactionTag;
}

} // namespace crpropa
