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

// differential cross section for pair production for x = Epositron/Egamma, compare Lee 96 arXiv:9604098
double PPSecondariesEnergyDistribution::dSigmadE_PPx(double x, double beta) {
	double A = (x / (1. - x) + (1. - x) / x );
	double B =  (1. / x + 1. / (1. - x) );
	double y = (1 - beta * beta);
	return A + y * B - y * y / 4 * B * B;
}

// Hold an data array to interpolate the energy distribution on
PPSecondariesEnergyDistribution::PPSecondariesEnergyDistribution() {
	N = 1000;
	size_t Ns = 1000;
	double s_min = 4 * mec2 * mec2;
	double s_max = 1e23 * eV * eV;
	double dls = log(s_max / s_min) / Ns;
	data = std::vector< std::vector<double> >(Ns, std::vector<double>(N));
	tab_s = std::vector<double>(Ns + 1);
	
	for (size_t i = 0; i < Ns + 1; ++i)
		tab_s[i] = s_min * exp(i*dls); // tabulate s bin borders
	tab_s[Ns] = s_max; // exact upper border (exp may round it)
	
	for (size_t i = 0; i < Ns; i++) {
		double s = s_min * exp(i*dls + 0.5*dls);
		double beta = sqrt(1 - s_min/s);
		double x0 = (1 - beta) / 2;
		double dx = log((1 + beta) / (1 - beta)) / N;
		
		// cumulative midpoint integration
		std::vector<double> data_i(1000);
		data_i[0] = dSigmadE_PPx(x0, beta) * expm1(dx);
		for (size_t j = 1; j < N; j++) {
			double x = x0 * exp(j*dx + 0.5*dx);
			double binWidth = exp((j+1)*dx)-exp(j*dx);
			data_i[j] = dSigmadE_PPx(x, beta) * binWidth + data_i[j-1];
		}
		data[i] = data_i;
	}
}

// sample positron energy from cdf(E, s_kin)
double PPSecondariesEnergyDistribution::sample(double E0, double s) {
	// get distribution for given s: bin idx holds tab_s[idx] <= s < tab_s[idx + 1], and s = s_max the last bin
	// (data has one bin less than tab_s has borders)
	if (!(s >= tab_s.front() && s <= tab_s.back()))
		return NAN;
	size_t idx = std::upper_bound(tab_s.begin(), tab_s.end(), s) - tab_s.begin() - 1;
	idx = std::min(idx, data.size() - 1);

	std::vector<double> s0 = data[idx];
	
	// draw random bin
	Random &random = Random::instance();
	size_t j = random.randBin(s0) + 1;
	
	double s_min = 4. * mec2 * mec2;
	double beta = sqrtl(1. - s_min / s);
	double x0 = (1. - beta) / 2.;
	double dx = log((1 + beta) / (1 - beta)) / N;
	double binWidth = x0 * (exp(j*dx) - exp((j-1)*dx));
	if (random.rand() < 0.5)
		return E0 * (x0 * exp((j-1) * dx) + binWidth);
	else
		return E0 * (1 - (x0 * exp((j-1) * dx) + binWidth));
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
	
	// sample electron / positron energy
	static PPSecondariesEnergyDistribution interpolation;
	double Ee = interpolation.sample(E, s);
	double Ep = E - Ee;
	double f = Ep / E;
	
	// for some backgrounds Ee=nan due to precision limitations.
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
