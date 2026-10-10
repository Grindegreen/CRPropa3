#ifndef CRPROPA_EMPAIRPRODUCTION_H
#define CRPROPA_EMPAIRPRODUCTION_H

#include <string>
#include <vector>

#include "crpropa/Module.h"
#include "crpropa/PhotonBackground.h"
#include "crpropa/Geometry.h"
#include "crpropa/InteractionRates.h"


namespace crpropa {
/**
 * \addtogroup EnergyLosses
 * @{
 */
 
/**
 @class EMPairProduction
 @brief Electron-pair production of photons with background photons.

 This module simulates electron-pair production of cosmic ray photons with background photons:
 gamma + gamma_b --> e+ + e- (Breit-Wheeler process).
 The resulting electron positron pair is optionally created (default = false).
 The module limits the propagation step size to a fraction of the mean free path (default = 0.1).
 Thinning is available. A thinning of 0 means that all particles are tracked. 
 For the maximum thinning of 1, only a few representative particles are added to the list of secondaries.
 Note that for thinning>0 the output must contain the column "weights", which should be included in the post-processing.
 The surface is defined to include the nodes of the grid contained within.
 */
class EMPairProduction: public Module {
private:
	ref_ptr<PhotonField> photonField; 	// target photon field
	bool haveElectrons;					// add secondary electrons to simulation
	double limit;						// limit the step to a fraction of the mean free path
	double thinning;					// factor of the thinning (0: no thinning, 1: maximum thinning)
	ref_ptr<Surface> surface; // surface that includes the nodes in the photonField grid to be included
	std::string interactionTag = "EMPP";
	ref_ptr<InteractionRates> interactionRates;
public:
	/** Constructor
	 The object used to load, store and access to the interaction rates of the process is the interactionRates pointer.
	 @param photonField		target photon field
	 @param haveElectrons	if true, add secondary electrons as candidates
	 @param thinning		weighted sampling of secondaries (0: all particles are tracked; 1: maximum thinning)
	 @param limit			step size limit as fraction of mean free path
	 @param surface  		Grid will be confined to `surface->distance(pos)<0`, so for example inside a closed surface
	 */
	EMPairProduction(ref_ptr<PhotonField> photonField, bool haveElectrons = false, double thinning = 0, double limit = 0.1, ref_ptr<Surface> surface = nullptr);

	/** Constructor with given interaction rates
	 No tables are loaded: the rates and cumulative rates come from interactionRates, e.g. computed
	 in memory for a photon field without tabulated rates (setPhotonField() would load tables again).
	 @param photonField			target photon field; gives the name and the redshift and time scaling
	 @param interactionRates	rates and cumulative rates of the process for this field; must not be null
	 @param haveElectrons		if true, add secondary electrons as candidates
	 @param thinning			weighted sampling of secondaries (0: all particles are tracked; 1: maximum thinning)
	 @param limit				step size limit as fraction of mean free path
	 */
	EMPairProduction(ref_ptr<PhotonField> photonField, ref_ptr<InteractionRates> interactionRates, bool haveElectrons = false, double thinning = 0, double limit = 0.1);

	// set the target photon field
	void setPhotonField(ref_ptr<PhotonField> photonField);

	// decide if secondary electrons are added to the simulation
	void setHaveElectrons(bool haveElectrons);

	/** Limit the propagation step to a fraction of the mean free path
	 * @param limit fraction of the mean free path
	 */
	void setLimit(double limit);

	/** Apply thinning with a given thinning factor
	 * @param thinning factor of thinning (0: no thinning, 1: maximum thinning)
	 */
	void setThinning(double thinning);

	/** Apply a surface that confine the position dependent photon field region.
	 * The rates are initialized only for distances smaller then 0, so `surface->distance(pos)<0`
	 * @param surface Grid will be confined to `surface->distance(pos)<0`, so for example inside a closed surface
	 */
	void setSurface(ref_ptr<Surface> surface);
	ref_ptr<Surface> getSurface() const;
		
	/** set a custom interaction tag to trace back this interaction
	 * @param tag string that will be added to the candidate and output
	 */
	void setInteractionTag(std::string tag);
	std::string getInteractionTag() const;

	/** set a custom interaction rate
	 * With this function you can change the type of interaction rate,
	 * if you would for example like to change from a homogeneous to a position
	 * dependent interaction rate.
	 * @param intRates ref_ptr to a InteractionRates class
	 */
	void setInteractionRates(ref_ptr<InteractionRates> intRates);
	ref_ptr<InteractionRates> getInteractionRates() const;
	
	/** Loads the interaction rate
	 * This function loads the interaction rate from a given file/folder.
	 * @param path The name of the file/folder containing the interaction rates
	 */
	void initRate(std::string path);
	
	/** Loads the cumulative interaction rate
	 * This function loads the interaction rate from a given file/folder.
	 * @param path The name of the file/folder containing the interaction rates
	 */
	void initCumulativeRate(std::string path);

	/**
	 * Get the interaction rate for a given energy, position, and redshift.
	 * For now, this function uses a redshift scaling factor for the interaction rate.
	 * Future releases will include a more accurate treatment of the redshift evolution of the photon field.
	 * @param E Energy of the primary particle
	 * @param position Position of the primary particle
	 * @param z Redshift of the primary particle
	 */
	double getRate(double E, const Vector3d &position = Vector3d(0.), double z = 0) const;
		
	void performInteraction(Candidate *candidate) const;
	void process(Candidate *candidate) const;
};

/**
 @class PPSecondariesEnergyDistribution
 @brief Energy distribution of the electron and positron from pair production, used by EMPairProduction.

 Holds the cumulative differential cross section in x = E_positron / E_photon (N = 1000 bins) for bins in the
 squared center-of-mass energy s, log-spaced from s_min = (2 m_e c^2)^2 with the width of 1000 bins up to 1e23 eV^2,
 continued up to at least s_max = 1e30 eV^2 (1638 bins). Beta and the bounds of x are computed from s_min / s without
 cancellation, so the table holds also for s >> s_min.
 */
class PPSecondariesEnergyDistribution {
private:
	std::vector<double> tab_s;					// bin borders in s (Ns + 1 values)
	std::vector< std::vector<double> > data;	// cumulative distribution in x for each bin in s
	size_t N;									// number of bins in x
	/** dsigma/dx without the prefactor, for y = 1 - beta^2 = s_min / s (exact also for s >> s_min) */
	double differentialCrossSection(double x, double y) const;
public:
	PPSecondariesEnergyDistribution();

	/** Differential cross section for x = E_positron / E_photon, compare Lee 96 arXiv:9604098 */
	double dSigmadE_PPx(double x, double beta);

	/** Sample the energy fraction of the softer lepton of the pair, min(x, 1 - x) of x = E_positron / E_photon
	 (the distribution is symmetric in x and 1 - x). Computed without cancellation: 1 - x as
	 x0 + (1 - x0) (1 - exp(-(N - j) dx)), since the upper edge of the last bin in x is 1 - x0.
	 @param s	squared center-of-mass energy
	 @returns	fraction in [x0, 1/2]; NaN for s outside the tabulated range [s_min, s_max]
	 */
	double sampleSoftFraction(double s);

	/** Sample the energy of one lepton of the pair (positron or electron at random).
	 @param E0	energy of the photon
	 @param s	squared center-of-mass energy
	 @returns	energy of the lepton; NaN for s outside the tabulated range [s_min, s_max]
	 */
	double sample(double E0, double s);
};

} // namespace crpropa

#endif // CRPROPA_EMPAIRPRODUCTION_H
