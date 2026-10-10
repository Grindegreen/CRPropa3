#ifndef CRPROPA_EMINVERSECOMPTONSCATTERING_H
#define CRPROPA_EMINVERSECOMPTONSCATTERING_H

#include <string>
#include <vector>

#include "crpropa/Module.h"
#include "crpropa/PhotonBackground.h"
#include "crpropa/InteractionRates.h"
#include "crpropa/Geometry.h"

namespace crpropa {
/**
 * \addtogroup EnergyLosses
 * @{
 */

/**
 @class EMInverseComptonScattering
 @brief Inverse Compton scattering of electrons with background photons.

 This module simulates inverse Compton scattering of electrons with background photons for several photon fields.
 The upscattered photons are optionally created as secondary particles (default = false).
 The module limits the propagation step size to a fraction of the mean free path (default = 0.1).
 Thinning is available. A thinning of 0 means that all particles are tracked. 
 For the maximum thinning of 1, only a few representative particles are added to the list of secondaries.
 Note that for thinning>0 the output must contain the column "weights", which should be included in the post-processing.
 The surface is defined to include the nodes of the grid contained within.
*/
class EMInverseComptonScattering: public Module {
private:
	ref_ptr<PhotonField> photonField;
	bool havePhotons;
	double limit;
	double thinning;
	ref_ptr<Surface> surface;
	std::string interactionTag = "EMIC";
	ref_ptr<InteractionRates> interactionRates;	
public:
	/** Constructor
	 The object used to load, store and access to the interaction rates of the process is the interactionRates pointer.
	 @param photonField		target photon field
	 @param havePhotons		if true, add secondary photons as candidates
	 @param thinning		weighted sampling of secondaries (0: all particles are tracked; 1: maximum thinning)
	 @param limit			step size limit as fraction of mean free path
	 @param surface  		Grid will be confined to `surface->distance(pos)<0`, so for example inside a closed surface
	 */
	EMInverseComptonScattering(ref_ptr<PhotonField> photonField, bool havePhotons = false, double thinning = 0, double limit = 0.1, ref_ptr<Surface> surface = nullptr);

	/** Constructor with given interaction rates
	 No tables are loaded: the rates and cumulative rates come from interactionRates, e.g. computed in memory for a photon
	 field without tabulated rates (setPhotonField() would load tables again).
	 @param photonField			target photon field; gives the name and the redshift and time scaling
	 @param interactionRates	rates and cumulative rates of the process for this field; must not be null
	 @param havePhotons		if true, add secondary photons as candidates
	 @param thinning			weighted sampling of secondaries (0: all particles are tracked; 1: maximum thinning)
	 @param limit				step size limit as fraction of mean free path
	 */
	EMInverseComptonScattering(ref_ptr<PhotonField> photonField, ref_ptr<InteractionRates> interactionRates, bool havePhotons = false, double thinning = 0, double limit = 0.1);
	
	// set the target photon field
	void setPhotonField(ref_ptr<PhotonField> photonField);
	
	// decide if secondary photons are added to the simulation
	void setHavePhotons(bool havePhotons);
	
	/** limit the step to a fraction of the mean free path
	 * @param limit	fraction of the mean free path, should be between 0 and 1
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

	void process(Candidate *candidate) const;
	void performInteraction(Candidate *candidate) const;
};

/**
 @class ICSSecondariesEnergyDistribution
 @brief Energy distribution of the electron after inverse Compton scattering, used by EMInverseComptonScattering.

 Holds the cumulative differential cross section in x = E_electron' / E_electron (Nrer = 1000 bins, log-spaced from
 x_0 = m_e^2 c^4 / s to 1) for bins in the squared center-of-mass energy s, log-spaced from s_min = m_e^2 c^4 with
 the width of 1000 bins up to 2e23 eV^2, continued up to at least s_max = 1e30 eV^2. The kinematics are computed
 from r = s_min / s (x_0 = r, beta = (1 - r) / (1 + r)) without cancellation, so the table holds also for s >> s_min.
 */
class ICSSecondariesEnergyDistribution {
private:
	std::vector< std::vector<double> > data;	// cumulative distribution in x for each bin in s
	std::vector<double> s_values;				// bin borders in s (Ns + 1 values)
	size_t Ns;									// number of bins in s
	size_t Nrer;								// number of bins in x
	double s_min;								// m_e^2 c^4
	/** dsigma/dx without the prefactor, for r = s_min / s (exact also for s >> s_min) */
	double differentialCrossSection(double x, double r) const;
public:
	ICSSecondariesEnergyDistribution();

	/** Differential cross-section, see Lee '96 (arXiv:9604098), eq. 23 for x = Ee'/Ee */
	double dSigmadE(double x, double beta);

	/** Sample the energy of the electron after the scattering.
	 @param Ee	energy of the electron before the scattering
	 @param s	squared center-of-mass energy; s >= s_min (s above s_max uses the distribution of the last bin)
	 @returns	energy of the electron after the scattering, between Ee s_min / s and Ee
	 */
	double sample(double Ee, double s);
};
/** @}*/

} // namespace crpropa

#endif // CRPROPA_EMINVERSECOMPTONSCATTERING_H
