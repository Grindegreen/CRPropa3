#include "crpropa/ParticleState.h"
#include "crpropa/Units.h"
#include "crpropa/Common.h"
#include "crpropa/ParticleID.h"
#include "crpropa/ParticleMass.h"

#include "HepPID/ParticleIDMethods.hh"

#include <cstdlib>
#include <sstream>

namespace crpropa {

ParticleState::ParticleState(int id, double E, Vector3d pos, Vector3d dir): id(0), energy(0.), position(0.), direction(0.), pmass(0.), charge(0.)
{
	setId(id);
	setEnergy(E);
	setPosition(pos);
	setDirection(dir);
}

void ParticleState::setPosition(const Vector3d &pos) {
	position = pos;
}

const Vector3d &ParticleState::getPosition() const {
	return position;
}

void ParticleState::setDirection(const Vector3d &dir) {
	direction = dir / dir.getR();
}

const Vector3d &ParticleState::getDirection() const {
	return direction;
}

void ParticleState::setEnergy(double newEnergy) {
	energy = std::max(0., newEnergy); // prevent negative energies
}

double ParticleState::getEnergy() const {
	return energy;
}

double ParticleState::getRigidity() const {
	return fabs(energy / charge);
}

void ParticleState::setId(int newId) {
	id = newId;
	pmass = particleMass(id);
	if (isNucleus(id)) {
		charge = chargeNumber(id) * eplus;
		if (id < 0)
			charge *= -1; // anti-nucleus
	} else {
		charge = HepPID::charge(id) * eplus;
	}
}

int ParticleState::getId() const {
	return id;
}

double ParticleState::getMass() const {
	return pmass;
}

double ParticleState::getCharge() const {
	return charge;
}

double ParticleState::getLorentzFactor() const {
	return energy / (pmass * c_squared) + 1;
}

void ParticleState::setLorentzFactor(double lf) {
	lf = std::max(1., lf); // prevent negative Lorentz factors
	energy = (lf - 1) * pmass * c_squared;
}

Vector3d ParticleState::getVelocity() const {
	return direction * getAbsoluteVelocity();
}

Vector3d ParticleState::getMomentum() const {
	// energy is the kinetic energy $K$, so $p c = \sqrt{K (K + 2 m c^{2})}$.
	// ($K / (\beta c)$, used before, is low by a factor $\gamma / (\gamma + 1)$.)
	// Massless particles: $p = E / c$, returned directly to avoid rounding.
	if (pmass == 0)
		return direction * (energy / c_light);
	return direction * (std::sqrt(energy * (energy + 2 * pmass * c_squared)) / c_light);
}

double ParticleState::getAbsoluteVelocity() const {
	// Massless particles move at $c$. Branching here avoids $0/0$ for $E = m = 0$
	// (e.g. a default-constructed candidate) and $\gamma = \infty$ for photons.
	if (pmass == 0)
		return c_light;
	// $v = p c^{2} / (K + m c^{2}) = c \sqrt{K (K + 2 m c^{2})} / (K + m c^{2})$.
	// Same result as $c \sqrt{1 - \gamma^{-2}}$, but without cancellation for $K \ll m c^{2}$.
	double restEnergy = pmass * c_squared;
	return c_light * std::sqrt(energy * (energy + 2 * restEnergy)) / (energy + restEnergy);
}

std::string ParticleState::getDescription() const {
	std::stringstream ss;
	ss << "Particle " << id << ", ";
	ss << "E = " << energy / EeV << " EeV, ";
	ss << "x = " << position / Mpc << " Mpc, ";
	ss << "p = " << direction;
	return ss.str();
}

} // namespace crpropa
