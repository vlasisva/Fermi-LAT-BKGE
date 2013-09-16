/**
 * @file Psf.cxx
 * @brief Wrapper class for IPsf for use from ROOT
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/GRBAnalysis-scons/irfs/rootIrfLoader/src/rootIrfLoader/Psf.cxx,v 1.1.1.1 2005/08/13 15:24:08 jchiang Exp $
 */

#include "irfInterface/Irfs.h"

#include "Psf.h"

namespace rootIrfLoader {

Psf::Psf(const std::string & irfs) : IrfBase(irfs), m_psf(0) {
   setIrfs(m_currentIrfs);
}

Psf::Psf(const Psf & rhs) : IrfBase(rhs) {
   m_psf = rhs.m_psf->clone();
}

Psf & Psf::operator=(const Psf & rhs) {
   if (*this != rhs) {
      IrfBase::operator=(rhs);
      delete m_psf;
      m_psf = rhs.m_psf->clone();
   }
   return *this;
}

Psf::~Psf() {
   delete m_psf;
}

bool Psf::operator==(const Psf & rhs) const {
   return this->currentIrfs() == rhs.currentIrfs();
}

void Psf::getRef(irfInterface::Irfs * irfs) {
   delete m_psf;
   m_psf = irfs->psf()->clone();
}

double Psf::operator()(double separation, double energy, double theta,
                       double phi) const {
   return m_psf->value(separation, energy, theta, phi);
}

double Psf::angularIntegral(double energy, double theta, double phi,
                            double radius) const {
   return m_psf->angularIntegral(energy, theta, phi, radius);
}

} // namespace rootIrfLoader
