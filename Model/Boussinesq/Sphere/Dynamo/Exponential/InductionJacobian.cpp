/**
 * @file InductionJacobian.cpp
 * @brief Source of the implementation of the vector induction equation in the
 * Boussinesq thermal convection dynamo in a sphere model
 */

// System includes
//

// Project includes
//
#include "Model/Boussinesq/Sphere/Dynamo/Exponential/InductionJacobian.hpp"
#include "Model/Boussinesq/Sphere/Dynamo/Exponential/InductionJacobianKernel.hpp"
#include "QuICC/PhysicalNames/JacobianMagnetic.hpp"
#include "QuICC/PhysicalNames/JacobianVelocity.hpp"
#include "QuICC/PhysicalNames/Magnetic.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"
#include "QuICC/SolveTiming/Prognostic.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/Transform/Path/CurlCurlNl.hpp"
#include "QuICC/Transform/Path/CurlNl.hpp"

namespace QuICC {

namespace Equations {

namespace Boussinesq {

namespace Sphere {

namespace Dynamo {

namespace Exponential {

InductionJacobian::InductionJacobian(SharedEquationParameters spEqParams,
   SpatialScheme::SharedCISpatialScheme spScheme,
   std::shared_ptr<Model::IModelBackend> spBackend,
   std::shared_ptr<EquationOptions> spOptions) :
    IVectorEquation(spEqParams, spScheme, spBackend)
{
   // Set the variable requirements
   this->setRequirements();
}

void InductionJacobian::setCoupling()
{
   int start;
   if (this->ss().has(SpatialScheme::Feature::SpectralOrdering132))
   {
      start = 1;
   }
   else if (this->ss().has(SpatialScheme::Feature::SpectralOrdering123))
   {
      start = 0;
   }
   else
   {
      throw std::logic_error(
         "Unknown spatial scheme was used to setup equations!");
   }

   auto features = defaultCouplingFeature();
   features.at(CouplingFeature::Nonlinear) = true;

   this->defineCoupling(FieldComponents::Spectral::TOR,
      CouplingInformation::PROGNOSTIC, start, features);

   this->defineCoupling(FieldComponents::Spectral::POL,
      CouplingInformation::PROGNOSTIC, start, features);
}

void InductionJacobian::setNLComponents()
{
   this->addNLComponent(FieldComponents::Spectral::POL,
      Transform::Path::CurlNl::id());

   this->addNLComponent(FieldComponents::Spectral::TOR,
      Transform::Path::CurlCurlNl::id());
}

void InductionJacobian::initNLKernel(const bool force)
{
   // Initialize if empty or forced
   if (force || !this->mspNLKernel)
   {
      // Initialize the physical kernel
      auto spNLKernel = std::make_shared<Physical::Kernel::InductionJacobianKernel>();
      spNLKernel->setJacobianMagnetic(this->name(), this->spUnknown());
      spNLKernel->setMagnetic(PhysicalNames::Magnetic::id(), this->spVector(PhysicalNames::Magnetic::id()));
      spNLKernel->setJacobianVelocity(PhysicalNames::JacobianVelocity::id(),
         this->spVector(PhysicalNames::JacobianVelocity::id()));
      spNLKernel->setVelocity(PhysicalNames::Velocity::id(),
         this->spVector(PhysicalNames::Velocity::id()));
      MHDFloat sgn = 1;
      if(!this->options().nonlinearIsLhs)
      {
         sgn = -1;
      }
      spNLKernel->init(1.0 * sgn);
      this->mspNLKernel = spNLKernel;
   }
}

void InductionJacobian::setRequirements()
{
   // Set velocity as equation unknown
   this->setName(PhysicalNames::JacobianMagnetic::id());

   // Set solver timing
   this->setSolveTiming(SolveTiming::Prognostic::id());

   // Forward transform generates nonlinear RHS
   this->setForwardPathsType(FWD_IS_NONLINEAR);

   // Get reference to spatial scheme
  const auto& ss = this->ss();
 
   // Add Magnetic to requirements
   auto& jmagReq = this->mRequirements.addField(PhysicalNames::JacobianMagnetic::id(),
      FieldRequirement(false, ss.spectral(), ss.physical()));
   jmagReq.enableSpectral();
   jmagReq.enablePhysical();

   // Add Magnetic to requirements
   auto& magReq = this->mRequirements.addField(PhysicalNames::Magnetic::id(),
      FieldRequirement(false, ss.spectral(), ss.physical()));
   magReq.enableSpectral();
   magReq.enablePhysical();

   // Add velocity to requirements: is scalar?
   auto& jvelReq = this->mRequirements.addField(PhysicalNames::JacobianVelocity::id(),
      FieldRequirement(false, ss.spectral(), ss.physical()));
   jvelReq.enableSpectral();
   jvelReq.enablePhysical();

   // Add velocity to requirements: is scalar?
   auto& velReq = this->mRequirements.addField(PhysicalNames::Velocity::id(),
      FieldRequirement(false, ss.spectral(), ss.physical()));
   velReq.enableSpectral();
   velReq.enablePhysical();
}

} // namespace Exponential
} // namespace Dynamo
} // namespace Sphere
} // namespace Boussinesq
} // namespace Equations
} // namespace QuICC
