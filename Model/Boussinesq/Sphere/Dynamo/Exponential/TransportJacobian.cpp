/**
 * @file TransportJacobian.cpp
 * @brief Source of the implementation of the transport equation in the
 * Boussinesq thermal convection dynamo in a sphere
 */

// System includes
//

// Project includes
//
#include "Model/Boussinesq/Sphere/Dynamo/Exponential/TransportJacobian.hpp"
#include "Model/Boussinesq/Sphere/Dynamo/Exponential/TransportJacobianKernel.hpp"
#include "QuICC/PhysicalNames/JacobianTemperature.hpp"
#include "QuICC/PhysicalNames/JacobianVelocity.hpp"
#include "QuICC/PhysicalNames/Temperature.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"
#include "QuICC/SolveTiming/Prognostic.hpp"
#include "QuICC/Transform/Path/ScalarNl.hpp"

namespace QuICC {

namespace Equations {

namespace Boussinesq {

namespace Sphere {

namespace Dynamo {

namespace Exponential {

TransportJacobian::TransportJacobian(SharedEquationParameters spEqParams,
   SpatialScheme::SharedCISpatialScheme spScheme,
   std::shared_ptr<Model::IModelBackend> spBackend,
   std::shared_ptr<EquationOptions> spOptions) :
    IScalarEquation(spEqParams, spScheme, spBackend, spOptions)
{
   // Set the variable requirements
   this->setRequirements();
}

void TransportJacobian::setCoupling()
{
   auto features = defaultCouplingFeature();
   features.at(CouplingFeature::Nonlinear) = true;

   this->defineCoupling(FieldComponents::Spectral::SCALAR,
      CouplingInformation::PROGNOSTIC, 0, features);
}

void TransportJacobian::setNLComponents()
{
   if(this->options().transformHasQi)
   {
      throw std::logic_error("Equation not setup with QI in transform stage");
   }
   else
   {
      this->addNLComponent(FieldComponents::Spectral::SCALAR,
         Transform::Path::ScalarNl::id());
   }
}

void TransportJacobian::initNLKernel(const bool force)
{
   // Initialize if empty or forced
   if (force || !this->mspNLKernel)
   {
      // Initialize the physical kernel
      auto spNLKernel = std::make_shared<Physical::Kernel::TransportJacobianKernel>();
      spNLKernel->setJacobianScalar(this->name(), this->spUnknown());
      spNLKernel->setScalar(PhysicalNames::Temperature::id(), this->spScalar(PhysicalNames::Temperature::id()));
      spNLKernel->setJacobianVector(PhysicalNames::JacobianVelocity::id(),
         this->spVector(PhysicalNames::JacobianVelocity::id()));
      spNLKernel->setVector(PhysicalNames::Velocity::id(),
         this->spVector(PhysicalNames::Velocity::id()));
      MHDFloat sgn = 1;
      if(!this->options().nonlinearIsLhs)
      {
         sgn = -1;
      }
      spNLKernel->init(1.0*sgn);
      this->mspNLKernel = spNLKernel;
   }
}

void TransportJacobian::setRequirements()
{
   // Set temperatur as equation unknown
   this->setName(PhysicalNames::JacobianTemperature::id());

   // Set solver timing
   this->setSolveTiming(SolveTiming::Prognostic::id());

   // Forward transform generates nonlinear RHS
   this->setForwardPathsType(FWD_IS_NONLINEAR);

   // Get reference to spatial scheme
   const auto& ss = this->ss();

   // Add temperature to requirements: is scalar?, need spectral?, need
   // physical?, need diff?
   auto& jtempReq =
      this->mRequirements.addField(PhysicalNames::JacobianTemperature::id(),
         FieldRequirement(true, ss.spectral(), ss.physical()));
   jtempReq.enableSpectral();
   jtempReq.enableGradient();

   // Add temperature to requirements: is scalar?, need spectral?, need
   // physical?, need diff?
   auto& tempReq =
      this->mRequirements.addField(PhysicalNames::Temperature::id(),
         FieldRequirement(true, ss.spectral(), ss.physical()));
   tempReq.enableSpectral();
   tempReq.enableGradient();

   // Add velocity to requirements: is scalar?, need spectral?, need physical?,
   // need diff?(, need curl?)
   auto& jvelReq = this->mRequirements.addField(PhysicalNames::JacobianVelocity::id(),
      FieldRequirement(false, ss.spectral(), ss.physical()));
   jvelReq.enableSpectral();
   jvelReq.enablePhysical();

   // Add velocity to requirements: is scalar?, need spectral?, need physical?,
   // need diff?(, need curl?)
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
