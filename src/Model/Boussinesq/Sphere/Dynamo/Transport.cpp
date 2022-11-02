/** 
 * @file Transport.cpp
 * @brief Source of the implementation of the transport equation in the Boussinesq thermal convection dynamo in a sphere
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Model/Boussinesq/Sphere/Dynamo/Transport.hpp"

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Math/Constants.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"
#include "QuICC/PhysicalNames/Temperature.hpp"
#include "QuICC/SolveTiming/Prognostic.hpp"
#include "QuICC/Transform/Path/I2ScalarNl.hpp"
#include "QuICC/Model/Boussinesq/Sphere/Dynamo/TransportKernel.hpp"

namespace QuICC {

namespace Equations {

namespace Boussinesq {

namespace Sphere {

namespace Dynamo {

   Transport::Transport(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend)
      : IScalarEquation(spEqParams, spScheme, spBackend)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   Transport::~Transport()
   {
   }

   void Transport::setCoupling()
   {
      auto features = defaultCouplingFeature();
      features.at(CouplingFeature::Nonlinear) = true;

      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::PROGNOSTIC, 0, features);
   }

   void Transport::setNLComponents()
   {
      this->addNLComponent(FieldComponents::Spectral::SCALAR, Transform::Path::I2ScalarNl::id());
   }

   void Transport::initNLKernel(const bool force)
   {
      // Initialize if empty or forced
      if(force || !this->mspNLKernel)
      {
         // Initialize the physical kernel
         auto spNLKernel = std::make_shared<Physical::Kernel::TransportKernel>();
         spNLKernel->setScalar(this->name(), this->spUnknown());
         spNLKernel->setVector(PhysicalNames::Velocity::id(), this->spVector(PhysicalNames::Velocity::id()));
         spNLKernel->init(1.0);
         this->mspNLKernel = spNLKernel;
      }
   }

   void Transport::setRequirements()
   {
      // Set temperatur as equation unknown
      this->setName(PhysicalNames::Temperature::id());

      // Set solver timing
      this->setSolveTiming(SolveTiming::Prognostic::id());

      // Forward transform generates nonlinear RHS
      this->setForwardPathsType(FWD_IS_NONLINEAR);

      // Get reference to spatial scheme
      const auto& ss = this->ss();

      // Add temperature to requirements: is scalar?, need spectral?, need physical?, need diff?
      auto& tempReq = this->mRequirements.addField(PhysicalNames::Temperature::id(), FieldRequirement(true, ss.spectral(), ss.physical()));
      tempReq.enableSpectral();
      tempReq.enableGradient();

      // Add velocity to requirements: is scalar?, need spectral?, need physical?, need diff?(, need curl?)
      auto& velReq = this->mRequirements.addField(PhysicalNames::Velocity::id(), FieldRequirement(false, ss.spectral(), ss.physical()));
      velReq.enableSpectral();
      velReq.enablePhysical();
   }

}
}
}
}
}
