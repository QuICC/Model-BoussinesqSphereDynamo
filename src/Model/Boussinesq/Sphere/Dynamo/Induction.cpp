/** 
 * @file Induction.cpp
 * @brief Source of the implementation of the vector induction equation in the Boussinesq thermal convection dynamo in a sphere model
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Model/Boussinesq/Sphere/Dynamo/Induction.hpp"

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Math/Constants.hpp"
#include "QuICC/PhysicalNames/Magnetic.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/SolveTiming/Prognostic.hpp"
#include "QuICC/Transform/Path/I2CurlNL.hpp"
#include "QuICC/Transform/Path/I2CurlCurlNL.hpp"
#include "QuICC/Model/Boussinesq/Sphere/Dynamo/InductionKernel.hpp"

namespace QuICC {

namespace Equations {

namespace Boussinesq {

namespace Sphere {

namespace Dynamo {

   Induction::Induction(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend)
      : IVectorEquation(spEqParams, spScheme, spBackend)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   Induction::~Induction()
   {
   }

   void Induction::setCoupling()
   {
      int start;
      if(this->ss().has(SpatialScheme::Feature::SpectralOrdering132))
      {
         start = 1;
      } else if(this->ss().has(SpatialScheme::Feature::SpectralOrdering123))
      {
         start = 0;
      } else
      {
         throw std::logic_error("Unknown spatial scheme was used to setup equations!");
      }

      auto features = defaultCouplingFeature();
      features.at(CouplingFeature::Nonlinear) = true;

      this->defineCoupling(FieldComponents::Spectral::TOR, CouplingInformation::PROGNOSTIC, start, features);

      this->defineCoupling(FieldComponents::Spectral::POL, CouplingInformation::PROGNOSTIC, start, features);
   }

   void Induction::setNLComponents()
   {
      this->addNLComponent(FieldComponents::Spectral::POL, Transform::Path::I2CurlNL::id());

      this->addNLComponent(FieldComponents::Spectral::TOR, Transform::Path::I2CurlCurlNL::id());
   }

   void Induction::initNLKernel(const bool force)
   {
      // Initialize if empty or forced
      if(force || !this->mspNLKernel)
      {
         // Initialize the physical kernel
         auto spNLKernel = std::make_shared<Physical::Kernel::InductionKernel>();
         spNLKernel->setMagnetic(this->name(), this->spUnknown());
         spNLKernel->setVelocity(PhysicalNames::Velocity::id(), this->spVector(PhysicalNames::Velocity::id()));
         spNLKernel->init(1.0);
         this->mspNLKernel = spNLKernel;
      }
   }

   void Induction::setRequirements()
   {
      // Set velocity as equation unknown
      this->setName(PhysicalNames::Magnetic::id());

      // Set solver timing
      this->setSolveTiming(SolveTiming::Prognostic::id());

      // Forward transform generates nonlinear RHS
      this->setForwardPathsType(FWD_IS_NONLINEAR);

      // Get reference to spatial scheme
      const auto& ss = this->ss();

      // Add Magnetic to requirements
      auto& magReq = this->mRequirements.addField(PhysicalNames::Magnetic::id(), FieldRequirement(false, ss.spectral(), ss.physical()));
      magReq.enableSpectral();
      magReq.enablePhysical();

      // Add velocity to requirements: is scalar?
      auto& velReq = this->mRequirements.addField(PhysicalNames::Velocity::id(), FieldRequirement(false, ss.spectral(), ss.physical()));
      velReq.enableSpectral();
      velReq.enablePhysical();
   }

}
}
}
}
}
