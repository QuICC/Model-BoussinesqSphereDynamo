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
#include "QuICC/SpatialScheme/3D/WLFl.hpp"
#include "QuICC/SpatialScheme/3D/WLFm.hpp"
#include "QuICC/Model/Boussinesq/Sphere/Dynamo/InductionKernel.hpp"

namespace QuICC {

namespace Equations {

namespace Boussinesq {

namespace Sphere {

namespace Dynamo {

   Induction::Induction(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme)
      : IVectorEquation(spEqParams,spScheme)
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
      if(this->ss().id() == SpatialScheme::WLFl::sId)
      {
         start = 1;
      } else if(this->ss().id() == SpatialScheme::WLFm::sId)
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
      this->addNLComponent(FieldComponents::Spectral::POL,0);

      this->addNLComponent(FieldComponents::Spectral::TOR,1);
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
      this->setSolveTiming(SolveTiming::PROGNOSTIC);

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
