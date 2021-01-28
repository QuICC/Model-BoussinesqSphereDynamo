/** 
 * @file Momentum.cpp
 * @brief Source of the implementation of the vector Navier-Stokes equation in the Boussinesq thermal convection dynamo in a sphere model
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Model/Boussinesq/Sphere/Dynamo/Momentum.hpp"

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Math/Constants.hpp"
#include "QuICC/NonDimensional/Ekman.hpp"
#include "QuICC/NonDimensional/MagPrandtl.hpp"
#include "QuICC/PhysicalNames/Magnetic.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"
#include "QuICC/SpatialScheme/3D/WLFl.hpp"
#include "QuICC/SpatialScheme/3D/WLFm.hpp"
#include "QuICC/Model/Boussinesq/Sphere/Dynamo/MomentumKernel.hpp"

namespace QuICC {

namespace Equations {

namespace Boussinesq {

namespace Sphere {

namespace Dynamo {

   Momentum::Momentum(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme)
      : IVectorEquation(spEqParams,spScheme)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   Momentum::~Momentum()
   {
   }

   void Momentum::setCoupling()
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

      this->defineCoupling(FieldComponents::Spectral::TOR, CouplingInformation::PROGNOSTIC, start, true, false);

      this->defineCoupling(FieldComponents::Spectral::POL, CouplingInformation::PROGNOSTIC, start, true, false);
   }

   void Momentum::setNLComponents()
   {
      this->addNLComponent(FieldComponents::Spectral::TOR, 0);

      this->addNLComponent(FieldComponents::Spectral::POL, 0);
   }

   void Momentum::initNLKernel(const bool force)
   {
      // Initialize the physical kernel
      MHDFloat T = 1.0/this->eqParams().nd(NonDimensional::Ekman::id());
      MHDFloat Pm = this->eqParams().nd(NonDimensional::MagPrandtl::id());
      auto spNLKernel = std::make_shared<Physical::Kernel::MomentumKernel>();
      spNLKernel->setVelocity(this->name(), this->spUnknown());
      spNLKernel->setMagnetic(PhysicalNames::Magnetic::id(), this->spVector(PhysicalNames::Magnetic::id()));
      spNLKernel->init(1.0, T*Pm, T*Pm);
      this->mspNLKernel = spNLKernel;
   }

   void Momentum::setRequirements()
   {
      // Set velocity as equation unknown
      this->setName(PhysicalNames::Velocity::id());

      // Set solver timing
      this->setSolveTiming(SolveTiming::PROGNOSTIC);

      // Forward transform generates nonlinear RHS
      this->setForwardPathsType(FWD_IS_NONLINEAR);

      // Get reference to spatial scheme
      const auto& ss = this->ss();

      // Add velocity to requirements: is scalar?
      auto& velReq = this->mRequirements.addField(PhysicalNames::Velocity::id(), FieldRequirement(false, ss.spectral(), ss.physical()));
      velReq.enableSpectral();
      velReq.enablePhysical();
      velReq.enableCurl();

      // Add magnetic to requirements: is scalar?
      auto& magReq = this->mRequirements.addField(PhysicalNames::Magnetic::id(), FieldRequirement(false, ss.spectral(), ss.physical()));
      magReq.enableSpectral();
      magReq.enablePhysical();
      magReq.enableCurl();
      assert(this->mRequirements.field(PhysicalNames::Velocity::id()).needSpectral());
      assert(this->mRequirements.field(PhysicalNames::Velocity::id()).needPhysical());
      assert(this->mRequirements.field(PhysicalNames::Velocity::id()).needPhysicalCurl());
      assert(this->mRequirements.field(PhysicalNames::Magnetic::id()).needSpectral());
      assert(this->mRequirements.field(PhysicalNames::Magnetic::id()).needPhysical());
      assert(this->mRequirements.field(PhysicalNames::Magnetic::id()).needPhysicalCurl());
   }

}
}
}
}
}
