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
#include "QuICC/NonDimensional/Rayleigh.hpp"
#include "QuICC/NonDimensional/Prandtl.hpp"
#include "QuICC/NonDimensional/Ekman.hpp"
#include "QuICC/NonDimensional/MagneticPrandtl.hpp"
#include "QuICC/PhysicalNames/Temperature.hpp"
#include "QuICC/PhysicalNames/Magnetic.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"
#include "QuICC/Bc/Name/StressFree.hpp"
#include "QuICC/SolveTiming/Prognostic.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/SpectralKernels/Sphere/ConserveAngularMomentum.hpp"
#include "QuICC/Transform/Path/I2CurlNl.hpp"
#include "QuICC/Transform/Path/NegI2CurlCurlNl.hpp"
#include "QuICC/Transform/Path/NegI4CurlCurlNl.hpp"
#include "QuICC/Model/Boussinesq/Sphere/Dynamo/MomentumKernel.hpp"

namespace QuICC {

namespace Equations {

namespace Boussinesq {

namespace Sphere {

namespace Dynamo {

   Momentum::Momentum(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend)
      : IVectorEquation(spEqParams, spScheme, spBackend)
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

   void Momentum::setNLComponents()
   {
      this->addNLComponent(FieldComponents::Spectral::TOR, Transform::Path::I2CurlNl::id());

      if(this->couplingInfo(FieldComponents::Spectral::POL).isSplitEquation())
      {
         this->addNLComponent(FieldComponents::Spectral::POL, Transform::Path::NegI2CurlCurlNl::id());
      }
      else
      {
         this->addNLComponent(FieldComponents::Spectral::POL, Transform::Path::NegI4CurlCurlNl::id());
      }
   }

   void Momentum::initNLKernel(const bool force)
   {
      // Initialize if empty or forced
      if(force || !this->mspNLKernel)
      {
         // Initialize the physical kernel
         MHDFloat T = 1.0/this->eqParams().nd(NonDimensional::Ekman::id());
         MHDFloat Ra = this->eqParams().nd(NonDimensional::Rayleigh::id());
         MHDFloat Pr = this->eqParams().nd(NonDimensional::Prandtl::id());
         MHDFloat Pm = this->eqParams().nd(NonDimensional::MagneticPrandtl::id());
         auto spNLKernel = std::make_shared<Physical::Kernel::MomentumKernel>();
         spNLKernel->setVelocity(this->name(), this->spUnknown());
         spNLKernel->setTemperature(PhysicalNames::Temperature::id(), this->spScalar(PhysicalNames::Temperature::id()));
         spNLKernel->setMagnetic(PhysicalNames::Magnetic::id(), this->spVector(PhysicalNames::Magnetic::id()));
         spNLKernel->init(1.0, T*Pm, Pm*Pm*Ra*T/Pr, T*Pm);
         this->mspNLKernel = spNLKernel;
      }
   }

   void Momentum::initConstraintKernel(const std::shared_ptr<std::vector<Array> >)
   {
      if(this->bcIds().bcId(this->name()) == Bc::Name::StressFree::id())
      {
         // Initialize the physical kernel
         auto spConstraint = std::make_shared<Spectral::Kernel::Sphere::ConserveAngularMomentum>(this->ss().has(SpatialScheme::Feature::ComplexSpectrum));
         spConstraint->setField(this->name(), this->spUnknown());
         spConstraint->setResolution(this->spRes());
         spConstraint->init(this->ss().has(SpatialScheme::Feature::SpectralOrdering123));
         this->setConstraintKernel(FieldComponents::Spectral::TOR, spConstraint);
      }
   }

   void Momentum::setRequirements()
   {
      // Set velocity as equation unknown
      this->setName(PhysicalNames::Velocity::id());

      // Set solver timing
      this->setSolveTiming(SolveTiming::Prognostic::id());

      // Forward transform generates nonlinear RHS
      this->setForwardPathsType(FWD_IS_NONLINEAR);

      // Get reference to spatial scheme
      const auto& ss = this->ss();

      // Add velocity to requirements: is scalar?
      auto& velReq = this->mRequirements.addField(PhysicalNames::Velocity::id(), FieldRequirement(false, ss.spectral(), ss.physical()));
      velReq.enableSpectral();
      velReq.enablePhysical();
      velReq.enableCurl();

      // Add temperature to requirements: is scalar?
      auto& tempReq = this->mRequirements.addField(PhysicalNames::Temperature::id(), FieldRequirement(true, ss.spectral(), ss.physical()));
      tempReq.enableSpectral();
      tempReq.enablePhysical();

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
