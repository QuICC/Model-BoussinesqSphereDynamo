/** 
 * @file Momentum.cpp
 * @brief Source of the implementation of the vector Navier-Stokes equation in the Boussinesq thermal convection dynamo in a sphere model
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//
#include "QuICC/TypeSelectors/TransformSelector.hpp"

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Model/Boussinesq/Sphere/Dynamo/Momentum.hpp"

// Project includes
//
#include "QuICC/Base/Typedefs.hpp"
#include "QuICC/Base/MathConstants.hpp"
#include "QuICC/Enums/NonDimensional.hpp"
#include "QuICC/PhysicalOperators/Cross.hpp"
#include "QuICC/PhysicalOperators/SphericalCoriolis.hpp"

namespace QuICC {

namespace Equations {

namespace Boussinesq {

namespace Sphere {

namespace Dynamo {

   Momentum::Momentum(SharedEquationParameters spEqParams)
      : IVectorEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   Momentum::~Momentum()
   {
   }

   void Momentum::setCoupling()
   {
      #if defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_WLFL
         int start = 1;
      #elif defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
         int start = 0;
      #endif //defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_WLFL

      this->defineCoupling(FieldComponents::Spectral::TOR, CouplingInformation::PROGNOSTIC, start, true, false);

      this->defineCoupling(FieldComponents::Spectral::POL, CouplingInformation::PROGNOSTIC, start, true, false);

      #if defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_WLFL
         // Create cos(theta) and sin(theta) data for Coriolis term
         int nTh = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
         Array thGrid = Transform::TransformSelector<Dimensions::Transform::TRA2D>::Type::generateGrid(nTh);
         this->mCosTheta = thGrid.array().cos();
         this->mSinTheta = thGrid.array().sin();
      #endif //defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_WLFL
   }

   void Momentum::setNLComponents()
   {
      this->addNLComponent(FieldComponents::Spectral::TOR, 0);

      this->addNLComponent(FieldComponents::Spectral::POL, 0);
   }

   void Momentum::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId) const
   {
      // Get square root of Taylor number
      MHDFloat T = 1.0/this->eqParams().nd(NonDimensional::EKMAN);
      MHDFloat Pm = this->eqParams().nd(NonDimensional::MAGPRANDTL);

      ///
      /// Compute \f$\vec u\wedge\left(\nabla\wedge\vec u\right) + \left(\nabla\wedge\vec B\right)\wedge\vec B\f$
      ///
      switch(compId)
      {
         case(FieldComponents::Physical::R):
            Physical::Cross<FieldComponents::Physical::THETA,FieldComponents::Physical::PHI>::set(rNLComp, this->unknown().dom(0).curl(), this->unknown().dom(0).phys(), 1.0);
            Physical::Cross<FieldComponents::Physical::THETA,FieldComponents::Physical::PHI>::add(rNLComp, this->vector(PhysicalNames::MAGNETIC).dom(0).phys(), this->vector(PhysicalNames::MAGNETIC).dom(0).curl(), T*Pm);
            break;
         case(FieldComponents::Physical::THETA):
            Physical::Cross<FieldComponents::Physical::PHI,FieldComponents::Physical::R>::set(rNLComp, this->unknown().dom(0).curl(), this->unknown().dom(0).phys(), 1.0);
            Physical::Cross<FieldComponents::Physical::PHI,FieldComponents::Physical::R>::add(rNLComp, this->vector(PhysicalNames::MAGNETIC).dom(0).phys(), this->vector(PhysicalNames::MAGNETIC).dom(0).curl(), T*Pm);
            break;
         case(FieldComponents::Physical::PHI):
            Physical::Cross<FieldComponents::Physical::R,FieldComponents::Physical::THETA>::set(rNLComp, this->unknown().dom(0).curl(), this->unknown().dom(0).phys(), 1.0);
            Physical::Cross<FieldComponents::Physical::R,FieldComponents::Physical::THETA>::add(rNLComp, this->vector(PhysicalNames::MAGNETIC).dom(0).phys(), this->vector(PhysicalNames::MAGNETIC).dom(0).curl(), T*Pm);
            break;
         default:
            assert(false);
            break;
      }

      #if defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_WLFL
         ///
         /// Compute Coriolis term
         ///
         Physical::SphericalCoriolis::add(rNLComp, compId, this->unknown().dom(0).spRes(), this->mCosTheta, this->mSinTheta, this->unknown().dom(0).phys(), T*Pm);
      #endif //defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_WLFL
   }

   void Momentum::setRequirements()
   {
      // Set velocity as equation unknown
      this->setName(PhysicalNames::VELOCITY);

      // Set solver timing
      this->setSolveTiming(SolveTiming::PROGNOSTIC);

      // Add velocity to requirements: is scalar?, need spectral?, need physical?, need diff?(, need curl?)
      this->mRequirements.addField(PhysicalNames::VELOCITY, FieldRequirement(false, true, true, false, true));

      // Add magnetic to requirements: is scalar?, need spectral?, need physical?, need diff?(, need curl?)
      this->mRequirements.addField(PhysicalNames::MAGNETIC, FieldRequirement(false, true, true, false, true));
   }

}
}
}
}
}
