/** 
 * @file Transport.cpp
 * @brief Source of the implementation of the transport equation in the Boussinesq thermal convection dynamo in a sphere
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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
#include "QuICC/Base/Typedefs.hpp"
#include "QuICC/Base/MathConstants.hpp"
#include "QuICC/Enums/NonDimensional.hpp"
#include "QuICC/PhysicalOperators/VelocityAdvection.hpp"

namespace QuICC {

namespace Equations {

namespace Boussinesq {

namespace Sphere {

namespace Dynamo {

   Transport::Transport(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   Transport::~Transport()
   {
   }

   void Transport::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::PROGNOSTIC, 0, true, false);
   }

   void Transport::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      /// 
      /// Computation of the advection:
      ///   \f$ \left(\vec u\cdot\nabla\right)\theta\f$
      ///
      Physical::VelocityAdvection<FieldComponents::Physical::R,FieldComponents::Physical::THETA,FieldComponents::Physical::PHI>::set(rNLComp, this->vector(PhysicalNames::VELOCITY).dom(0).phys(), this->unknown().dom(0).grad(), 1.0);
   }

   void Transport::setRequirements()
   {
      // Set temperatur as equation unknown
      this->setName(PhysicalNames::TEMPERATURE);

      // Set solver timing
      this->setSolveTiming(SolveTiming::PROGNOSTIC);

      // Add temperature to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::TEMPERATURE, FieldRequirement(true, true, false, true));

      // Add velocity to requirements: is scalar?, need spectral?, need physical?, need diff?(, need curl?)
      this->mRequirements.addField(PhysicalNames::VELOCITY, FieldRequirement(false, true, true, false, false));
   }

}
}
}
}
}
