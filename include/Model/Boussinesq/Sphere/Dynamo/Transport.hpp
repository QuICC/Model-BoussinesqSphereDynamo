/**
 * @file Transport.hpp
 * @brief Implementation of the transport equation for the Boussinesq thermal convection dynamo in a sphere
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_TRANSPORT_HPP
#define QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_TRANSPORT_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "TypeSelectors/ScalarSelector.hpp"
#include "Equations/IScalarEquation.hpp"

namespace QuICC {

namespace Equations {

namespace Boussinesq {

namespace Sphere {

namespace Dynamo {

   /**
    * @brief Implementation of the transport equation for the Boussinesq thermal convection dynamo in a sphere 
    */
   class Transport: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams  Shared equation parameters
          */
         Transport(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~Transport();

         /**
          * @brief Compute the nonlinear interaction term
          *
          * @param rNLComp Nonlinear term component
          * @param id      ID of the component (allows for a more general implementation)
          */
         virtual void computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const;

      protected:
         /**
          * @brief Set variable requirements
          */
         virtual void setRequirements();

         /**
          * @brief Set the equation coupling information
          */
         virtual void setCoupling();

      private:
   };

}
}
}
}
}

#endif // QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_TRANSPORT_HPP
