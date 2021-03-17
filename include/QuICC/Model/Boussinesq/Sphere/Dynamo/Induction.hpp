/**
 * @file Induction.hpp
 * @brief Implementation of the vector induction equation for the Boussinesq thermal convection dynamo sphere
 */

#ifndef QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_INDUCTION_HPP
#define QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_INDUCTION_HPP

// Configuration includes
//

// System includes
//
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Equations/IVectorEquation.hpp"

namespace QuICC {

namespace Equations {

namespace Boussinesq {

namespace Sphere {

namespace Dynamo {

   /**
    * @brief Implementation of the vector induction equation for the Boussinesq thermal convection dynamo in a sphere
    */
   class Induction: public IVectorEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams  Shared equation parameters
          */
         Induction(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend);

         /**
          * @brief Simple empty destructor
          */
         virtual ~Induction();

         /**
          * @brief Initialize nonlinear interaction kernel
          */
         virtual void initNLKernel(const bool force = false);

      protected:
         /**
          * @brief Set variable requirements
          */
         virtual void setRequirements();

         /**
          * @brief Set the equation coupling information
          */
         virtual void setCoupling();

         /**
          * @brief Set the nonlinear integration components
          */
         virtual void setNLComponents();

      private:
   };

}
}
}
}
}

#endif // QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_INDUCTION_HPP
