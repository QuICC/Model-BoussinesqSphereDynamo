/**
 * @file Momentum.hpp
 * @brief Implementation of the vector Navier-Stokes equation for the Boussinesq thermal convection dynamo sphere
 */

#ifndef QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_MOMENTUM_HPP
#define QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_MOMENTUM_HPP

// Configuration includes
//

// System includes
//
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Equations/IVectorEquation.hpp"

namespace QuICC {

namespace Equations {

namespace Boussinesq {

namespace Sphere {

namespace Dynamo {

   /**
    * @brief Implementation of the vector Navier-Stokes equation for the Boussinesq thermal convection dynamo in a sphere
    */
   class Momentum: public IVectorEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams  Shared equation parameters
          */
         Momentum(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend);

         /**
          * @brief Simple empty destructor
          */
         virtual ~Momentum();

         /**
          * @brief Initialize constraint kernel
          *
          * @param spMesh  Physical space mesh
          */
         void initConstraintKernel(const std::shared_ptr<std::vector<Array> > spMesh) final;

         /**
          * @brief Initialize nonlinear interaction kernel
          */
         void initNLKernel(const bool force = false) final;

      protected:
         /**
          * @brief Set variable requirements
          */
         void setRequirements() final;

         /**
          * @brief Set the equation coupling information
          */
         void setCoupling() final;

         /**
          * @brief Set the nonlinear integration components
          */
         void setNLComponents() final;

      private:
   };

}
}
}
}
}

#endif // QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_MOMENTUM_HPP
