/**
 * @file PhysicalModel.hpp
 * @brief Implementation of the Boussinesq rotating thermal dynamo in a sphere (Toroidal/Poloidal formulation)
 */

#ifndef QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_IMPLICIT_PHYSICALMODEL_HPP
#define QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_IMPLICIT_PHYSICALMODEL_HPP

// System includes
//
#include <string>

// Project includes
//
#include "Model/Boussinesq/Sphere/Dynamo/IDynamoModel.hpp"
#include "QuICC/SpatialScheme/3D/WLFm.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Sphere {

namespace Dynamo {

namespace Implicit {

   /**
    * @brief Implementation of the Boussinesq rotating thermal dynamo sphere model (Toroidal/Poloidal formulation)
    */
   class PhysicalModel: public IDynamoModel
   {
      public:
         /// Typedef for the spatial scheme used
         typedef SpatialScheme::WLFm SchemeType;

         /**
          * @brief Constructor
          */
         PhysicalModel() = default;

         /**
          * @brief Destructor
          */
         virtual ~PhysicalModel() = default;

         /// Python script/module name
         virtual std::string PYMODULE() override;

         /**
          * @brief Initialize specialized backend
          */
         void init() final;

      protected:

      private:
   };

} // Implicit
} // Dynamo
} // Sphere
} // Boussinesq
} // Model
} // QuICC

#endif // QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_IMPLICIT_PHYSICALMODEL_HPP
