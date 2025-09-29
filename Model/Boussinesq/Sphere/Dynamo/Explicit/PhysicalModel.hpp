/**
 * @file PhysicalModel.hpp
 * @brief Implementation of the Boussinesq thermal convection dynamo in a sphere
 * model (Toroidal/Poloidal formulation) without coupled solve (standard
 * implementation)
 */

#ifndef QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_EXPLICIT_PHYSICALMODEL_HPP
#define QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_EXPLICIT_PHYSICALMODEL_HPP

// System includes
//
#include <string>

// Project includes
//
#include "QuICC/SpatialScheme/3D/WLFl.hpp"
#include "QuICC/Model/PyModelBackend.hpp"
#include "Model/Boussinesq/Sphere/Dynamo/Explicit/ModelBackend.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Sphere {

namespace Dynamo {

namespace Explicit {

/**
 * @brief Implementation of the Boussinesq thermal convection dynamo in a sphere
 * model (Toroidal/Poloidal formulation) without coupled solve (standard
 * implementation)
 */
template <typename TBuilder> class PhysicalModel : public TBuilder
{
public:
   /// Typedef for the spatial scheme used
   typedef SpatialScheme::WLFl SchemeType;

   /**
    * @brief Constructor
    */
   PhysicalModel() = default;

   /**
    * @brief Destructor
    */
   virtual ~PhysicalModel() = default;

   /**
    * @brief Initialize specialized backend
    */
   void init() final;

protected:
private:
};

template <typename TBuilder> void PhysicalModel<TBuilder>::init()
{
   TBuilder::init();
#ifdef QUICC_MODEL_BOUSSINESQSPHERERTC_EXPLICIT_BACKEND_CPP

   this->mpBackend = std::make_shared<ModelBackend>();
#else
   std::string pyModule = "boussinesq.sphere.dynamo.explicit.physical_model";
   std::string pyClass = "PhysicalModel";

   this->mpBackend =
      std::make_shared<PyModelBackend>(pyModule, pyClass);
#endif
}

} // namespace Explicit
} // namespace Dynamo
} // namespace Sphere
} // namespace Boussinesq
} // namespace Model
} // namespace QuICC

#endif // QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_EXPLICIT_PHYSICALMODEL_HPP
