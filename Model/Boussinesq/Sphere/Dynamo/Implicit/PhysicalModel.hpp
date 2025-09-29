/**
 * @file PhysicalModel.hpp
 * @brief Implementation of the Boussinesq rotating thermal dynamo in a sphere
 * (Toroidal/Poloidal formulation)
 */

#ifndef QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_IMPLICIT_PHYSICALMODEL_HPP
#define QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_IMPLICIT_PHYSICALMODEL_HPP

// System includes
//
#include <string>

// Project includes
//
#include "QuICC/SpatialScheme/3D/WLFm.hpp"
#include "QuICC/Model/PyModelBackend.hpp"
#include "Model/Boussinesq/Sphere/Dynamo/Implicit/ModelBackend.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Sphere {

namespace Dynamo {

namespace Implicit {

/**
 * @brief Implementation of the Boussinesq rotating thermal dynamo sphere model
 * (Toroidal/Poloidal formulation)
 */
template <typename TBuilder> class PhysicalModel : public TBuilder
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

#ifdef QUICC_MODEL_BOUSSINESQSPHERERTC_IMPLICIT_BACKEND_CPP
   this->mpBackend = std::make_shared<ModelBackend>();
#else
   std::string pyModule = "boussinesq.sphere.dynamo.implicit.physical_model";
   std::string pyClass = "PhysicalModel";

   this->mpBackend =
      std::make_shared<PyModelBackend>(pyModule, pyClass);
#endif
}

} // namespace Implicit
} // namespace Dynamo
} // namespace Sphere
} // namespace Boussinesq
} // namespace Model
} // namespace QuICC

#endif // QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_IMPLICIT_PHYSICALMODEL_HPP
