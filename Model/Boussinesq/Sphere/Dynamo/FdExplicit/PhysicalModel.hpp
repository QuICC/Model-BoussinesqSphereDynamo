/**
 * @file PhysicalModel.hpp
 * @brief Implementation of the Boussinesq thermal convection dynamo in a sphere
 * model (Toroidal/Poloidal formulation) without coupled solve (standard
 * implementation)
 */

#ifndef QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_FDEXPLICIT_PHYSICALMODEL_HPP
#define QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_FDEXPLICIT_PHYSICALMODEL_HPP

// System includes
//
#include <string>

// Project includes
//
#include "Model/Boussinesq/Sphere/Dynamo/IDynamoModel.hpp"
#include "QuICC/SpatialScheme/3D/FdWLFl.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Sphere {

namespace Dynamo {

namespace FdExplicit {

/**
 * @brief Implementation of the Boussinesq thermal convection dynamo in a sphere
 * model (Toroidal/Poloidal formulation) without coupled solve (standard
 * implementation)
 */
class PhysicalModel : public IDynamoModel
{
public:
   /// Typedef for the spatial scheme used
   typedef SpatialScheme::FdWLFl SchemeType;

   /**
    * @brief Constructor
    */
   PhysicalModel() = default;

   /**
    * @brief Destructor
    */
   virtual ~PhysicalModel() = default;

   /// Python script/module name
   std::string PYMODULE() final;

   /**
    * @brief Initialize specialized backend
    */
   void init() final;

protected:
private:
};

} // namespace FdExplicit
} // namespace Dynamo
} // namespace Sphere
} // namespace Boussinesq
} // namespace Model
} // namespace QuICC

#endif // QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_FDEXPLICIT_PHYSICALMODEL_HPP
