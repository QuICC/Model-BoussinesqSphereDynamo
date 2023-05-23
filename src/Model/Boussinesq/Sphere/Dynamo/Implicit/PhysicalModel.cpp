/**
 * @file PhysicalModel.cpp
 * @brief Source of the Boussinesq thermal convection dynamo in a sphere (Toroidal/Poloidal formulation)
 */

// System includes
//

// Project includes
//
#include "QuICC/Model/Boussinesq/Sphere/Dynamo/Implicit/PhysicalModel.hpp"
#include "QuICC/Model/Boussinesq/Sphere/Dynamo/Implicit/ModelBackend.hpp"
#include "QuICC/Model/PyModelBackend.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Sphere {

namespace Dynamo {

namespace Implicit {

   std::string PhysicalModel::PYMODULE()
   {
      return "boussinesq.sphere.dynamo.implicit.physical_model";
   }

   void PhysicalModel::init()
   {
#ifdef QUICC_MODEL_BOUSSINESQSPHEREDYNAMO_IMPLICIT_BACKEND_CPP
      IPhysicalModel<Simulation,StateGenerator,VisualizationGenerator>::init();

      this->mpBackend = std::make_shared<ModelBackend>();
#else
      IPhysicalPyModel<Simulation,StateGenerator,VisualizationGenerator>::init();

      this->mpBackend = std::make_shared<PyModelBackend>(this->PYMODULE(), this->PYCLASS());
#endif
   }

} // Implicit
} // Dynamo
} // Sphere
} // Boussinesq
} // Model
} // QuICC
