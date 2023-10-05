/**
 * @file PhysicalModel.cpp
 * @brief Source of the Boussinesq thermal convection dynamo in a sphere (Toroidal/Poloidal formulation) without coupled solve (standard implementation)
 */

// System includes
//

// Project include
//
#include "Model/Boussinesq//Sphere/Dynamo/Explicit/PhysicalModel.hpp"
#include "Model/Boussinesq//Sphere/Dynamo/Explicit/ModelBackend.hpp"
#include "QuICC/Model/PyModelBackend.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Sphere {

namespace Dynamo {

namespace Explicit {

   std::string PhysicalModel::PYMODULE()
   {
      return "boussinesq.sphere.dynamo.explicit.physical_model";
   }

   void PhysicalModel::init()
   {
#ifdef QUICC_MODEL_BOUSSINESQSPHEREDYNAMO_EXPLICIT_BACKEND_CPP
      IPhysicalModel<Simulation,StateGenerator,VisualizationGenerator>::init();

      this->mpBackend = std::make_shared<ModelBackend>();
#else
      IPhysicalPyModel<Simulation,StateGenerator,VisualizationGenerator>::init();

      this->mpBackend = std::make_shared<PyModelBackend>(this->PYMODULE(), this->PYCLASS());
#endif
   }

} // Explicit
} // Dynamo
} // Sphere
} // Boussinesq
} // Model
} // QuICC
