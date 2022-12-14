/**
 * @file PhysicalModel.cpp
 * @brief Source of the Boussinesq thermal convection dynamo in a sphere (Toroidal/Poloidal formulation) without coupled solve (standard implementation)
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Model/Boussinesq/Sphere/Dynamo/Explicit/PhysicalModel.hpp"
#include "QuICC/Model/Boussinesq/Sphere/Dynamo/Explicit/ModelBackend.hpp"
#include "QuICC/Model/PyModelBackend.hpp"

// Project includes
//

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
#if 0
      IPhysicalPyModel<Simulation,StateGenerator,VisualizationGenerator>::init();

      this->mpBackend = std::make_shared<PyModelBackend>(this->PYMODULE(), this->PYCLASS());
#else
      IPhysicalModel<Simulation,StateGenerator,VisualizationGenerator>::init();

      this->mpBackend = std::make_shared<ModelBackend>();
#endif
   }

}
}
}
}
}
}
