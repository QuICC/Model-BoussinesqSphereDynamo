/**
 * @file IExponentialModel.cpp
 * @brief Source of the Boussinesq thermal convection dynamo in a sphere
 * (Toroidal/Poloidal formulation)
 */

// System includes
//

// Project includes
//
#include "Model/Boussinesq/Sphere/Dynamo/Exponential/IExponentialModel.hpp"
#include "Model/Boussinesq/Sphere/Dynamo/Induction.hpp"
#include "Model/Boussinesq/Sphere/Dynamo/Momentum.hpp"
#include "Model/Boussinesq/Sphere/Dynamo/Transport.hpp"
#include "Model/Boussinesq/Sphere/Dynamo/Exponential/InductionJacobian.hpp"
#include "Model/Boussinesq/Sphere/Dynamo/Exponential/MomentumJacobian.hpp"
#include "Model/Boussinesq/Sphere/Dynamo/Exponential/TransportJacobian.hpp"
#include "Model/Boussinesq/Sphere/Dynamo/gitHash.hpp"
#include "QuICC/PhysicalNames/JacobianMagnetic.hpp"
#include "QuICC/PhysicalNames/JacobianTemperature.hpp"
#include "QuICC/PhysicalNames/JacobianVelocity.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Sphere {

namespace Dynamo {

namespace Exponential {

std::vector<std::size_t> IExponentialModel::excludedFieldIds() const
{
   std::vector<std::size_t> fields = {
      PhysicalNames::JacobianMagnetic::id(),
      PhysicalNames::JacobianVelocity::id(),
      PhysicalNames::JacobianTemperature::id(),
   };

   return fields;
}

void IExponentialModel::addEquations(SharedSimulation spSim)
{
   auto optZero = std::make_shared<Equations::EquationOptions>(0, false, false);

   // Add transport equation
   spSim->addEquation<Equations::Boussinesq::Sphere::Dynamo::Transport>(
      this->spBackend(), optZero);

   // Add Navier-Stokes equation
   spSim->addEquation<Equations::Boussinesq::Sphere::Dynamo::Momentum>(
      this->spBackend(), optZero);

   // Add induction equation
   spSim->addEquation<Equations::Boussinesq::Sphere::Dynamo::Induction>(
      this->spBackend(), optZero);

   auto optOne = std::make_shared<Equations::EquationOptions>(1, false, false, false);

   // Add transport jacobian equation
   spSim->addEquation<Equations::Boussinesq::Sphere::Dynamo::Exponential::TransportJacobian>(
      this->spBackend(), optOne);

   // Add Navier-Stokes jacobian equation
   spSim->addEquation<Equations::Boussinesq::Sphere::Dynamo::Exponential::MomentumJacobian>(
      this->spBackend(), optOne);

   // Add induction equation
   spSim->addEquation<Equations::Boussinesq::Sphere::Dynamo::Exponential::InductionJacobian>(
      this->spBackend(), optOne);

   #ifdef QUICC_USE_MLIR_GRAPH
      #error "MLIR Graph setup not implemented for exponential Dynamo model"
   #endif
}

} // namespace Exponential
} // namespace Dynamo
} // namespace Sphere
} // namespace Boussinesq
} // namespace Model
} // namespace QuICC
