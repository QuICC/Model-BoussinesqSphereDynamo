/**
 * @file IExponentialDynamoModel.hpp
 * @brief Implementation of the Boussinesq thermal convection dynamo in a sphere
 * model (Toroidal/Poloidal formulation)
 */

#ifndef QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_EXPONENTIAL_IEXPONENTIALDYNAMOMODEL_HPP
#define QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_EXPONENTIAL_IEXPONENTIALDYNAMOMODEL_HPP

// System includes
//
#include <string>

// Project includes
//
#include "QuICC/Simulation/Simulation.hpp"
#include "Model/Boussinesq/Sphere/Dynamo/IDynamoModel.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Sphere {

namespace Dynamo {

namespace Exponential {

/**
 * @brief Implementation of the Boussinesq thermal convection dynamo in a sphere
 * model (Toroidal/Poloidal formulation)
 */
class IExponentialDynamoModel : public IDynamoModel
{
public:
   /**
    * @brief Constructor
    */
   IExponentialDynamoModel() = default;

   /**
    * @brief Destructor
    */
   virtual ~IExponentialDynamoModel() = default;

   /**
    * @brief Exclude fields from initial state
    */
   virtual std::vector<std::size_t> excludedFieldIds() const override;

   /**
    * @brief Add the required equations
    *
    * @param spSim   Shared simulation object
    */
   virtual void addEquations(SharedSimulation spSim) override;

protected:
private:
};

} // namespace Exponential
} // namespace Dynamo
} // namespace Sphere
} // namespace Boussinesq
} // namespace Model
} // namespace QuICC

#endif // QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_EXPONENTIAL_IEXPONENTIALDYNAMOMODEL_HPP
