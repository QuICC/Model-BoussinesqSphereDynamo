/**
 * @file IDynamoState.hpp
 * @brief Implementation of the Boussinesq thermal convection dynamo in a sphere
 * model (Toroidal/Poloidal formulation)
 */

#ifndef QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_IDYNAMOSTATE_HPP
#define QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_IDYNAMOSTATE_HPP

// System includes
//
#include <string>

// Project includes
//
#include "QuICC/Generator/StateGenerator.hpp"
#include "QuICC/Model/IStateGeneratorBuilder.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Sphere {

namespace Dynamo {

/**
 * @brief Implementation of the Boussinesq thermal convection dynamo in a sphere
 * model (Toroidal/Poloidal formulation)
 */
class IDynamoState : public IStateGeneratorBuilder<StateGenerator>
{
public:
   /**
    * @brief Constructor
    */
   IDynamoState() = default;

   /**
    * @brief Destructor
    */
   virtual ~IDynamoState() = default;

   /// Formulation used for vector fields
   virtual VectorFormulation::Id SchemeFormulation() override;

   /**
    * @brief Version string
    */
   std::string version() const final;

   /**
    * @brief Add the initial state generation equations
    *
    * @param spGen   Shared generator object
    */
   virtual void addStates(SharedStateGenerator spGen) override;

protected:
private:
};

} // namespace Dynamo
} // namespace Sphere
} // namespace Boussinesq
} // namespace Model
} // namespace QuICC

#endif // QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_IDYNAMOSTATE_HPP
