/**
 * @file TransportJacobian.hpp
 * @brief Implementation of the transport equation for the Boussinesq thermal
 * convection dynamo in a sphere
 */

#ifndef QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_EXPONENTIAL_TRANSPORTJACOBIAN_HPP
#define QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_EXPONENTIAL_TRANSPORTJACOBIAN_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Equations/IScalarEquation.hpp"

namespace QuICC {

namespace Equations {

namespace Boussinesq {

namespace Sphere {

namespace Dynamo {

namespace Exponential {

/**
 * @brief Implementation of the transport equation for the Boussinesq thermal
 * convection dynamo in a sphere
 */
class TransportJacobian : public IScalarEquation
{
public:
   /**
    * @brief Simple constructor
    *
    * @param spEqParams  Shared equation parameters
    */
   TransportJacobian(SharedEquationParameters spEqParams,
      SpatialScheme::SharedCISpatialScheme spScheme,
      std::shared_ptr<Model::IModelBackend> spBackend,
      std::shared_ptr<EquationOptions> spOptions);

   /**
    * @brief Simple empty destructor
    */
   virtual ~TransportJacobian() = default;

   /**
    * @brief Initialize nonlinear interaction kernel
    */
   void initNLKernel(const bool force = false) final;

protected:
   /**
    * @brief Set nonlinear component path
    */
   void setNLComponents() final;

   /**
    * @brief Set variable requirements
    */
   void setRequirements() final;

   /**
    * @brief Set the equation coupling information
    */
   void setCoupling() final;

private:
};

} // namespace Exponential
} // namespace Dynamo
} // namespace Sphere
} // namespace Boussinesq
} // namespace Equations
} // namespace QuICC

#endif // QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_EXPONENTIAL_TRANSPORTJACOBIAN_HPP
