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
   virtual void initNLKernel(const bool force = false) override;

protected:
   /**
    * @brief Set nonlinear component path
    */
   virtual void setNLComponents() override;

   /**
    * @brief Set variable requirements
    */
   virtual void setRequirements() override;

   /**
    * @brief Set the equation coupling information
    */
   virtual void setCoupling() override;

private:
};

} // namespace Exponential
} // namespace Dynamo
} // namespace Sphere
} // namespace Boussinesq
} // namespace Equations
} // namespace QuICC

#endif // QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_EXPONENTIAL_TRANSPORTJACOBIAN_HPP
