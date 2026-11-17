/**
 * @file InductionJacobian.hpp
 * @brief Implementation of the vector induction equation for the Boussinesq
 * thermal convection dynamo sphere
 */

#ifndef QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_EXPONENTIAL_INDUCTIONJACOBIAN_HPP
#define QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_EXPONENTIAL_INDUCTIONJACOBIAN_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Equations/IVectorEquation.hpp"

namespace QuICC {

namespace Equations {

namespace Boussinesq {

namespace Sphere {

namespace Dynamo {

namespace Exponential {

/**
 * @brief Implementation of the vector induction equation for the Boussinesq
 * thermal convection dynamo in a sphere
 */
class InductionJacobian : public IVectorEquation
{
public:
   /**
    * @brief Simple constructor
    *
    * @param spEqParams  Shared equation parameters
    */
   InductionJacobian(SharedEquationParameters spEqParams,
      SpatialScheme::SharedCISpatialScheme spScheme,
      std::shared_ptr<Model::IModelBackend> spBackend,
      std::shared_ptr<EquationOptions> spOptions);

   /**
    * @brief Simple empty destructor
    */
   virtual ~InductionJacobian() = default;

   /**
    * @brief Initialize nonlinear interaction kernel
    */
   void initNLKernel(const bool force = false) final;

protected:
   /**
    * @brief Set variable requirements
    */
   void setRequirements() final;

   /**
    * @brief Set the equation coupling information
    */
   void setCoupling() final;

   /**
    * @brief Set the nonlinear integration components
    */
   void setNLComponents() final;

private:
};

} // namespace Exponential
} // namespace Dynamo
} // namespace Sphere
} // namespace Boussinesq
} // namespace Equations
} // namespace QuICC

#endif // QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_EXPONENTIAL_INDUCTIONJACOBIAN_HPP
