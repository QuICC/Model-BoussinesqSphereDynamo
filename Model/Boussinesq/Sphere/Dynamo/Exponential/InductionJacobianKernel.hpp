/**
 * @file InductionJacobianKernel.hpp
 * @brief Physical kernel for the InductionJacobian nonlinear kernel
 */

#ifndef QUICC_PHYSICAL_INDUCTIONJACOBIANKERNEL_HPP
#define QUICC_PHYSICAL_INDUCTIONJACOBIANKERNEL_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/PhysicalKernels/IPhysicalKernel.hpp"

namespace QuICC {

namespace Physical {

namespace Kernel {

/**
 * @brief Physical kernel for the induction nonlinear kernel
 */
class InductionJacobianKernel : public IPhysicalKernel
{
public:
   /**
    * @brief Simple constructor
    */
   InductionJacobianKernel() = default;

   /**
    * @brief Simple empty destructor
    */
   virtual ~InductionJacobianKernel() = default;

   /**
    * @brief Set the smart pointer to the velocity field
    *
    * \param name Name of the field
    * \param spField Shared pointer to the vector field
    */
   void setJacobianVelocity(std::size_t,
      Framework::Selector::VariantSharedVectorVariable spField);

   /**
    * @brief Set the smart pointer to the velocity field
    *
    * \param name Name of the field
    * \param spField Shared pointer to the vector field
    */
   void setVelocity(std::size_t,
      Framework::Selector::VariantSharedVectorVariable spField);

   /**
    * @brief Set the smart pointer to the magnetic field
    *
    * \param name Name of the field
    * \param spField Shared pointer to the vector field
    */
   void setJacobianMagnetic(std::size_t,
      Framework::Selector::VariantSharedVectorVariable spField);

   /**
    * @brief Set the smart pointer to the magnetic field
    *
    * \param name Name of the field
    * \param spField Shared pointer to the vector field
    */
   void setMagnetic(std::size_t,
      Framework::Selector::VariantSharedVectorVariable spField);

   /**
    * @brief Initialize kernel
    */
   void init(const MHDFloat induction);

   /**
    * @brief Compute the physical kernel
    *
    * @param rNLComp Nonlinear term component
    * @param id      ID of the component (allows for a more general
    * implementation)
    */
   void compute(Framework::Selector::PhysicalScalarField& rNLComp,
      FieldComponents::Physical::Id id) const final;

protected:

private:
   /**
    * @brief Name ID of the jacobian magnetic field
    */
   std::size_t mJMName;

   /**
    * @brief Name ID of the magnetic field
    */
   std::size_t mMName;

   /**
    * @brief Name ID of the jacobian velocity field
    */
   std::size_t mJVName;

   /**
    * @brief Name ID of the velocity field
    */
   std::size_t mVName;

   /**
    * @brief Scaling constant for inertial term
    */
   MHDFloat mInduction;
};

/// Typedef for a smart InductionJacobianKernel
typedef std::shared_ptr<InductionJacobianKernel> SharedInductionJacobianKernel;

} // namespace Kernel
} // namespace Physical
} // namespace QuICC

#endif // QUICC_PHYSICAL_INDUCTIONJACOBIANKERNEL_HPP
