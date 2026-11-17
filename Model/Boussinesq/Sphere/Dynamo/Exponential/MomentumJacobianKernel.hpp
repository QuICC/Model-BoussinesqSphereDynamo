/**
 * @file MomentumJacobianKernel.hpp
 * @brief Physical kernel for the MomentumJacobian nonlinear kernel
 */

#ifndef QUICC_PHYSICAL_MOMENTUMJACOBIANKERNEL_HPP
#define QUICC_PHYSICAL_MOMENTUMJACOBIANKERNEL_HPP

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
 * @brief Physical kernel for the MomentumJacobian nonlinear kernel
 */
class MomentumJacobianKernel : public IPhysicalKernel
{
public:
   /**
    * @brief Simple constructor
    */
   MomentumJacobianKernel() = default;

   /**
    * @brief Simple empty destructor
    */
   virtual ~MomentumJacobianKernel() = default;

   /**
    * @brief Set the physical mesh on which kernel is working
    */
   virtual void setMesh(std::shared_ptr<std::vector<Array>> spMesh) override;

   /**
    * @brief Set the smart pointer to the temperature field
    *
    * \param name Name of the field
    * \param spField Shared pointer to the scalar field
    */
   void setJacobianTemperature(std::size_t name,
      Framework::Selector::VariantSharedScalarVariable spField);

   /**
    * @brief Set the smart pointer to the temperature field
    *
    * \param name Name of the field
    * \param spField Shared pointer to the scalar field
    */
   void setTemperature(std::size_t name,
      Framework::Selector::VariantSharedScalarVariable spField);

   /**
    * @brief Set the smart pointer to the velocity field
    *
    * \param name Name of the field
    * \param spField Shared pointer to the vector field
    */
   void setJacobianVelocity(std::size_t name,
      Framework::Selector::VariantSharedVectorVariable spField);

   /**
    * @brief Set the smart pointer to the velocity field
    *
    * \param name Name of the field
    * \param spField Shared pointer to the vector field
    */
   void setVelocity(std::size_t name,
      Framework::Selector::VariantSharedVectorVariable spField);

   /**
    * @brief Set the smart pointer to the magnetic field
    *
    * \param name Name of the field
    * \param spField Shared pointer to the vector field
    */
   void setJacobianMagnetic(std::size_t name,
      Framework::Selector::VariantSharedVectorVariable spField);

   /**
    * @brief Set the smart pointer to the magnetic field
    *
    * \param name Name of the field
    * \param spField Shared pointer to the vector field
    */
   void setMagnetic(std::size_t name,
      Framework::Selector::VariantSharedVectorVariable spField);

   /**
    * @brief Initialize kernel
    */
   void init(const MHDFloat inertia, const MHDFloat coriolis,
      const MHDFloat buoyancy, const MHDFloat lorentz);

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
    * @brief Name ID of the jacobian velocity field
    */
   std::size_t mJVName;

   /**
    * @brief Name ID of the velocity field
    */
   std::size_t mVName;

   /**
    * @brief Name ID of the jacobian temperature field
    */
   std::size_t mJTName;

   /**
    * @brief Name ID of the temperature field
    */
   std::size_t mTName;

   /**
    * @brief Name ID of the jacobian magnetic field
    */
   std::size_t mJMName;

   /**
    * @brief Name ID of the magnetic field
    */
   std::size_t mMName;

   /**
    * @brief Scaling constant for inertial term
    */
   MHDFloat mInertia;

   /**
    * @brief Scaling constant for Coriolis term
    */
   MHDFloat mCoriolis;

   /**
    * @brief Scaling constant for Buoyancy term
    */
   MHDFloat mBuoyancy;

   /**
    * @brief Scaling constant for Lorentz term
    */
   MHDFloat mLorentz;

   /**
    * @brief Storage for the radial grid values (if required)
    */
   Array mRadius;

   /**
    * @brief Storage for the cos(theta) grid values (if required)
    */
   Array mCosTheta;

   /**
    * @brief Storage for the sin(theta) grid values (if required)
    */
   Array mSinTheta;
};

/// Typedef for a smart MomentumJacobianKernel
typedef std::shared_ptr<MomentumJacobianKernel> SharedMomentumJacobianKernel;

} // namespace Kernel
} // namespace Physical
} // namespace QuICC

#endif // QUICC_PHYSICAL_MOMENTUMJACOBIANKERNEL_HPP
