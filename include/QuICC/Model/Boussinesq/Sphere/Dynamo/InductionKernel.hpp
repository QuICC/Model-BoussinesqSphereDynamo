/**
 * @file InductionKernel.hpp
 * @brief Physical kernel for the Induction nonlinear kernel
 */

#ifndef QUICC_PHYSICAL_INDUCTIONKERNEL_HPP
#define QUICC_PHYSICAL_INDUCTIONKERNEL_HPP

// First include
//

// Configuration includes
//
#include <memory>

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/PhysicalKernels/IPhysicalKernel.hpp"

namespace QuICC {

namespace Physical {

namespace Kernel {

   /**
    * @brief Physical kernel for the induction nonlinear kernel
    */
   class InductionKernel: public IPhysicalKernel
   {
      public:
         /**
          * @brief Simple constructor
          */
         explicit InductionKernel();

         /**
          * @brief Simple empty destructor
          */
         virtual ~InductionKernel();

         /**
          * @brief Set the smart pointer to the velocity field
          *
          * \param name Name of the field
          * \param spField Shared pointer to the vector field
          */
         virtual void setVelocity(std::size_t, Framework::Selector::VariantSharedVectorVariable spField);

         /**
          * @brief Set the smart pointer to the magnetic field
          *
          * \param name Name of the field
          * \param spField Shared pointer to the vector field
          */
         virtual void setMagnetic(std::size_t, Framework::Selector::VariantSharedVectorVariable spField);

         /**
          * @brief Initialize kernel
          */
         void init(const MHDFloat induction);

         /**
          * @brief Compute the physical kernel
          *
          * @param rNLComp Nonlinear term component
          * @param id      ID of the component (allows for a more general implementation)
          */
         virtual void compute(Framework::Selector::PhysicalScalarField& rNLComp, FieldComponents::Physical::Id id) const;

      protected:
         /**
          * @brief Get name ID of the unknown
          */
         std::size_t name() const;

      private:
         /**
          * @brief Name ID of the magnetic field
          */
         std::size_t mName;

         /**
          * @brief Name ID of the velocity field
          */
         std::size_t mVelName;

         /**
          * @brief Scaling constant for inertial term
          */
         MHDFloat mInduction;

   };

   /// Typedef for a smart InductionKernel
   typedef std::shared_ptr<InductionKernel> SharedInductionKernel;

}
}
}

#endif // QUICC_PHYSICAL_INDUCTIONKERNEL_HPP
