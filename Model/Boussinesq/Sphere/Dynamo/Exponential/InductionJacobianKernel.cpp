/**
 * @file InductionJacobianKernel.cpp
 * @brief Source of physical space kernel for the InductionJacobian equation
 */

// System includes
//

// Project includes
//
#include "Model/Boussinesq/Sphere/Dynamo/Exponential/InductionJacobianKernel.hpp"
#include "QuICC/PhysicalOperators/Cross.hpp"

namespace QuICC {

namespace Physical {

namespace Kernel {

void InductionJacobianKernel::setJacobianMagnetic(std::size_t name,
   Framework::Selector::VariantSharedVectorVariable spField)
{
   // Safety assertion
   assert(this->mScalars.count(name) + this->mVectors.count(name) == 0);

   this->mJMName = name;

   this->setField(name, spField);
}

void InductionJacobianKernel::setMagnetic(std::size_t name,
   Framework::Selector::VariantSharedVectorVariable spField)
{
   // Safety assertion
   assert(this->mScalars.count(name) + this->mVectors.count(name) == 0);

   this->mMName = name;

   this->setField(name, spField);
}

void InductionJacobianKernel::setJacobianVelocity(std::size_t name,
   Framework::Selector::VariantSharedVectorVariable spField)
{
   // Safety assertion
   assert(this->mScalars.count(name) + this->mVectors.count(name) == 0);

   this->mJVName = name;

   this->setField(name, spField);
}

void InductionJacobianKernel::setVelocity(std::size_t name,
   Framework::Selector::VariantSharedVectorVariable spField)
{
   // Safety assertion
   assert(this->mScalars.count(name) + this->mVectors.count(name) == 0);

   this->mVName = name;

   this->setField(name, spField);
}

void InductionJacobianKernel::init(const MHDFloat induction)
{
   // Set scaling constants
   this->mInduction= induction;
}

void InductionJacobianKernel::compute(Framework::Selector::PhysicalScalarField& rNLComp,
   FieldComponents::Physical::Id id) const
{
   ///
   /// Compute \f$\left(\vec u\wedge\vec B\right)\f$
   ///
   std::visit(
      [&](auto&& m, auto&& jm, auto&& v, auto&& jv)
      {
         switch (id)
         {
         case (FieldComponents::Physical::R):
            Physical::Cross<FieldComponents::Physical::THETA,
               FieldComponents::Physical::PHI>::set(rNLComp, jm->dom(0).phys(),
               v->dom(0).phys(), this->mInduction);
            Physical::Cross<FieldComponents::Physical::THETA,
               FieldComponents::Physical::PHI>::add(rNLComp, m->dom(0).phys(),
               jv->dom(0).phys(), this->mInduction);
            break;
         case (FieldComponents::Physical::THETA):
            Physical::Cross<FieldComponents::Physical::PHI,
               FieldComponents::Physical::R>::set(rNLComp, jm->dom(0).phys(),
               v->dom(0).phys(), this->mInduction);
            Physical::Cross<FieldComponents::Physical::PHI,
               FieldComponents::Physical::R>::add(rNLComp, m->dom(0).phys(),
               jv->dom(0).phys(), this->mInduction);
            break;
         case (FieldComponents::Physical::PHI):
            Physical::Cross<FieldComponents::Physical::R,
               FieldComponents::Physical::THETA>::set(rNLComp, jm->dom(0).phys(),
               v->dom(0).phys(), this->mInduction);
            Physical::Cross<FieldComponents::Physical::R,
               FieldComponents::Physical::THETA>::add(rNLComp, m->dom(0).phys(),
               jv->dom(0).phys(), this->mInduction);
            break;
         default:
            assert(false);
            break;
         }
      },
      this->vector(this->mMName), this->vector(this->mJMName), this->vector(this->mVName), this->vector(this->mJVName));
}

} // namespace Kernel
} // namespace Physical
} // namespace QuICC
