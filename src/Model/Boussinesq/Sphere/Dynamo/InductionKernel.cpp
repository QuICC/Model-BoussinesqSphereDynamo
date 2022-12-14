/**
 * @file InductionKernel.cpp
 * @brief Source of physical space kernel for the Induction equation
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Model/Boussinesq/Sphere/Dynamo/InductionKernel.hpp"

// Project includes
//
#include "QuICC/PhysicalOperators/Cross.hpp"

namespace QuICC {

namespace Physical {

namespace Kernel {

   InductionKernel::InductionKernel()
      : IPhysicalKernel()
   {
   }

   InductionKernel::~InductionKernel()
   {
   }

   std::size_t InductionKernel::name() const
   {
      return this->mName;
   }

   void InductionKernel::setMagnetic(std::size_t name, Framework::Selector::VariantSharedVectorVariable spField)
   {
      // Safety assertion
      assert(this->mScalars.count(name) + this->mVectors.count(name) == 0);

      this->mName = name;

      this->setField(name, spField);
   }

   void InductionKernel::setVelocity(std::size_t name, Framework::Selector::VariantSharedVectorVariable spField)
   {
      // Safety assertion
      assert(this->mScalars.count(name) + this->mVectors.count(name) == 0);

      this->mVelName = name;

      this->setField(name, spField);
   }

   void InductionKernel::init(const MHDFloat induction)
   {
      // Set scaling constants
      this->mInduction = induction;
   }

   void InductionKernel::compute(Framework::Selector::PhysicalScalarField& rNLComp, FieldComponents::Physical::Id id) const
   {
      ///
      /// Compute \f$\left(\vec u\wedge\vec B\right)\f$
      ///
      switch(id)
      {
         case(FieldComponents::Physical::R):
            std::visit([&](auto&& m, auto&& v){Physical::Cross<FieldComponents::Physical::THETA,FieldComponents::Physical::PHI>::set(rNLComp, m->dom(0).phys(), v->dom(0).phys(), this->mInduction);}, this->vector(this->name()), this->vector(this->mVelName));
            break;
         case(FieldComponents::Physical::THETA):
            std::visit([&](auto&& m, auto&& v){Physical::Cross<FieldComponents::Physical::PHI,FieldComponents::Physical::R>::set(rNLComp, m->dom(0).phys(), v->dom(0).phys(), this->mInduction);}, this->vector(this->name()), this->vector(this->mVelName));
            break;
         case(FieldComponents::Physical::PHI):
            std::visit([&](auto&& m, auto&& v){Physical::Cross<FieldComponents::Physical::R,FieldComponents::Physical::THETA>::set(rNLComp, m->dom(0).phys(), v->dom(0).phys(), this->mInduction);}, this->vector(this->name()), this->vector(this->mVelName));
            break;
         default:
            assert(false);
            break;
      }
   }

}
}
}
