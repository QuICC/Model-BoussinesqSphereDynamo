/** 
 * @file TransportKernelHeat.cpp
 * @brief Source of physical space kernel for the Transport equation
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Model/Boussinesq/Sphere/Dynamo/TransportKernelHeat.hpp"

// Project includes
//
#include "QuICC/PhysicalOperators/SphericalHeatAdvection.hpp"

namespace QuICC {

namespace Physical {

namespace Kernel {

   TransportKernelHeat::TransportKernelHeat()
      : IPhysicalKernel()
   {
   }

   TransportKernelHeat::~TransportKernelHeat()
   {
   }

   std::size_t TransportKernelHeat::name() const
   {
      return this->mName;
   }

   void TransportKernelHeat::setScalar(std::size_t name, Framework::Selector::VariantSharedScalarVariable spField)
   {
      // Safety assertion
      assert(this->mScalars.count(name) + this->mVectors.count(name) == 0);

      this->mName = name;

      this->setField(name, spField);
   }

   void TransportKernelHeat::setVector(std::size_t name, Framework::Selector::VariantSharedVectorVariable spField)
   {
      // Safety assertion
      assert(this->mScalars.count(name) + this->mVectors.count(name) == 0);

      this->mVName = name;

      this->setField(name, spField);
   }

   void TransportKernelHeat::init(const MHDFloat transport)
   {
      this->mTransport = transport;
   }

   void TransportKernelHeat::setMesh(std::shared_ptr<std::vector<Array> > spMesh)
   {
      IPhysicalKernel::setMesh(spMesh);
      this->mRadius = this->mspMesh->at(0);
   }

   void TransportKernelHeat::compute(Framework::Selector::PhysicalScalarField& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      /// 
      /// Computation of the advection:
      ///   \f$ \left(\vec u\cdot\nabla\right)\theta\f$
      ///
      std::visit([&](auto&& v, auto&& t){Physical::SphericalHeatAdvection<FieldComponents::Physical::R,FieldComponents::Physical::THETA,FieldComponents::Physical::PHI>::set(rNLComp, v->dom(0).res(), this->mRadius, v->dom(0).phys(), t->dom(0).grad(), this->mTransport);}, this->vector(this->mVName), this->scalar(this->name()));
   }

}
}
}
