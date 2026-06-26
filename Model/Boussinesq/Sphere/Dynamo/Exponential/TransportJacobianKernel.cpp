/**
 * @file TransportJacobianKernel.cpp
 * @brief Source of physical space kernel for the TransportJacobian equation
 */

// System includes
//

// Project includes
//
#include "Model/Boussinesq/Sphere/Dynamo/Exponential/TransportJacobianKernel.hpp"
#include "QuICC/PhysicalOperators/VelocityAdvection.hpp"
#include "QuICC/PhysicalOperators/SphericalHeatAdvection.hpp"

namespace QuICC {

namespace Physical {

namespace Kernel {

void TransportJacobianKernel::setJacobianScalar(std::size_t name,
   Framework::Selector::VariantSharedScalarVariable spField)
{
   // Safety assertion
   assert(this->mScalars.count(name) + this->mVectors.count(name) == 0);

   this->mJSName = name;

   this->setField(name, spField);
}

void TransportJacobianKernel::setScalar(std::size_t name,
   Framework::Selector::VariantSharedScalarVariable spField)
{
   // Safety assertion
   assert(this->mScalars.count(name) + this->mVectors.count(name) == 0);

   this->mSName = name;

   this->setField(name, spField);
}

void TransportJacobianKernel::setJacobianVector(std::size_t name,
   Framework::Selector::VariantSharedVectorVariable spField)
{
   // Safety assertion
   assert(this->mScalars.count(name) + this->mVectors.count(name) == 0);

   this->mJVName = name;

   this->setField(name, spField);
}

void TransportJacobianKernel::setVector(std::size_t name,
   Framework::Selector::VariantSharedVectorVariable spField)
{
   // Safety assertion
   assert(this->mScalars.count(name) + this->mVectors.count(name) == 0);

   this->mVName = name;

   this->setField(name, spField);
}

void TransportJacobianKernel::init(const MHDFloat transport)
{
   this->mTransport= transport;
}

void TransportJacobianKernel::setMesh(std::shared_ptr<std::vector<Array>> spMesh)
{
   IPhysicalKernel::setMesh(spMesh);

   this->mRadius = spMesh->at(0);
}

void TransportJacobianKernel::compute(Framework::Selector::PhysicalScalarField& rNLComp,
   FieldComponents::Physical::Id id) const
{
   // Assert on scalar component is used
   assert(id == FieldComponents::Physical::SCALAR);

   ///
   /// Computation of the advection:
   ///   \f$ \left(\vec u\cdot\nabla\right)\theta\f$
   ///
   std::visit(
      [&](auto&& v, auto&& jv, auto&& t, auto&& jt)
      {
         Physical::VelocityAdvection<FieldComponents::Physical::R,
            FieldComponents::Physical::THETA,
            FieldComponents::Physical::PHI>::set(rNLComp, v->dom(0).phys(), jt->dom(0).grad(),
            this->mTransport);

         Physical::SphericalHeatAdvection<FieldComponents::Physical::R,
            FieldComponents::Physical::THETA,
            FieldComponents::Physical::PHI>::add(rNLComp, jv->dom(0).res(),
            this->mRadius, jv->dom(0).phys(), t->dom(0).grad(),
            this->mTransport);
      },
      this->vector(this->mVName), this->vector(this->mJVName), this->scalar(this->mSName), this->scalar(this->mJSName));
}

} // namespace Kernel
} // namespace Physical
} // namespace QuICC
