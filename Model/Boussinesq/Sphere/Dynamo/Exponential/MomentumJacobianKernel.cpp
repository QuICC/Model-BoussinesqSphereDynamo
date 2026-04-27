/**
 * @file MomentumJacobianKernel.cpp
 * @brief Source of physical space kernel for the MomentumJacobian equation
 */

// System includes
//

// Project includes
//
#include "Model/Boussinesq/Sphere/Dynamo/Exponential/MomentumJacobianKernel.hpp"
#include "QuICC/PhysicalOperators/Cross.hpp"
#include "QuICC/PhysicalOperators/SphericalBuoyancy.hpp"
#include "QuICC/PhysicalOperators/SphericalCoriolis.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"

namespace QuICC {

namespace Physical {

namespace Kernel {

void MomentumJacobianKernel::setJacobianVelocity(std::size_t name,
   Framework::Selector::VariantSharedVectorVariable spField)
{
   // Safety assertion
   assert(this->mScalars.count(name) + this->mVectors.count(name) == 0);

   this->mJVName = name;

   this->setField(name, spField);
}

void MomentumJacobianKernel::setVelocity(std::size_t name,
   Framework::Selector::VariantSharedVectorVariable spField)
{
   // Safety assertion
   assert(this->mScalars.count(name) + this->mVectors.count(name) == 0);

   this->mVName = name;

   this->setField(name, spField);
}

void MomentumJacobianKernel::setJacobianMagnetic(std::size_t name,
   Framework::Selector::VariantSharedVectorVariable spField)
{
   // Safety assertion
   assert(this->mScalars.count(name) + this->mVectors.count(name) == 0);

   this->mJMName = name;

   this->setField(name, spField);
}

void MomentumJacobianKernel::setMagnetic(std::size_t name,
   Framework::Selector::VariantSharedVectorVariable spField)
{
   // Safety assertion
   assert(this->mScalars.count(name) + this->mVectors.count(name) == 0);

   this->mMName = name;

   this->setField(name, spField);
}

void MomentumJacobianKernel::setJacobianTemperature(std::size_t name,
   Framework::Selector::VariantSharedScalarVariable spField)
{
   // Safety assertion
   assert(this->mScalars.count(name) + this->mVectors.count(name) == 0);

   this->mJTName = name;

   this->setField(name, spField);
}

void MomentumJacobianKernel::setTemperature(std::size_t name,
   Framework::Selector::VariantSharedScalarVariable spField)
{
   // Safety assertion
   assert(this->mScalars.count(name) + this->mVectors.count(name) == 0);

   this->mTName = name;

   this->setField(name, spField);
}

void MomentumJacobianKernel::init(const MHDFloat inertia, const MHDFloat coriolis,
   const MHDFloat buoyancy, const MHDFloat lorentz)
{
   // Set scaling constants
   this->mInertia = inertia;
   this->mCoriolis = coriolis;
   this->mLorentz = lorentz;
   this->mBuoyancy = buoyancy;
}

void MomentumJacobianKernel::setMesh(std::shared_ptr<std::vector<Array>> spMesh)
{
   IPhysicalKernel::setMesh(spMesh);

   this->mRadius = this->mspMesh->at(0);

   if (std::visit(
          [&](auto&& v) -> bool
          {
             return (v->dom(0).res().sim().ss().has(
                SpatialScheme::Feature::SpectralOrdering132));
          },
          this->vector(this->mJVName)))
   {
      this->mCosTheta = this->mspMesh->at(1).array().cos();
      this->mSinTheta = this->mspMesh->at(1).array().sin();
   }
}

void MomentumJacobianKernel::compute(Framework::Selector::PhysicalScalarField& rNLComp,
   FieldComponents::Physical::Id id) const
{
   ///
   /// Compute \f$\left(\nabla\wedge\vec u\right)\wedge\vec u\f$
   ///
   std::visit(
      [&](auto&& v, auto&& jv, auto&& t, auto&& jt, auto&& m, auto&& jm)
      {
         switch (id)
         {
         case (FieldComponents::Physical::R):
            Physical::Cross<FieldComponents::Physical::THETA,
               FieldComponents::Physical::PHI>::set(rNLComp, v->dom(0).curl(),
               jv->dom(0).phys(), this->mInertia);
            Physical::Cross<FieldComponents::Physical::THETA,
               FieldComponents::Physical::PHI>::add(rNLComp, jv->dom(0).curl(),
               v->dom(0).phys(), this->mInertia);
            Physical::Cross<FieldComponents::Physical::THETA,
               FieldComponents::Physical::PHI>::add(rNLComp, m->dom(0).phys(),
               jm->dom(0).curl(), this->mLorentz);
            Physical::Cross<FieldComponents::Physical::THETA,
               FieldComponents::Physical::PHI>::add(rNLComp, jm->dom(0).phys(),
               m->dom(0).curl(), this->mLorentz);
            break;
         case (FieldComponents::Physical::THETA):
            Physical::Cross<FieldComponents::Physical::PHI,
               FieldComponents::Physical::R>::set(rNLComp, v->dom(0).curl(),
               jv->dom(0).phys(), this->mInertia);
            Physical::Cross<FieldComponents::Physical::PHI,
               FieldComponents::Physical::R>::add(rNLComp, jv->dom(0).curl(),
               v->dom(0).phys(), this->mInertia);
            Physical::Cross<FieldComponents::Physical::PHI,
               FieldComponents::Physical::R>::add(rNLComp, m->dom(0).phys(),
               jm->dom(0).curl(), this->mLorentz);
            Physical::Cross<FieldComponents::Physical::PHI,
               FieldComponents::Physical::R>::add(rNLComp, jm->dom(0).phys(),
               m->dom(0).curl(), this->mLorentz);
      break;
         case (FieldComponents::Physical::PHI):
            Physical::Cross<FieldComponents::Physical::R,
               FieldComponents::Physical::THETA>::set(rNLComp, v->dom(0).curl(),
               jv->dom(0).phys(), this->mInertia);
            Physical::Cross<FieldComponents::Physical::R,
               FieldComponents::Physical::THETA>::add(rNLComp, jv->dom(0).curl(),
               v->dom(0).phys(), this->mInertia);
            Physical::Cross<FieldComponents::Physical::R,
               FieldComponents::Physical::THETA>::add(rNLComp, m->dom(0).phys(),
               jm->dom(0).curl(), this->mLorentz);
            Physical::Cross<FieldComponents::Physical::R,
               FieldComponents::Physical::THETA>::add(rNLComp, jm->dom(0).phys(),
               m->dom(0).curl(), this->mLorentz);
            break;
         default:
            assert(false);
            break;
         }

         Physical::SphericalBuoyancy::sub(rNLComp, id, jt->dom(0).res(),
            this->mRadius, jt->dom(0).phys(), this->mBuoyancy);

         if (v->dom(0).res().sim().ss().has(
                      SpatialScheme::Feature::SpectralOrdering132))
         {
            ///
            /// Compute Coriolis term
            ///
            Physical::SphericalCoriolis::add(rNLComp, id, jv->dom(0).res(),
               this->mCosTheta, this->mSinTheta, jv->dom(0).phys(),
               this->mCoriolis);
         }
      },
      this->vector(this->mVName), this->vector(this->mJVName), this->scalar(this->mTName), this->scalar(this->mJTName), this->vector(this->mMName), this->vector(this->mJMName));
}

} // namespace Kernel
} // namespace Physical
} // namespace QuICC
