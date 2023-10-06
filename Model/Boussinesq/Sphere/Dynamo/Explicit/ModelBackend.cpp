/**
 * @file ModelBackend.cpp
 * @brief Source of the interface for model backend
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "Model/Boussinesq/Sphere/Dynamo/Explicit/ModelBackend.hpp"
#include "QuICC/ModelOperator/Time.hpp"
#include "QuICC/ModelOperator/ImplicitLinear.hpp"
#include "QuICC/ModelOperator/ExplicitLinear.hpp"
#include "QuICC/ModelOperator/ExplicitNonlinear.hpp"
#include "QuICC/ModelOperator/ExplicitNextstep.hpp"
#include "QuICC/ModelOperator/Stencil.hpp"
#include "QuICC/ModelOperator/Boundary.hpp"
#include "QuICC/ModelOperator/SplitImplicitLinear.hpp"
#include "QuICC/ModelOperator/SplitBoundary.hpp"
#include "QuICC/ModelOperator/SplitBoundaryValue.hpp"
#include "QuICC/ModelOperatorBoundary/FieldToRhs.hpp"
#include "QuICC/ModelOperatorBoundary/SolverHasBc.hpp"
#include "QuICC/ModelOperatorBoundary/SolverNoTau.hpp"
#include "QuICC/ModelOperatorBoundary/Stencil.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Bc/Name/FixedTemperature.hpp"
#include "QuICC/Bc/Name/FixedFlux.hpp"
#include "QuICC/Bc/Name/StressFree.hpp"
#include "QuICC/Bc/Name/NoSlip.hpp"
#include "QuICC/Bc/Name/Insulating.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"
#include "QuICC/PhysicalNames/Magnetic.hpp"
#include "QuICC/PhysicalNames/Temperature.hpp"
#include "QuICC/NonDimensional/Prandtl.hpp"
#include "QuICC/NonDimensional/MagneticPrandtl.hpp"
#include "QuICC/NonDimensional/Ekman.hpp"
#include "QuICC/NonDimensional/Rayleigh.hpp"
#include "QuICC/NonDimensional/CflInertial.hpp"
#include "QuICC/NonDimensional/CflTorsional.hpp"
#include "QuICC/NonDimensional/CflAlfvenScale.hpp"
#include "QuICC/NonDimensional/CflAlfvenDamping.hpp"
#include "QuICC/Tools/IdToHuman.hpp"
#include "QuICC/Resolutions/Tools/IndexCounter.hpp"
#include "QuICC/SparseSM/Worland/Id.hpp"
#include "QuICC/SparseSM/Worland/I2.hpp"
#include "QuICC/SparseSM/Worland/I4.hpp"
#include "QuICC/SparseSM/Worland/I2Lapl.hpp"
#include "QuICC/SparseSM/Worland/I4Lapl.hpp"
#include "QuICC/SparseSM/Worland/I4Lapl2.hpp"
#include "QuICC/SparseSM/Worland/Boundary/Value.hpp"
#include "QuICC/SparseSM/Worland/Boundary/D1.hpp"
#include "QuICC/SparseSM/Worland/Boundary/D2.hpp"
#include "QuICC/SparseSM/Worland/Boundary/R1D1DivR1.hpp"
#include "QuICC/SparseSM/Worland/Boundary/InsulatingSphere.hpp"
#include "QuICC/SparseSM/Worland/Boundary/Operator.hpp"
#include "QuICC/SparseSM/Worland/Stencil/Value.hpp"
#include "QuICC/SparseSM/Worland/Stencil/D1.hpp"
#include "QuICC/SparseSM/Worland/Stencil/R1D1DivR1.hpp"
#include "QuICC/SparseSM/Worland/Stencil/InsulatingSphere.hpp"
#include "QuICC/SparseSM/Worland/Stencil/ValueD1.hpp"
#include "QuICC/SparseSM/Worland/Stencil/ValueD2.hpp"

#include "QuICC/Polynomial/Worland/WorlandBase.hpp"
#include "QuICC/Equations/CouplingIndexType.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Sphere {

namespace Dynamo {

namespace Explicit {

   ModelBackend::ModelBackend()
      : IDynamoBackend(),
#ifdef QUICC_TRANSFORM_WORLAND_TRUNCATE_QI
      mcTruncateQI(true)
#else
      mcTruncateQI(false)
#endif // QUICC_TRANSFORM_WORLAND_TRUNCATE_QI
   {
   }

   ModelBackend::SpectralFieldIds ModelBackend::implicitFields(const SpectralFieldId& fId) const
   {
      SpectralFieldIds fields = {fId};

      return fields;
   }

   void ModelBackend::equationInfo(EquationInfo& info, const SpectralFieldId& fId, const Resolution& res) const
   {
      // Operators are real
      info.isComplex = false;

      // Splitting 4th poloidal equation into two systems
      if(fId == std::make_pair(PhysicalNames::Velocity::id(), FieldComponents::Spectral::POL))
      {
         info.isSplitEquation = this->useSplitEquation();
      }
      else
      {
         info.isSplitEquation = false;
      }

      // Implicit coupled fields
      info.im = this->implicitFields(fId);

      // Explicit linear terms
      info.exL.clear();

      // Explicit nonlinear terms
      info.exNL.clear();

      // Explicit nextstep terms
      info.exNS.clear();

      // Index mode
      info.indexMode = static_cast<int>(Equations::CouplingIndexType::SLOWEST_MULTI_RHS);
   }

   void ModelBackend::blockSize(int& tN, int& gN, ArrayI& shift, int& rhs, const SpectralFieldId& fId, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs) const
   {
      assert(eigs.size() == 1);
      int l = eigs.at(0);
      this->blockInfo(tN, gN, shift, rhs, fId, res, l,bcs);
   }

   void ModelBackend::operatorInfo(OperatorInfo& info, const SpectralFieldId& fId, const Resolution& res, const Equations::Tools::ICoupling& coupling, const BcMap& bcs) const
   {
      // Loop overall matrices/eigs
      for(int idx = 0; idx < info.tauN.size(); ++idx)
      {
         auto eigs = coupling.getIndexes(res, idx);

         int tN, gN, rhs;
         ArrayI shift(3);

         this->blockSize(tN, gN, shift, rhs, fId, res, eigs, bcs);

         info.tauN(idx) = tN;
         info.galN(idx) = gN;
         info.galShift.row(idx) = shift;
         info.rhsCols(idx) = rhs;

         // Compute system size
         int sN = 0;
         for(auto f: this->implicitFields(fId))
         {
            this->blockSize(tN, gN, shift, rhs, f, res, eigs, bcs);
            sN += gN;
         }

         if(sN == 0)
         {
            sN = info.galN(idx);
         }

         info.sysN(idx) = sN;
      }
   }

   void ModelBackend::implicitBlock(DecoupledZSparse& decMat, const SpectralFieldId& rowId, const SpectralFieldId& colId, const int matIdx, const Resolution& res, const std::vector<MHDFloat>& eigs, const NonDimensional::NdMap& nds, const bool isSplitOperator) const
   {
      assert(eigs.size() == 1);
      int l = eigs.at(0);

      auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);

      auto a = Polynomial::Worland::WorlandBase::ALPHA_CHEBYSHEV;
      auto b = Polynomial::Worland::WorlandBase::DBETA_CHEBYSHEV;

      auto Pr = nds.find(NonDimensional::Prandtl::id())->second->value();
      auto Pm = nds.find(NonDimensional::MagneticPrandtl::id())->second->value();

      if(rowId == std::make_pair(PhysicalNames::Velocity::id(),FieldComponents::Spectral::TOR) && rowId == colId)
      {
         SparseSM::Worland::I2Lapl spasm(nN, nN, a, b, l, 1*this->mcTruncateQI);
         decMat.real() = Pm*spasm.mat();
      }
      else if(rowId == std::make_pair(PhysicalNames::Velocity::id(),FieldComponents::Spectral::POL) && rowId == colId)
      {
         // Split fourth order system
         if(this->useSplitEquation())
         {
            if(isSplitOperator)
            {
               SparseSM::Worland::I2Lapl spasm(nN, nN, a, b, l, 1*this->mcTruncateQI);
               decMat.real() = Pm*spasm.mat();
            }
            else
            {
               SparseSM::Worland::I2Lapl spasm(nN, nN, a, b, l, 1*this->mcTruncateQI);
               decMat.real() = spasm.mat();
            }
         }
         else
         {
            SparseSM::Worland::I4Lapl2 spasm(nN, nN, a, b, l, 2*this->mcTruncateQI);
            decMat.real() = Pm*spasm.mat();
         }
      }
      else if(rowId == std::make_pair(PhysicalNames::Magnetic::id(),FieldComponents::Spectral::TOR) && rowId == colId)
      {
         SparseSM::Worland::I2Lapl spasm(nN, nN, a, b, l, 1*this->mcTruncateQI);
         decMat.real() = spasm.mat();
      }
      else if(rowId == std::make_pair(PhysicalNames::Magnetic::id(),FieldComponents::Spectral::POL) && rowId == colId)
      {
         SparseSM::Worland::I2Lapl spasm(nN, nN, a, b, l, 1*this->mcTruncateQI);
         decMat.real() = spasm.mat();
      }
      else if(rowId == std::make_pair(PhysicalNames::Temperature::id(), FieldComponents::Spectral::SCALAR) && rowId == colId)
      {
         SparseSM::Worland::I2Lapl spasm(nN, nN, a, b, l, 1*this->mcTruncateQI);
         decMat.real() = (Pm/Pr)*spasm.mat();
      }
      else
      {
         throw std::logic_error("Equations are not setup properly");
      }
   }

   void ModelBackend::timeBlock(DecoupledZSparse& decMat, const SpectralFieldId& fieldId, const int matIdx, const Resolution& res, const std::vector<MHDFloat>& eigs, const NonDimensional::NdMap& nds) const
   {
      assert(eigs.size() == 1);
      int l = eigs.at(0);

      auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);

      auto a = Polynomial::Worland::WorlandBase::ALPHA_CHEBYSHEV;
      auto b = Polynomial::Worland::WorlandBase::DBETA_CHEBYSHEV;

      if(fieldId == std::make_pair(PhysicalNames::Velocity::id(),FieldComponents::Spectral::TOR))
      {
         SparseSM::Worland::I2 spasm(nN, nN, a, b, l, 1*this->mcTruncateQI);
         decMat.real() = spasm.mat();
      }
      else if(fieldId == std::make_pair(PhysicalNames::Velocity::id(),FieldComponents::Spectral::POL))
      {
         if(this->useSplitEquation())
         {
            SparseSM::Worland::I2 spasm(nN, nN, a, b, l, 1*this->mcTruncateQI);
            decMat.real() = spasm.mat();
         }
         else
         {
            SparseSM::Worland::I4Lapl spasm(nN, nN, a, b, l, 2*this->mcTruncateQI);
            decMat.real() = spasm.mat();
         }
      }
      else if(fieldId == std::make_pair(PhysicalNames::Magnetic::id(),FieldComponents::Spectral::TOR))
      {
         SparseSM::Worland::I2 spasm(nN, nN, a, b, l, 1*this->mcTruncateQI);
         decMat.real() = spasm.mat();
      }
      else if(fieldId == std::make_pair(PhysicalNames::Magnetic::id(),FieldComponents::Spectral::POL))
      {
         SparseSM::Worland::I2 spasm(nN, nN, a, b, l, 1*this->mcTruncateQI);
         decMat.real() = spasm.mat();
      }
      else if(fieldId == std::make_pair(PhysicalNames::Temperature::id(), FieldComponents::Spectral::SCALAR))
      {
         SparseSM::Worland::I2 spasm(nN, nN, a, b, l, 1*this->mcTruncateQI);
         decMat.real() = spasm.mat();
      }
   }

   void ModelBackend::splitBoundaryValueBlock(DecoupledZSparse& decMat, const SpectralFieldId& fieldId, const int matIdx, const Resolution& res, const std::vector<MHDFloat>& eigs, const NonDimensional::NdMap& nds) const
   {
      assert(eigs.size() == 1);

      int l = eigs.at(0);
      auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);

      if(fieldId == std::make_pair(PhysicalNames::Velocity::id(),FieldComponents::Spectral::POL))
      {
         decMat.real().resize(nN, 1);
         decMat.imag().resize(nN, 1);

         Eigen::Triplet<MHDFloat> val = {0, 0, 1.0};
         std::vector<Eigen::Triplet<MHDFloat> > triplets = {val};
         decMat.real().setFromTriplets(triplets.begin(), triplets.end());
         decMat.imag().setFromTriplets(triplets.begin(), triplets.end());
      }
   }

   void ModelBackend::modelMatrix(DecoupledZSparse& rModelMatrix, const std::size_t opId, const Equations::CouplingInformation::FieldId_range imRange, const int matIdx, const std::size_t bcType, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds) const
   {
      assert(eigs.size() == 1);
      int l = eigs.at(0);

      // Time operator
      if(opId == ModelOperator::Time::id())
      {
         bool needStencil = (this->useGalerkin() && bcType == ModelOperatorBoundary::SolverNoTau::id());
         bool needTau = bcType == ModelOperatorBoundary::SolverHasBc::id();

         for(auto pRowId = imRange.first; pRowId != imRange.second; pRowId++)
         {
            this->timeBlock(rModelMatrix, *pRowId, matIdx, res, eigs, nds);

            // Apply boundary condition
            if(needStencil)
            {
               this->applyGalerkinStencil(rModelMatrix.real(), *pRowId, *pRowId, l, res, bcs, nds);
            }
            else if(needTau)
            {
               this->applyTau(rModelMatrix.real(), *pRowId, *pRowId, l, res, bcs, nds, false);
            }
         }
      }
      // Linear operator
      else if(opId == ModelOperator::ImplicitLinear::id() || opId == ModelOperator::SplitImplicitLinear::id())
      {
         bool isSplit = (opId == ModelOperator::SplitImplicitLinear::id());
         bool needStencil = (this->useGalerkin() && bcType == ModelOperatorBoundary::SolverNoTau::id());
         bool needTau = bcType == ModelOperatorBoundary::SolverHasBc::id();

         for(auto pRowId = imRange.first; pRowId != imRange.second; pRowId++)
         {
            for(auto pColId = imRange.first; pColId != imRange.second; pColId++)
            {
               this->implicitBlock(rModelMatrix, *pRowId, *pColId, matIdx, res, eigs, nds, isSplit);

               // Apply boundary condition
               if(needStencil)
               {
                  this->applyGalerkinStencil(rModelMatrix.real(), *pRowId, *pColId, l, res, bcs, nds);
               }
               else if(needTau)
               {
                  this->applyTau(rModelMatrix.real(), *pRowId, *pColId, l, res, bcs, nds, isSplit);
               }
            }
         }
      }
      // Boundary operator
      else if(opId == ModelOperator::Boundary::id() || opId == ModelOperator::SplitBoundary::id())
      {
         bool isSplit = (opId == ModelOperator::SplitBoundary::id());
         bool needStencil = this->useGalerkin();
         bool needTau = (bcType == ModelOperatorBoundary::SolverHasBc::id());

         auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);

         for(auto pRowId = imRange.first; pRowId != imRange.second; pRowId++)
         {
            for(auto pColId = imRange.first; pColId != imRange.second; pColId++)
            {
               rModelMatrix.real().resize(nN, nN);

               // Apply boundary condition
               if(needStencil)
               {
                  this->applyGalerkinStencil(rModelMatrix.real(), *pRowId, *pColId, l, res, bcs, nds);
               }
               else if(needTau)
               {
                  this->applyTau(rModelMatrix.real(), *pRowId, *pColId, l, res, bcs, nds, isSplit);
               }
            }
         }
      }
      // Split equation boundary value
      else if(opId == ModelOperator::SplitBoundaryValue::id())
      {
         for(auto pRowId = imRange.first; pRowId != imRange.second; pRowId++)
         {
            this->splitBoundaryValueBlock(rModelMatrix, *pRowId, matIdx, res, eigs, nds);
         }
      }
      else
      {
         throw std::logic_error("Requested operator type is not implemented");
      }
   }

   void ModelBackend::galerkinStencil(SparseMatrix& mat, const SpectralFieldId& fieldId, const int matIdx, const Resolution& res, const std::vector<MHDFloat>& eigs, const bool makeSquare, const BcMap& bcs, const NonDimensional::NdMap& nds) const
   {
      assert(eigs.size() == 1);
      int l = eigs.at(0);
      this->stencil(mat, fieldId, l, res, makeSquare, bcs, nds);
   }

   void ModelBackend::explicitBlock(DecoupledZSparse& mat, const SpectralFieldId& fId, const std::size_t opId,  const SpectralFieldId fieldId, const int matIdx, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds) const
   {
      // Explicit linear operator
      if(opId == ModelOperator::ExplicitLinear::id())
      {
         // Nothing to be done
         throw std::logic_error("There are no explicit linear operators");
      }
      // Explicit nonlinear operator
      else if(opId == ModelOperator::ExplicitNonlinear::id())
      {
         throw std::logic_error("There are no explicit nonlinear operators");
      }
      // Explicit nextstep operator
      else if(opId == ModelOperator::ExplicitNextstep::id())
      {
         throw std::logic_error("There are no explicit nextstep operators");
      }
   }

} // Explicit
} // Dynamo
} // Sphere
} // Boussinesq
} // Model
} // QuICC
