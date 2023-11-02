/**
 * @file ModelBackend.cpp
 * @brief Source of the interface for model backend
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "Model/Boussinesq/Sphere/Dynamo/Implicit/ModelBackend.hpp"
#include "QuICC/Bc/Name/FixedFlux.hpp"
#include "QuICC/Bc/Name/FixedTemperature.hpp"
#include "QuICC/Bc/Name/NoSlip.hpp"
#include "QuICC/Bc/Name/StressFree.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Equations/CouplingIndexType.hpp"
#include "QuICC/ModelOperator/Boundary.hpp"
#include "QuICC/ModelOperator/ExplicitLinear.hpp"
#include "QuICC/ModelOperator/ExplicitNextstep.hpp"
#include "QuICC/ModelOperator/ExplicitNonlinear.hpp"
#include "QuICC/ModelOperator/ImplicitLinear.hpp"
#include "QuICC/ModelOperator/SplitBoundary.hpp"
#include "QuICC/ModelOperator/SplitBoundaryValue.hpp"
#include "QuICC/ModelOperator/SplitImplicitLinear.hpp"
#include "QuICC/ModelOperator/Stencil.hpp"
#include "QuICC/ModelOperator/Time.hpp"
#include "QuICC/ModelOperatorBoundary/FieldToRhs.hpp"
#include "QuICC/ModelOperatorBoundary/SolverHasBc.hpp"
#include "QuICC/ModelOperatorBoundary/SolverNoTau.hpp"
#include "QuICC/ModelOperatorBoundary/Stencil.hpp"
#include "QuICC/NonDimensional/CflInertial.hpp"
#include "QuICC/NonDimensional/Ekman.hpp"
#include "QuICC/NonDimensional/MagneticPrandtl.hpp"
#include "QuICC/NonDimensional/Prandtl.hpp"
#include "QuICC/NonDimensional/Rayleigh.hpp"
#include "QuICC/PhysicalNames/Magnetic.hpp"
#include "QuICC/PhysicalNames/Temperature.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"
#include "QuICC/Polynomial/Worland/WorlandBase.hpp"
#include "QuICC/Resolutions/Tools/IndexCounter.hpp"
#include "QuICC/SparseSM/Worland/Boundary/D1.hpp"
#include "QuICC/SparseSM/Worland/Boundary/D2.hpp"
#include "QuICC/SparseSM/Worland/Boundary/Operator.hpp"
#include "QuICC/SparseSM/Worland/Boundary/R1D1DivR1.hpp"
#include "QuICC/SparseSM/Worland/Boundary/Value.hpp"
#include "QuICC/SparseSM/Worland/I2.hpp"
#include "QuICC/SparseSM/Worland/I2Lapl.hpp"
#include "QuICC/SparseSM/Worland/I2Qm.hpp"
#include "QuICC/SparseSM/Worland/I2Qp.hpp"
#include "QuICC/SparseSM/Worland/I4.hpp"
#include "QuICC/SparseSM/Worland/I4Lapl.hpp"
#include "QuICC/SparseSM/Worland/I4Lapl2.hpp"
#include "QuICC/SparseSM/Worland/I4Qm.hpp"
#include "QuICC/SparseSM/Worland/I4Qp.hpp"
#include "QuICC/SparseSM/Worland/Id.hpp"
#include "QuICC/SparseSM/Worland/Stencil/D1.hpp"
#include "QuICC/SparseSM/Worland/Stencil/R1D1DivR1.hpp"
#include "QuICC/SparseSM/Worland/Stencil/Value.hpp"
#include "QuICC/SparseSM/Worland/Stencil/ValueD1.hpp"
#include "QuICC/SparseSM/Worland/Stencil/ValueD2.hpp"
#include "QuICC/Tools/IdToHuman.hpp"
#include "Types/Internal/Math.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Sphere {

namespace Dynamo {

namespace Implicit {

namespace implDetails {

/**
 * @brief Specific options for current model
 */
struct BlockOptionsImpl : public details::BlockOptions
{
   /**
    * @brief default ctor
    */
   BlockOptionsImpl() = default;

   /**
    * @brief default dtor
    */
   virtual ~BlockOptionsImpl() = default;

   /// Jones-Worland alpha
   Internal::MHDFloat a;
   /// Jones-Worland beta
   Internal::MHDFloat b;
   /// Harmonic order m
   int m;
   /// Use truncated quasi-inverse?
   bool truncateQI;
   /// Boundary condition
   std::size_t bcId;
   /// Split operator for influence matrix?
   bool isSplitOperator;
};
} // namespace implDetails

ModelBackend::ModelBackend() :
    IDynamoBackend(),
#ifdef QUICC_TRANSFORM_WORLAND_TRUNCATE_QI
    mcTruncateQI(true)
#else
    mcTruncateQI(false)
#endif // QUICC_TRANSFORM_WORLAND_TRUNCATE_QI
{}

void ModelBackend::enableSplitEquation(const bool flag)
{
   if (flag)
   {
      throw std::logic_error(
         "Split equation for implicit model is not implemented");
   }
   else
   {
      IDynamoBackend::enableSplitEquation(flag);
   }
}

bool ModelBackend::isComplex(const SpectralFieldId& fId) const
{
   return true;
}

ModelBackend::SpectralFieldIds ModelBackend::implicitFields(
   const SpectralFieldId& fId) const
{
   SpectralFieldId velTor = std::make_pair(PhysicalNames::Velocity::id(),
      FieldComponents::Spectral::TOR);
   SpectralFieldId velPol = std::make_pair(PhysicalNames::Velocity::id(),
      FieldComponents::Spectral::POL);
   SpectralFieldId temp = std::make_pair(PhysicalNames::Temperature::id(),
      FieldComponents::Spectral::SCALAR);

   SpectralFieldIds fields;
   if (fId == velTor || fId == velPol || fId == temp)
   {
      fields = {velTor, velPol, temp};
   }
   else
   {
      fields = {fId};
   }

   // Make sure fields are sorted
   std::sort(fields.begin(), fields.end());
   return fields;
}

void ModelBackend::equationInfo(EquationInfo& info, const SpectralFieldId& fId,
   const Resolution& res) const
{
   // Operators are real
   info.isComplex = this->isComplex(fId);

   // Operators are real
   if (fId == std::make_pair(PhysicalNames::Velocity::id(),
                 FieldComponents::Spectral::POL))
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
   info.indexMode =
      static_cast<int>(Equations::CouplingIndexType::SLOWEST_SINGLE_RHS);
}

void ModelBackend::operatorInfo(OperatorInfo& info, const SpectralFieldId& fId,
   const Resolution& res, const Equations::Tools::ICoupling& coupling,
   const BcMap& bcs) const
{
   // Loop overall matrices/eigs
   for (int idx = 0; idx < info.tauN.size(); ++idx)
   {
      auto eigs = coupling.getIndexes(res, idx);

      int tN, gN, rhs;
      ArrayI shift(3);

      this->blockInfo(tN, gN, shift, rhs, fId, res, eigs.at(0), bcs);

      info.tauN(idx) = tN;
      info.galN(idx) = gN;
      info.galShift.row(idx) = shift;
      info.rhsCols(idx) = rhs;

      // Compute system size
      int sN = 0;
      for (auto f: this->implicitFields(fId))
      {
         this->blockInfo(tN, gN, shift, rhs, f, res, eigs.at(0), bcs);
         sN += gN;
      }

      if (sN == 0)
      {
         sN = info.galN(idx);
      }

      info.sysN(idx) = sN;
   }
}

std::vector<details::BlockDescription> ModelBackend::implicitBlockBuilder(
   const SpectralFieldId& rowId, const SpectralFieldId& colId,
   const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs,
   const NonDimensional::NdMap& nds, const bool isSplitOperator) const
{
   std::vector<details::BlockDescription> descr;

   // Create description with common options
   auto getDescription = [&]() -> details::BlockDescription&
   {
      descr.push_back({});
      auto& d = descr.back();
      auto opts = std::make_shared<implDetails::BlockOptionsImpl>();
      opts->a = Polynomial::Worland::WorlandBase::ALPHA_CHEBYSHEV;
      opts->b = Polynomial::Worland::WorlandBase::DBETA_CHEBYSHEV;
      opts->m = eigs.at(0);
      opts->bcId = bcs.find(colId.first)->second;
      opts->truncateQI = this->mcTruncateQI;
      opts->isSplitOperator = isSplitOperator;
      d.opts = opts;

      return d;
   };

   if (rowId == std::make_pair(PhysicalNames::Velocity::id(),
                   FieldComponents::Spectral::TOR))
   {
      if (rowId == colId)
      {
         // Real part of operator
         auto realOp = [](const int nNr, const int nNc, const int l,
                          std::shared_ptr<details::BlockOptions> opts,
                          const NonDimensional::NdMap& nds)
         {
            SparseMatrix bMat(nNr, nNc);

            if (l > 0)
            {
               auto& o =
                  *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(
                     opts);
               const auto Pm = nds.find(NonDimensional::MagneticPrandtl::id())
                                  ->second->value();
               SparseSM::Worland::I2Lapl i2lapl(nNr, nNc, o.a, o.b, l,
                  1 * o.truncateQI);
               bMat = Pm * i2lapl.mat();
            }

            return bMat;
         };

         // Imaginary part of operator
         auto imagOp = [](const int nNr, const int nNc, const int l,
                          std::shared_ptr<details::BlockOptions> opts,
                          const NonDimensional::NdMap& nds)
         {
            SparseMatrix bMat(nNr, nNc);

            if (l > 0)
            {
               auto& o =
                  *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(
                     opts);

               const auto T =
                  1.0 / nds.find(NonDimensional::Ekman::id())->second->value();
               const auto Pm = nds.find(NonDimensional::MagneticPrandtl::id())
                                  ->second->value();
               const auto dl = static_cast<MHDFloat>(l);
               const auto invlapl = 1.0 / (dl * (dl + 1.0));

               SparseSM::Worland::I2 i2(nNr, nNc, o.a, o.b, l,
                  1 * o.truncateQI);
               bMat = o.m * T * Pm * invlapl * i2.mat();
            }

            return bMat;
         };

         // Create block diagonal operator
         auto& d = getDescription();
         d.nRowShift = 0;
         d.nColShift = 0;
         d.realOp = realOp;
         d.imagOp = imagOp;
      }
      else if (colId == std::make_pair(PhysicalNames::Velocity::id(),
                           FieldComponents::Spectral::POL))
      {
         // Real part of first lower diagonal
         auto realOpLower = [](const int nNr, const int nNc, const int l,
                               std::shared_ptr<details::BlockOptions> opts,
                               const NonDimensional::NdMap& nds)
         {
            SparseMatrix bMat(nNr, nNc);

            if (l > 0)
            {
               auto& o =
                  *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(
                     opts);

               auto coriolis = [](const int l, const int m)
               {
                  return (l - MHD_MP(1.0)) * (l + MHD_MP(1.0)) *
                         Internal::Math::sqrt(
                            ((l - m) * (l + m)) /
                            ((MHD_MP(2.0) * l - MHD_MP(1.0)) *
                               (MHD_MP(2.0) * l + MHD_MP(1.0))));
               };

               const auto T =
                  1.0 / nds.find(NonDimensional::Ekman::id())->second->value();
               const auto Pm = nds.find(NonDimensional::MagneticPrandtl::id())
                                  ->second->value();
               const auto dl = static_cast<Internal::MHDFloat>(l);
               const auto invlapl = 1.0 / (dl * (dl + 1.0));

               SparseSM::Worland::I2Qm corQm(nNr, nNc, o.a, o.b, l,
                  1 * o.truncateQI);
               auto norm = coriolis(l, o.m);
               bMat =
                  -static_cast<MHDFloat>(norm * T * Pm * invlapl) * corQm.mat();
            }

            return bMat;
         };

         // Create first lower diagonal operator
         auto& dLow = getDescription();
         dLow.nRowShift = 1;
         dLow.nColShift = 0;
         dLow.realOp = realOpLower;
         dLow.imagOp = nullptr;

         // Real part of first upper diagonal
         auto realOpUpper = [](const int nNr, const int nNc, const int l,
                               std::shared_ptr<details::BlockOptions> opts,
                               const NonDimensional::NdMap& nds)
         {
            SparseMatrix bMat(nNr, nNc);

            if (l > 0)
            {
               auto& o =
                  *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(
                     opts);

               auto coriolis = [](const int l, const int m)
               {
                  return (l - MHD_MP(1.0)) * (l + MHD_MP(1.0)) *
                         Internal::Math::sqrt(
                            ((l - m) * (l + m)) /
                            ((MHD_MP(2.0) * l - MHD_MP(1.0)) *
                               (MHD_MP(2.0) * l + MHD_MP(1.0))));
               };

               const auto T =
                  1.0 / nds.find(NonDimensional::Ekman::id())->second->value();
               const auto Pm = nds.find(NonDimensional::MagneticPrandtl::id())
                                  ->second->value();
               const auto dl = static_cast<Internal::MHDFloat>(l);
               const auto invlapl = MHD_MP(1.0) / (dl * (dl + MHD_MP(1.0)));
               SparseSM::Worland::I2Qp corQp(nNr, nNc, o.a, o.b, l,
                  1 * o.truncateQI);
               auto norm = -coriolis(l + 1, o.m);
               bMat =
                  -static_cast<MHDFloat>(norm * T * Pm * invlapl) * corQp.mat();
            }

            return bMat;
         };

         // Create first upper diagonal operator
         auto& dUp = getDescription();
         dUp.nRowShift = 0;
         dUp.nColShift = 1;
         dUp.realOp = realOpUpper;
         dUp.imagOp = nullptr;
      }
   }
   else if (rowId == std::make_pair(PhysicalNames::Velocity::id(),
                        FieldComponents::Spectral::POL))
   {
      if (rowId == colId)
      {
         // Real part of block
         auto realOp = [](const int nNr, const int nNc, const int l,
                          std::shared_ptr<details::BlockOptions> opts,
                          const NonDimensional::NdMap& nds)
         {
            SparseMatrix bMat(nNr, nNc);

            if (l > 0)
            {
               auto& o =
                  *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(
                     opts);
               const auto Pm = nds.find(NonDimensional::MagneticPrandtl::id())
                                  ->second->value();

               SparseSM::Worland::I4Lapl2 i4lapl2(nNr, nNc, o.a, o.b, l,
                  2 * o.truncateQI);
               bMat = Pm * i4lapl2.mat();
            }

            return bMat;
         };

         // Imaginary part of block
         auto imagOp = [](const int nNr, const int nNc, const int l,
                          std::shared_ptr<details::BlockOptions> opts,
                          const NonDimensional::NdMap& nds)
         {
            SparseMatrix bMat(nNr, nNc);

            if (l > 0)
            {
               auto& o =
                  *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(
                     opts);

               const auto dl = static_cast<Internal::MHDFloat>(l);
               const auto invlapl = MHD_MP(1.0) / (dl * (dl + MHD_MP(1.0)));
               const auto T =
                  1.0 / nds.find(NonDimensional::Ekman::id())->second->value();
               const auto Pm = nds.find(NonDimensional::MagneticPrandtl::id())
                                  ->second->value();
               SparseSM::Worland::I4Lapl coriolis(nNr, nNc, o.a, o.b, l,
                  2 * o.truncateQI);

               // Correct Laplacian for 4th order system according to:
               // McFadden,Murray,Boisvert,
               // Elimination of Spurious Eigenvalues in the
               // Chebyshev Tau Spectral Method,
               // JCP 91, 228-239 (1990)
               // We simply drop the last column
               if (o.bcId == Bc::Name::NoSlip::id())
               {
                  SparseSM::Worland::Id qid(nNr, nNc, o.a, o.b, l, -1);
                  bMat = static_cast<MHDFloat>(o.m * T * Pm * invlapl) *
                         coriolis.mat() * qid.mat();
               }
               else
               {
                  bMat = static_cast<MHDFloat>(o.m * T * Pm * invlapl) *
                         coriolis.mat();
               }
            }

            return bMat;
         };

         // Create diagonal block
         auto& d = getDescription();
         d.nRowShift = 0;
         d.nColShift = 0;
         d.realOp = realOp;
         d.imagOp = imagOp;
      }
      else if (colId == std::make_pair(PhysicalNames::Velocity::id(),
                           FieldComponents::Spectral::TOR))
      {
         // Create real part of block
         auto realOpLower = [](const int nNr, const int nNc, const int l,
                               std::shared_ptr<details::BlockOptions> opts,
                               const NonDimensional::NdMap& nds)
         {
            SparseMatrix bMat(nNr, nNc);

            if (l > 0)
            {
               auto& o =
                  *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(
                     opts);

               auto coriolis = [](const int l, const int m)
               {
                  return (l - MHD_MP(1.0)) * (l + MHD_MP(1.0)) *
                         Internal::Math::sqrt(
                            ((l - m) * (l + m)) /
                            ((MHD_MP(2.0) * l - MHD_MP(1.0)) *
                               (MHD_MP(2.0) * l + MHD_MP(1.0))));
               };

               const auto T =
                  1.0 / nds.find(NonDimensional::Ekman::id())->second->value();
               const auto Pm = nds.find(NonDimensional::MagneticPrandtl::id())
                                  ->second->value();

               const auto dl = static_cast<MHDFloat>(l);
               const auto invlapl = 1.0 / (dl * (dl + 1.0));
               SparseSM::Worland::I4Qm corQm(nNr, nNc, o.a, o.b, l,
                  2 * o.truncateQI);
               auto norm = coriolis(l, o.m);
               bMat =
                  static_cast<MHDFloat>(norm * T * Pm * invlapl) * corQm.mat();
            }

            return bMat;
         };

         // Create first lower diagonal operator
         auto& dLow = getDescription();
         dLow.nRowShift = 1;
         dLow.nColShift = 0;
         dLow.realOp = realOpLower;
         dLow.imagOp = nullptr;

         // Create real part of block
         auto realOpUpper = [](const int nNr, const int nNc, const int l,
                               std::shared_ptr<details::BlockOptions> opts,
                               const NonDimensional::NdMap& nds)
         {
            SparseMatrix bMat(nNr, nNc);

            if (l > 0)
            {
               auto& o =
                  *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(
                     opts);

               auto coriolis = [](const int l, const int m)
               {
                  return (l - MHD_MP(1.0)) * (l + MHD_MP(1.0)) *
                         Internal::Math::sqrt(
                            ((l - m) * (l + m)) /
                            ((MHD_MP(2.0) * l - MHD_MP(1.0)) *
                               (MHD_MP(2.0) * l + MHD_MP(1.0))));
               };

               const auto T =
                  1.0 / nds.find(NonDimensional::Ekman::id())->second->value();
               const auto Pm = nds.find(NonDimensional::MagneticPrandtl::id())
                                  ->second->value();

               const auto dl = static_cast<MHDFloat>(l);
               const auto invlapl = 1.0 / (dl * (dl + 1.0));
               SparseSM::Worland::I4Qp corQp(nNr, nNc, o.a, o.b, l,
                  2 * o.truncateQI);
               auto norm = -coriolis(l + 1, o.m);
               bMat =
                  static_cast<MHDFloat>(norm * T * Pm * invlapl) * corQp.mat();
            }

            return bMat;
         };

         // Create first upper diagonal operator
         auto& dUp = getDescription();
         dUp.nRowShift = 0;
         dUp.nColShift = 1;
         dUp.realOp = realOpUpper;
         dUp.imagOp = nullptr;
      }
   }
   else if (rowId == std::make_pair(PhysicalNames::Magnetic::id(),
                        FieldComponents::Spectral::TOR))
   {
      if (rowId == colId)
      {
         // Creat real part of block
         auto realOp = [](const int nNr, const int nNc, const int l,
                          std::shared_ptr<details::BlockOptions> opts,
                          const NonDimensional::NdMap& nds)
         {
            auto& o =
               *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(opts);

            SparseSM::Worland::I2Lapl i2lapl(nNr, nNc, o.a, o.b, l,
               1 * o.truncateQI);
            SparseMatrix bMat = i2lapl.mat();

            return bMat;
         };

         // Create diagonal block
         auto& d = getDescription();
         d.nRowShift = 0;
         d.nColShift = 0;
         d.realOp = realOp;
         d.imagOp = nullptr;
      }
   }
   else if (rowId == std::make_pair(PhysicalNames::Magnetic::id(),
                        FieldComponents::Spectral::POL))
   {
      if (rowId == colId)
      {
         // Creat real part of block
         auto realOp = [](const int nNr, const int nNc, const int l,
                          std::shared_ptr<details::BlockOptions> opts,
                          const NonDimensional::NdMap& nds)
         {
            auto& o =
               *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(opts);

            SparseSM::Worland::I2Lapl i2lapl(nNr, nNc, o.a, o.b, l,
               1 * o.truncateQI);
            SparseMatrix bMat = i2lapl.mat();

            return bMat;
         };

         // Create diagonal block
         auto& d = getDescription();
         d.nRowShift = 0;
         d.nColShift = 0;
         d.realOp = realOp;
         d.imagOp = nullptr;
      }
   }
   else if (rowId == std::make_pair(PhysicalNames::Temperature::id(),
                        FieldComponents::Spectral::SCALAR))
   {
      if (rowId == colId)
      {
         // Creat real part of block
         auto realOp = [](const int nNr, const int nNc, const int l,
                          std::shared_ptr<details::BlockOptions> opts,
                          const NonDimensional::NdMap& nds)
         {
            auto& o =
               *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(opts);

            const auto Pr =
               nds.find(NonDimensional::Prandtl::id())->second->value();
            const auto Pm =
               nds.find(NonDimensional::MagneticPrandtl::id())->second->value();

            SparseSM::Worland::I2Lapl i2lapl(nNr, nNc, o.a, o.b, l,
               1 * o.truncateQI);
            SparseMatrix bMat = (Pm / Pr) * i2lapl.mat();

            return bMat;
         };

         // Create diagonal block
         auto& d = getDescription();
         d.nRowShift = 0;
         d.nColShift = 0;
         d.realOp = realOp;
         d.imagOp = nullptr;
      }
   }
   else
   {
      throw std::logic_error("Equations are not setup properly");
   }

   return descr;
}

std::vector<details::BlockDescription> ModelBackend::timeBlockBuilder(
   const SpectralFieldId& rowId, const SpectralFieldId& colId,
   const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs,
   const NonDimensional::NdMap& nds) const
{
   assert(rowId == colId);
   auto fieldId = rowId;

   std::vector<details::BlockDescription> descr;

   // Create description with common options
   auto getDescription = [&]() -> details::BlockDescription&
   {
      descr.push_back({});
      auto& d = descr.back();
      auto opts = std::make_shared<implDetails::BlockOptionsImpl>();
      opts->a = Polynomial::Worland::WorlandBase::ALPHA_CHEBYSHEV;
      opts->b = Polynomial::Worland::WorlandBase::DBETA_CHEBYSHEV;
      opts->m = eigs.at(0);
      opts->bcId = bcs.find(colId.first)->second;
      opts->truncateQI = this->mcTruncateQI;
      opts->isSplitOperator = false;
      d.opts = opts;

      return d;
   };

   if (fieldId == std::make_pair(PhysicalNames::Velocity::id(),
                     FieldComponents::Spectral::TOR))
   {
      // Real part of operator
      auto realOp = [](const int nNr, const int nNc, const int l,
                       std::shared_ptr<details::BlockOptions> opts,
                       const NonDimensional::NdMap& nds)
      {
         assert(nNr == nNc);

         SparseMatrix bMat;
         auto& o =
            *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(opts);

         if (l > 0)
         {
            SparseSM::Worland::I2 spasm(nNr, nNc, o.a, o.b, l,
               1 * o.truncateQI);
            bMat = spasm.mat();
         }
         else
         {
            SparseSM::Worland::Id qid(nNr, nNc, o.a, o.b, l);
            bMat = qid.mat();
         }

         return bMat;
      };

      // Create block diagonal operator
      auto& d = getDescription();
      d.nRowShift = 0;
      d.nColShift = 0;
      d.realOp = realOp;
      d.imagOp = nullptr;
   }
   else if (fieldId == std::make_pair(PhysicalNames::Velocity::id(),
                          FieldComponents::Spectral::POL))
   {
      // Real part of operator
      auto realOp = [](const int nNr, const int nNc, const int l,
                       std::shared_ptr<details::BlockOptions> opts,
                       const NonDimensional::NdMap& nds)
      {
         assert(nNr == nNc);

         SparseMatrix bMat;
         auto& o =
            *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(opts);

         if (l > 0)
         {
            SparseSM::Worland::I4Lapl spasm(nNr, nNc, o.a, o.b, l,
               2 * o.truncateQI);

            // Correct Laplacian for 4th order system according to:
            // McFadden,Murray,Boisvert,
            // Elimination of Spurious Eigenvalues in the
            // Chebyshev Tau Spectral Method,
            // JCP 91, 228-239 (1990)
            // We simply drop the last column
            if (o.bcId == Bc::Name::NoSlip::id())
            {
               SparseSM::Worland::Id qid(nNr, nNc, o.a, o.b, l, -1);
               bMat = spasm.mat() * qid.mat();
            }
            else
            {
               bMat = spasm.mat();
            }
         }
         else
         {
            SparseSM::Worland::Id qid(nNr, nNc, o.a, o.b, l);
            bMat = qid.mat();
         }

         return bMat;
      };

      // Create block diagonal operator
      auto& d = getDescription();
      d.nRowShift = 0;
      d.nColShift = 0;
      d.realOp = realOp;
      d.imagOp = nullptr;
   }
   else if (fieldId == std::make_pair(PhysicalNames::Magnetic::id(),
                          FieldComponents::Spectral::TOR))
   {
      // Real part of operator
      auto realOp = [](const int nNr, const int nNc, const int l,
                       std::shared_ptr<details::BlockOptions> opts,
                       const NonDimensional::NdMap& nds)
      {
         assert(nNr == nNc);

         SparseMatrix bMat;
         auto& o =
            *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(opts);

         if (l > 0)
         {
            SparseSM::Worland::I2 spasm(nNr, nNc, o.a, o.b, l,
               1 * o.truncateQI);
            bMat = spasm.mat();
         }
         else
         {
            SparseSM::Worland::Id qid(nNr, nNc, o.a, o.b, l);
            bMat = qid.mat();
         }

         return bMat;
      };

      // Create block diagonal operator
      auto& d = getDescription();
      d.nRowShift = 0;
      d.nColShift = 0;
      d.realOp = realOp;
      d.imagOp = nullptr;
   }
   else if (fieldId == std::make_pair(PhysicalNames::Magnetic::id(),
                          FieldComponents::Spectral::POL))
   {
      // Real part of operator
      auto realOp = [](const int nNr, const int nNc, const int l,
                       std::shared_ptr<details::BlockOptions> opts,
                       const NonDimensional::NdMap& nds)
      {
         assert(nNr == nNc);

         SparseMatrix bMat;
         auto& o =
            *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(opts);

         if (l > 0)
         {
            SparseSM::Worland::I2 spasm(nNr, nNc, o.a, o.b, l,
               1 * o.truncateQI);
            bMat = spasm.mat();
         }
         else
         {
            SparseSM::Worland::Id qid(nNr, nNc, o.a, o.b, l);
            bMat = qid.mat();
         }

         return bMat;
      };

      // Create block diagonal operator
      auto& d = getDescription();
      d.nRowShift = 0;
      d.nColShift = 0;
      d.realOp = realOp;
      d.imagOp = nullptr;
   }
   else if (fieldId == std::make_pair(PhysicalNames::Temperature::id(),
                          FieldComponents::Spectral::SCALAR))
   {
      // Real part of operator
      auto realOp = [](const int nNr, const int nNc, const int l,
                       std::shared_ptr<details::BlockOptions> opts,
                       const NonDimensional::NdMap& nds)
      {
         auto& o =
            *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(opts);

         SparseSM::Worland::I2 spasm(nNr, nNc, o.a, o.b, l, 1 * o.truncateQI);
         SparseMatrix bMat = spasm.mat();

         return bMat;
      };

      // Create block diagonal operator
      auto& d = getDescription();
      d.nRowShift = 0;
      d.nColShift = 0;
      d.realOp = realOp;
      d.imagOp = nullptr;
   }

   return descr;
}

void ModelBackend::modelMatrix(DecoupledZSparse& rModelMatrix,
   const std::size_t opId,
   const Equations::CouplingInformation::FieldId_range imRange,
   const int matIdx, const std::size_t bcType, const Resolution& res,
   const std::vector<MHDFloat>& eigs, const BcMap& bcs,
   const NonDimensional::NdMap& nds) const
{
   assert(eigs.size() == 1);
   int m = eigs.at(0);
   auto maxL = res.counter().dim(Dimensions::Simulation::SIM2D,
                  Dimensions::Space::SPECTRAL, m) -
               1;

   // Time operator
   if (opId == ModelOperator::Time::id())
   {
      for (auto pRowId = imRange.first; pRowId != imRange.second; pRowId++)
      {
         auto rowId = *pRowId;
         auto colId = rowId;
         const auto& fields = this->implicitFields(rowId);
         auto descr = timeBlockBuilder(rowId, colId, res, eigs, bcs, nds);
         buildBlock(rModelMatrix, descr, rowId, colId, fields, matIdx, bcType,
            res, m, maxL, bcs, nds, false);
      }
   }
   // Linear operator
   else if (opId == ModelOperator::ImplicitLinear::id() ||
            opId == ModelOperator::SplitImplicitLinear::id())
   {
      bool isSplit = (opId == ModelOperator::SplitImplicitLinear::id());

      for (auto pRowId = imRange.first; pRowId != imRange.second; pRowId++)
      {
         auto rowId = *pRowId;
         const auto& fields = this->implicitFields(rowId);
         for (auto pColId = imRange.first; pColId != imRange.second; pColId++)
         {
            auto colId = *pColId;
            auto descr =
               implicitBlockBuilder(rowId, colId, res, eigs, bcs, nds, isSplit);
            buildBlock(rModelMatrix, descr, rowId, colId, fields, matIdx,
               bcType, res, m, maxL, bcs, nds, isSplit);
         }
      }
   }
   // Boundary operator
   else if (opId == ModelOperator::Boundary::id() ||
            opId == ModelOperator::SplitBoundary::id())
   {
      bool isSplit = (opId == ModelOperator::SplitBoundary::id());

      for (auto pRowId = imRange.first; pRowId != imRange.second; pRowId++)
      {
         auto rowId = *pRowId;
         const auto& fields = this->implicitFields(rowId);
         for (auto pColId = imRange.first; pColId != imRange.second; pColId++)
         {
            auto colId = *pColId;
            auto descr =
               boundaryBlockBuilder(rowId, colId, res, eigs, bcs, nds, isSplit);
            buildBlock(rModelMatrix, descr, rowId, colId, fields, matIdx,
               bcType, res, m, maxL, bcs, nds, isSplit);
         }
      }
   }
   else
   {
      throw std::logic_error("Requested operator type is not implemented");
   }
}

std::vector<details::BlockDescription> ModelBackend::boundaryBlockBuilder(
   const SpectralFieldId& rowId, const SpectralFieldId& colId,
   const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs,
   const NonDimensional::NdMap& nds, const bool isSplit) const
{
   std::vector<details::BlockDescription> descr;

   // Create description with common options
   auto getDescription = [&]() -> details::BlockDescription&
   {
      descr.push_back({});
      auto& d = descr.back();
      auto opts = std::make_shared<implDetails::BlockOptionsImpl>();
      opts->a = Polynomial::Worland::WorlandBase::ALPHA_CHEBYSHEV;
      opts->b = Polynomial::Worland::WorlandBase::DBETA_CHEBYSHEV;
      opts->m = eigs.at(0);
      opts->bcId = bcs.find(colId.first)->second;
      opts->truncateQI = this->mcTruncateQI;
      opts->isSplitOperator = isSplit;
      d.opts = opts;

      return d;
   };

   if (rowId == colId)
   {
      // Real part of operator
      auto realOp = [](const int nNr, const int nNc, const int l,
                       std::shared_ptr<details::BlockOptions> opts,
                       const NonDimensional::NdMap& nds)
      {
         SparseMatrix bMat(nNr, nNc);

         return bMat;
      };

      // Create block diagonal operator
      auto& d = getDescription();
      d.nRowShift = 0;
      d.nColShift = 0;
      d.realOp = realOp;
      d.imagOp = nullptr;
   }
   else
   {
      // Create block diagonal operator
      auto& d = getDescription();
      d.nRowShift = 0;
      d.nColShift = 0;
      d.realOp = nullptr;
      d.imagOp = nullptr;
   }

   return descr;
}

void ModelBackend::galerkinStencil(SparseMatrix& mat,
   const SpectralFieldId& fieldId, const int matIdx, const Resolution& res,
   const std::vector<MHDFloat>& eigs, const bool makeSquare, const BcMap& bcs,
   const NonDimensional::NdMap& nds) const
{
   assert(this->useGalerkin());
   assert(eigs.size() == 1);
   int m = eigs.at(0);
   auto maxL = res.counter().dim(Dimensions::Simulation::SIM2D,
                  Dimensions::Space::SPECTRAL, m) -
               1;

   // Compute system size
   const auto& fields = this->implicitFields(fieldId);
   const auto sysRows =
      systemInfo(fieldId, fieldId, fields, m, maxL, res, bcs, makeSquare, false)
         .blockRows;
   const auto sysCols =
      systemInfo(fieldId, fieldId, fields, m, maxL, res, bcs, true, false)
         .blockCols;

   auto nL = res.counter().dim(Dimensions::Simulation::SIM2D,
      Dimensions::Space::SPECTRAL, m);

   if (mat.size() == 0)
   {
      mat.resize(sysRows, sysCols);
   }
   assert(mat.rows() == sysRows);
   assert(mat.cols() == sysCols);

   assert(eigs.size() == 1);

   int rowShift = 0;
   int colShift = 0;
   for (int l = m; l < nL; l++)
   {
      SparseMatrix S;
      this->stencil(S, fieldId, l, res, makeSquare, bcs, nds);
      this->addBlock(mat, S, rowShift, colShift);

      rowShift += S.rows();
      colShift += S.cols();
   }
}

void ModelBackend::explicitBlock(DecoupledZSparse& mat,
   const SpectralFieldId& fId, const std::size_t opId,
   const SpectralFieldId fieldId, const int matIdx, const Resolution& res,
   const std::vector<MHDFloat>& eigs, const BcMap& bcs,
   const NonDimensional::NdMap& nds) const
{
   // Explicit linear operator
   if (opId == ModelOperator::ExplicitLinear::id())
   {
      // Nothing to be done
      throw std::logic_error("There are no explicit linear operators");
   }
   // Explicit nonlinear operator
   else if (opId == ModelOperator::ExplicitNonlinear::id())
   {
      throw std::logic_error("There are no explicit nonlinear operators");
   }
   // Explicit nextstep operator
   else if (opId == ModelOperator::ExplicitNextstep::id())
   {
      throw std::logic_error("There are no explicit nextstep operators");
   }
}

} // namespace Implicit
} // namespace Dynamo
} // namespace Sphere
} // namespace Boussinesq
} // namespace Model
} // namespace QuICC
