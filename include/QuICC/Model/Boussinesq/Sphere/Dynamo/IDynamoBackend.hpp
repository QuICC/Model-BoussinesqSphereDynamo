/**
 * @file IDynamoBackend.hpp
 * @brief Base model backend for dynamo model
 */

#ifndef QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_IDYNAMOBACKEND_HPP
#define QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_IDYNAMOBACKEND_HPP

// System includes
//
#include <string>
#include <vector>
#include <map>
#include <memory>

// Project includes
//
#include "QuICC/Model/IModelBackend.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Sphere {

namespace Dynamo {

   /**
    * @brief Base model backend for Dynamo model
    */
   class IDynamoBackend: public IModelBackend
   {
      public:
         /**
          * @brief Constructor
          */
         IDynamoBackend();

         /**
          * @brief Destructor
          */
         virtual ~IDynamoBackend() = default;

         /**
          * @brief Get vector of names for the physical fields
          */
         virtual std::vector<std::string> fieldNames() const override;

         /**
          * @brief Get vector of names for the nondimensional parameters
          */
         virtual std::vector<std::string> paramNames() const override;

         /**
          * @brief Get vector of bools about periodic box
          */
         virtual std::vector<bool> isPeriodicBox() const override;

         /**
          * @brief Enable galerkin basis
          */
         virtual void enableGalerkin(const bool flag) override;

         /**
          * @brief Get auotmatically computed parameters based on input parameters
          *
          * @param cfg  Input parameters
          */
         virtual std::map<std::string,MHDFloat> automaticParameters(const std::map<std::string,MHDFloat>& cfg) const override;

      protected:
         /**
          * @brief Use Galerkin basis?
          */
         bool useGalerkin() const;

         /**
          * @brief Split high order equations?
          */
         bool useSplitEquation() const;

         /**
          * @brief Enable split equation
          */
         void enableSplitEquation(const bool flag);

         /**
          * @brief Number of boundary conditions
          */
         int nBc(const SpectralFieldId& fId) const;

         /**
          * @brief Get operator block information
          */
         void blockInfo(int& tN, int& gN, ArrayI& shift, int& rhs, const SpectralFieldId& fId, const Resolution& res, const MHDFloat l, const BcMap& bcs) const;

         /**
          * @brief Apply tau line for boundary condition
          */
         void applyTau(SparseMatrix& mat, const SpectralFieldId& rowId, const SpectralFieldId& colId, const int l, const Resolution& res, const BcMap& bcs, const NonDimensional::NdMap& nds, const bool isSplitOperator) const;

         /**
          * @brief Boundary condition stencil
          */
         virtual void stencil(SparseMatrix& mat, const SpectralFieldId& fId, const int l, const Resolution& res, const bool makeSquare, const BcMap& bcs, const NonDimensional::NdMap& nds) const;

         /**
          * @brief Apply galerkin stencil for boundary condition
          */
         void applyGalerkinStencil(SparseMatrix& decMat, const SpectralFieldId& rowId, const SpectralFieldId& colId, const int l, const Resolution& res, const BcMap& bcs, const NonDimensional::NdMap& nds) const;

      private:
         /**
          * @brief Use Galerkin basis?
          */
         bool mUseGalerkin;

         /**
          * @brief Split high order equation?
          */
         bool mUseSplitEquation;
   };

} // Dynamo
} // Sphere
} // Boussinesq
} // Model
} // QuICC

#endif // QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_IDYNAMOBACKEND_HPP
