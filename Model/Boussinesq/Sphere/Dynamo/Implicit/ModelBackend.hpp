/**
 * @file ModelBackend.hpp
 * @brief Model backend
 */

#ifndef QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_IMPLICIT_MODELBACKEND_HPP
#define QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_IMPLICIT_MODELBACKEND_HPP

// System includes
//
#include <string>
#include <vector>
#include <map>
#include <memory>

// Project includes
//
#include "Model/Boussinesq/Sphere/Dynamo/IDynamoBackend.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Sphere {

namespace Dynamo {

namespace Implicit {

   namespace internal {
      struct SystemInfo
      {
         int systemSize;
         int blockRows;
         int blockCols;
         int startRow;
         int startCol;

         SystemInfo(const int size, const int rows, const int cols, const int row, const int col)
            : systemSize(size), blockRows(rows), blockCols(cols), startRow(row), startCol(col)
         {
         };
      };
   }

   /**
    * @brief Interface for model backend
    */
   class ModelBackend: public IDynamoBackend
   {
      public:
         /**
          * @brief Constructor
          */
         ModelBackend();

         /**
          * @brief Destructor
          */
         virtual ~ModelBackend() = default;

         /**
          * @brief Enable splitting 4th equation into two 2nd order?
          *
          * @param flag True/False to enable option
          */
         virtual void enableSplitEquation(const bool flag) override;

         /**
          * @brief Get equation information
          *
          * @param info Equation information
          * @param fId  Field ID
          * @param res  Resolution object
          */
         virtual void equationInfo(EquationInfo& info, const SpectralFieldId& fId, const Resolution& res) const override;

         /**
          * @brief Get operator information
          *
          * @param info       Equation information
          * @param fId        Field ID
          * @param res        Resolution object
          * @param coupling   Equation/Field coupling information
          * @param bcs        Boundary conditions
          */
         virtual void operatorInfo(OperatorInfo& info, const SpectralFieldId& fId, const Resolution& res, const Equations::Tools::ICoupling& coupling, const BcMap& bcs) const override;

         /**
          * @brief Build model matrix
          *
          * @param rModelMatrix  Input/Output matrix to fill with operators
          * @param opId          Type of model matrix
          * @param imRange       Coupled fields
          * @param matIdx        Matrix index
          * @param bcType        Boundary condition scheme (Tau vs Galerkin)
          * @param res           Resolution object
          * @param eigs          Indexes of other dimensions
          * @param bcs           Boundary conditions
          * @param nds           Nondimensional parameters
          */
         virtual void modelMatrix(DecoupledZSparse& rModelMatrix, const std::size_t opId, const Equations::CouplingInformation::FieldId_range imRange, const int matIdx, const std::size_t bcType, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds) const override;

         /**
          * @brief Build galerkin stencil
          *
          * @param mat     Input/Output matrix to fill with stencil
          * @param fId     Field ID
          * @param matIdx  Matrix index
          * @param res     Resolution object
          * @param eigs          Indexes of other dimensions
          * @param makeSquare Truncate stencil to obtain square matrix?
          * @param bcs           Boundary conditions
          * @param nds           Nondimensional parameters
          */
         virtual void galerkinStencil(SparseMatrix& mat, const SpectralFieldId& fId, const int matIdx, const Resolution& res, const std::vector<MHDFloat>& eigs, const bool makeSquare, const BcMap& bcs, const NonDimensional::NdMap& nds) const override;

         /**
          * @brief Build explicit block
          *
          * @param rModelMatrix  Input/Output matrix to fill with operators
          * @param fId           Equation field ID
          * @param opId          Type of explicit operator
          * @param fieldId       Coupled field ID
          * @param matIdx        Matrix index
          * @param res           Resolution object
          * @param eigs          Indexes of other dimensions
          * @param bcs           Boundary conditions
          * @param nds           Nondimensional parameters
          */
         virtual void explicitBlock(DecoupledZSparse& mat, const SpectralFieldId& fId, const std::size_t opId,  const SpectralFieldId fieldId, const int matIdx, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds) const override;

      protected:
         /**
          * @brief Build model matrix
          *
          * @param fId  Field ID
          */
         SpectralFieldIds implicitFields(const SpectralFieldId& fId) const;

         /**
          * @brief Build implicit matrix block
          *
          * @param decMat  Ouput matrix
          * @param rowId   Field ID of block matrix row
          * @param colId   Field ID of block matrix column
          * @param matIdx  Matrix ID
          * @param res     Resolution object
          * @param eigs    Slow indexes
          * @param nds     Nondimension parameters
          * @param isSplitOperator  Set operator of split system
          */
         void implicitBlock(DecoupledZSparse& decMat, const SpectralFieldId& rowId, const SpectralFieldId& colId, const int matIdx, const std::size_t bcType, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds, const bool isSplitOperator) const;

         /**
          * @brief Build time matrix block
          *
          * @param decMat  Input/Output matrix to fill with operators
          * @param fieldId   ID of field 
          * @param matIdx        Matrix index
          * @param bcType        Boundary condition scheme (Tau vs Galerkin)
          * @param res           Resolution object
          * @param eigs          Indexes of other dimensions
          * @param bcs           Boundary conditions
          * @param nds           Nondimensional parameters
          */
         void timeBlock(DecoupledZSparse& decMat, const SpectralFieldId& fieldId, const int matIdx, const std::size_t bcType, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds) const;

         /**
          * @brief Build boundary matrix block
          *
          * @param decMat  Input/Output matrix to fill with operators
          * @param rowId   ID of field of equation
          * @param colId   ID of field 
          * @param matIdx        Matrix index
          * @param bcType        Boundary condition scheme (Tau vs Galerkin)
          * @param res           Resolution object
          * @param eigs          Indexes of other dimensions
          * @param bcs           Boundary conditions
          * @param nds           Nondimensional parameters
          * @param isSplitEquation Is second operator of split 4th order equation
          */
         void boundaryBlock(DecoupledZSparse& decMat, const SpectralFieldId& rowId, const SpectralFieldId& colId, const int matIdx, const std::size_t bcType, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds, const bool isSplitEquation) const;

      private:
         /**
          * @brief Add block matrix to full system matrix
          *
          * @param mat Input/Output matrix to add matrix block to
          * @param block matrix block to add
          * @param rowShift   Start row of matrix block
          * @param colShift   Start column of matrix block
          * @param coeff      Scaling coefficient of matrix block
          */
         void addBlock(SparseMatrix& mat, const SparseMatrix& block, const int rowShift, const int colShift, const MHDFloat coeff = 1.0) const;

         /**
          * @brief Get operator information
          *
          * @param fId  Field ID
          * @param m    Harmonic order m
          * @param res  Resolution object
          * @param bcs  Boundary conditions
          * @param isGalerkin Use Galerkin scheme?
          */
         int blockSize(const SpectralFieldId& fId, const int m, const Resolution& res, const BcMap& bcs, const bool isGalerkin) const;

         /**
          * @brief Get operator block shape
          * 
          * @param rowId   Equation Field ID
          * @param colId    Field Id
          * @param m       Harmonic order m
          * @param res  Resolution object
          * @param bcs  Boundary conditions
          * @param isGalerkin Use Galerkin scheme?
          * @param dropRows?  Drop Tau line rows
          */
         std::pair<int,int> blockShape(const SpectralFieldId& rowId, const SpectralFieldId& colId, const int m, const Resolution& res, const BcMap& bcs, const bool isGalerkin, const bool dropRows) const;

         /**
          * @brief Compute size information of full system
          *
          * @param rowId   Equation Field ID
          * @param colId    Field Id
          * @param m       Harmonic order m
          * @param res  Resolution object
          * @param bcs  Boundary conditions
          * @param isGalerkin Use Galerkin scheme?
          * @param dropRows?  Drop Tau line rows
          */
         internal::SystemInfo systemInfo(const SpectralFieldId& colId, const SpectralFieldId& rowId, const int m, const Resolution& res, const BcMap& bcs, const bool isGalerkin, const bool dropRows) const;

         /**
          * @brief Truncate quasi-inverse operators?
          */
         const bool mcTruncateQI;
   };

} // Implicit
} // Dynamo
} // Sphere
} // Boussines
} // Model
} // QuICC

#endif // QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_IMPLICIT_MODELBACKEND_HPP
