/**
 * @file IExponentialBackend.hpp
 * @brief Base model backend for dynamo model
 */

#ifndef QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_EXPONENTIAL_IEXPONENTIALBACKEND_HPP
#define QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_EXPONENTIAL_IEXPONENTIALBACKEND_HPP

// System includes
//
#include <map>
#include <memory>
#include <string>
#include <vector>

// Project includes
//
#include "Model/Boussinesq/Sphere/Dynamo/IDynamoBackend.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Sphere {

namespace Dynamo {

namespace Exponential {

/**
 * @brief Base model backend for Dynamo model
 */
class IExponentialBackend : public IDynamoBackend
{
public:
   /**
    * @brief Constructor
    */
   IExponentialBackend() = default;

   /**
    * @brief Destructor
    */
   virtual ~IExponentialBackend() = default;

   /**
    * @brief Get vector of names for the physical fields
    */
   virtual std::vector<std::string> fieldNames() const override;

protected:
   /**
    * @brief Number of boundary conditions
    *
    * @fId  Field ID for which to get number of BC
    */
   int nBc(const SpectralFieldId& fId) const override;

   /**
    * @brief Apply tau line for boundary condition
    *
    * @param mat     Input/Output matrix to apply tau line to
    * @param rowId   ID of field of equation
    * @param colId   ID of field
    * @param l       Harmonic degree
    * @param opts    Options
    * @param nN      1D dimension
    * @param bcs     Boundary conditions
    * @param nds     Nondimensional parameters
    * @param isSplitOperator  Is second operator of split 4th order system?
    */
   void applyTau(SparseMatrix& mat, const SpectralFieldId& rowId,
      const SpectralFieldId& colId, const int l,
      std::shared_ptr<details::BlockOptions> opts, const int nN,
      const BcMap& bcs, const NonDimensional::NdMap& nds,
      const bool isSplitOperator) const override;

   /**
    * @brief Boundary condition stencil
    *
    * @param mat        Input/Output matrix to store galerkin stencil
    * @param fID        Field ID
    * @param l          Harmonic degree
    * @param nN         1D dimension
    * @param makeSquare Truncate operator to make square
    * @param bcs        Boundary conditions
    * @param nds        Nondimensional parameters
    */
   void stencil(SparseMatrix& mat, const SpectralFieldId& fId, const int l,
      const int nN, const bool makeSquare, const BcMap& bcs,
      const NonDimensional::NdMap& nds) const;

private:
};

} // namespace Exponential
} // namespace Dynamo
} // namespace Sphere
} // namespace Boussinesq
} // namespace Model
} // namespace QuICC

#endif // QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_EXPONENTIAL_IEXPONENTIALBACKEND_HPP
