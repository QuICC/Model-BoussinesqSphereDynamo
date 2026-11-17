target_sources(${QUICC_CURRENT_MODEL_LIB}_exponential ${QUICC_CMAKE_SRC_VISIBILITY}
  IExponentialBackend.cpp
  IExponentialModel.cpp
  InductionJacobian.cpp
  InductionJacobianKernel.cpp
  ModelBackend.cpp
  MomentumJacobian.cpp
  MomentumJacobianKernel.cpp
  TransportJacobian.cpp
  TransportJacobianKernel.cpp
  )
