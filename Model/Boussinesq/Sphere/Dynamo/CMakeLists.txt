quicc_target_sources(${QUICC_CURRENT_MODEL_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
  IDynamoModel.cpp
  IDynamoBackend.cpp
  Induction.cpp
  InductionKernel.cpp
  Momentum.cpp
  MomentumKernel.cpp
  Transport.cpp
  TransportKernel.cpp
  )

add_subdirectory(Explicit)
add_subdirectory(Implicit)
