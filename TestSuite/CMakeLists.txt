include(ConfigureBenchmarks)
quicc_add_benchmark(Explicit
  MODEL "BoussinesqSphereDynamo"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  )

quicc_add_benchmark(Explicit
  MODEL "BoussinesqSphereDynamo"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  VARIANTS "framework/parallel/algorithm:serial" "framework/timestepping/scheme:ImEx_PC2" "setup/model/split_equation:On"
  FILTER "algorithm" "split_equation"
  )

quicc_add_benchmark(Implicit
    MODEL "BoussinesqSphereDynamo"
    WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
    ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
    TIMEOUT 500
  )

if(QUICC_USE_MPI)
  quicc_add_benchmark(Implicit
    MODEL "BoussinesqSphereDynamo"
    WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
    ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
    VARIANTS "framework/parallel/algorithm:single2d"
    TIMEOUT 500
  )
endif()
