include(ConfigureBenchmarks)
quicc_add_benchmark(Explicit
  MODEL "BoussinesqSphereDynamo"
  TYPE "Config"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  )

quicc_add_benchmark(Explicit
  MODEL "BoussinesqSphereDynamo"
  TYPE "State"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  )

quicc_add_benchmark(Explicit
  MODEL "BoussinesqSphereDynamo"
  TYPE "Visualization"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  )

if(NOT QUICC_USE_MPI)
  quicc_add_benchmark(Explicit
    MODEL "BoussinesqSphereDynamo"
    TYPE "Model"
    WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
    ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
    VARIANTS "framework/parallel/algorithm:serial" "framework/timestepping/scheme:ImEx_RK_CB2"
    FILTER "algorithm" "scheme"
    DATAFILTER "timestepping/scheme"
    )

  quicc_add_benchmark(Explicit
    MODEL "BoussinesqSphereDynamo"
    TYPE "Model"
    WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
    ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
    VARIANTS "framework/parallel/algorithm:serial" "framework/timestepping/scheme:ImEx_RK_CB2" "setup/boundary/scheme:Galerkin"
    FILTER "algorithm" "scheme"
    VALIDATIONFILTER "boundary/scheme"
    DATAFILTER "timestepping/scheme"
    )

  quicc_add_benchmark(Explicit
    MODEL "BoussinesqSphereDynamo"
    TYPE "Model"
    WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
    ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
    VARIANTS "framework/parallel/algorithm:serial" "framework/timestepping/scheme:ImEx_PC2"
    FILTER "algorithm" "scheme"
    DATAFILTER "timestepping/scheme"
    )

  quicc_add_benchmark(Explicit
    MODEL "BoussinesqSphereDynamo"
    TYPE "Model"
    WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
    ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
    VARIANTS "framework/parallel/algorithm:serial" "framework/timestepping/scheme:ImEx_PC2" "setup/boundary/scheme:Galerkin"
    FILTER "algorithm" "scheme"
    DATAFILTER "timestepping/scheme"
    VALIDATIONFILTER "boundary/scheme"
    )

  quicc_add_benchmark(Explicit
    MODEL "BoussinesqSphereDynamo"
    TYPE "Model"
    WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
    ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
    VARIANTS "framework/parallel/algorithm:serial" "framework/timestepping/scheme:ImEx_PC2" "setup/model/split_equation:On"
    FILTER "algorithm" "timestepping/scheme" "split_equation"
    DATAFILTER "timestepping/scheme"
    VALIDATIONFILTER "split_equation"
    )
else()
  quicc_add_benchmark(Explicit
    MODEL "BoussinesqSphereDynamo"
    TYPE "Model"
    WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
    ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
    VARIANTS "framework/parallel/algorithm:tubular" "framework/timestepping/scheme:ImEx_PC2"
    FILTER "algorithm" "scheme"
    DATAFILTER "timestepping/scheme"
    MPIRANKS 6
    )

  quicc_add_benchmark(Explicit
    MODEL "BoussinesqSphereDynamo"
    TYPE "Model"
    WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
    ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
    VARIANTS "framework/parallel/algorithm:tubular" "framework/timestepping/scheme:ImEx_PC2" "setup/boundary/scheme:Galerkin"
    FILTER "algorithm" "scheme"
    DATAFILTER "timestepping/scheme"
    VALIDATIONFILTER "boundary/scheme"
    MPIRANKS 6
    )

  quicc_add_benchmark(Explicit
    MODEL "BoussinesqSphereDynamo"
    TYPE "Model"
    WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
    ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
    VARIANTS "framework/parallel/algorithm:tubular" "framework/timestepping/scheme:ImEx_PC2" "setup/model/split_equation:On"
    FILTER "algorithm" "timestepping/scheme" "split_equation"
    DATAFILTER "timestepping/scheme"
    VALIDATIONFILTER "split_equation"
    MPIRANKS 6
    )
endif()

quicc_add_benchmark(Implicit
    MODEL "BoussinesqSphereDynamo"
    TYPE "Config"
    WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
    ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
    TIMEOUT 500
  )

quicc_add_benchmark(Implicit
    MODEL "BoussinesqSphereDynamo"
    TYPE "State"
    WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
    ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
    TIMEOUT 500
  )

quicc_add_benchmark(Implicit
    MODEL "BoussinesqSphereDynamo"
    TYPE "Visualization"
    WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
    ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
    TIMEOUT 500
  )

quicc_add_benchmark(Implicit
    MODEL "BoussinesqSphereDynamo"
    TYPE "Model"
    WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
    ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
    TIMEOUT 500
  )

if(QUICC_USE_MPI)
  quicc_add_benchmark(Implicit
    MODEL "BoussinesqSphereDynamo"
    TYPE "Model"
    WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
    ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
    VARIANTS "framework/parallel/algorithm:single2d"
    TIMEOUT 500
  )
endif()
