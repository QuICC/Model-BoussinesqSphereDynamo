list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake.d")

# Set Model directory name
get_filename_component(QUICC_CURRENT_MODEL_DIRNAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(QUICC_CURRENT_MODEL_DIR "Models/${QUICC_CURRENT_MODEL_DIRNAME}")
string(TOLOWER "quicc_${QUICC_CURRENT_MODEL_DIRNAME}" QUICC_CURRENT_MODEL_LIB)

# Model library
include(ConfigureModel)
quicc_add_model(Boussinesq/Sphere/Dynamo
  TYPES Explicit Implicit
  MODEL_DIRNAME ${QUICC_CURRENT_MODEL_DIRNAME})

add_subdirectory(Executables)

include (setup/TestSuite)
