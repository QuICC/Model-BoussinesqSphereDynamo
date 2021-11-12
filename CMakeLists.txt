# Set Model name and library name
get_filename_component(modName ${CMAKE_CURRENT_SOURCE_DIR} NAME)
string(TOLOWER "quicc_${modName}" QUICC_CURRENT_MODEL_LIB)
set(modImpls "Explicit" "Implicit")

set(_broken "Tubular" "Single1D")
list(FIND _broken ${QUICC_MPIALGO} _pos)
if(NOT _pos EQUAL -1)
  list(REMOVE_ITEM modImpls "Implicit")
endif()

function(quicc_model_id impl id)
  set(ModelId
    "Boussinesq/Sphere/Dynamo/${impl}")
  set(${id} "${ModelId}" PARENT_SCOPE)
endfunction()

# Set library visibility
set(QUICC_CMAKE_SRC_VISIBILITY PRIVATE)

# Setup generic model library
add_library(${QUICC_CURRENT_MODEL_LIB} "")
set_target_properties(${QUICC_CURRENT_MODEL_LIB} PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(${QUICC_CURRENT_MODEL_LIB} PUBLIC
  include/
  )
target_link_libraries(${QUICC_CURRENT_MODEL_LIB} PUBLIC
  quicc_framework
  )

# Create model implementations libraries
foreach(impl ${modImpls})
  string(TOLOWER "${QUICC_CURRENT_MODEL_LIB}_${impl}" modLib)
  add_library(${modLib} "")
  set_target_properties(${modLib} PROPERTIES LINKER_LANGUAGE CXX)
  target_include_directories(${modLib} PUBLIC
    include/
    )
  target_link_libraries(${modLib} PUBLIC
    ${QUICC_CURRENT_MODEL_LIB}
    )
endforeach()

# Update python files
add_custom_target(${QUICC_CURRENT_MODEL_LIB}_updatepy)
add_custom_command(TARGET ${QUICC_CURRENT_MODEL_LIB}_updatepy POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_directory
  "${CMAKE_CURRENT_SOURCE_DIR}/Python"
  "${CMAKE_BINARY_DIR}/Python"
  COMMENT "Copying Python files for ${modName}"
  VERBATIM
  )
add_dependencies(${QUICC_CURRENT_MODEL_LIB} ${QUICC_CURRENT_MODEL_LIB}_updatepy)

add_subdirectory(src)

# Create model executables
foreach(impl ${modImpls})
  quicc_model_id(${impl} ModelId)
  string(TOLOWER "${QUICC_CURRENT_MODEL_LIB}_${impl}" modLib)
  quicc_create_config_exe(${ModelId} "${modLib}")
  quicc_create_state_exe(${ModelId} "${modLib}")
  quicc_create_model_exe(${ModelId} "${modLib}")
  quicc_create_visu_exe(${ModelId} "${modLib}")
endforeach()
unset(QUICC_CURRENT_MODEL_LIB)
unset(modName)
unset(modImpls)
unset(modLib)
