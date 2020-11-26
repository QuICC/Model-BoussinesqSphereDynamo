# Set Model name and library name
get_filename_component(modName ${CMAKE_CURRENT_SOURCE_DIR} NAME)
string(TOLOWER "quicc_${modName}" libName)
set(modImpls "Explicit" "Implicit")

# Setup generic model library
quicc_add_library(${libName} "")
set_target_properties(${libName} PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(${libName} PUBLIC include/)
target_link_libraries(${libName}
  quicc_framework
)
# Create model implementations libraries
foreach(impl ${modImpls})
  string(TOLOWER "${libName}_${impl}" modLib)
  quicc_add_library(${modLib} "")
  set_target_properties(${modLib} PROPERTIES LINKER_LANGUAGE CXX)
  target_include_directories(${modLib} PUBLIC include/)
  target_link_libraries(${modLib}
    ${libName}
  )
endforeach()

# Update python files
add_custom_target(${libName}_updatepy)
add_custom_command(TARGET ${libName}_updatepy POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_directory
   "${CMAKE_CURRENT_SOURCE_DIR}/Python"
   "${CMAKE_BINARY_DIR}/Python"
   COMMENT "Copying Python files for ${modName}"
  VERBATIM
)
add_dependencies(${libName} ${libName}_updatepy)

add_subdirectory(src)

# Create model executables
foreach(impl ${modImpls})
  set(ModelId
    "Boussinesq/Sphere/Dynamo/${impl}")
  set(ModelFile
    "include/QuICC/Model/${ModelId}/PhysicalModel.hpp")
  string(TOLOWER "${libName}_${impl}" modLib)
  quicc_create_config_exe(${ModelId} "${modLib}")
  quicc_create_state_exe(${ModelId} "${modLib}")
  quicc_create_model_exe(${ModelId} "${modLib}")
  quicc_create_visu_exe(${ModelId} "${modLib}")
endforeach()
unset(libName)
unset(modName)
unset(modImpls)
unset(modLib)
