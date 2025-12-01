if(TARGET ${QUICC_CURRENT_MODEL_LIB}_implicit)
  target_sources(${QUICC_CURRENT_MODEL_LIB}_implicit ${QUICC_CMAKE_SRC_VISIBILITY}
    ModelBackend.cpp
    )
endif()
