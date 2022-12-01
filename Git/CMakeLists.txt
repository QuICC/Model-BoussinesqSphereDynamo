include(gitUtils/AddGitHashLib)
AddGitHashLib("Model/${target}"
  PARENT ${QUICC_CURRENT_MODEL_LIB}
)
