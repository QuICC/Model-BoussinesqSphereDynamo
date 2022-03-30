include(ConfigureModel)
quicc_add_model(Boussinesq/Sphere/Dynamo
  MODIMPLS Explicit Implicit
  BROKENMPI Tubular Single1D)
