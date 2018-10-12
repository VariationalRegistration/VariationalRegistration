# the top-level README is used for describing this module, just
# re-used it for documentation here
get_filename_component(MY_CURRENT_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
file(READ "${MY_CURRENT_DIR}/README.rst" DOCUMENTATION)

# itk_module() defines the module dependencies in VariationalRegistration
# The testing module in VariationalRegistration depends on ITKTestKernel
# By convention those modules outside of ITK are not prefixed with
# ITK.

# define the dependencies of the include module and the tests
itk_module(VariationalRegistration
  DEPENDS
    ITKCommon
    ITKIOImageBase
    ITKImageFilterBase
    ITKSmoothing
    ITKFFT
    ITKFiniteDifference
    ITKDisplacementField
    ITKRegistrationCommon
    ITKMathematicalMorphology
    ITKBinaryMathematicalMorphology
    ITKTestKernel #necessary to handle IO in src
  TEST_DEPENDS
    ITKTestKernel
 EXCLUDE_FROM_DEFAULT
 DESCRIPTION
   "${DOCUMENTATION}"
)


