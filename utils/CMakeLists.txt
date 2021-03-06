find_package(ROOT REQUIRED COMPONENTS GenVector MathCore)

#ROOT_GENERATE_DICTIONARY(G__NPDetGeoCad
#  include/TGeoToStep.h
#  OPTIONS -I${OpenCASCADE_INCLUDE_DIR} 
#  OPTIONS -I${CMAKE_CURRENT_SOURCE_DIR} 
#  OPTIONS -I${CMAKE_CURRENT_SOURCE_DIR}/include 
#  OPTIONS -I${CMAKE_CURRENT_SOURCE_DIR}/src 
#  LINKDEF include/LinkDef.h
#  )
#add_custom_target(G__NPDetGeoCad_ROOTDICT DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/G__NPDetGeoCad.cxx)

add_library(eicd_utils SHARED
  src/utils.cxx
  )

target_compile_features(eicd_utils
    PUBLIC cxx_auto_type
    PUBLIC cxx_trailing_return_types
    PUBLIC cxx_std_17
    PRIVATE cxx_variadic_templates
    )

target_compile_options(eicd_utils PRIVATE 
  -Wno-extra
  -Wno-ignored-qualifiers
  -Wno-overloaded-virtual
  -Wno-shadow)

target_link_libraries(eicd_utils 
  PUBLIC eicd
  PUBLIC ROOT::GenVector ROOT::MathCore)

target_include_directories(eicd_utils 
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  PUBLIC $<INSTALL_INTERFACE:include>
  )

install(FILES 
  include/eicd/analysis_utils.h
  include/eicd/vector_utils.h
  include/eicd/vector_utils_legacy.h
  DESTINATION include/eicd
  )
#install(FILES
#  "${CMAKE_CURRENT_BINARY_DIR}/libNPDetGeoCad.rootmap"
#  "${CMAKE_CURRENT_BINARY_DIR}/libNPDetGeoCad_rdict.pcm"
#  DESTINATION lib)


install(TARGETS eicd_utils
  EXPORT ${PROJECT_NAME}Targets
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES DESTINATION include
  )
