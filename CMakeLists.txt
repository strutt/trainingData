# CMakeLists.txt for event package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 2.8.10)
project (decimatedDistributions)

set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_COLOR_MAKEFILE ON)





# If ROOTSYS is defined in the usual way then includes are at $ROOTSYS/include, same for libs.
if(DEFINED $ENV{ROOTSYS})
  list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
# If not then we can access the variables via the ROOT CMake stuff directly
else()
  execute_process(COMMAND root-config --etcdir OUTPUT_VARIABLE ROOT_ETC_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${ROOT_ETC_DIR}/cmake)
endif()

if(DEFINED ENV{ANITA_UTIL_INSTALL_DIR})
    message("ANITA_UTIL_INSTALL_DIR is set to $ENV{ANITA_UTIL_INSTALL_DIR}")
    set(UTIL_LIB_DIR $ENV{ANITA_UTIL_INSTALL_DIR}/lib)
    set(UTIL_INC_DIR $ENV{ANITA_UTIL_INSTALL_DIR}/include)
    set(UTIL_BIN_DIR $ENV{ANITA_UTIL_INSTALL_DIR}/bin)
    set(UTIL_SHARE_DIR $ENV{ANITA_UTIL_INSTALL_DIR}/share)
    set(LD_UTIL $ENV{ANITA_UTIL_INSTALL_DIR}/lib)
    set(INC_UTIL $ENV{ANITA_UTIL_INSTALL_DIR}/include)
elseif(DEFINED ENV{ARA_UTIL_INSTALL_DIR})
    message("ARA_UTIL_INSTALL_DIR is set to $ENV{ANITA_UTIL_INSTALL_DIR}")
    set(UTIL_LIB_DIR $ENV{ARA_UTIL_INSTALL_DIR}/lib)
    set(UTIL_INC_DIR $ENV{ARA_UTIL_INSTALL_DIR}/include)
    set(UTIL_BIN_DIR $ENV{ARA_UTIL_INSTALL_DIR}/bin)
    set(UTIL_SHARE_DIR $ENV{ARA_UTIL_INSTALL_DIR}/share)
    set(LD_UTIL $ENV{ARA_UTIL_INSTALL_DIR}/lib)
    set(INC_UTIL $ENV{ARA_UTIL_INSTALL_DIR}/include)
else()
    message("Will install to /usr/local  ANITA_UTIL_INSTALL_DIR not set")
    set(UTIL_LIB_DIR /usr/local/lib)
    set(UTIL_INC_DIR /usr/local/include)
    set(UTIL_BIN_DIR /usr/local/bin)
    set(UTIL_SHARE_DIR /usr/local/share)
    set(LD_UTIL /usr/local/lib)
    set(INC_UTIL /usr/local/include)
endif()


#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${UTIL_SHARE_DIR}/cmake/modules/")

find_package(ROOT REQUIRED COMPONENTS MathMore MathCore RIO Hist Tree Net Minuit Minuit2 TreePlayer)
message("ROOT_INCLUDE_DIRS is set to ${ROOT_INCLUDE_DIRS}")

include_directories(${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/include ${ZLIB_INCLUDES} ${UTIL_INC_DIR})
add_definitions(${ROOT_CXX_FLAGS} -DFFTW_64_BIT)
link_directories(${UTIL_LIB_DIR} ./)


set(libname "AnitaCluster")
set(DICTNAME G__${libname})
set(SOURCE_FILES AnitaClusterer.cxx BaseList.cxx)
set(HEADER_FILES AnitaClusterer.h BaseList.h)
ROOT_GENERATE_DICTIONARY(${DICTNAME} ${HEADER_FILES} LINKDEF LinkDef.h)

##---Create a shared library with geneated dictionary
add_library(${libname} SHARED ${SOURCE_FILES} ${DICTNAME}.cxx)
#target_link_libraries(${libname} ${ZLIB_LIBRARIES} AnitaEvent AnitaAnalysisTools RootFftwWrapper ${ROOT_LIBRARIES} ${FFTW_LIBRARIES})
target_link_libraries(${libname} ${ZLIB_LIBRARIES} AnitaEvent AnitaCorrelator AnitaAnalysisTools ${ROOT_LIBRARIES})
#target_link_libraries(${libname} ${ZLIB_LIBRARIES})

#set(ANITA_LIBS AnitaEvent AnitaAnalysis AnitaAnalysisTools AnitaCorrelator AnalysisCuts)
set(ANITA_LIBS AnitaEvent AnitaAnalysis AnitaAnalysisTools AnitaCorrelator RootFftwWrapper)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wshadow -Werror")


# removed powerSpectra imagePeakHilbertPeak getWaisHeadings plotInitialDistributions project plotMinBiasDistributions makeDataQualityTrees makeMinBiasDataQualityTrees plotDataQuality templateSelfTriggeredBlast generateAngularResolutionTree reconstructPeaks waisDistributions)
set(BINARIES reconstructDecimated reconstructMinBias reconstructWais reconstructMonteCarlo plotReconstructedWais plotReconstructedMinBias plotReconstructedDecimated makeIndexedHeaderTrees examineCoherentlySummedWaveforms makeDecimatedDataQualityTrees makeSlimDecimatedDataQualityTrees makeWaisDataQualityTrees makeMinBiasDataQualityTrees makeThermalBackgroundTree makeSignalTree projection projectionDecimated dummyProjection getRidOfTheseEvents drawProjectionPlots clusteringBeta plotDataQuality slimTree)

FOREACH(binary ${BINARIES})
  MESSAGE(STATUS "Process file: ${binary}")
  add_executable(${binary} ${binary}.cxx)
  target_link_libraries(${binary} ${ZLIB_LIBRARIES} ${ANITA_LIBS} ${ROOT_LIBRARIES} ${FFTW_LIBRARIES} ${libname})
ENDFOREACH(binary)
