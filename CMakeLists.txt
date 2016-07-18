# CMakeLists.txt for event package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 2.8.10)
project (decimatedDistributions)

set(CMAKE_DISABLE_SOURCE_CHANGES  ON)
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_COLOR_MAKEFILE ON)

if(DEFINED $ENV{ROOTSYS})
  list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
else()
  execute_process(COMMAND root-config --etcdir OUTPUT_VARIABLE ROOT_ETC_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
  #set(ROOT_ETC_DIR "${ROOT_ETC_DIR}/cmake")
  #message("is this root-config --etcdir? ${ROOT_ETC_DIR}")
endif()
#message("CMAKE_PREFIX_PATH = ${CMAKE_PREFIX_PATH}")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${ROOT_ETC_DIR}/cmake)

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



find_package(ROOT REQUIRED COMPONENTS MathMore MathCore RIO Hist Tree Net Minuit Minuit2)
#find_package(FFTW REQUIRED)
message("ROOT_INCLUDE_DIRS is set to ${ROOT_INCLUDE_DIRS}")

include_directories(${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/include ${FFTW_INCLUDES} ${ZLIB_INCLUDES} ${UTIL_INC_DIR})
add_definitions(${ROOT_CXX_FLAGS})
link_directories(${UTIL_LIB_DIR})

set(ANITA_LIBS AnitaEvent AnitaAnalysis AnitaAnalysisTools AnitaCorrelator RootFftwWrapper)

# removed imagePeakHilbertPeak getWaisHeadings
set(BINARIES powerSpectra initialDistributions minBiasDistributions plotInitialDistributions project plotMinBiasDistributions makeDataQualityTrees makeMinBiasDataQualityTrees plotDataQuality templateSelfTriggeredBlast generateAngularResolutionTree reconstructPeaks waisDistributions)

FOREACH(binary ${BINARIES})
  MESSAGE(STATUS "Process file: ${binary}")
  add_executable(${binary} ${binary}.cxx)
  target_link_libraries(${binary} ${ZLIB_LIBRARIES} ${ANITA_LIBS} ${ROOT_LIBRARIES} ${FFTW_LIBRARIES})  
ENDFOREACH(binary)



