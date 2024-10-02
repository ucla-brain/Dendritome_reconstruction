# set eigen related variables
set(EIGEN_INCLUDE_PATH ${3RD_PARTY_PATH}/eigen)
if(NOT EXISTS ${EIGEN_INCLUDE_PATH})
    message(FATAL_ERROR "${EIGEN_INCLUDE_PATH} does not exist. ensure  eigen is built")
endif()