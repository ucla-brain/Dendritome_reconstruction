# set gtest related variables
set(GTEST_PATH ${3RD_PARTY_PATH}/gtest)
set(GTEST_INCLUDE_PATH
    ${GTEST_PATH}/googletest/include
    ${GTEST_PATH}/googlemock/include)
foreach(gtest_include_path ${GTEST_INCLUDE_PATH})
    if(NOT EXISTS ${gtest_include_path})
        message(FATAL_ERROR "${gtest_include_path} does not exist. ensure gtest is built")
    endif()
endforeach()
set(GTEST_LIBRARIES
    ${GTEST_PATH}/lib/libgmock.a
    ${GTEST_PATH}/lib/libgtest.a)
foreach(gtest_library ${GTEST_LIBRARIES})
    if(NOT EXISTS ${gtest_library})
        message(FATAL_ERROR "${gtest_library} does not exist. ensure gtest is built")
    endif()
endforeach()