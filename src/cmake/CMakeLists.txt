include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

write_basic_package_version_file(
    PyHMMERConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
    PyHMMERConfig.cmake.in
    PyHMMERConfig.cmake
    INSTALL_DESTINATION ${PYHMMER_INSTALL_LIBS_DIR}/cmake/PyHMMER
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/PyHMMERConfigVersion.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/PyHMMERConfig.cmake
    DESTINATION ${PYHMMER_INSTALL_LIBS_DIR}/cmake/PyHMMER
)
